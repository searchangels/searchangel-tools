"""birthfamily.py - offline DNA-match analysis for unknown-parentage cases. Stdlib only.
Never contacts any testing company; it only reads files you already have on disk.

  python birthfamily.py start                 (what '1 Capture a case' runs: asks a case name, opens the guide)
  python birthfamily.py menu                  (what '2 Work the case' runs: numbered menu of the commands below)
  python birthfamily.py guide CASE_FOLDER [--floor 90]   (step-by-step window: the normal way to capture a case)
  python birthfamily.py watch [CASE_FOLDER]        (files every page you copy; the easy way)
  python birthfamily.py paste clip matches.csv                    (after Ctrl+A, Ctrl+C on the match list)
  python birthfamily.py paste clip shared.csv --with "Jane Doe"   (same, on Jane Doe's shared-matches page)
       ('clip' reads the clipboard; a .txt filename works too. Rows are appended.)
  python birthfamily.py cluster matches.csv shared.csv [--trees DIR] [--min 40] [--max 400]
  python birthfamily.py find CASE_FOLDER Surname [Surname ...]   (does a spouse's family appear in the matches?)
  python birthfamily.py adoptee CASE_FOLDER 1965 California   (birth year and place: turns on ages and generation levels)
  python birthfamily.py skeleton CASE_FOLDER Adoptee Name   (starting tree: adoptee, blank parents, closest match per side)
  python birthfamily.py walk CASE_FOLDER      (step down the generations from a common-ancestor couple)
  python birthfamily.py check 497 1C1R        (does this cM fit the relationship a tree claims?)
  python birthfamily.py odds case.txt [--trials 1000]
  python birthfamily.py gedcom case.txt HYPOTHESIS out.ged
  python birthfamily.py selftest

cluster inputs
  matches.csv  columns: name,cM[,tree]   (tree = file inside --trees DIR: a .ged, or a .txt of a pasted tree page)
  shared.csv   columns: a,b              (one row per pair of matches who also match each other)

odds input (case.txt) - the adoptee is always called TARGET, names have no spaces:
  Ann = Gpa Gma                      known pedigree: child = parent parent
  match Cara 850                     cM that TARGET shares with Cara
  hyp child_of_Bob: TARGET = Bob U1  a hypothesis; ';' separates extra lines, new names = unknown people
"""
import csv, math, os, random, re, sys
from collections import Counter

# Approximate sex-averaged genetic lengths (cM), chr1..22. Calibration knob: testing companies
# report slightly different totals (Ancestry parent/child is ~3,475 cM).
CHROMS = [267, 251, 219, 203, 198, 187, 181, 164, 160, 171, 153,
          163, 124, 115, 118, 128, 126, 118, 104, 104, 61, 70]


# ---------- cluster ----------

def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return [{k.strip().lower(): (v or "").strip() for k, v in row.items() if k} for row in csv.DictReader(f)]


def cluster(cm, pairs, lo, hi):
    # ponytail: connected components. Matches above `hi` are left out because close relatives
    # bridge every cluster (Leeds method). Swap in real community detection if clusters smear together.
    keep = {n for n, c in cm.items() if lo <= c <= hi}
    adj = {n: set() for n in keep}
    for a, b in pairs:
        if a in keep and b in keep and a != b:
            adj[a].add(b); adj[b].add(a)
    seen, out = set(), []
    for n in sorted(keep, key=cm.get, reverse=True):
        if n in seen:
            continue
        comp, stack = [], [n]
        seen.add(n)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x] - seen:
                seen.add(y); stack.append(y)
        out.append(sorted(comp, key=cm.get, reverse=True))
    return out


def gedcom_people(path):
    """Set of (first given name, surname, birth year) for everyone in a GEDCOM."""
    # ponytail: exact name + exact year match. Add fuzzy spelling / +-2 years when real trees miss overlaps.
    people, name, in_birt = set(), None, False
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            p = line.strip().split(" ", 2)
            if len(p) < 2:
                continue
            lvl, tag, val = p[0], p[1], p[2] if len(p) > 2 else ""
            if lvl == "0":
                name = None
            if lvl == "1":
                in_birt = tag == "BIRT"
                if tag == "NAME":
                    m = re.match(r"\s*([^\s/]+)[^/]*/([^/]+)/", val)
                    name = (m.group(1).lower(), m.group(2).strip().lower()) if m else None
            if lvl == "2" and tag == "DATE" and in_birt and name:
                y = re.search(r"\b(1[5-9]\d\d|20\d\d)\b", val)
                if y:
                    people.add(name + (y.group(1),))
    return people


YEAR = re.compile(r"\b(1[5-9]\d\d|20[0-2]\d)\b")


def text_people(path):
    """Same as gedcom_people, for a tree page pasted as text: a name line followed by a line with a year."""
    # Layout checked against Ancestry's tree preview on a match page (Nov 2026), where the surname is
    # repeated as a link ("John Bernard KeelerKeeler", sometimes followed by " shared surname").
    people, prev = set(), ""
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            line = line.strip()
            y = YEAR.search(line)
            words = re.sub(r"\s*shared surname$", "", prev).split()
            dbl = words and re.fullmatch(r"([A-Z].*?)\1", words[-1])
            if dbl:
                words[-1] = dbl.group(1)
            if y and len(line) < 40 and len(words) >= 2 and not re.search(r"\d", prev):
                people.add((words[0].lower(), words[-1].lower(), y.group(1)))
            if line:
                prev = line
    return people


CM = re.compile(r"([\d,]+(?:\.\d+)?)\s*cM\b", re.I)
NOISE = re.compile(r"\b(cousins?|siblings?|sister|brother|parent|child|grand\w*|aunt|uncle|niece|nephew|family|"
                   r"relative|twin|side|managed by|trees?|people|shared|ancestors?|segments?|notes?|message|"
                   r"unviewed|new|groups?|view|learn)\b|^\W*$|\d", re.I)


def parse_matches(text):
    """[(name, cM, relationship, side, tree_info)] from a match-list page copied with Ctrl+A / Ctrl+C.
    Ancestry layout (checked against a real paste, Nov 2026):
      name / [Managed by ...] / blank / relationship / side / 'NNN cM | N% shared DNA' / tree lines / Add"""
    lines = [l.strip() for l in text.splitlines()]
    out = []
    for i, l in enumerate(lines):
        m = CM.search(l)
        if not m:
            continue
        blank = next((j for j in range(i - 1, max(i - 5, 0), -1) if not lines[j]), None)
        if blank:  # anchored on the blank line under the name
            j = blank - 1
            if lines[j].lower().startswith("managed by") and j > 0:
                j -= 1
            name, between = lines[j], lines[blank + 1:i]
        else:
            # ponytail: fallback for other page layouts - nearest line above that doesn't look like a label.
            # Skips names containing label words or digits; compare the count with the screen.
            j = next((j for j in range(i - 1, max(i - 8, -1), -1) if not NOISE.search(lines[j])), None)
            if j is None:
                continue
            name, between = lines[j], []
        rel = next((b for b in between if "side" not in b.lower()), "")
        side = next((b for b in between if "side" in b.lower()), "")
        info = []
        for t in lines[i + 1:i + 6]:
            if t in ("Add", "More options", ""):
                break
            if t not in info:
                info.append(t)
        if not out or out[-1][0] != name:
            out.append((name, float(m.group(1).replace(",", "")), rel, side, "; ".join(info)))
    return out


def save_rows(text, out, with_name):
    """Append parsed matches (or shared pairs) to a CSV, skipping rows already there. Returns new rows."""
    found = parse_matches(text)
    old = set()
    if os.path.exists(out):
        old = {(r["a"], r["b"]) if with_name else r["name"] for r in read_csv(out)}
    rows = [((with_name, r[0]) if with_name else (r[0], r[1], "", r[2], r[3], r[4])) for r in found]
    rows = [r for r in rows if (r if with_name else r[0]) not in old]
    with open(out, "a", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        if not os.path.getsize(out):
            w.writerow(["a", "b"] if with_name else ["name", "cM", "tree", "relationship", "side", "tree_info"])
        w.writerows(rows)
    return rows


def cmd_paste(src, out, with_name):
    if src == "clip":
        import tkinter
        tk = tkinter.Tk(); tk.withdraw()
        text = tk.clipboard_get(); tk.destroy()
    else:
        with open(src, encoding="utf-8-sig", errors="replace") as f:
            text = f.read()
    rows = save_rows(text, out, with_name)
    print(f"added {len(rows)} new rows to {out}. Check this count against the screen. First few:")
    for r in rows[:5]:
        print("  " + "  ".join(str(x) for x in r[:2]))


def cmd_watch(folder):
    """Sit in the background; every page the volunteer copies (Ctrl+A, Ctrl+C) is filed automatically.
    Only reads the local clipboard. Never touches any website."""
    import time, tkinter
    pastes = os.path.join(folder, "pastes")
    os.makedirs(pastes, exist_ok=True)
    tk = tkinter.Tk(); tk.withdraw()

    def clip():
        try:
            tk.update()
            return tk.clipboard_get()
        except tkinter.TclError:
            return None

    last, n = clip(), len(os.listdir(pastes))
    print(f"Watching the clipboard. On each Ancestry page press Ctrl+A then Ctrl+C. Files go to {folder}. Ctrl+C here to stop.", flush=True)
    while True:
        time.sleep(0.5)
        text = clip()
        if not text or text == last or len(text) < 200:
            continue
        last, n = text, n + 1
        raw = os.path.join(pastes, f"{n:03d}.txt")
        with open(raw, "w", encoding="utf-8") as f:
            f.write(text)
        # ponytail: only the main match list is recognised so far ("All matches" tab bar). Shared-match
        # and tree pages are kept as raw files until real samples show how to tell them apart.
        if "All matches" in text and CM.search(text):
            rows = save_rows(text, os.path.join(folder, "matches.csv"), None)
            print(f"match list: {len(rows)} new matches added to matches.csv", flush=True)
        else:
            print(f"saved {raw} (page type not recognised yet)", flush=True)


def make_groups(rows, pairs, lo, hi):
    """(close matches, [(label, [names])]) - shared-match clusters, or parent sides when there are no links."""
    cm = {r["name"]: float(r["cm"]) for r in rows}
    close = sorted((n for n in cm if cm[n] > hi), key=cm.get, reverse=True)
    if any(a != b for a, b in pairs):
        raw = cluster(cm, pairs, lo, hi)
        comps = [c for c in raw if len(c) > 1]
        # People linked only through a close match (left out of clustering above) would end up alone.
        # Group them under the closest such match instead. ponytail: each goes to one close match only.
        loose = [c[0] for c in raw if len(c) == 1]
        link = {frozenset(p) for p in pairs}
        named = []
        for c in close:
            got = [n for n in loose if frozenset((c, n)) in link]
            if got:
                named.append((f"Shared with close match {c}", got))
                loose = [n for n in loose if n not in got]
        comps = [(f"Cluster {i}", c) for i, c in enumerate(comps, 1)] + named
        if loose:
            comps.append(("Not grouped (no shared-match link captured)", loose))
    else:  # no shared-match data: fall back to Ancestry's own parent-side split
        sides = {}
        for r in sorted(rows, key=lambda r: -float(r["cm"])):
            if lo <= float(r["cm"]) <= hi:
                sides.setdefault(r.get("side") or "Unassigned", []).append(r["name"])
        comps = sorted(sides.items())
    return close, comps


def surname_hits(folder, sur):
    """[(match row, [people in their tree with that surname])] for matches carrying or listing a surname."""
    # ponytail: plain spelling match. Pass variants yourself (Smith Smyth) until fuzzy matching is added.
    pat = re.compile(rf"\b{re.escape(sur)}", re.I)
    hits = []
    for r in sorted(read_csv(os.path.join(folder, "matches.csv")), key=lambda r: -float(r["cm"])):
        people = []
        if r.get("tree") and os.path.exists(os.path.join(folder, "trees", r["tree"])):
            with open(os.path.join(folder, "trees", r["tree"]), encoding="utf-8-sig", errors="replace") as f:
                for line in f:
                    line = re.sub(r"\s*shared surname$", "", line.strip())
                    line = re.sub(r"([A-Z][\w'-]+)\1$", r"\1", line)   # "KeelerKeeler" -> "Keeler"
                    if pat.search(line) and not YEAR.search(line) and line not in people:
                        people.append(line)
        if pat.search(r["name"]) or people:
            hits.append((r, people))
    return hits


def cmd_find(folder, surnames):
    """The keystone check: does a spouse's family show up among the adoptee's matches?
    For each surname, list matches who carry it or have it in their collected tree, with side and group."""
    rows = read_csv(os.path.join(folder, "matches.csv"))
    s = os.path.join(folder, "shared.csv")
    pairs = [(r["a"], r["b"]) for r in read_csv(s)] if os.path.exists(s) else []
    close, comps = make_groups(rows, pairs, 40, 400)
    group = {n: label for label, comp in comps for n in comp}
    group.update({n: "close match" for n in close})
    for sur in surnames:
        hits = surname_hits(folder, sur)
        print(f"\n{sur}: {len(hits)} matches" + ("" if hits else "  (nothing in the data collected so far)"))
        for r, people in hits:
            where = "in their tree: " + "; ".join(people[:3]) if people else "it is the match's own name"
            print(f"  {r['name']}  {float(r['cm']):.0f} cM  [{r.get('side') or 'Unassigned'}, {group.get(r['name'], 'below cut-off')}]  {where}")
        sides = Counter(r.get("side") or "Unassigned" for r, _ in hits)
        if hits:
            print("  => " + ", ".join(f"{c} on {s}" for s, c in sides.most_common()))


GEN = {0: "the adoptee's own", 1: "parent", 2: "grandparent", 3: "great-grandparent",
       4: "2x-great-grandparent", 5: "3x-great-grandparent"}


def adoptee(folder):
    """Known facts about the adoptee for this case: {'birth_year': '1965', 'birth_place': 'California'}."""
    p = os.path.join(folder, "adoptee.txt")
    if not os.path.exists(p):
        return {}
    with open(p, encoding="utf-8") as f:
        return dict(l.strip().split("=", 1) for l in f if "=" in l)


def cmd_adoptee(folder, year, place):
    # Deliberately no Social Security number here: only the state it points to is useful, and a typed
    # state carries none of the risk.
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "adoptee.txt"), "w", encoding="utf-8") as f:
        f.write(f"birth_year={int(year)}\nbirth_place={place}\n")
    print(f"Saved: born {int(year)}" + (f" in {place}" if place else ""))


def generation(year, folder):
    """'great-grandparent level' for someone born in `year`, judged from the adoptee's birth year."""
    # ponytail: assumes 28 years per generation. Real gaps run 15 to 50, so this can be one level off.
    by = adoptee(folder).get("birth_year")
    if not by or not year:
        return ""
    n = round((int(by) - int(year)) / 28)
    return "younger than the adoptee" if n < 0 else GEN.get(n, f"{n} generations up") + " level"


def parent_age(year, folder):
    """Could someone born in `year` be the adoptee's birth parent, by age alone?"""
    by = adoptee(folder).get("birth_year")
    if not by or not year:
        return ""
    age = int(by) - int(year)
    if age < 14:
        return f"too young to be a birth parent (age {age} at the adoptee's birth)" if age >= 0 else "born after the adoptee"
    # Birth parents of adoptees were most often young. A tendency, not a rule, so it only changes the
    # wording and the order of candidates; nobody is ruled out by it. ponytail: 25 is a round-number cut.
    if age <= 25:
        return f"age {age} at the adoptee's birth: COULD BE A BIRTH PARENT, and the typical age for one"
    if age <= 55:
        return f"age {age} at the adoptee's birth: COULD BE A BIRTH PARENT, though older than most"
    return f"age {age} at the adoptee's birth: more likely a grandparent or earlier"


def records_hint(folder):
    a = adoptee(folder)
    if not a.get("birth_year"):
        return ("Tip: save the adoptee's birth year and place with the 'adoptee' command to get ages and "
                "generation levels here.")
    place = a.get("birth_place", "")
    out = f"Records to check next: births in {place or 'the birth state'} in {a['birth_year']}, and the candidates' marriage and death records."
    if place.strip().lower() in ("ca", "california") and 1905 <= int(a["birth_year"]) <= 1995:
        out += ("\nCalifornia: the California Birth Index (1905-1995) lists each child's name, birth date, county "
                "and mother's maiden name. Search it for the birth date and for the candidate mother's maiden name.")
    return out


def cmd_walk(folder):
    """The pendulum: start at a common-ancestor couple, list their children and who each married, see
    which spouse's family also shows up in the adoptee's matches, move the focus to that child, repeat."""
    us = lambda s: "_".join(s.split())
    evidence = lambda sur: (lambda h: (sum(float(r["cm"]) for r, _ in h), len(h),
                                       Counter(r.get("side") or "Unassigned" for r, _ in h).most_common(1)[0][0] if h else ""))(surname_hits(folder, sur))
    first = [x.strip() for x in input("Starting couple, as: husband, wife, birth year of either (year optional)\n> ").split(",")]
    if len(first) < 2:
        sys.exit("need two names separated by a comma")
    couple, year = [us(first[0]), us(first[1])], first[2] if len(first) > 2 and first[2].isdigit() else ""
    line_side = evidence(couple[0].split("_")[-1])[2] or evidence(couple[1].split("_")[-1])[2]
    print(f"This line's matches are mostly on: {line_side or 'unknown (no matches carry these surnames)'}")
    # one tree file per parent's side, so the two family trees stay separate until they join
    tree_path = os.path.join(folder, "tree_" + (re.sub(r"\W+", "_", line_side.lower()).strip("_") or "unknown") + ".txt")
    while True:
        gen = generation(year, folder)
        print(f"\nFOCUS: {couple[0]} + {couple[1]}" + (f"   (born about {year}: {gen})" if gen else ""))
        print(f"  Where to find their children: public trees hide anyone still living as 'Private', so the\n"
              f"  younger generations are usually missing. Look up the obituaries of {couple[0].replace('_', ' ')} and\n"
              f"  {couple[1].replace('_', ' ')} (newspapers.com, Find a Grave): they normally list surviving children,\n"
              f"  the children's spouses, and who died before them.")
        print("Type each child, who they married and the child's birth year, as: child, spouse, year\n"
              "(spouse and year can be blank). Empty line when done.")
        kids, born = [], {}
        while True:
            line = input("> ").strip()
            if not line:
                break
            parts = [p.strip() for p in line.split(",")] + ["", ""]
            kids.append((us(parts[0]), us(parts[1])))
            born[us(parts[0])] = parts[2] if parts[2].isdigit() else ""
        if not kids:
            break
        # ponytail: the tree file lists the focus couple as father, mother in the order typed; a married-in
        # husband typed second would be exported with the wrong sex. Ask for sex if GEDCOM export needs it.
        with open(tree_path, "a", encoding="utf-8") as f:
            f.writelines(f"{c} = {couple[0]} {couple[1]}\n" for c, _ in kids)
        scored = sorted(((evidence(s.split("_")[-1]) if s else (0.0, 0, "")) + (c, s) for c, s in kids), reverse=True)
        for total, n, side, c, s in scored:
            print(f"  {c} + {s or '(spouse unknown)'}: " + (f"{n} matches from the {s.split('_')[-1]} family, "
                  f"{total:.0f} cM in total, mostly {side}" if n else "no matches from the spouse's family")
                  + ("   [WEAK: rests on one match. Confirm with records before relying on it]" if n == 1
                     else "   [moderate: two matches]" if n == 2 else "   [strong]" if n else ""))
            if parent_age(born[c], folder):
                print(f"      born {born[c]}, {parent_age(born[c], folder)}")
        total, n, side, c, s = scored[0]
        if not n:
            print("\nNo spouse's family shows up in the matches. Either that side's trees were not collected, or the"
                  "\nbirth parent is in this generation: one of the children just listed.")
            able = sorted((k for k, _ in kids if "COULD BE" in parent_age(born[k], folder)), key=lambda k: -int(born[k]))
            if able:
                print("By age, the possible birth parents here are (youngest first, the more common pattern): " + ", ".join(able))
            print(records_hint(folder))
            break
        if line_side and side and side != line_side:
            print(f"\n{s.split('_')[-1]} matches are on the OTHER parent's side ({side}). That is the sign that"
                  f"\n{c} and {s} are the two birth parents: each brings one half of the adoptee's matches.")
            age = parent_age(born[c], folder)
            if age and "COULD BE" not in age:
                print(f"BUT the age does not fit: {c} was {age}. Re-check this step.")
            print(records_hint(folder))
            break
        ans = input(f"\nStrongest: {c} + {s}. Press Enter to move the focus there, or type another child's name: ").strip()
        pick = next(((kc, ks) for _, _, _, kc, ks in scored if kc == us(ans)), (c, s)) if ans else (c, s)
        couple, year = [pick[0], pick[1] or f"Unknown_spouse_of_{pick[0]}"], born[pick[0]]
    print(f"\nTree so far is saved in {tree_path}")


def level(cm):
    """Rough reading of how far up the shared ancestors are. Rule of thumb from typical cM averages
    (1C ~870, 2C ~230, 3C ~75); ranges overlap, so `check` and records decide, not this label."""
    if cm >= 1300:
        return "VERY CLOSE: half-sibling, aunt/uncle, grandparent or closer. One step from a birth parent"
    if cm >= 550:
        return ("first-cousin level. If near the adoptee's age: a first cousin, so their parent is a sibling of the "
                "birth parent. If a generation or two older: a great-aunt/uncle")
    if cm >= 200:
        return "great-grandparent level (second cousin, or first cousin once removed)"
    if cm >= 90:
        return "2x-great-grandparent level (second to third cousin)"
    return "distant"


def cmd_cluster(matches_csv, shared_csv, trees_dir, lo, hi):
    rows = read_csv(matches_csv)
    cm = {r["name"]: float(r["cm"]) for r in rows}
    tree = {r["name"]: r.get("tree") for r in rows if r.get("tree")}
    pairs = [(r["a"], r["b"]) for r in read_csv(shared_csv)] if os.path.exists(shared_csv) else []
    close, comps = make_groups(rows, pairs, lo, hi)
    side_of = {r["name"]: r.get("side") or "Unassigned" for r in rows}
    # split the leftover pile by side, then file every group under the parent side most of its members are on
    split = []
    for label, comp in comps:
        if label.startswith("Not grouped"):
            for s in sorted({side_of[n] for n in comp}):
                split.append((label, [n for n in comp if side_of[n] == s]))
        else:
            split.append((label, comp))
    majority = lambda comp: Counter(side_of[n] for n in comp).most_common(1)[0][0]
    sides = sorted({majority(c) for _, c in split} | {side_of[n] for n in close},
                   key=lambda s: ("parent" not in s.lower(), s))
    print("Each parent's side is its own family tree. Work them separately; they only join at the birth parents.")
    for side in sides:
        mine = [(label, comp) for label, comp in split if majority(comp) == side]
        everyone = sorted([n for n in close if side_of[n] == side] + [n for _, c in mine for n in c], key=cm.get, reverse=True)
        print(f"\n{'=' * 70}\n{side.upper()}: {len(everyone)} matches\n{'=' * 70}")
        if everyone:
            top = everyone[0]
            print(f"Closest match, start here: {top}  {cm[top]:.0f} cM\n  {level(cm[top])}.")
        for n in [n for n in close if side_of[n] == side]:
            print(f"  close match: {n}  {cm[n]:.0f} cM")
        families = [label for label, _ in mine if not label.startswith("Not grouped")]
        if len(families) > 1:
            print("Families on this side that should join by marriage on the way down: " + "  +  ".join(families))
        for label, comp in mine:
            print(f"\n{label}  ({len(comp)} matches)")
            for n in comp:
                print(f"  {n}  {cm[n]:.0f} cM")
            seen, surnames, ntrees = Counter(), Counter(), 0
            for n in comp:
                if trees_dir and n in tree and os.path.exists(os.path.join(trees_dir, tree[n])):
                    p = os.path.join(trees_dir, tree[n])
                    people = gedcom_people(p) if p.lower().endswith(".ged") else text_people(p)
                    seen.update(people)
                    surnames.update({sur for _, sur, _ in people})   # count each surname once per tree
                    ntrees += 1
            top = [(s, c) for s, c in surnames.most_common(8) if c >= 2]
            if top:
                print(f"  Surnames in 2+ of the {ntrees} trees collected: "
                      + ", ".join(f"{s.title()} ({c})" for s, c in top))
            # most recent first: the latest-born shared ancestor is the one to focus on and build down from
            common = sorted(((k, c) for k, c in seen.items() if c >= 2), key=lambda kc: (-int(kc[0][2]), -kc[1]))
            if common:
                print("  People in 2+ of these matches' trees, most recent first (start with the top one):")
                for (given, sur, year), c in common[:15]:
                    gen = generation(year, os.path.dirname(matches_csv))
                    print(f"    {given.title()} {sur.title()} b.{year}  in {c} trees" + (f"   ({gen})" if gen else ""))
    comps = split

    # Endogamy check: cousin marriages and close-knit communities make matches look closer than they are.
    # ponytail: three rule-of-thumb signs with round-number thresholds, not a statistical test. Tune on real
    # endogamous cases. Ancestry's longest-segment figure would be a better sign but isn't in the copied pages.
    signs = []
    both = sum("both" in s.lower() for s in side_of.values())
    if both / max(len(side_of), 1) >= 0.15:
        signs.append(f"{both} of {len(side_of)} matches are on BOTH parents' sides")
    for label, comp in comps:
        mix = Counter(side_of[n] for n in comp if "parent" in side_of[n].lower())
        if label.startswith("Cluster") and len(mix) > 1 and min(mix.values()) / sum(mix.values()) >= 0.25:
            signs.append(f"{label} mixes both parents' sides ({dict(mix)})")
    seg, pastes = [], os.path.join(os.path.dirname(matches_csv), "pastes")
    if os.path.isdir(pastes):
        for fn in os.listdir(pastes):
            with open(os.path.join(pastes, fn), encoding="utf-8", errors="replace") as f:
                m = re.search(r"shared DNA: ([\d,]+) cM across (\d+) segments", f.read())
            if m and float(m.group(1).replace(",", "")) <= hi:
                seg.append(float(m.group(1).replace(",", "")) / int(m.group(2)))
    if len(seg) >= 5:
        med = sorted(seg)[len(seg) // 2]
        if med < 12:
            signs.append(f"shared DNA is in small pieces (typically {med:.0f} cM per segment; unrelated-community cases run 15+)")
    print("\nEndogamy check (cousin marriages, isolated or close-knit communities):")
    if signs:
        for s in signs:
            print(f"  WARNING: {s}")
        print("  Matches are probably more distant than their cM suggests. Do not trust 'check', 'odds' or the\n"
              "  relationship labels at face value on this case; lean on records and the largest matches.")
    else:
        print("  No signs in the data collected"
              + (f" (typical segment {sorted(seg)[len(seg) // 2]:.0f} cM, {both} matches on both sides)." if len(seg) >= 5 else "."))


# ---------- guide (step-by-step window; reads only the local clipboard) ----------

class Guide:
    """Walks a volunteer through a case: match list pages -> each close match's page -> shared matches.
    The person opens every page; this only files what they copy."""

    HEAD = ["name", "cM", "tree", "relationship", "side", "tree_info"]

    def __init__(self, folder, floor=90.0, top=12):
        self.folder, self.floor, self.top = folder, floor, top
        for d in ("pastes", "trees"):
            os.makedirs(os.path.join(folder, d), exist_ok=True)
        self.m, self.s = os.path.join(folder, "matches.csv"), os.path.join(folder, "shared.csv")
        self.seen = os.path.join(folder, "pages_done.txt")
        self.stage, self.queue, self.total = "matches", [], 0
        self.n = len(os.listdir(os.path.join(folder, "pastes")))
        if self.rows() and min(float(r["cm"]) for r in self.rows()) < floor:
            self.next_stage()   # reopening a case whose match list is already complete

    def rows(self):
        rows = read_csv(self.m) if os.path.exists(self.m) else []
        return sorted(rows, key=lambda r: -float(r["cm"]))

    def done_pages(self):
        return set(open(self.seen, encoding="utf-8").read().splitlines()) if os.path.exists(self.seen) else set()

    @property
    def count(self):
        """(number, label) for the big counter in the window."""
        if self.stage == "matches":
            return len(self.rows()), "matches collected"
        if self.stage == "pages":
            return len(self.done_pages()), "match pages collected"
        return (len(read_csv(self.s)) if os.path.exists(self.s) else 0), "shared-match links collected"

    def where(self, name):
        rows = self.rows()
        i = next(k for k, r in enumerate(rows) if r["name"] == name)
        return f"{name}   ({float(rows[i]['cm']):.0f} cM)"

    @property
    def instruction(self):
        pos = f"({self.total - len(self.queue) + 1} of {self.total})"
        find = "\n\nCan't see the name? Use Search above the match list. Any order is fine."
        if self.stage == "matches":
            return ("STEP 1 of 3: MATCH LIST\n1. On Ancestry open DNA > Matches.\n"
                    f"2. Click Filter > Shared DNA > Custom centimorgan range. Type {self.floor:.0f} in Min, "
                    "leave Max empty, click Apply filters.\n3. Set 50 per page at the bottom of the list.\n"
                    "4. Press Ctrl+A, then Ctrl+C.\n5. Go to the next page and repeat until this window says the list is done.")
        if self.stage == "pages":
            return (f"STEP 2 of 3: PUBLIC TREES {pos}\nFirst, on the main match LIST (not a person's page; press Back if you are on one): "
                    "click Filter > Trees > Public linked trees > Apply filters. The list now shows only people with trees.\n\n"
                    f"1. Click the name:\n\n    {self.where(self.queue[0])}\n\n"
                    "2. Press Ctrl+A, then Ctrl+C.\n3. Press the browser Back button." + find)
        if self.stage == "shared":
            return (f"STEP 3 of 3: SHARED MATCHES {pos}\nFirst, on the main match LIST (not a person's page): click Filter > Trees > All > Apply filters.\n\n"
                    f"1. Click the name:\n\n    {self.where(self.queue[0])}\n\n"
                    "2. Click the \"Shared Matches\" tab on that page.\n3. Press Ctrl+A, then Ctrl+C.\n"
                    "4. Press Back until you see the match list." + find)
        return f"FINISHED\nResults are in:\n{os.path.join(self.folder, 'clusters.txt')}\nYou can close this window."

    def next_stage(self):
        rows = [r for r in self.rows() if float(r["cm"]) >= self.floor]
        if self.stage == "matches":
            # matches above the cut-off with a public tree of 8+ people, in list order (matches the
            # volunteer's screen once Filter > Trees > Public linked trees is applied)
            self.stage, done = "pages", self.done_pages()
            big = lambda r: "public linked tree" in r.get("tree_info", "").lower() and int(
                (re.search(r"([\d,]+) people", r["tree_info"]) or [0, "0"])[1].replace(",", "")) >= 8
            self.queue = [r["name"] for i, r in enumerate(rows) if r["name"] not in done and big(r)]
        elif self.stage == "pages":
            self.stage = "shared"
            done = {r["a"] for r in read_csv(self.s)} if os.path.exists(self.s) else set()
            self.queue = [r["name"] for r in rows[:self.top] if r["name"] not in done]
        elif self.stage == "shared":
            import contextlib
            self.stage = "done"
            with open(os.path.join(self.folder, "clusters.txt"), "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
                cmd_cluster(self.m, self.s, os.path.join(self.folder, "trees"), 40, 400)
        self.total = len(self.queue)
        if self.stage in ("pages", "shared") and not self.queue:
            self.next_stage()

    def skip(self, name=None):
        if self.queue:
            self.queue.remove(name or self.queue[0])
        if not self.queue and self.stage != "done":
            self.next_stage()

    def finish_step(self):
        if self.stage != "done":
            self.queue = []
            self.next_stage()

    def set_tree(self, name, fn):
        rows = self.rows()
        with open(self.m, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            w.writerow(self.HEAD)
            w.writerows([r["name"], r["cm"], fn if r["name"] == name else r.get("tree", ""),
                         r.get("relationship", ""), r.get("side", ""), r.get("tree_info", "")] for r in rows)

    def on_copy(self, text):
        """Handle one copied page. Returns (ok, message)."""
        self.n += 1
        with open(os.path.join(self.folder, "pastes", f"{self.n:03d}.txt"), "w", encoding="utf-8") as f:
            f.write(text)
        main_list = "All matches" in text and "By cluster" in text
        if self.stage == "matches":
            found = parse_matches(text)
            if not main_list or not found:
                return False, "That doesn't look like the DNA match list. Nothing was saved."
            new, low = save_rows(text, self.m, None), min(r[1] for r in found)
            msg = f"Received {len(new)} new matches ({len(self.rows())} total). Lowest on this page: {low:.0f} cM."
            # footer is "Items per page / 3 / of 3"; a list that fits on one page shows no "of N" at all
            pg = re.search(r"Items per page\n+([\d,]+)\n(?:of ([\d,]+)\n)?", text)
            if pg and pg.group(2) in (None, pg.group(1)):
                self.next_stage()
                return True, msg + "\nThat was the last page, so the match list is done."
            if low < self.floor:
                self.next_stage()
                return True, msg + f"\nThat is below the {self.floor:.0f} cM cut-off, so the match list is done."
            return True, msg + "\nReady for the next 50."
        if self.stage == "done":
            return False, "This case is finished."
        # a single match's page: header line is "You and <name>"
        who = re.search(r"^You and (.+)$", text, re.M)
        name = who and who.group(1).strip()
        if main_list:
            # a (filtered) match list copied mid-way: keep any new people and refresh the worklist
            new = save_rows(text, self.m, None)
            if self.stage == "pages":
                self.stage = "matches"
                self.next_stage()
            return True, f"Match list received ({len(new)} new people added). Now click the name shown above."
        if not name:
            return False, "That doesn't look like a match's page. Nothing was saved."
        if name not in self.queue:
            return False, f"That is {name}'s page, which isn't on the worklist for this step (or is already done)."
        trees_tab = "How are you and" in text
        if self.stage == "pages":
            if not trees_tab:
                return False, "That looks like a different tab. Use the page you land on right after clicking the name."
            # tree preview sits between "...'s Linked Tree43 People" and the "Surnames" heading
            m = re.search(r"Linked Tree[\d,]* People\n(.*?)\nSurnames\n", text, re.S)
            n = 0
            if m:
                fn = re.sub(r"[^\w.-]+", "_", name) + ".txt"
                path = os.path.join(self.folder, "trees", fn)
                with open(path, "w", encoding="utf-8") as f:
                    f.write(m.group(1))
                n = len(text_people(path))
                self.set_tree(name, fn)
            with open(self.seen, "a", encoding="utf-8") as f:
                f.write(name + "\n")
            self.skip(name)
            hidden = len(re.findall(r"^Private$", m.group(1), re.M)) if m else 0
            return True, (f"Received {name}: " + (f"tree preview with {n} dated people." if m else "no public tree on this page.")
                          + (f" {hidden} living people hidden as Private." if hidden else ""))
        # shared stage
        if trees_tab:
            return False, "That is the first tab. Click \"Shared Matches\" on that page, then copy again."
        # ponytail: layout-independent - any known match whose name appears on its own line counts as shared.
        # Unverified against a real Shared Matches tab; misses shared matches beyond the captured match list.
        lines = {l.strip() for l in text.splitlines()}
        pairs = [(name, r["name"]) for r in self.rows() if r["name"] in lines and r["name"] != name]
        new = not os.path.exists(self.s)
        with open(self.s, "a", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["a", "b"])
            w.writerows(pairs or [(name, name)])   # self-pair marks "done, none found"
        self.skip(name)
        return True, f"Received {len(pairs)} shared matches for {name}."


def cmd_guide(folder, floor):
    import tkinter
    g = Guide(folder, floor)
    root = tkinter.Tk()
    root.title("Birth family: case capture")
    root.attributes("-topmost", True)
    root.geometry("480x540")
    title = tkinter.Label(root, justify="left", anchor="w", wraplength=450, font=("Segoe UI", 15, "bold"))
    title.pack(fill="x", padx=14, pady=(12, 0))
    instr = tkinter.Label(root, justify="left", anchor="nw", wraplength=450, font=("Segoe UI", 11))
    instr.pack(fill="x", padx=14, pady=(4, 6))
    big = tkinter.Label(root, text="0", font=("Segoe UI", 44, "bold"), fg="#1a7f37")
    big.pack()
    biglabel = tkinter.Label(root, font=("Segoe UI", 10))
    biglabel.pack(pady=(0, 6))
    shown = [0]

    def tick():   # count up one at a time so the volunteer sees rows arriving
        target, label = g.count
        biglabel.config(text=label)
        if shown[0] != target:
            shown[0] = shown[0] + 1 if shown[0] < target else target
            big.config(text=str(shown[0]))
        root.after(40, tick)
    state = tkinter.Label(root, text="WAITING for a copy...", fg="white", bg="#555", font=("Segoe UI", 12, "bold"))
    state.pack(fill="x", padx=14)
    msg = tkinter.Label(root, justify="left", anchor="nw", wraplength=430, font=("Segoe UI", 10))
    msg.pack(fill="both", expand=True, padx=14, pady=6)

    def refresh():
        head, _, body = g.instruction.partition("\n")
        title.config(text=head)
        instr.config(text=body)

    def skip():
        g.skip(); refresh()
        msg.config(text="Skipped.")

    def finish():
        g.finish_step(); refresh()
        msg.config(text="Moved on to the next step.")

    bar = tkinter.Frame(root)
    bar.pack(pady=(0, 10))
    tkinter.Button(bar, text="Skip this person", command=skip).pack(side="left", padx=6)
    tkinter.Button(bar, text="Finish this step", command=finish).pack(side="left", padx=6)

    def clip():
        try:
            return root.clipboard_get()
        except tkinter.TclError:
            return None

    last = [clip()]

    def poll():
        text = clip()
        if text and text != last[0] and len(text) >= 100 and g.stage != "done":
            last[0] = text
            state.config(text="PROCESSING...", bg="#b36b00")
            root.update()
            ok, m = g.on_copy(text)
            msg.config(text=m)
            done = g.stage == "done"
            state.config(text="FINISHED" if done else ("RECEIVED. Waiting for the next copy..." if ok else "NOT RECOGNISED. Waiting..."),
                         bg="#1a7f37" if ok else "#b3261e")
            refresh()
        root.after(400, poll)

    refresh()
    tick()
    poll()
    root.mainloop()


# ---------- odds (pedigree simulation, same idea as WATO / BanyanDNA) ----------

def poisson(lam):
    k, p, limit = 0, 1.0, math.exp(-lam)
    while True:
        p *= random.random()
        if p <= limit:
            return k
        k += 1


def gamete(a, b, L):
    """One recombined chromosome from a parent's two copies. A copy = [(start_cM, founder_id), ...]."""
    # ponytail: no crossover interference, sex-averaged map. Fine for ranking hypotheses, not for forensics.
    cuts = sorted(random.uniform(0, L) for _ in range(poisson(L / 100))) + [L]
    src = (a, b) if random.random() < 0.5 else (b, a)
    out, lo = [], 0.0
    for i, hi in enumerate(cuts):
        h = src[i % 2]
        for j, (s, fid) in enumerate(h):
            e = h[j + 1][0] if j + 1 < len(h) else L
            if e > lo and s < hi and not (out and out[-1][1] == fid):
                out.append((max(s, lo), fid))
        lo = hi
    return out


def ibd(x, y, L):
    pts = sorted({s for s, _ in x} | {s for s, _ in y}) + [L]
    out, i, j = [], 0, 0
    for k in range(len(pts) - 1):
        s = pts[k]
        while i + 1 < len(x) and x[i + 1][0] <= s:
            i += 1
        while j + 1 < len(y) and y[j + 1][0] <= s:
            j += 1
        if x[i][1] == y[j][1]:
            out.append((s, pts[k + 1]))
    return out


def shared_cm(p, q, min_seg):
    total = 0.0
    for c, L in enumerate(CHROMS):
        ivs = sorted(iv for hp in p for hq in q for iv in ibd(hp[c], hq[c], L))
        cur_s = cur_e = None
        for s, e in ivs + [(float("inf"), 0)]:
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                if cur_e - cur_s >= min_seg:
                    total += cur_e - cur_s
                cur_s, cur_e = s, e
    return total


def simulate(ped, a, others, trials, min_seg=8.0):
    """ped: {child: (parent, parent)}. Returns {other: [shared cM with a, one per trial]}."""
    res = {o: [] for o in others}
    for _ in range(trials):
        genomes, nid = {}, [0]

        def get(n, depth=0):
            if n not in genomes:
                if depth > 60:
                    sys.exit(f"pedigree loop at {n}")
                if n in ped:
                    genomes[n] = tuple([gamete(g[0][c], g[1][c], L) for c, L in enumerate(CHROMS)]
                                       for g in (get(par, depth + 1) for par in ped[n]))
                else:  # founder: two unique chromosomes sets
                    nid[0] += 2
                    genomes[n] = tuple([[(0.0, nid[0] + k)] for _ in CHROMS] for k in (0, 1))
            return genomes[n]

        ga = get(a)
        for o in others:
            res[o].append(shared_cm(ga, get(o), min_seg))
    return res


def cmd_check(cm, rel, trials=2000):
    """Does a shared-cM amount fit a relationship a tree claims? A guard against flawed public trees.
    rel examples: sibling, aunt (also uncle/niece/nephew), 1C, 1C1R, 2C, 3C2R; prefix 'half' for half."""
    r = rel.lower().replace(" ", "").replace("-", "")
    half = r.startswith("half")
    r = r[4:] if half else r
    m = re.fullmatch(r"(\d)c(?:(\d)r)?", r)
    if m:
        up, down = int(m.group(1)) + 1, int(m.group(1)) + 1 + int(m.group(2) or 0)
    elif r in ("sibling", "sister", "brother"):
        up, down = 1, 1
    elif r in ("aunt", "uncle", "niece", "nephew"):
        up, down = 1, 2
    elif r in ("greataunt", "greatuncle", "grandaunt", "granduncle", "grandniece", "grandnephew"):
        up, down = 1, 3
    else:
        sys.exit(f"don't know the relationship '{rel}'. Try: sibling, aunt, granduncle, 1C, 1C1R, 2C, half 1C")
    # two lines of descent from a couple (or from one shared ancestor, for half relationships)
    ped = {"A1": ("G1", "G2"), "B1": ("G1", "G3" if half else "G2")}
    for side, n in (("A", up), ("B", down)):
        for k in range(2, n + 1):
            ped[f"{side}{k}"] = (f"{side}{k - 1}", f"S{side}{k}")
    v = sorted(simulate(ped, f"A{up}", [f"B{down}"], trials)[f"B{down}"])
    pct = sum(x <= cm for x in v) / trials
    lo, med, hi = v[trials // 100], v[trials // 2], v[-1 - trials // 100]
    verdict = ("FITS" if 0.05 <= pct <= 0.95 else "UNUSUAL but possible" if lo <= cm <= hi
               else "DOES NOT FIT: the tree is probably wrong about this relationship")
    print(f"{cm:.0f} cM as {rel}: {verdict}")
    print(f"  simulated range for that relationship: {lo:.0f} to {hi:.0f} cM (typical {med:.0f})")
    return verdict


def parse_rel(text):
    child, parents = text.split("=")
    parents = parents.split()
    if len(parents) != 2:
        sys.exit(f"need exactly two parents: {text.strip()}")
    return child.strip(), tuple(parents)


def parse_case(path):
    ped, obs, hyps, auto = {}, {}, {}, False
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            if line == "hyp auto":
                auto = True
            elif line.startswith("match "):
                _, n, c = line.split()
                obs[n] = float(c)
            elif line.startswith("hyp "):
                label, body = line[4:].split(":", 1)
                hyps[label.strip()] = dict(parse_rel(x) for x in body.split(";") if x.strip())
            else:
                c, p = parse_rel(line)
                ped[c] = p
    if auto:
        # walk every person in the tree and try TARGET as their child, grandchild and great-grandchild
        for p in sorted({q for c, ps in ped.items() for q in (c, *ps)} - {"TARGET"}):
            hyps[f"child_of_{p}"] = {"TARGET": (p, f"U1_{p}")}
            hyps[f"grandchild_of_{p}"] = {f"X1_{p}": (p, f"U1_{p}"), "TARGET": (f"X1_{p}", f"U2_{p}")}
            hyps[f"great_grandchild_of_{p}"] = {f"X1_{p}": (p, f"U1_{p}"), f"X2_{p}": (f"X1_{p}", f"U2_{p}"),
                                                "TARGET": (f"X2_{p}", f"U3_{p}")}
    return ped, obs, hyps


def write_gedcom(ped, out, unknown=(), born=None):
    """ped: {child: (parent, parent)}. Names in `unknown` are written as sexless 'Unknown' placeholders."""
    people = sorted({p for c, ps in ped.items() for p in (c, *ps)})
    pid = {n: f"@I{i}@" for i, n in enumerate(people, 1)}
    fams = {}
    for c, ps in ped.items():
        fams.setdefault(ps, []).append(c)
    fid = {ps: f"@F{i}@" for i, ps in enumerate(fams, 1)}
    sex = {p: s for ps in fams for p, s in zip(ps, "MF")}
    lines = ["0 HEAD", "1 SOUR birthfamily", "1 GEDC", "2 VERS 5.5.1", "2 FORM LINEAGE-LINKED", "1 CHAR UTF-8"]
    for n in people:
        parts = n.replace("_", " ").split()
        name = "Unknown //" if n in unknown else f"{' '.join(parts[:-1]) or parts[0]} /{parts[-1] if len(parts) > 1 else ''}/"
        lines += [f"0 {pid[n]} INDI", f"1 NAME {name}"]
        if n in sex and n not in unknown:
            lines.append(f"1 SEX {sex[n]}")
        if born and n in born:
            lines += ["1 BIRT", f"2 DATE {born[n]}"]
        if n in ped:
            lines.append(f"1 FAMC {fid[ped[n]]}")
        lines += [f"1 FAMS {fid[ps]}" for ps in fams if n in ps]
    for ps, kids in fams.items():
        lines += [f"0 {fid[ps]} FAM", f"1 HUSB {pid[ps[0]]}", f"1 WIFE {pid[ps[1]]}"]
        lines += [f"1 CHIL {pid[k]}" for k in kids]
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines + ["0 TRLR"]) + "\n")
    print(f"wrote {out}: {len(people)} people, {len(fams)} families")


def cmd_gedcom(path, hyp, out):
    """Write the known pedigree plus one chosen hypothesis as a GEDCOM for upload as a private tree.
    Convention: in 'child = father mother' the father is listed first. Underscores become spaces."""
    ped, _, hyps = parse_case(path)
    if hyp not in hyps:
        sys.exit(f"no hypothesis '{hyp}'. Available: {', '.join(hyps)}")
    write_gedcom({**ped, **hyps[hyp]}, out)


def rel_shape(rel):
    """(adoptee's generations below the shared ancestors, match's generations below, half?) or None.
    Takes the first option of an Ancestry label like '1st cousin 1x removed or half granduncle'."""
    # ponytail: for 'removed' cousins the match is assumed to be the OLDER generation. Ancestry's label
    # doesn't say which way round; the searcher moves the match if their age says otherwise.
    r = rel.lower().split(" or ")[0].strip()
    half = r.startswith("half")
    r = r[4:].strip() if half else r
    m = re.fullmatch(r"(\d)(?:st|nd|rd|th) cousin(?: (\d)x removed)?", r)
    if m:
        n, k = int(m.group(1)), int(m.group(2) or 0)
        return n + 1 + k, n + 1, half
    r = r.replace(" ", "").replace("-", "")
    table = {"sibling": (1, 1), "sister": (1, 1), "brother": (1, 1), "aunt": (2, 1), "uncle": (2, 1),
             "niece": (1, 2), "nephew": (1, 2), "granduncle": (3, 1), "grandaunt": (3, 1),
             "greatuncle": (3, 1), "greataunt": (3, 1), "grandniece": (1, 3), "grandnephew": (1, 3)}
    return table[r] + (half,) if r in table else None


def cmd_skeleton(folder, name):
    """The blank canvas a searcher starts with: the adoptee, two unknown parents, and on each parent's
    side the closest DNA match hung from unknown, sexless placeholder ancestors up to the shared couple."""
    rows = sorted(read_csv(os.path.join(folder, "matches.csv")), key=lambda r: -float(r["cm"]))
    me = "_".join(name.split())
    ped, unknown, placed = {}, set(), []

    def blank(tag):
        unknown.add(tag)
        return tag

    parents = {s: blank(f"?{s}") for s in ("Parent 1", "Parent 2")}
    ped[me] = (parents["Parent 1"], parents["Parent 2"])
    for s, p in parents.items():
        top = next((r for r in rows if (r.get("side") or "").startswith(s) and rel_shape(r.get("relationship", ""))), None)
        if not top:
            print(f"{s}: no match with a usable relationship label; left as a blank parent.")
            continue
        down, mdown, half = rel_shape(top["relationship"])
        if down == 1:     # sibling or niece/nephew level: the shared couple are the adoptee's own parents
            couple = ped[me] if not half else (p, blank(f"?{s} other partner"))
        else:
            couple = (blank(f"?{s} ancestor A"), blank(f"?{s} ancestor B"))
            line = [p] + [blank(f"?{s} line gen {i}") for i in range(2, down)]
            for i, person in enumerate(line):
                ped[person] = couple if i == len(line) - 1 else (line[i + 1], blank(f"?{s} spouse of gen {i + 2}"))
        match = "_".join(top["name"].split())
        mline = [match] + [blank(f"?{s} match line gen {i}") for i in range(1, mdown)]
        for i, person in enumerate(mline):
            if i < len(mline) - 1:
                ped[person] = (mline[i + 1], blank(f"?{s} match line spouse {i + 1}"))
            else:
                ped[person] = couple if not half else (couple[0], blank(f"?{s} half partner"))
        placed.append(f"{s}: {top['name']} ({float(top['cm']):.0f} cM) placed as '{top['relationship'].split(' or ')[0]}'"
                      f", shared couple {down - 1} generation(s) above the parent")
    by = adoptee(folder).get("birth_year")
    out = os.path.join(folder, "skeleton.ged")
    write_gedcom(ped, out, unknown, {me: by} if by else None)
    for line in placed:
        print("  " + line)
    print("Next: on Ancestry create a new tree by uploading this file, keep it private, and link the adoptee's\n"
          "DNA test to the adoptee in it. Replace each Unknown as the walk identifies them.")


def cmd_odds(path, trials):
    ped, obs, hyps = parse_case(path)
    if not obs or not hyps:
        sys.exit("case file needs at least one 'match' line and one 'hyp' line")
    scores, detail = {}, {}
    for label, extra in hyps.items():
        sims = simulate({**ped, **extra}, "TARGET", list(obs), trials)
        # ponytail: likelihood = product over matches of P(simulated cM lands near observed), treating
        # matches as independent (WATO does the same). Joint density needs far more trials.
        like, rows = 1.0, []
        for n, o in obs.items():
            w = max(15.0, 0.12 * o)
            v = sorted(sims[n])
            hit = sum(abs(x - o) <= w for x in v)
            like *= max(hit, 0.5) / trials
            rows.append((n, o, v[len(v) // 20], v[len(v) // 2], v[-1 - len(v) // 20], hit == 0))
        scores[label], detail[label] = like, rows
    tot = sum(scores.values())
    for rank, label in enumerate(sorted(scores, key=scores.get, reverse=True)):
        if rank >= 6 and len(scores) > 8:   # long auto lists: detail for the leaders only
            print(f"{label}: {100 * scores[label] / tot:.1f}%")
            continue
        print(f"\n{label}: {100 * scores[label] / tot:.1f}%")
        for n, o, p5, med, p95, miss in detail[label]:
            print(f"  {n}: observed {o:.0f} cM, simulated {p5:.0f}-{p95:.0f} (median {med:.0f})"
                  + ("   <-- never happened in simulation" if miss else ""))
    print("\nPercentages compare only the hypotheses listed, with equal priors. Not proof: confirm with records.")
    return scores


def case_folder(name):
    """Cases live in a 'cases' folder next to this file, so volunteers never type a path."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases", "_".join(name.split()))


def write_report(folder):
    import contextlib
    out = os.path.join(folder, "clusters.txt")
    with open(out, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        cmd_cluster(os.path.join(folder, "matches.csv"), os.path.join(folder, "shared.csv"),
                    os.path.join(folder, "trees"), 40, 400)
    return out


def show_report(folder):
    out = write_report(folder)
    print(f"The report is saved as {out}")
    if hasattr(os, "startfile"):
        os.startfile(out)   # opens in Notepad on Windows


def ask(prompt):
    try:
        return input(prompt).strip()
    except EOFError:
        sys.exit()


VERSION = "0.1.0"
REPO = "searchangels/searchangel-tools"            # where updates come from
REPO_FILE = "BirthFamilyTool/birthfamily.py"


def check_update():
    """Offer a newer copy of this file from GitHub. Download only: nothing about any case is sent.
    A private repo needs a read-only token saved next to this file as update_token.txt; without one
    (or offline) the check is skipped silently and the tool runs as it is."""
    import shutil, urllib.request
    here = os.path.abspath(__file__)
    tok = os.path.join(os.path.dirname(here), "update_token.txt")
    ver = lambda s: tuple(int(x) for x in s.split("."))
    try:
        if os.environ.get("BIRTHFAMILY_UPDATE_URL"):          # for testing the updater itself
            req = urllib.request.Request(os.environ["BIRTHFAMILY_UPDATE_URL"])
        elif os.path.exists(tok):
            with open(tok, encoding="utf-8") as f:
                token = f.read().strip()
            req = urllib.request.Request(f"https://api.github.com/repos/{REPO}/contents/{REPO_FILE}",
                                         headers={"Authorization": f"Bearer {token}", "Accept": "application/vnd.github.raw"})
        else:
            req = urllib.request.Request(f"https://raw.githubusercontent.com/{REPO}/main/{REPO_FILE}")
        with urllib.request.urlopen(req, timeout=6) as r:
            new = r.read().decode("utf-8")
        m = re.search(r'^VERSION = "(\d+(?:\.\d+)*)"', new, re.M)
        if not m or ver(m.group(1)) <= ver(VERSION):
            return
        compile(new, "birthfamily.py", "exec")               # never install a download that isn't valid code
    except Exception:
        return
    ans = ask(f"\nVersion {m.group(1)} of this tool is available (you have {VERSION}).\n"
              "Press Enter to update now, or type n to keep this version: ")
    if ans.lower().startswith("n"):
        return
    shutil.copyfile(here, here + ".previous")                # one step back, in case the new one misbehaves
    with open(here, "w", encoding="utf-8", newline="\n") as f:
        f.write(new)
    print(f"Updated to version {m.group(1)}. Close this window and start the tool again.")
    sys.exit()


def cmd_start():
    """Double-click entry point 1: name the case, open the capture window, then show the report."""
    check_update()
    print(f"Birth Family Tool, version {VERSION}\n")
    print("Give this case a short name, for example: smith")
    print("Use the same name every time you come back to this case.\n")
    name = ask("Case name: ")
    if not name:
        return
    folder = case_folder(name)
    print("\nA small window is opening. Keep it next to your browser and follow what it says.")
    print("Close the window when you are finished or want a break. Your progress is saved.\n")
    cmd_guide(folder, 90.0)
    if os.path.exists(os.path.join(folder, "matches.csv")):
        show_report(folder)
    else:
        print("Nothing was captured, so there is no report yet.")


def cmd_menu(open_files=True):
    """Double-click entry point 2: a numbered menu over the commands a searcher uses after capturing."""
    check_update()
    print(f"Birth Family Tool, version {VERSION}\n")
    name = ask("Case name, the same one you used when capturing: ")
    folder = case_folder(name)
    if not name or not os.path.exists(os.path.join(folder, "matches.csv")):
        print(f'\nNothing has been captured for "{name}" yet. Run "1 Capture a case" first.')
        return
    while True:
        print("\n---------------------------------------------------------\n"
              "  1  Save the adoptee's birth year and place\n"
              "  2  Show the report\n"
              "  3  Make the starting tree file to upload to Ancestry\n"
              "  4  Walk down the generations\n"
              "  5  Does a surname show up in the matches?\n"
              "  6  Does a cM amount fit a relationship?\n"
              "  7  Quit\n"
              "---------------------------------------------------------")
        ch = ask("Type a number and press Enter: ")
        try:
            if ch == "1":
                year = ask("Adoptee's birth year, for example 1970: ")
                place = ask("State or place of birth, for example California: ")
                if year.isdigit():
                    cmd_adoptee(folder, year, place)
                else:
                    print("That isn't a year. Nothing saved.")
            elif ch == "2":
                show_report(folder) if open_files else print(f"The report is saved as {write_report(folder)}")
            elif ch == "3":
                who = ask("Adoptee's name as it should appear in the tree: ")
                if who:
                    cmd_skeleton(folder, who)
            elif ch == "4":
                cmd_walk(folder)
            elif ch == "5":
                names = ask("Surname or surnames, separated by spaces: ").split()
                if names:
                    cmd_find(folder, names)
            elif ch == "6":
                cm = ask("Shared cM, for example 497: ").replace(",", "")
                rel = ask("Relationship, for example 1C1R, 2C, half 1C, aunt, granduncle: ")
                if cm.replace(".", "", 1).isdigit() and rel:
                    cmd_check(float(cm), rel)
                else:
                    print("Need a number and a relationship.")
            elif ch == "7":
                return
        except SystemExit as e:   # a command refused its input: show why and return to the menu
            if e.code:
                print(e.code)


# ---------- self-check ----------

def selftest():
    random.seed(1)
    ped = {"A": ("G1", "G2"), "B": ("G1", "G2"), "T": ("A", "X"), "C": ("B", "Y")}
    s = simulate(ped, "T", ["A", "B", "C"], 300)
    mean = lambda v: sum(v) / len(v)
    assert abs(mean(s["A"]) - sum(CHROMS)) < 1, mean(s["A"])            # parent/child = whole genome
    # reported cM is half-identical coverage of one genome copy, so aunt ~50% (~1,750), 1C ~25% (~870)
    assert 0.45 < mean(s["B"]) / sum(CHROMS) < 0.55, mean(s["B"])
    assert 0.21 < mean(s["C"]) / sum(CHROMS) < 0.28, mean(s["C"])
    cm = {"m1": 300, "m2": 200, "m3": 150, "m4": 90, "sis": 2500}
    cl = cluster(cm, [("m1", "m2"), ("m3", "m4"), ("sis", "m1"), ("sis", "m3")], 40, 400)
    assert cl == [["m1", "m2"], ["m3", "m4"]], cl
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        g = os.path.join(d, "t.ged")
        open(g, "w").write("0 @I1@ INDI\n1 NAME John Q /Smith/\n1 BIRT\n2 DATE 12 MAR 1880\n1 DEAT\n2 DATE 1950\n")
        assert gedcom_people(g) == {("john", "smith", "1880")}
        case = os.path.join(d, "case.txt")
        open(case, "w").write("A = G1 G2\nB = G1 G2\nC = B Y\nmatch C 850\n"
                              "hyp first_cousin: TARGET = A U1\n"
                              "hyp first_cousin_once_removed: P = A U1; TARGET = P U2\n")
        sc = cmd_odds(case, 300)
        assert sc["first_cousin"] > 5 * sc["first_cousin_once_removed"], sc
        assert cmd_check(870, "1C", 300) == "FITS" and cmd_check(870, "3C", 300).startswith("DOES NOT")
        cmd_adoptee(os.path.join(d, "case"), "1970", "California")
        cd = os.path.join(d, "case")
        assert generation("1885", cd) == "great-grandparent level" and generation("1945", cd) == "parent level"
        assert "COULD BE" in parent_age("1945", cd) and "too young" in parent_age("1960", cd) and "grandparent" in parent_age("1900", cd)
        assert "California Birth Index" in records_hint(cd)
        assert rel_shape("granduncle") == (3, 1, False) and rel_shape("1st cousin 1x removed or half granduncle") == (3, 2, False)
        assert rel_shape("half 2nd cousin or 2nd cousin 1x removed") == (3, 3, True) and rel_shape("close family") is None
        auto_case = os.path.join(d, "auto.txt")
        open(auto_case, "w").write("A = G1 G2\nB = G1 G2\nC = B Y\nD = A Z\nmatch C 850\nmatch D 1700\nhyp auto\n")
        sc = cmd_odds(auto_case, 150)   # truth: TARGET is D's half or full sibling, i.e. a child of A
        best = max(sc, key=sc.get)
        assert best in ("child_of_A", "child_of_Z"), best
        t = os.path.join(d, "t.txt")
        open(t, "w").write("Mary Ann Jones\n1880–1950\n\nHints\nJohn Smith\nBorn 1875\n")
        assert text_people(t) == {("mary", "jones", "1880"), ("john", "smith", "1875")}, text_people(t)
        paste = ("Jane Doe\nManaged by Bob Roe\n1st–2nd Cousin\nPaternal side\n412 cM | 6% shared DNA\n"
                 "Public linked tree\n1,204 People\nAdd note\nSam O'Neil\nClose family\n1,690 cM | 24% shared DNA\n")
        got = parse_matches(paste)   # no blank line under the name: fallback path
        assert [r[:2] for r in got] == [("Jane Doe", 412.0), ("Sam O'Neil", 1690.0)], got
        real = ("More options\nnote - x\nP\nsarah62\nManaged by Al Roe\n\n2nd cousin\nParent 2's side\n"
                "222 cM | 3% shared DNA\nPublic linked tree\nPublic linked tree\n4 people\nCommon ancestor\nAdd\n")
        assert parse_matches(real) == [("sarah62", 222.0, "2nd cousin", "Parent 2's side",
                                        "Public linked tree; 4 people; Common ancestor")], parse_matches(real)
        row = lambda n, cm, tree: f"{n}\n\n2nd cousin\nParent 1's side\n{cm} cM | 3% shared DNA\n{tree}\nAdd\nMore options\n"
        g = Guide(os.path.join(d, "case"))
        assert g.on_copy("random text")[0] is False
        ok, m = g.on_copy("All matches\nBy cluster\n" + row("Ann Lee", 300, "Public linked tree\n43 people")
                          + row("Bo Wu", 250, "No trees") + row("Cy Low", 80, "No trees"))
        assert ok and g.stage == "pages" and g.queue == ["Ann Lee"], (m, g.stage, g.queue)
        ok, m = g.on_copy("All matches\nBy cluster\n" + row("Dee Fox", 120, "Public linked tree\n9 people"))
        assert ok and g.queue == ["Ann Lee", "Dee Fox"], (m, g.queue)                         # list copied mid-way extends worklist
        g.skip("Dee Fox")
        page = lambda n, extra: f"Match list{n}\nYou and {n}\nHow are you and {n} related?\n{extra}"
        g3 = Guide(os.path.join(d, "case3"))                                                  # filtered list: ends on last page
        ok, m = g3.on_copy("All matches\nBy cluster\n" + row("Ann Lee", 300, "No trees") + "50\nItems per page\n\n2\nof 2\n")
        assert ok and "last page" in m and g3.stage == "shared", (m, g3.stage)
        g4 = Guide(os.path.join(d, "case4"))                                                  # one-page list: no "of N"
        assert "last page" in g4.on_copy("All matches\nBy cluster\n" + row("Ann Lee", 300, "No trees")
                                         + "20\nItems per page\n1\nSupport Center\n")[1]
        g5 = Guide(os.path.join(d, "case5"))                                                  # page 1 of 3: keep going
        assert "next 50" in g5.on_copy("All matches\nBy cluster\n" + row("Ann Lee", 300, "No trees")
                                       + "50\nItems per page\n\n1\nof 3\n")[1]
        assert g.on_copy(page("Bo Wu", "No trees\n"))[0] is False                             # no tree: not on the worklist
        assert g.on_copy(page("Zed Nobody", ""))[0] is False                                  # not on the worklist
        ok, m = g.on_copy(page("Ann Lee", "Ann's Linked Tree43 People\nPrivate\n\nMary Ann JonesJones shared surname\n"
                                          "1880-1950\n\nJohn SmithSmith\nSurnames\nShared Surnames\n"))
        assert ok and "1 dated" in m and g.stage == "shared" and g.queue == ["Ann Lee", "Bo Wu", "Dee Fox"], (m, g.stage, g.queue)
        assert text_people(os.path.join(d, "case", "trees", "Ann_Lee.txt")) == {("mary", "jones", "1880")}
        assert g.on_copy(page("Ann Lee", ""))[0] is False                                     # wrong tab refused
        assert g.on_copy("You and Ann Lee\nShared Matches\nBo Wu\n250 cM\n")[0]
        g.finish_step()                                                                       # leave Bo Wu, Dee Fox
        assert g.stage == "done", g.stage
        rep = open(os.path.join(d, "case", "clusters.txt"), encoding="utf-8").read()
        assert "Cluster 1" in rep and "Ann Lee" in rep and "Bo Wu" in rep, rep
        assert g.rows()[0]["tree"] == "Ann_Lee.txt"
        cmd_skeleton(cd, "Test Adoptee")   # case has Ann Lee, '2nd cousin', Parent 1's side
        sk = open(os.path.join(cd, "skeleton.ged")).read()
        assert "Test /Adoptee/" in sk and "Ann /Lee/" in sk and "Unknown //" in sk and "DATE 1970" in sk, sk
        g2 = Guide(os.path.join(d, "case"))                                                   # reopening resumes
        assert g2.stage == "pages" and g2.queue == ["Dee Fox"], (g2.stage, g2.queue)                  # skipped one is re-offered
        out = os.path.join(d, "o.ged")
        cmd_gedcom(case, "first_cousin", out)
        t = open(out).read()
        assert t.count(" INDI") == 8 and t.count(" FAM\n") == 3 and "1 CHIL" in t, t
    print("\nselftest OK")


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda flag, default: a[a.index(flag) + 1] if flag in a else default
    if a[:1] == ["selftest"]:
        selftest()
    elif a[:1] == ["cluster"] and len(a) >= 3:
        cmd_cluster(a[1], a[2], opt("--trees", None), float(opt("--min", 40)), float(opt("--max", 400)))
    elif a[:1] == ["paste"] and len(a) >= 3:
        cmd_paste(a[1], a[2], opt("--with", None))
    elif a[:1] == ["find"] and len(a) >= 3:
        cmd_find(a[1], a[2:])
    elif a[:1] == ["walk"] and len(a) >= 2:
        cmd_walk(a[1])
    elif a[:1] == ["check"] and len(a) >= 3:
        cmd_check(float(a[1]), " ".join(a[2:]))
    elif a[:1] == ["adoptee"] and len(a) >= 3:
        cmd_adoptee(a[1], a[2], " ".join(a[3:]))
    elif a[:1] == ["skeleton"] and len(a) >= 3:
        cmd_skeleton(a[1], " ".join(a[2:]))
    elif a[:1] == ["start"]:
        cmd_start()
    elif a[:1] == ["menu"]:
        cmd_menu()
    elif a[:1] == ["guide"] and len(a) >= 2:
        cmd_guide(a[1], float(opt("--floor", 90)))
    elif a[:1] == ["watch"]:
        cmd_watch(a[1] if len(a) > 1 else ".")
    elif a[:1] == ["gedcom"] and len(a) >= 4:
        cmd_gedcom(a[1], a[2], a[3])
    elif a[:1] == ["odds"] and len(a) >= 2:
        cmd_odds(a[1], int(opt("--trials", 1000)))
    else:
        print(__doc__)
