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
    if not os.path.exists(path):   # a case that has just been started has nothing copied into it yet
        return []
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
        if not m or "across" in l or "shared DNA:" in l:   # "5% shared DNA: 352 cM across 20 segments" heads a match's own page; it is not a row
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
            f.writelines(f"# born {c} {y}\n" for c, y in born.items() if y)
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
            with open(tree_path, "a", encoding="utf-8") as f:   # lets the tree page hang the adoptee under them
                f.write(f"@parents|{line_side}|{c}|{s}\n")
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
                    "and in Max the highest number Ancestry accepts (try 3490). Click Apply filters.\n"
                    "3. Set 50 per page at the bottom of the list.\n"
                    "4. Press Ctrl+A, then Ctrl+C.\n5. Go to the next page and repeat until this window says the list is done.")
        if self.stage == "pages":
            return (f"STEP 2 of 3: PUBLIC TREES {pos}\nFirst, on the match list you just copied: "
                    "click Filter > Trees > Public linked trees > Apply filters. This only makes the names below easier to find.\n\n"
                    f"1. Click the name:\n\n    {self.where(self.queue[0])}\n\n"
                    "2. Press Ctrl+A, then Ctrl+C.\n3. Press the browser Back button." + find)
        if self.stage == "shared":
            return (f"STEP 3 of 3: SHARED MATCHES {pos}\nFirst, back on the match list: click Filter > Trees > All > Apply filters.\n\n"
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

    # Labels kept readable in layout_report.txt; every other word is masked so no names leave the computer.
    KEEP = set("you and shared dna cm across segments how are related linked unlinked tree trees people surnames "
               "matches match all by cluster private public common ancestor ancestors ethnicity relationship "
               "maternal paternal parent both sides unassigned items per of managed options more".split())
    EARLY = ("That page had not finished loading when it was copied, so nothing was saved. Wait until the list "
             "has stopped filling in and you can see \"Items per page\" at the bottom, then press Ctrl+A and Ctrl+C again.")
    UNKNOWN = ("\nA copy of the page's layout with all names blanked out was saved as layout_report.txt in the case "
               "folder. Please email that file to support@searchangels.org so the tool can be fixed.")

    def page_owner(self, text):
        """Name of the match whose page was copied. On your own test the header is "You and <name>"; on a
        test someone shared with you it names the test taker instead, so look for a worklist name too."""
        who = re.search(r"^\s*You and\s+(.+?)\s*$", text, re.M)
        if who:
            return who.group(1)
        # ponytail: assumed header "<test taker> and <name>"; unverified against a real shared-test page
        return next((q for q in self.queue
                     if re.search(r"\band\s+" + re.escape(q) + r"(\s+related\?)?\s*$", text, re.M)), None)

    def on_copy(self, text):
        """Handle one copied page. Returns (ok, message). A refused page leaves a name-free layout report."""
        ok, msg = self.take(text)
        if not ok:
            mask = lambda m: m.group() if m.group().lower() in self.KEEP else "x" * len(m.group())
            with open(os.path.join(self.folder, "layout_report.txt"), "w", encoding="utf-8") as f:
                f.write(f"Birth Family Tool {VERSION}, step: {self.stage}\n{msg}\n\n" + re.sub(r"[^\W\d_]+", mask, text))
        return ok, msg

    def take(self, text):
        self.n += 1
        with open(os.path.join(self.folder, "pastes", f"{self.n:03d}.txt"), "w", encoding="utf-8") as f:
            f.write(text)
        main_list = "All matches" in text and "By cluster" in text
        if self.stage == "matches":
            found = parse_matches(text)
            if not main_list or not found:
                return False, "That doesn't look like the DNA match list. Nothing was saved." + self.UNKNOWN
            new, low = save_rows(text, self.m, None), min(r[1] for r in found)
            msg = f"Received {len(new)} new matches ({len(self.rows())} total). Lowest on this page: {low:.0f} cM."
            if "Items per page" not in text:   # kept, because every row read is real, but the page was probably still loading
                msg += ("\nCHECK: this page may not have finished loading when you copied it (no \"Items per page\" line at the "
                        "bottom). If the list on your screen is longer, wait a moment and copy the same page again.")
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
        name = self.page_owner(text)
        if main_list:
            # a (filtered) match list copied mid-way: keep any new people and refresh the worklist
            new = save_rows(text, self.m, None)
            if self.stage == "pages":
                self.stage = "matches"
                self.next_stage()
            return True, f"Match list received ({len(new)} new people added). Now click the name shown above."
        if not name:
            return False, ("That isn't the page of anyone on the worklist for this step, or the tool can't read "
                           "this page. Nothing was saved." + self.UNKNOWN)
        if name not in self.queue:
            return False, f"That is {name}'s page, which isn't on the worklist for this step (or is already done)."
        trees_tab = re.search(r"How are .+ related\?", text)
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
        # The Shared Matches tab lists people in the same layout as the main match list and ends with an
        # "Items per page" line (seen on real copied pages, Nov 2026). Without that line the list had not finished
        # loading when it was copied, and saving it would record a partial list as if it were complete.
        if "Items per page" not in text:
            return False, self.EARLY
        listed = [r for r in parse_matches(text) if r[0] != name]
        if listed:
            save_rows(text, self.m, None)   # shared matches are the adoptee's matches too: keep any not seen before
            pairs = [(name, r[0]) for r in listed]
        else:   # a layout this reader doesn't know: fall back to known match names appearing on their own line
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
    ped, unknown, placed, hang = {}, set(), [], {}   # hang: how many generations above each parent the shared couple sits

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
        hang[s] = down - 1
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
    return write_tree_page(folder, me, ped, unknown, rows, hang)


def write_tree_page(folder, me, ped, unknown, rows, hang=None):
    """tree.html: the starting tree joined with everything the walk found, drawn as a tree the searcher
    can edit in the browser and download as a GEDCOM file. One offline file, nothing is sent anywhere."""
    import glob, json, time
    born, walk, parents = {}, {}, None
    if adoptee(folder).get("birth_year"):
        born[me] = adoptee(folder)["birth_year"]
    for fn in sorted(glob.glob(os.path.join(folder, "tree_*.txt"))):
        with open(fn, encoding="utf-8") as f:
            for line in f:
                w = line.split()
                if line.startswith("# born") and len(w) == 4:
                    born[w[2]] = w[3]
                elif line.startswith("@parents|"):
                    _, side, c, s = line.strip().split("|")
                    # ponytail: only "Parent 2" swaps the order; Maternal/Paternal labels land as typed
                    parents = (s, c) if side.startswith("Parent 2") else (c, s)
                elif len(w) == 4 and w[1] == "=":
                    walk[w[0]] = (w[2], w[3])
    if parents:   # the walk named the birth parents: they take the place of the two blank parents
        swap = dict(zip(ped[me], parents))
        ped = {swap.get(c, c): tuple(swap.get(p, p) for p in ps) for c, ps in ped.items()}
        unknown = unknown - set(swap)
    ped = {**ped, **walk}   # where both describe someone's parents, what the searcher typed in the walk wins
    while True:   # a blank placeholder left with no children after the join says nothing: drop it
        dead = [c for c in ped if c in unknown and not any(c in ps for ps in ped.values())]
        if not dead:
            break
        for c in dead:
            del ped[c]
    names = sorted({p for c, ps in ped.items() for p in (c, *ps)})
    pid = {n: i for i, n in enumerate(names, 1)}
    cm = {"_".join(r["name"].split()): r for r in rows}
    fams = {}
    for c, ps in ped.items():
        fams.setdefault(ps, []).append(pid[c])
    # Sex is only set for a top couple typed as "husband, wife". Everyone else stays blank (gender neutral)
    # until the searcher sets it, because the walk doesn't know which of a later couple is which.
    sex = {p: s for ps in fams if not any(q in ped or q in unknown for q in ps) for p, s in zip(ps, "MF")}
    people = [{"id": pid[n], "name": "Unknown" if n in unknown else n.replace("_", " "), "sex": sex.get(n, ""),
               "b": born.get(n, ""), "d": "", "root": n == me, "ph": n if n in unknown else "",   # ph: which placeholder this is
               "note": (f"DNA match, {float(cm[n]['cm']):.0f} cM, {cm[n].get('relationship', '')}" if n in cm
                        else "The adoptee" if n == me else "")} for n in names]
    data = {"case": os.path.basename(os.path.abspath(folder)), "stamp": time.strftime("%Y-%m-%d %H:%M:%S"), "people": people,
            "couple_depth": hang or {},
            "fams": [{"a": pid[a], "b": pid[b], "kids": k} for (a, b), k in fams.items()]}
    out = os.path.join(folder, "tree.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(TREE_PAGE.replace("/*DATA*/", json.dumps(data).replace("</", "<\\/")))
    print(f"wrote {out}: {len(people)} people. It opens in your browser; edit it there and use its\n"
          "Download GEDCOM button to get the file to upload to Ancestry as a new private tree.")
    return out


def record_files():
    """{'numident.db': path, ...} for the record files found on this computer. They are large and optional.
    Looked for in: the BIRTHFAMILY_RECORDS setting, folders listed in records_path.txt beside this file, and
    a 'records' folder beside this file."""
    dirs = [d for d in os.environ.get("BIRTHFAMILY_RECORDS", "").split(os.pathsep) if d]
    for here in dict.fromkeys((data_dir(), program_dir())):
        listing = os.path.join(here, "records_path.txt")
        if os.path.exists(listing):
            with open(listing, encoding="utf-8") as f:
                dirs += [l.strip() for l in f if l.strip()]
        dirs.append(os.path.join(here, "records"))
    return {n: os.path.join(d, n) for d in reversed(dirs) for n in ("numident.db", "texas.db", "birls.db", "va_burials.db")
            if os.path.exists(os.path.join(d, n))}


def record_lookup(q, files=None):
    """Look one person up in whichever record files are present. q: first, last, maiden, b, d (their own),
    ff, fl (father's first and last name), ml (mother's maiden name). Read-only; nothing leaves this computer."""
    import sqlite3
    files = record_files() if files is None else files
    up = lambda k: re.sub(r"[^A-Z' -]", "", (q.get(k) or "").upper()).strip()
    first, last, maiden = up("first"), up("last"), up("maiden")
    year = lambda k: int(m.group()) if (m := re.search(r"\d{4}", q.get(k) or "")) else None
    b = year("b")
    nice = lambda *parts: " ".join(p for p in parts if p).title()
    out = {"have": sorted(files), "self": [], "children": [], "siblings": [], "marriages": [], "veterans": [], "burials": []}
    if not first or not last:
        return out

    def con(name):
        db = sqlite3.connect("file:" + files[name].replace("\\", "/") + "?mode=ro", uri=True)
        db.execute("PRAGMA case_sensitive_like=ON")   # lets 'first LIKE ...' use the (last, first) index
        return db

    if "numident.db" in files:
        db = con("numident.db")
        cols = "pid, first, middle, last, dob, pob_city, pob_st, ffirst, fmid, flast, mfirst, mmid, mlast"

        def people(where, args, other_parent=None, limit=400):
            """Rows grouped into people: one person can have several rows (married names, corrections)."""
            got = {}
            for r in db.execute(f"SELECT {cols} FROM rec WHERE {where} LIMIT {limit}", args):
                p = got.setdefault(r[0], {"name": nice(r[1], r[2], r[3]), "also": [], "b": r[4][:4], "pob": "", "father": "", "mother": "", "other": ""})
                if r[3].title() not in p["name"].split() and r[3].title() not in p["also"]:
                    p["also"].append(r[3].title())
                p["pob"] = p["pob"] or nice(r[5]) + (", " + r[6] if r[6] else "")
                p["father"] = p["father"] or nice(r[7], r[8], r[9])
                p["mother"] = p["mother"] or nice(r[10], r[11], r[12])
                p["b"] = p["b"] or r[4][:4]
            for pid, p in got.items():
                d = db.execute("SELECT dod FROM death WHERE pid = ? LIMIT 1", (pid,)).fetchone()
                p["d"] = d[0][:4] if d else ""
                if other_parent:
                    p["other"] = p[other_parent]
            return sorted(got.values(), key=lambda p: p["b"] or "9999")[:60]

        span = lambda lo, hi: (f" AND dob BETWEEN '{lo}' AND '{hi + 1}'" if b else "")
        names = [n for n in dict.fromkeys((last, maiden)) if n]
        marks = ",".join("?" * len(names))
        out["self"] = people(f"last IN ({marks}) AND first LIKE ?" + span(b - 2, b + 2) if b else f"last IN ({marks}) AND first LIKE ?",
                             (*names, first + "%"))
        kid_span = span(b + 14, b + 60) if b else ""
        as_father = people("flast = ? AND ffirst LIKE ?" + kid_span, (last, first + "%"), "mother")
        as_mother = people(f"mlast IN ({marks}) AND mfirst LIKE ?" + kid_span, (*names, first + "%"), "father")
        # ponytail: both are shown; a common name will list other people's children too. The researcher picks by the other parent.
        out["children"] = (as_father + as_mother)[:80]
        if up("ff") and up("fl"):
            where, args = "flast = ? AND ffirst LIKE ?", [up("fl"), up("ff") + "%"]
            if up("ml"):
                where, args = where + " AND mlast = ?", args + [up("ml")]
            me = {p["name"] for p in out["self"]}
            out["siblings"] = [p for p in people(where, args) if p["name"] not in me]
    if "texas.db" in files:
        db = con("texas.db")
        for kind in ("marriage", "divorce"):
            for mine, theirs in (("1", "2"), ("2", "1")):
                for r in db.execute(f"SELECT name{mine}, age{mine}, name{theirs}, age{theirs}, year, county, last{theirs}, first{theirs} FROM {kind} "
                                    f"WHERE last{mine} = ? AND first{mine} LIKE ? LIMIT 80", (last, first + "%")):
                    if b and r[1].isdigit() and abs((r[4] - int(r[1])) - b) > 3:
                        continue   # the age in the index puts this person's birth too far from the tree's
                    out["marriages"].append({"kind": kind, "year": r[4], "county": r[5].title(), "self": r[0].title(), "age": r[1],
                                             "spouse": nice(r[7], r[6]) if r[6] else r[2].title(), "spouse_age": r[3]})
        out["marriages"] = sorted(out["marriages"], key=lambda m: m["year"])[:60]
    born = (f" AND dob BETWEEN '{b - 2}' AND '{b + 3}'" if b else "")
    if "birls.db" in files:
        out["veterans"] = [{"name": nice(r[1], r[2], r[0]), "b": r[3][:4], "d": r[4][:4], "bfull": r[3], "dfull": r[4], "branch": r[5]}
                           for r in con("birls.db").execute("SELECT last, first, middle, dob, dod, branch FROM veteran WHERE last = ? AND first LIKE ?"
                                                            + born + " LIMIT 40", (last, first + "%"))]
    if "va_burials.db" in files:
        db = con("va_burials.db")
        row = lambda r: {"name": nice(r[1], r[2], r[0]), "b": r[3][:4], "d": r[4][:4], "cemetery": r[5], "state": r[6],
                         "relationship": r[7], "veteran": nice(r[9], r[8])}
        sel = "SELECT last, first, middle, dob, dod, cemetery, state, relationship, vet_last, vet_first FROM burial WHERE "
        own = [row(r) for r in db.execute(sel + "last = ? AND first LIKE ?" + born + " LIMIT 30", (last, first + "%"))]
        family = [row(r) for r in db.execute(sel + "vet_last = ? AND vet_first LIKE ? AND relationship <> 'Veteran (Self)' LIMIT 40", (last, first + "%"))]
        out["burials"] = own + [f for f in family if f not in own]
    return out


def auto_fill(tree, files=None):
    """What the Social Security records settle about people already on the tree, with no judgement call:
      - PARENTS, when exactly one person in the whole file has that first name, surname and birth year (and a
        middle initial that does not disagree), and that record names both parents;
      - CHILDREN of a couple, when a record names BOTH of them as father and mother by first name and surname.
    Anything less certain is left for the researcher's own lookup. Returns [{pid, items}] for the page to add."""
    import sqlite3
    files = record_files() if files is None else files
    if "numident.db" not in files or not isinstance(tree, dict):
        return []
    db = sqlite3.connect("file:" + files["numident.db"].replace("\\", "/") + "?mode=ro", uri=True)
    people = {p.get("id"): p for p in tree.get("people", []) if isinstance(p, dict)}
    fams = [f for f in tree.get("fams", []) if isinstance(f, dict)]
    child_of = {k: f for f in fams for k in f.get("kids", [])}
    words = lambda s: re.sub(r"[^A-Za-z' -]", "", str(s or "")).upper().split()
    known = lambda p: bool(p) and not re.match(r"(?i)unknown", p.get("name") or "") and len(words(p.get("name"))) >= 2
    maiden = lambda p: (m.group(1).upper() if (m := re.search(r"maiden name probably (\S+)", p.get("note") or "")) else words(p["name"])[-1])
    nice = lambda *parts: " ".join(x for x in parts if x).title()
    src, out = "Suggested by Social Security records, check", []
    yr = lambda q: (int(m.group()) if q and (m := re.search(r"\d{4}", str(q.get("b") or ""))) else 0)
    dod = lambda pid: (d[0][:4] if (d := db.execute("SELECT dod FROM death WHERE pid = ? LIMIT 1", (pid,)).fetchone()) else "")
    for p in people.values():
        if len(out) >= 60:
            break
        if not known(p) or p.get("root") or str(p.get("note") or "").startswith("DNA match"):
            continue
        w, y, f = words(p["name"]), yr(p), child_of.get(p.get("id"))
        has_parents = bool(f and any(known(people.get(x)) for x in (f.get("a"), f.get("b"))))
        if has_parents and y:
            continue
        if y:
            lo = hi = y
        else:
            # No birth year on the tree (typical for a parent named on a child's record). Work out when they must have
            # been born from their family: 15 to 60 years before their eldest child, and within 15 years of their spouse.
            mine = [g for g in fams if p["id"] in (g.get("a"), g.get("b"))]
            kids = [k for g in mine for k in (yr(people.get(x)) for x in g.get("kids", [])) if k]
            mates = [m for g in mine for m in (yr(people.get(g.get("b") if g.get("a") == p["id"] else g.get("a"))),) if m]
            if not kids and not mates:
                continue
            lo = max(([min(kids) - 60] if kids else []) + [m - 15 for m in mates])
            hi = min(([min(kids) - 15] if kids else []) + [m + 15 for m in mates])
            if lo > hi:
                continue
        names = tuple(dict.fromkeys((w[-1], maiden(p))))
        marks = ",".join("?" * len(names))
        rows = db.execute(f"SELECT pid, middle, ffirst, fmid, flast, mfirst, mmid, mlast, dob FROM rec WHERE last IN ({marks}) AND first = ? "
                          "AND substr(dob, 1, 4) BETWEEN ? AND ? LIMIT 300", (*names, w[0], str(lo), str(hi))).fetchall()
        mid = w[1] if len(w) > 2 else ""
        if mid:                                                  # a middle name or initial on the tree narrows it
            rows = [r for r in rows if not r[1] or r[1][0] == mid[0]]
        if len({r[0] for r in rows}) != 1:                       # nobody, or more than one person it could be
            continue
        items = [] if y else [{"self": True, "b": rows[0][8][:4], "d": dod(rows[0][0])}]
        full = next((r for r in rows if r[2] and r[4] and r[5] and r[7]), None)
        if full and not has_parents:
            items += [{"role": "parent", "name": nice(full[2], full[3], full[4]), "src": src},
                      {"role": "parent", "name": nice(full[5], full[6], full[7]), "src": src}]
        if items:
            out.append({"pid": p["id"], "items": items})
    # MARRIAGES. A wife's record is filed under her maiden name and again under her married name, as one person.
    # So for a couple on the tree: her children's records give her maiden name (when they all agree), and the one
    # woman with her first name who appears under BOTH surnames is her. That gives her birth year and her parents.
    named1 = lambda p: bool(p) and not re.match(r"(?i)unknown", p.get("name") or "") and len(words(p.get("name"))) >= 1
    for f in fams:
        if len(out) >= 90:
            break
        a, b = people.get(f.get("a")), people.get(f.get("b"))
        for man, wife in ((a, b), (b, a)):
            if not (known(man) and named1(wife)):
                continue
            mw, ww = words(man["name"]), words(wife["name"])
            noted = re.search(r"maiden name probably (\S+)", wife.get("note") or "")
            born_as = noted.group(1).upper() if noted else (ww[-1] if len(ww) >= 2 and ww[-1] != mw[-1] else "")
            if not born_as:
                seen = {r[0] for r in db.execute("SELECT DISTINCT mlast FROM rec WHERE flast = ? AND ffirst = ? AND mfirst = ? AND mlast <> '' LIMIT 40",
                                                 (mw[-1], mw[0], ww[0]))} - {mw[-1]}
                if len(seen) != 1:
                    continue
                born_as = seen.pop()
            wf = child_of.get(wife.get("id"))
            if wf and any(known(people.get(x)) for x in (wf.get("a"), wf.get("b"))) and wife.get("b"):
                continue                                         # already has her parents and birth year
            ids = [{r[0] for r in db.execute("SELECT pid FROM rec WHERE last = ? AND first = ? LIMIT 3000", (sur, ww[0]))} for sur in (born_as, mw[-1])]
            both = ids[0] & ids[1]
            if len(both) != 1 or max(len(ids[0]), len(ids[1])) >= 3000:   # not exactly one woman, or too common a name to be sure
                continue
            pid = both.pop()
            rows = db.execute("SELECT dob, ffirst, fmid, flast, mfirst, mmid, mlast FROM rec WHERE pid = ?", (pid,)).fetchall()
            full = next((r for r in rows if r[1] and r[3] and r[4] and r[6]), None)
            items = [{"self": True, "b": next((r[0][:4] for r in rows if r[0]), ""), "d": dod(pid), "maiden": born_as.title()}]
            if full:
                items += [{"role": "parent", "name": nice(full[1], full[2], full[3]), "src": src},
                          {"role": "parent", "name": nice(full[4], full[5], full[6]), "src": src}]
            out.append({"pid": wife["id"], "items": items})
            break
    for f in fams:
        if len(out) >= 120:
            break
        a, b = people.get(f.get("a")), people.get(f.get("b"))
        if not (known(a) and known(b)):
            continue
        for dad, mom in ((a, b), (b, a)):
            dw, mw = words(dad["name"]), words(mom["name"])
            got = {}
            for r in db.execute("SELECT pid, first, middle, last, dob FROM rec WHERE flast = ? AND ffirst = ? AND mlast = ? AND mfirst = ? LIMIT 120",
                                (dw[-1], dw[0], maiden(mom), mw[0])):
                if r[0] not in got or r[3] == dw[-1]:            # prefer the row under the birth surname
                    got[r[0]] = r
            items = [{"role": "child", "partner": mom["id"], "name": nice(r[1], r[2], r[3]), "b": r[4][:4], "d": dod(pid), "src": src}
                     for pid, r in got.items()]
            if items:
                out.append({"pid": dad["id"], "items": items})
                break
    return out


def paste_page(folder, text):
    """An Ancestry page pasted into the tree page instead of the collecting window: file it in the case the same way.
    Returns (ok, message)."""
    g = Guide(folder)
    if "All matches" in text and "By cluster" in text:
        g.stage = "matches"
    else:   # a single match's page: the first tab holds their tree, the Shared Matches tab their shared matches
        g.queue = [r["name"] for r in g.rows()]
        g.total = len(g.queue)
        g.stage = "pages" if re.search(r"How are .+ related\?", text) else "shared"
    ok, msg = g.on_copy(text)
    if ok and g.stage != "matches" and not re.search(r"How are .+ related\?", text):
        # their Shared Matches tab was copied. If their tree is still wanted, say so plainly instead of asking again in silence
        msg += " That was the Shared Matches tab. If you were asked for this person's tree, click the Trees tab on their page and copy that too."
    if ok and os.path.exists(os.path.join(folder, "matches.csv")):
        write_report(folder)
    return ok, msg


def save_problem_report(source, what, text, mask=True):
    """A file the researcher can email when a pasted page was not read. With names blanked, only the words the
    tool looks for stay readable, which is enough to see the page's layout and wording."""
    import time
    keep = Guide.KEEP | set("survived survivors include including wife husband spouse spouses son sons daughter daughters children child "
                            "brother brothers sister sisters siblings sibling parents mother father preceded predeceased death died born "
                            "married the late his her their mr mrs dr in on at is was were he she grandchildren great family members "
                            "memorial birth burial flowers inscription gravesite details added photos bio a an to with one two three four "
                            "five six seven eight nine ten half step law name names year years age aged".split())
    body = re.sub(r"[^\W\d_]+", lambda m: m.group() if m.group().lower() in keep else "x" * len(m.group()), text) if mask else text
    folder = os.path.join(os.path.expanduser("~"), "Downloads")
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "BirthFamilyTool-report-" + time.strftime("%Y%m%d-%H%M%S") + ".txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"Birth Family Tool {VERSION} problem report\nCopied from: {source}\nWhat went wrong: {what}\n"
                f"Names blanked: {'yes' if mask else 'no'}\n{'-' * 60}\n{body}\n")
    return path


def ancestry_couples(folder):
    """The shared ancestors Ancestry itself names on a match's page, read from the pages already copied:

        <match> could be your 1st cousin 1x removed through:
        <name>
        great-grandfather
        1868–1939
        View Relationship
        <name>
        great-grandmother
        ...

    Each gives a person's exact relationship to the adoptee, so the generation is known, not estimated.
    Layout taken from real copied pages (Nov 2026). Returns one entry per distinct couple."""
    side_of = {r["name"]: r.get("side") or "" for r in read_csv(os.path.join(folder, "matches.csv"))}
    cm_of = {r["name"]: float(r["cm"]) for r in read_csv(os.path.join(folder, "matches.csv"))}
    entry = re.compile(r"^([^\n]+)\n((?:(\d+)(?:st|nd|rd|th) )?((?:great-)*)(grand)?(father|mother))\n\s*\n(\d{4})?\s*[–—-]\s*(\d{4})?", re.M)
    found = {}
    pastes = os.path.join(folder, "pastes")
    for fn in sorted(os.listdir(pastes)) if os.path.isdir(pastes) else []:
        with open(os.path.join(pastes, fn), encoding="utf-8", errors="replace") as f:
            text = f.read()
        who, at = re.search(r"^You and (.+)$", text, re.M), text.find("could be your")
        if not who or at < 0:
            continue
        rel = re.search(r"could be your (.+?) through:", text[at:at + 200])
        block = text[at:at + 1500].split("Linked Tree")[0]
        people = []
        for m in entry.finditer(block):
            # father 1, grandfather 2, great-grandfather 3, "2nd great-grandfather" 4, and so on
            gen = (1 if not m.group(5) else 2 + (int(m.group(3)) if m.group(3) else m.group(4).count("great-")))
            people.append({"name": m.group(1).strip(), "rel": m.group(2), "gen": gen, "sex": "M" if m.group(6) == "father" else "F",
                           "b": m.group(7) or "", "d": m.group(8) or ""})
        match = who.group(1).strip()
        if not people or match not in side_of:
            continue
        key = tuple(sorted((p["name"], p["b"]) for p in people))
        hit = found.setdefault(key, {"people": people, "side": side_of[match], "matches": []})
        if match not in [x[0] for x in hit["matches"]]:
            hit["matches"].append((match, cm_of.get(match, 0), rel.group(1) if rel else ""))
    out = []
    for hit in found.values():
        hit["matches"].sort(key=lambda x: -x[1])
        hit["note"] = "Common ancestor named by Ancestry, shared with " + ", ".join(f"{n} ({cm:.0f} cM, {r})" if r else f"{n} ({cm:.0f} cM)" for n, cm, r in hit["matches"][:4])
        hit["matches"] = [n for n, _, _ in hit["matches"]]
        out.append(hit)
    return sorted(out, key=lambda h: (h["people"][0]["gen"], -len(h["matches"])))


_census_at = [0.0]


def census_1950(first, last, state=""):
    """Look a person up in the National Archives' 1950 census name index (1950census.archives.gov), the way its own
    search page does. The site's rules for programs allow every visitor and ask for a one-second pause between
    requests, which is kept here; each request says plainly that it comes from this tool. Only the one name being
    looked up is sent. Returns the census sheets on which that name appears, each with the names written on it:
    volunteer-corrected where someone has done that, otherwise as the handwriting software read them."""
    import json, time, urllib.parse, urllib.request
    time.sleep(max(0.0, 1.1 - (time.time() - _census_at[0])))
    q = {"page": 1, "size": 40, "name": f"{first} {last}".strip()}
    if re.fullmatch(r"[A-Za-z]{2}", state or ""):
        q["state"] = state.upper()
    req = urllib.request.Request("https://1950census.archives.gov/api/search?" + urllib.parse.urlencode(q),
                                 headers={"User-Agent": f"BirthFamilyTool/{VERSION} (genealogy research; one lookup at a time)"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            data = json.load(r)
    finally:
        _census_at[0] = time.time()
    want = {first.lower(), last.lower()} - {""}
    out = []
    for r in data.get("results", []):
        checked = [n for n in r.get("transcribe") or [] if n.get("name")]
        lines = checked or [n for n in r.get("names") or [] if n.get("name")]
        words = lambda s: set(re.sub(r"[^a-z ]", " ", s.lower()).split())
        hit = next((i for i, n in enumerate(lines) if want <= words(n["name"])), None)
        if hit is None:
            continue   # the index matches loosely; keep only sheets where this exact name is written
        out.append({"state": r.get("state") or "", "county": r.get("county") or "", "ed": r.get("ed") or "", "checked": bool(checked),
                    "at": hit, "lines": [{"row": n.get("row"), "name": n["name"]} for n in lines],
                    "view": "https://1950census.archives.gov/search/?" + urllib.parse.urlencode(
                        {"county": r.get("county") or "", "ed": r.get("ed") or "", "name": q["name"], "page": 1, "state": r.get("abbr") or ""})})
    return {"results": out[:12], "searched": q["name"], "state": q.get("state", "")}


def preview_pedigrees(folder):
    """The family trees inside the copied match pages, with who is whose parent.

    A match's tree preview lists up to five generations in a fixed order: a person, then the father's whole
    line, then the mother's (worked out from real copied previews, Nov 2026: in every full one the surnames run
    self, father, father's father, ... exactly as that order predicts). Only a FULL preview of 31 people can be
    read this way; in a shorter one somebody is missing and the positions can no longer be trusted, so it is
    left as a plain list. Living people show as "Private" and are skipped by the tree."""
    father, mother = {}, {}

    def lay(i, depth):
        if depth < 4:
            father[i], mother[i] = i + 1, i + 2 ** (4 - depth)
            lay(father[i], depth + 1)
            lay(mother[i], depth + 1)
    lay(0, 0)
    out = []
    for r in sorted(read_csv(os.path.join(folder, "matches.csv")), key=lambda r: -float(r["cm"])):
        p = os.path.join(folder, "trees", r.get("tree") or "")
        if not r.get("tree") or not os.path.exists(p) or p.lower().endswith(".ged"):
            continue
        with open(p, encoding="utf-8", errors="replace") as f:
            blocks = [b.strip() for b in re.split(r"\n\s*\n", f.read()) if b.strip()]
        people = []
        for b in blocks:
            lines = b.split("\n")
            if " is a preview of " in b or lines[0].lower().startswith("surnames"):
                continue
            name = re.sub(r"([A-Z][\w'-]+)\1$", r"\1", re.sub(r"\s*shared surname$", "", lines[0].strip()))
            years = re.search(r"(\d{4})?\s*[–—-]\s*(\d{4})?", lines[1]) if len(lines) > 1 else None
            people.append({"name": "" if name == "Private" else name, "b": (years and years.group(1)) or "", "d": (years and years.group(2)) or ""})
        entry = {"match": r["name"], "cm": round(float(r["cm"])), "side": r.get("side") or "", "people": people}
        if len(people) == 31:
            out.append({**entry, "father": father, "mother": mother})
        elif len(people) > 1:
            # A shorter preview: positions can't be trusted, but one link still can. In this listing a father comes
            # straight after his child, so a next entry with the same surname, 15 to 60 years older where both years
            # are known, is that person's father. Mothers can't be told apart here and are left out.
            sur = lambda p: p["name"].split()[-1].lower() if p["name"] else ""
            ok = lambda c, f: not (c["b"] and f["b"]) or 15 <= int(c["b"]) - int(f["b"]) <= 60
            near = {i: i + 1 for i in range(len(people) - 1)
                    if sur(people[i]) and sur(people[i]) == sur(people[i + 1]) and ok(people[i], people[i + 1])}
            if near:
                out.append({**entry, "father": near, "mother": {}})
    return out


def shared_ancestors(folder):
    """People who appear in the copied trees of two or more DNA matches, most recent first. The most recent
    such couple on a parent's side is the likely shared ancestor that side's matches descend from."""
    seen = {}
    for r in read_csv(os.path.join(folder, "matches.csv")):
        p = os.path.join(folder, "trees", r.get("tree") or "")
        if r.get("tree") and os.path.exists(p):
            for key in (gedcom_people(p) if p.lower().endswith(".ged") else text_people(p)):
                seen.setdefault(key, []).append(r)
    out = []
    for (first, sur, year), rows in seen.items():
        if len(rows) < 2 or not str(year).isdigit():
            continue
        sides = Counter(r.get("side") or "side not assigned" for r in rows)
        who = ", ".join(f"{r['name']} ({float(r['cm']):.0f} cM)" for r in sorted(rows, key=lambda r: -float(r["cm"]))[:5])
        out.append({"name": f"{first} {sur}".title(), "b": str(year), "n": len(rows), "side": sides.most_common(1)[0][0],
                    "matches": [r["name"] for r in rows],
                    "note": f"Shared ancestor: in {len(rows)} matches' trees, {sides.most_common(1)[0][0]} ({who})"})
    return sorted(out, key=lambda c: (-int(c["b"]), -c["n"]))[:40]


def pages_to_copy(folder):
    """Which Ancestry pages would add most to this case: matches with a public tree not yet copied, and the
    closest matches whose shared-match list has not been copied."""
    rows = sorted(read_csv(os.path.join(folder, "matches.csv")), key=lambda r: -float(r["cm"]))
    label = lambda r: f"{r['name']} ({float(r['cm']):.0f} cM)"
    seen = os.path.join(folder, "pages_done.txt")
    visited = set(open(seen, encoding="utf-8").read().splitlines()) if os.path.exists(seen) else set()   # pages already copied, tree or no tree
    # only trees big enough to be worth a copy (8 people or more), the same rule the collecting window uses, so the
    # tool does not come back later asking for small trees it chose to skip the first time
    size = lambda r: int((re.search(r"([\d,]+) people", r.get("tree_info") or "") or [0, "0"])[1].replace(",", ""))
    trees = [label(r) for r in rows if "public linked tree" in (r.get("tree_info") or "").lower() and size(r) >= 8
             and not r.get("tree") and r["name"] not in visited]
    s = os.path.join(folder, "shared.csv")
    done = {r["a"] for r in read_csv(s)} if os.path.exists(s) else set()
    shared = [label(r) for r in rows[:15] if r["name"] not in done]
    sides = Counter((r.get("side") or "").split("'")[0] for r in rows if float(r["cm"]) >= 200)
    thin = [s for s in ("Parent 1", "Parent 2") if not sides.get(s)]
    return {"trees": trees[:8], "more_trees": max(0, len(trees) - 8), "shared": shared[:8], "thin_sides": thin, "matches": len(rows)}


def dna_surnames(folder, surnames):
    """{SURNAME: {n, cm, side, who}}: how each surname in the tree shows up among the adoptee's DNA matches."""
    out = {}
    for s in surnames[:300]:
        s = re.sub(r"[^A-Za-z' -]", "", str(s)).strip()
        hits = surname_hits(folder, s) if len(s) > 2 else []
        sides = Counter(r.get("side") or "" for r, _ in hits)
        out[s.upper()] = {"n": len(hits), "cm": round(sum(float(r["cm"]) for r, _ in hits)),
                          "side": sides.most_common(1)[0][0] if hits else "",
                          "who": [f"{r['name']} ({float(r['cm']):.0f} cM)" for r, _ in hits[:6]]}
    return out


TREE_RULES = 2   # raise this (here and RULES in the tree page) when the rules for building a tree change


def tree_server(folder=None):
    """A page server reachable only from this computer and only with the key in the address it opens.
    Returns (server, address). It serves the tree page for the open case, saves the tree in the case
    folder, and answers the page's questions about DNA matches and the record files. Started without a
    case it also serves the home page, where the researcher picks or starts a case (see cmd_app)."""
    import contextlib, http.server, io, json, queue, secrets, time, urllib.parse
    key, at = secrets.token_urlsafe(16), {"folder": folder}
    saved = lambda: os.path.join(at["folder"], "tree_saved.json")
    case_of = lambda name: case_folder(re.sub(r"[^A-Za-z0-9 _-]", "", str(name)).strip() or "case")

    def cases():
        root, out = os.path.dirname(case_folder("x")), []
        for d in sorted(os.listdir(root)) if os.path.isdir(root) else []:
            f = os.path.join(root, d)
            if os.path.exists(os.path.join(f, "matches.csv")) or os.path.exists(os.path.join(f, "tree.html")):
                out.append({"name": d, "matches": len(read_csv(os.path.join(f, "matches.csv"))),
                            "tree": os.path.exists(os.path.join(f, "tree.html")), "born": adoptee(f).get("birth_year", ""),
                            "trees": len(os.listdir(os.path.join(f, "trees"))) if os.path.isdir(os.path.join(f, "trees")) else 0,
                            "shared": len(read_csv(os.path.join(f, "shared.csv"))),
                            "when": max(os.path.getmtime(os.path.join(f, n)) for n in os.listdir(f))})
        return sorted(out, key=lambda c: -c["when"])   # the case worked on most recently comes first

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def reply(self, body, kind="application/json", code=200):
            body = body if isinstance(body, bytes) else json.dumps(body).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", kind)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def allowed(self, u):
            # the key travels in the page address and in a header no other website can add to a request
            host = (self.headers.get("Host") or "").split(":")[0]
            page_key = dict(urllib.parse.parse_qsl(u.query)).get("k")
            return host in ("127.0.0.1", "localhost") and key in (page_key if u.path in ("/", "/home") else self.headers.get("X-Key"),)

        def do_GET(self):
            u = urllib.parse.urlparse(self.path)
            if u.path == "/mark.png":   # the faint logo behind the empty window; nothing private, so no key is asked for
                with open(os.path.join(getattr(sys, "_MEIPASS", program_dir()), "searchangels_mark.png"), "rb") as f:
                    return self.reply(f.read(), "image/png")
            if not self.allowed(u):
                return self.reply(b"Open the tree from the tool's menu.", "text/plain", 403)
            q = dict(urllib.parse.parse_qsl(u.query))
            srv.seen = time.time()
            if u.path == "/home" or (u.path == "/" and not at["folder"]):
                return self.reply(HOME_PAGE.encode("utf-8"), "text/html; charset=utf-8")
            if u.path == "/":
                # Always the tree screen of THIS version of the tool, with the case's data put into it. The copy of
                # the screen saved in the case folder dates from when the tree was first built; serving that file
                # as it stood meant an updated tool kept showing an old screen for every existing case.
                with open(os.path.join(at["folder"], "tree.html"), encoding="utf-8") as f:
                    data = f.read().split("const DATA = ", 1)[1].split(";\nconst KEY", 1)[0]
                return self.reply(TREE_PAGE.replace("/*DATA*/", data).encode("utf-8"), "text/html; charset=utf-8")
            if u.path == "/api/ping":
                return self.reply({"ok": True})
            if u.path == "/api/cases":
                return self.reply({"version": VERSION, "update": srv.update[0] if srv.update else "", "busy": srv.busy,
                                   "records": sorted(record_files()), "cases": cases()})
            if u.path == "/api/report":
                with open(write_report(case_of(q.get("name"))), encoding="utf-8") as f:
                    return self.reply({"text": f.read()})
            if u.path == "/api/tree":
                if at["folder"] and os.path.exists(saved()):
                    with open(saved(), "rb") as f:
                        raw = f.read()
                    try:
                        t = json.loads(raw)
                    except ValueError:
                        t = None
                    # A tree saved under older, looser rules (people put on twice, census households that did not
                    # belong) is set aside as a dated backup and the page builds it again from what was collected.
                    if isinstance(t, dict) and t.get("rules") != TREE_RULES:
                        import shutil
                        shutil.copyfile(saved(), saved()[:-5] + time.strftime(".before-fix-%Y%m%d-%H%M%S.json"))
                        t["outdated"] = True
                        return self.reply(t)
                    return self.reply(raw)
                return self.reply({})
            if u.path == "/api/records":
                return self.reply(record_lookup(q))
            if u.path == "/api/census":   # the 1950 census name index, one person at a time, when the researcher asks
                try:
                    return self.reply(census_1950(q.get("first", ""), q.get("last", ""), q.get("state", "")))
                except Exception as e:
                    return self.reply({"results": [], "error": f"The 1950 census site could not be reached ({type(e).__name__})."})
            if u.path == "/api/pulse" and at["folder"]:   # has anything new been copied into this case?
                stamp = 0
                for root_, _, names in os.walk(at["folder"]):
                    if os.path.basename(root_) in ("pastes",):
                        continue
                    stamp += sum(int(os.path.getmtime(os.path.join(root_, n))) for n in names
                                 if n.endswith(".csv") or os.path.basename(root_) == "trees")
                return self.reply({"busy": srv.busy, "stamp": stamp})
            if u.path == "/api/base" and at["folder"]:    # the starting tree worked out afresh from everything copied so far
                who = os.path.join(at["folder"], "adoptee_name.txt")
                name = open(who, encoding="utf-8").read().strip() if os.path.exists(who) else "Adoptee"
                with open(cmd_skeleton(at["folder"], name), encoding="utf-8") as f:
                    return self.reply(json.loads(f.read().split("const DATA = ", 1)[1].split(";\nconst KEY", 1)[0]))
            if u.path == "/api/leads" and at["folder"]:   # what the latest import adds: shared ancestors, and pages still worth copying
                top = sorted(read_csv(os.path.join(at["folder"], "matches.csv")), key=lambda r: -float(r["cm"]))
                closest = {s: next((r["name"] for r in top if (r.get("side") or "").startswith(s) and rel_shape(r.get("relationship", ""))), "")
                           for s in ("Parent 1", "Parent 2")}   # the match the starting tree hangs on each side
                return self.reply({"ancestors": shared_ancestors(at["folder"]), "pages": pages_to_copy(at["folder"]), "closest": closest,
                                   "couples": ancestry_couples(at["folder"]), "pedigrees": preview_pedigrees(at["folder"]),
                                   "name": os.path.basename(at["folder"])})
            self.reply({}, code=404)

        def do_POST(self):
            u = urllib.parse.urlparse(self.path)
            if not self.allowed(u):
                return self.reply({}, code=403)
            data = json.loads(self.rfile.read(int(self.headers.get("Content-Length") or 0)) or b"null")
            if u.path == "/api/close" and isinstance(data, dict):
                srv.next = "again" if data.get("next") == "again" else ""
                self.reply({"ok": True})
                return srv.done.set()
            if u.path == "/api/tree" and at["folder"] and isinstance(data, dict) and isinstance(data.get("people"), list):
                with open(saved() + ".part", "w", encoding="utf-8") as f:
                    json.dump(data, f)
                os.replace(saved() + ".part", saved())   # never leave a half-written tree behind
                return self.reply({"ok": True})
            if u.path == "/api/autofill":
                return self.reply(auto_fill(data))
            if u.path == "/api/paste" and at["folder"] and isinstance(data, dict):
                ok, msg = paste_page(at["folder"], str(data.get("text") or ""))
                return self.reply({"ok": ok, "msg": msg})
            if u.path == "/api/report" and isinstance(data, dict):
                return self.reply({"path": save_problem_report(str(data.get("source") or ""), str(data.get("what") or ""),
                                                               str(data.get("text") or ""), bool(data.get("mask", True)))})
            if u.path == "/api/dna" and at["folder"] and isinstance(data, list):
                return self.reply(dna_surnames(at["folder"], data))
            if u.path == "/api/open" and isinstance(data, dict):       # home page: open a case's tree, building it if asked
                f = case_of(data.get("name"))
                if not os.path.exists(os.path.join(f, "matches.csv")) and not data.get("adoptee"):
                    return self.reply({"need": "capture"})
                if str(data.get("year") or "").isdigit():
                    cmd_adoptee(f, data["year"], str(data.get("place") or ""))
                if data.get("adoptee"):   # a brand-new case starts as the adoptee under two unknown parents, before any match is copied
                    os.makedirs(f, exist_ok=True)
                    with open(os.path.join(f, "adoptee_name.txt"), "w", encoding="utf-8") as fh:
                        fh.write(str(data["adoptee"]).strip())
                    cmd_skeleton(f, str(data["adoptee"]))
                    if os.path.exists(os.path.join(f, "tree_saved.json")):
                        os.replace(os.path.join(f, "tree_saved.json"), os.path.join(f, "tree_saved.previous.json"))
                if not os.path.exists(os.path.join(f, "tree.html")):
                    # Matches were collected but no tree was ever built. Don't stop to ask questions: build it from
                    # everything already collected. The adoptee's name and birth year are filled in on the tree itself.
                    who = os.path.join(f, "adoptee_name.txt")
                    cmd_skeleton(f, open(who, encoding="utf-8").read().strip() if os.path.exists(who) else "Adoptee")
                at["folder"] = f
                return self.reply({"ok": True})
            if u.path == "/api/export" and isinstance(data, dict):     # tree page: save the GEDCOM where the researcher will find it
                folder = os.path.join(os.path.expanduser("~"), "Downloads")
                os.makedirs(folder, exist_ok=True)
                path = os.path.join(folder, re.sub(r"[^A-Za-z0-9 ._-]", "", str(data.get("name"))) or "birth-family-tree.ged")
                with open(path, "w", encoding="utf-8", newline="\n") as f:
                    f.write(str(data.get("text") or ""))
                return self.reply({"path": path})
            if u.path == "/api/capture" and isinstance(data, dict):    # home page: open the copy-and-paste window for a case
                srv.jobs.put(case_of(data.get("name")))
                return self.reply({"ok": True})
            if u.path == "/api/check" and isinstance(data, dict):      # home page: does a cM amount fit a relationship?
                text = io.StringIO()
                try:
                    with contextlib.redirect_stdout(text):
                        cmd_check(float(str(data.get("cm")).replace(",", "")), str(data.get("rel")))
                except (SystemExit, ValueError) as e:
                    text.write((str(e) if isinstance(e, SystemExit) and e.code else "")
                               or "Need a number of cM and a relationship such as 1C1R, 2C, half 1C, aunt.")
                return self.reply({"text": text.getvalue()})
            if u.path == "/api/update" and srv.update:
                install_update(srv.update[1])
                self.reply({"ok": True})
                return srv.jobs.put(None)
            if u.path == "/api/quit":
                self.reply({"ok": True})
                return srv.jobs.put(None)
            self.reply({}, code=404)

    import threading
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    srv.done, srv.next = threading.Event(), ""   # set when the page's "Save and close" is used
    srv.jobs, srv.seen, srv.busy, srv.update = queue.Queue(), time.time(), False, None   # used by cmd_app
    return srv, f"http://127.0.0.1:{srv.server_address[1]}/?k={key}"


def cmd_app():
    """The windowed way to run the tool: no command line. Opens the home page in the browser; the researcher
    picks or starts a case there. This thread only waits, and opens the copy-and-paste window when asked."""
    import queue, subprocess, threading, time, webbrowser
    srv, address = tree_server(None)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    threading.Thread(target=lambda: setattr(srv, "update", newer_version()), daemon=True).start()
    home = address.replace("/?k=", "/home?k=")
    try:
        import webview   # the app's own window (uses the WebView2 that ships with Windows 10 and 11)
    except ImportError:
        webview = None
    if webview and os.environ.get("BIRTHFAMILY_WINDOW") != "0":
        def jobs():
            while True:
                job = srv.jobs.get()
                if job is None:          # Quit on the home page, or an update: close the window
                    for w in list(webview.windows):
                        w.destroy()
                    return
                srv.busy = True
                try:                     # the copy-and-paste window runs as its own small program beside the browser
                    subprocess.run([sys.executable] + ([] if FROZEN else [os.path.abspath(__file__)]) + ["guide", job])
                    if os.path.exists(os.path.join(job, "matches.csv")):
                        write_report(job)
                finally:
                    srv.busy = False
        threading.Thread(target=jobs, daemon=True).start()
        webview.settings["ALLOW_DOWNLOADS"] = True
        webview.create_window("Birth Family Tool", home, width=1320, height=900, min_size=(940, 620))
        webview.start()                  # returns when the window is closed
        return srv.shutdown()
    webbrowser.open(home)                # no window library on this computer: fall back to the web browser
    while True:
        try:
            job = srv.jobs.get(timeout=30)
        except queue.Empty:
            # ponytail: pages check in every 20 seconds; none for 30 minutes means the browser was closed without Quit
            if time.time() - srv.seen > 1800:
                break
            continue
        if job is None:
            break
        srv.busy = True
        try:
            cmd_guide(job, 90.0)   # the capture window must run on this main thread
            if os.path.exists(os.path.join(job, "matches.csv")):
                write_report(job)
        finally:
            srv.busy, srv.seen = False, time.time()
    time.sleep(0.5)   # let the last reply reach the page
    srv.shutdown()


def show_tree(folder):
    """Open the tree page for this case in the browser and keep answering it until the researcher is done."""
    import threading, webbrowser
    srv, address = tree_server(folder)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    found = record_files()
    print("\nThe tree is opening in your browser. Your changes are saved in the case folder as you make them.")
    print("Record files found: " + (", ".join(sorted(found)) if found else "none (the tree still works; record lookups will be empty)"))
    webbrowser.open(address)
    print("LEAVE THIS WINDOW OPEN while you work on the tree.\n"
          "When you have finished, choose File, then Save and close, on the tree page (or press Enter here).")
    try:
        import msvcrt   # Windows: notice either the page's Save and close or the Enter key, whichever comes first
        while not srv.done.wait(0.2):
            if msvcrt.kbhit() and msvcrt.getwch() in "\r\n":
                break
    except ImportError:
        ask("")
    srv.shutdown()
    print("Tree saved and closed.")
    return srv.next


HOME_PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Birth Family Tool</title>
<style>
:root { --bg:#fff; --fg:#1b1b1b; --mute:#666; --card:#f4f6f8; --edge:#c9d1d9; --me:#dff5e3; --btn:#1a7f37; --warn:#fff7e0; --bar:#f3f3f3; --hover:#dbe7f5; }
@media (prefers-color-scheme: dark) { :root { --bg:#16181b; --fg:#e8e8e8; --mute:#a0a0a0; --card:#24282d; --edge:#444c56; --me:#1d3b25; --btn:#2ea043; --warn:#3d3417; --bar:#22252a; --hover:#2f3b4d; } }
html, body { height:100%; }
body { background:var(--bg); color:var(--fg); font:15px/1.45 "Segoe UI", system-ui, sans-serif; margin:0; display:flex; flex-direction:column; }
h2 { font-size:19px; margin:0 0 8px; } h3 { font-size:15px; margin:14px 0 4px; } .mute { color:var(--mute); }
button { font:inherit; padding:7px 12px; border-radius:6px; border:1px solid var(--edge); background:var(--card); color:var(--fg); cursor:pointer; }
button.go { background:var(--btn); color:#fff; border-color:var(--btn); font-weight:600; }
input { font:inherit; padding:6px; background:var(--bg); color:var(--fg); border:1px solid var(--edge); border-radius:4px; }
/* a plain menu bar across the top, like any Windows program */
.menubar { display:flex; align-items:center; background:var(--bar); border-bottom:1px solid var(--edge); padding:0 4px; user-select:none; }
.menu { position:relative; } .menu summary { list-style:none; cursor:pointer; padding:5px 11px; border-radius:4px; }
.menu summary::-webkit-details-marker { display:none; } .menu[open] summary, .menu summary:hover { background:var(--hover); }
.menu > div { position:absolute; z-index:5; top:100%; left:0; min-width:250px; background:var(--bg); border:1px solid var(--edge); border-radius:0 0 6px 6px; box-shadow:0 6px 18px rgba(0,0,0,.2); padding:4px; }
.menu > div button { display:block; width:100%; text-align:left; border:0; background:none; border-radius:4px; padding:6px 10px; }
.menu > div button:hover { background:var(--hover); } .menu hr { border:0; border-top:1px solid var(--edge); margin:4px 0; }
.menubar .ver { margin-left:auto; padding:0 10px; color:var(--mute); font-size:13px; }
/* the empty work area: a faint Search Angels mark until a case is open */
main { flex:1; display:flex; flex-direction:column; align-items:center; justify-content:center; text-align:center; padding:20px; }
main img { width:min(38vh, 300px); opacity:.16; } main p { color:var(--mute); max-width:520px; }
#banner { width:100%; max-width:760px; }
.note { background:var(--warn); border:1px solid var(--edge); border-radius:6px; padding:10px 12px; margin:8px 0; text-align:left; }
table { border-collapse:collapse; width:100%; } td, th { text-align:left; padding:8px 10px 8px 0; border-bottom:1px solid var(--edge); vertical-align:middle; }
.bar { display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin:8px 0; }
dialog { background:var(--bg); color:var(--fg); border:1px solid var(--edge); border-radius:8px; width:min(760px, 94vw); box-shadow:0 10px 40px rgba(0,0,0,.35); }
dialog::backdrop { background:rgba(0,0,0,.25); }
dialog label { display:block; margin:8px 0 2px; font-size:13px; color:var(--mute); } dialog label input { width:100%; box-sizing:border-box; }
#start { width:min(520px, 92vw); } #start button { display:block; width:100%; text-align:left; margin:8px 0; padding:11px 14px; } #start small { display:block; color:var(--mute); font-weight:400; }
#start button.go small { color:#e6f4ea; }
pre { white-space:pre-wrap; font:13px/1.4 Consolas, monospace; max-height:60vh; overflow:auto; text-align:left; }
details.q { margin:6px 0; } details.q summary { cursor:pointer; font-weight:600; }
</style></head><body>
<nav class="menubar">
  <details class="menu"><summary>File</summary><div>
    <button data-do="start">Start window...</button>
    <button data-do="cases">Open a case...</button>
    <button data-do="new">New case...</button>
    <hr><button data-do="quit">Exit</button>
  </div></details>
  <details class="menu"><summary>Tools</summary><div>
    <button data-do="cm">Does a cM amount fit a relationship?...</button>
    <button data-do="records">Record files on this computer...</button>
  </div></details>
  <details class="menu"><summary>Help</summary><div>
    <button data-do="how">How to use the tool</button>
    <button data-do="about">About and source credits</button>
  </div></details>
  <span class="ver" id="ver"></span>
</nav>
<main>
  <div id="banner"></div>
  <img src="/mark.png" alt="">
  <p>Open a case or start a new one from the <strong>File</strong> menu.</p>
</main>

<dialog id="start">
  <h2>Birth Family Tool</h2>
  <p class="mute">Where would you like to start?</p>
  <div id="start_rows"></div>
</dialog>
<dialog id="cases"><h2>Open a case</h2><div id="case_rows"></div><div class="bar"><button data-close>Close</button></div></dialog>
<dialog id="new"><form method="dialog">
  <h2>New case</h2>
  <p class="mute">The tree opens straight away with the adoptee under two unknown parents. A small window then tells you which Ancestry pages to open and copy, and the tree fills in as you go.</p>
  <label>Case name, for example the adoptee's surname <input id="new_name"></label>
  <label>Adoptee's name as it should appear in the tree <input id="new_who"></label>
  <label>Adoptee's birth year, if known <input id="new_year" inputmode="numeric"></label>
  <label>State or place of birth, if known <input id="new_place"></label>
  <div class="bar"><button class="go" value="ok">Open the tree and start collecting</button><button value="cancel">Cancel</button></div>
</form></dialog>
<dialog id="who"><form method="dialog">
  <h2 id="who_head"></h2>
  <label>Adoptee's name as it should appear in the tree <input id="who_name"></label>
  <label>Adoptee's birth year, if known <input id="who_year" inputmode="numeric"></label>
  <label>State or place of birth, if known <input id="who_place"></label>
  <div class="bar"><button class="go" value="ok">Build the tree</button><button value="cancel">Cancel</button></div>
</form></dialog>
<dialog id="cm"><h2>Does a cM amount fit a relationship?</h2>
  <div class="bar"><input id="cm_n" placeholder="Shared cM, e.g. 497" size="18"> <input id="cm_rel" placeholder="Relationship, e.g. 1C1R, 2C, half 1C, aunt" size="36"> <button class="go" id="cm_go">Check</button></div>
  <pre id="cm_out"></pre><div class="bar"><button data-close>Close</button></div></dialog>
<dialog id="rep"><h2 id="rep_head"></h2><pre id="rep_text"></pre><div class="bar"><button data-close>Close</button></div></dialog>
<dialog id="info"><h2 id="info_head"></h2><div id="info_body" style="text-align:left"></div><div class="bar"><button data-close>Close</button></div></dialog>

<script>
const KEY = new URLSearchParams(location.search).get("k"), $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const api = (path, body) => fetch(path, {method: body ? "POST" : "GET", headers: {"X-Key": KEY}, body: body ? JSON.stringify(body) : undefined}).then(r => r.json());
let pending = null, R = {cases: [], records: []}, greeted = false;
const HOW = [
  ["Start a case", "File, New case. Type a short name. A small window opens beside your browser and tells you which Ancestry pages to open and copy with Ctrl+A, Ctrl+C. Close it when you have finished or want a break; your progress is kept."],
  ["Open the tree", "File, Open a case, then Open the tree. The first time, it asks for the adoptee's name and birth year. The tree shows the adoptee under two parents with each parent's line above."],
  ["Add what you find", "On the tree choose Add, Paste a page I copied. Paste any Ancestry match list or match page, a Find a Grave memorial or an obituary. The tree grows from each one."],
  ["Follow the DNA", "Under each family the tree shows which DNA matches carry that surname. When a child's spouse's family has matches, press the green button to move the adoptee's line under that couple."],
  ["Send the tree to Ancestry", "On the tree choose File, Export GEDCOM file. The steps for uploading it and linking the adoptee's DNA are shown straight afterwards."],
  ["When a page is not read", "In the paste box open 'Send us a report'. It saves a file with names blanked. Email that file to support@searchangels.org and an update will follow."],
  ["Keep case data private", "Cases hold real names. They are kept on this computer only. Delete a case's folder when the case is closed."]];
const ABOUT = `<p>Birth Family Tool, for Search Angels researchers. It never visits Ancestry or any site you copy from: you browse and copy, and it reads what you paste. Suggestions and problems: support@searchangels.org</p>
<h3>Sources and credits</h3><ul>
<li><strong>Historic newspapers to 1963:</strong> Library of Congress, <em>Chronicling America</em>. Searched only when you ask.</li>
<li><strong>1950 census:</strong> National Archives and Records Administration, 1950census.archives.gov. Public domain. Looked up only when you ask.</li>
<li><strong>Social Security applications and claims</strong> (the file runs 1936 to 2007, but it is thin before the 1980s: about nine in ten of its death records are from 1980 or later, and almost none are before 1960): Wray, Anthony; U.S. Federal Security Agency, Social Security Board; and U.S. Social Security Administration. "Numerical Identification Files (NUMIDENT), 1936 to 2007." ICPSR [distributor], doi.org/10.3886/E207202V3. Creative Commons Attribution 4.0. Reduced to names, dates, places and parents' names; Social Security numbers removed. Also cite Hollingsworth, Karbownik, Thomasson and Wray (2024), <em>American Economic Review</em> 114 (7): 2201 to 38.</li>
<li><strong>Texas marriage and divorce indexes, 1966 to 2024:</strong> Texas Department of State Health Services, Vital Statistics Section.</li>
<li><strong>Deceased veterans (BIRLS death file):</strong> U.S. Department of Veterans Affairs, obtained and published by Reclaim The Records, birls.org. Public domain.</li>
<li><strong>Veterans cemetery burials:</strong> U.S. Department of Veterans Affairs, National Cemetery Administration, Gravesite Locator. CC0 public domain.</li>
<li><strong>Read from text you paste, and not stored:</strong> Ancestry match pages, Find a Grave memorial pages, and obituaries.</li></ul>`;
const show = id => { document.querySelectorAll("dialog[open]").forEach(d => d.close()); $(id).showModal(); };
const info = (head, html) => { $("info_head").textContent = head; $("info_body").innerHTML = html; show("info"); };
function caseRows() {
  $("case_rows").innerHTML = R.cases.length ? "<table><tr><th>Case</th><th>DNA matches</th><th></th></tr>" + R.cases.map(c =>
    `<tr><td><strong>${esc(c.name)}</strong>${c.born ? ` <span class="mute">born ${esc(c.born)}</span>` : ""}</td><td>${c.matches}</td><td><div class="bar">`
    + `<button class="go" data-do="open" data-name="${esc(c.name)}">${c.tree ? "Open the tree" : "Build the tree"}</button>`
    + (c.tree ? `<button data-do="rebuild" data-name="${esc(c.name)}">Rebuild</button>` : "")
    + `<button data-do="report" data-name="${esc(c.name)}">Report</button><button data-do="capture" data-name="${esc(c.name)}">Collect more matches</button></div></td></tr>`).join("") + "</table>"
    : '<p class="mute">No cases yet. Use File, New case.</p>';
}
function startRows() {
  const last = R.cases[0];   // cases arrive most recently worked first
  $("start_rows").innerHTML =
    (last ? `<button class="go" data-do="open" data-name="${esc(last.name)}">Continue with ${esc(last.name)}<small>Everything you collected is kept: ${last.matches} DNA matches, ${last.trees} match trees, ${last.shared} shared-match links. ${last.tree ? "Opens the tree where you left off." : "Builds the tree from them now."}</small></button>` : "")
    + (R.cases.length > 1 ? `<button data-do="cases">Open another case<small>${R.cases.length} cases on this computer</small></button>` : "")
    + `<button ${last ? "" : 'class="go" '}data-do="new">Start a new case<small>Collect an adoptee's DNA matches from Ancestry</small></button>`
    + `<button data-do="how">Show me how it works<small>A short guide to each step</small></button>`
    + `<p class="mute" style="text-align:right;margin:4px 0 0"><a href="#" data-close>Close this and use the menus</a></p>`;
}
function load() {
  return api("/api/cases").then(r => {
    R = r; $("ver").textContent = "Version " + r.version;
    $("banner").innerHTML = (r.update ? `<div class="note">Version ${esc(r.update)} is available. <button data-do="update">Update now</button> The tool closes and reopens.</div>` : "")
      + (r.busy ? '<div class="note"><strong>The collecting window is open.</strong> Follow what it says. When you close it, come back here and choose File, Open a case.</div>' : "");
    caseRows(); startRows();
    if (!greeted && !r.busy) { greeted = true; $("start").showModal(); }
  }, () => { $("banner").innerHTML = '<p class="note">The tool is no longer running. Close this window and start it again.</p>'; });
}
// Opening a case always lands on its tree. `collect` also opens the collecting window beside it.
let collectAfter = false;
function open_(name, extra, collect) {
  collectAfter = !!collect;
  api("/api/open", {name, ...extra}).then(r => {
    if (r.ok) return (collectAfter ? api("/api/capture", {name}) : Promise.resolve()).then(() => location = "/?app=1&k=" + KEY);
    ask(name, r.need === "adoptee" ? "Before the first tree for " + name : "Set up " + name);   // either way the tree needs to know who the adoptee is
  });
}
function ask(name, head) { pending = name; $("who_head").textContent = head; show("who"); }
$("who").addEventListener("close", () => {
  if ($("who").returnValue !== "ok" || !pending || !$("who_name").value.trim()) return;
  open_(pending, {adoptee: $("who_name").value.trim(), year: $("who_year").value.trim(), place: $("who_place").value.trim()}, collectAfter);
});
$("new").addEventListener("close", () => {
  const n = $("new_name").value.trim(), w = $("new_who").value.trim();
  if ($("new").returnValue !== "ok" || !n) return;
  open_(n, {adoptee: w || "Adoptee", year: $("new_year").value.trim(), place: $("new_place").value.trim()}, true);
});
function capture(name) { document.querySelectorAll("dialog[open]").forEach(d => d.close()); open_(name, {}, true); }
document.addEventListener("click", e => {
  document.querySelectorAll(".menu[open]").forEach(m => { if (!m.contains(e.target) || e.target.tagName === "BUTTON") m.open = false; });
  if (e.target.closest("[data-close]")) return e.target.closest("dialog").close();
  const b = e.target.closest("[data-do]"); if (!b) return;
  const act = b.dataset.do, name = b.dataset.name;
  if (act === "start") return load().then(() => show("start"));
  if (act === "cases") return load().then(() => show("cases"));
  if (act === "new") { ["new_name", "new_who", "new_year", "new_place"].forEach(i => $(i).value = ""); return show("new"); }
  if (act === "cm") return show("cm");
  if (act === "how") return info("How to use the tool", HOW.map(([q, a], i) => `<details class="q"${i ? "" : " open"}><summary>${q}</summary><p>${a}</p></details>`).join(""));
  if (act === "about") return info("About", ABOUT);
  if (act === "records") return info("Record files on this computer", R.records.length ? "<p>Found: " + R.records.map(esc).join(", ") + ".</p><p class='mute'>These are the public record files the tree looks people up in.</p>"
                                      : "<p>No record files are installed on this computer. The tree works without them; record lookups will be empty.</p>");
  if (act === "open") return open_(name, {});
  if (act === "rebuild") return ask(name, "Rebuild the tree for " + name + ". The tree you have now is kept as a backup file in the case folder.");
  if (act === "capture") return capture(name);
  if (act === "report") { $("rep_head").textContent = "Report for " + name; $("rep_text").textContent = "Working..."; show("rep"); return api("/api/report?name=" + encodeURIComponent(name)).then(r => $("rep_text").textContent = r.text); }
  if (act === "update") return api("/api/update", {}).then(() => document.body.innerHTML = "<main><h2>Updating</h2><p>The tool will close and reopen.</p></main>");
  if (act === "quit") return api("/api/quit", {});
});
$("cm_go").onclick = () => api("/api/check", {cm: $("cm_n").value, rel: $("cm_rel").value}).then(r => $("cm_out").textContent = r.text);
setInterval(() => api("/api/ping").catch(() => {}), 20000);   // tells the tool this window is still open
load();
</script></body></html>
"""


TREE_PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Birth family tree</title>
<style>
:root { --bg:#fff; --fg:#1b1b1b; --mute:#666; --line:#999; --card:#f4f6f8; --edge:#c9d1d9; --unk:#fff7e0; --me:#dff5e3; --hit:#e3eeff; --btn:#1a7f37; }
@media (prefers-color-scheme: dark) { :root { --bg:#16181b; --fg:#e8e8e8; --mute:#a0a0a0; --line:#777; --card:#24282d; --edge:#444c56; --unk:#3d3417; --me:#1d3b25; --hit:#1c2f4d; --btn:#2ea043; } }
body { background:var(--bg); color:var(--fg); font:15px/1.45 "Segoe UI", system-ui, sans-serif; margin:0; padding:16px; }
h1 { font-size:18px; margin:0 0 4px; } .goal summary { cursor:pointer; } p { margin:6px 0; max-width:900px; } .mute { color:var(--mute); }
button { font:inherit; padding:7px 12px; border-radius:6px; border:1px solid var(--edge); background:var(--card); color:var(--fg); cursor:pointer; }
button.go { background:var(--btn); color:#fff; border-color:var(--btn); font-weight:600; }
.bar { display:flex; flex-wrap:wrap; gap:8px; margin:12px 0; align-items:center; }
.key span { display:inline-block; padding:1px 8px; border-radius:4px; border:1px solid var(--edge); margin-right:6px; font-size:13px; }
.scroll { overflow-x:auto; padding-bottom:12px; }
.tree { display:inline-block; min-width:100%; margin:8px 0 20px; }
.tree ul { display:flex; justify-content:center; padding:20px 0 0; margin:0; position:relative; }
.tree li { list-style:none; position:relative; padding:20px 6px 0; display:flex; flex-direction:column; align-items:center; }
.tree li::before, .tree li::after { content:""; position:absolute; top:0; right:50%; width:50%; height:20px; border-top:1px solid var(--line); }
.tree li::after { right:auto; left:50%; border-left:1px solid var(--line); }
.tree li:only-child::before, .tree li:only-child::after { border-top:0; }
.tree li:first-child::before, .tree li:last-child::after { border-top:0; }
.tree li:last-child::before { border-right:1px solid var(--line); } .tree li:last-child::after { border-left:0; }
.tree ul ul::before { content:""; position:absolute; top:0; left:50%; height:20px; border-left:1px solid var(--line); }
.tree > ul { padding-top:0; } .tree > ul > li { padding-top:0; } .tree > ul > li::before, .tree > ul > li::after { display:none; }
.couple { display:flex; gap:4px; align-items:stretch; }
.p { background:var(--card); border:1px solid var(--edge); border-radius:6px; padding:5px 9px; min-width:96px; max-width:170px; text-align:center; cursor:pointer; }
.p:hover, .p:focus { outline:2px solid var(--btn); }
.p.unk { background:var(--unk); font-style:italic; } .p.me { background:var(--me); font-weight:700; } .p.hit { background:var(--hit); }
.p small { display:block; color:var(--mute); font-style:normal; font-weight:400; font-size:12px; }
dialog { background:var(--bg); color:var(--fg); border:1px solid var(--edge); border-radius:8px; width:min(420px, 92vw); }
.next { background:var(--unk); border:1px solid var(--edge); border-radius:6px; padding:10px 14px; max-width:900px; margin:8px 0; } .next ol { margin:6px 0 0; padding-left:20px; } .next li { margin:4px 0; }
#prog { display:flex; align-items:center; gap:12px; margin:0 0 10px; max-width:900px; } #prog[hidden] { display:none; }
#prog_bar { flex:0 0 260px; height:12px; border-radius:6px; background:var(--edge); overflow:hidden; } #prog_fill { height:100%; width:0; background:var(--btn); transition:width .5s; }
#prog_txt { font-weight:600; } #ob_next { background:var(--me); border:1px solid var(--edge); border-radius:6px; padding:8px 10px; margin:0 0 8px; } #ob_next:empty { display:none; }
/* the paste box floats over the tree: it can be dragged anywhere and resized, and the tree behind keeps working */
#ob { position:fixed; margin:0; left:auto; right:16px; top:56px; width:min(560px, 92vw); max-height:calc(100vh - 80px); overflow:auto; resize:both; z-index:8; padding:0 14px 12px; box-shadow:0 8px 30px rgba(0,0,0,.35); }
#ob_head { position:sticky; top:0; margin:0 -14px 8px; padding:8px 14px; background:var(--card); border-bottom:1px solid var(--edge); cursor:move; font-weight:600; user-select:none; touch-action:none; }
#ob_head .mute { font-weight:400; font-size:12px; margin-left:8px; } #ob_x { float:right; border:0; background:none; padding:0 4px; }
#mini { position:fixed; left:10px; bottom:10px; z-index:4; background:var(--bg); border:1px solid var(--edge); border-radius:8px; padding:6px 8px; box-shadow:0 2px 10px rgba(0,0,0,.25); font-size:12px; }
#mini:empty { display:none; } #mini line { stroke:var(--line); stroke-width:1; } #mini circle { cursor:pointer; stroke:var(--bg); stroke-width:1; }
#mini .me { fill:var(--btn); } #mini .dna { fill:#2f6fd0; } #mini .known { fill:var(--mute); } #mini .unk { fill:var(--edge); } #mini circle.here { stroke:#e0a800; stroke-width:3; }
.minihead { color:var(--mute); margin-bottom:2px; } .minikey i { display:inline-block; width:9px; height:9px; border-radius:50%; margin:0 3px 0 8px; } .minikey i.dna { background:#2f6fd0; margin-left:0; } .minikey i.known { background:var(--mute); } .minikey i.unk { background:var(--edge); }
body { padding-bottom:210px !important; }
@media (min-width: 1150px) { body { padding-left:270px !important; padding-bottom:16px !important; } }   /* the map gets its own left margin */
.pedigree { display:inline-block; min-width:100%; text-align:center; padding:8px 0; }
/* the sheet the tree is drawn on: fixed window, movable and zoomable contents */
.view { position:relative; height:62vh; min-height:340px; overflow:hidden; border:1px solid var(--edge); border-radius:8px; background:var(--card); cursor:grab; touch-action:none; user-select:none; }
.view:active { cursor:grabbing; } .view .pedigree { position:absolute; left:0; top:0; min-width:0; width:max-content; transform-origin:0 0; padding:12px; }
.view .p { background:var(--bg); } .view .p.unk { background:var(--unk); } .view .p.me { background:var(--me); } .view .p.hit { background:var(--hit); }
button:disabled { opacity:.45; cursor:default; } .todo { background:var(--hit); border-radius:6px; padding:6px 9px; } .done { background:var(--me); border-radius:6px; padding:6px 9px; }
.viewbar { margin:4px 0 6px; } .viewbar .mute { font-size:13px; }
.sibs { display:flex; align-items:flex-start; gap:8px; } .stem { display:flex; flex-direction:column; align-items:center; align-self:stretch; }
.stem::after { content:""; flex:1; min-height:10px; border-left:2px solid var(--line); }
/* every link is a drawn line: person down to the couple's bar, bar down to the children's bar, a tick to each child */
.par { padding-bottom:0; } .pedigree > .anc > .stem::after, .pedigree > .anc > .sibs > .stem::after { display:none; }
.anc > .sibs { border-top:2px solid var(--line); } .sibs > .desc::before, .sibs > .stem::before { content:""; height:9px; border-left:2px solid var(--line); }
.anc > .sibs, .sibs.below { margin-top:0; padding-top:0; } .par { margin-bottom:9px; } .desc > .couple::after { content:none; }
.desc:has(> .sibs.below) > .couple { margin-bottom:0; } .desc > .sibs.below { margin-top:9px; position:relative; }
.desc > .sibs.below::before { content:""; position:absolute; left:50%; top:-11px; height:9px; border-left:2px solid var(--line); }
.desc { display:flex; flex-direction:column; align-items:center; } .sibs.below { border-top:2px solid var(--line); margin-top:6px; padding-top:6px; }
.anc { display:inline-flex; flex-direction:column; align-items:center; vertical-align:bottom; padding:0 3px; }
.par { display:flex; align-items:flex-end; justify-content:center; gap:4px; padding-bottom:9px; margin-bottom:9px; border-bottom:2px solid var(--line); position:relative; }
.par::after { content:""; position:absolute; left:50%; bottom:-11px; height:9px; border-left:2px solid var(--line); }
.more { font-size:12px; color:var(--mute); cursor:pointer; margin-bottom:4px; text-decoration:underline; }
.p.focus { outline:3px solid var(--btn); } .p .dna { color:var(--btn); font-weight:600; }
.fam { border:1px solid var(--edge); border-radius:6px; padding:8px 10px; margin:8px 0; } .kids { display:flex; flex-wrap:wrap; gap:10px; }
.trail a { color:inherit; } h2 { font-size:17px; margin:18px 0 4px; } h3 { font-size:15px; margin:14px 0 4px; }
table { border-collapse:collapse; } td { padding:4px 10px 4px 0; vertical-align:top; border-bottom:1px solid var(--edge); }
/* the same plain menu bar as the start window, fixed across the top */
.menubar { position:fixed; top:0; left:0; right:0; z-index:6; display:flex; align-items:center; background:var(--card); border-bottom:1px solid var(--edge); padding:0 4px; user-select:none; }
body { padding-top:46px !important; }
.menu { position:relative; } .menu summary, .menubar > button { list-style:none; cursor:pointer; padding:5px 11px; border:0; border-radius:4px; background:none; }
.menu summary::-webkit-details-marker { display:none; } .menu[open] summary, .menu summary:hover, .menubar > button:hover { background:var(--edge); }
.menu > div { position:absolute; z-index:7; top:100%; left:0; min-width:300px; padding:4px; background:var(--bg); border:1px solid var(--edge); border-radius:0 0 6px 6px; box-shadow:0 6px 18px rgba(0,0,0,.2); }
.menu > div button { display:block; width:100%; text-align:left; border:0; background:none; border-radius:4px; padding:6px 10px; font-weight:400; color:var(--fg); }
.menu > div button:hover { background:var(--edge); } .menu > div span { display:block; }
.menubar #saved { margin-left:auto; padding:0 10px; font-size:13px; }
.navsep { width:1px; align-self:stretch; background:var(--edge); margin:4px 8px; } .menubar > button:disabled { opacity:.4; cursor:default; }
.goal { background:var(--me); border:1px solid var(--edge); border-radius:6px; padding:10px 12px; max-width:900px; margin:8px 0; }
#how details { margin:6px 0; } #how summary { cursor:pointer; font-weight:600; }
dialog label { display:block; margin:8px 0 2px; font-size:13px; color:var(--mute); }
dialog input, dialog select { font:inherit; width:100%; box-sizing:border-box; padding:6px; background:var(--bg); color:var(--fg); border:1px solid var(--edge); border-radius:4px; }
</style></head><body>
<nav class="menubar">
  <details class="menu"><summary>File</summary><div>
    <button class="go" id="dl">Export GEDCOM file</button>
    <button id="m_close">Save and close this case</button>
    <button id="m_other">Save, then open or start another case</button>
    <span class="mute" style="padding:0 6px">Your work also saves by itself after every change.</span>
    <button id="reset">Undo all my edits</button>
  </div></details>
  <details class="menu"><summary>Edit</summary><div>
    <button id="m_undo">Undo&nbsp;&nbsp;<span class="mute" style="display:inline">Ctrl+Z</span></button>
    <button id="m_redo">Redo&nbsp;&nbsp;<span class="mute" style="display:inline">Ctrl+Y</span></button>
  </div></details>
  <details class="menu"><summary>Add</summary><div>
    <button id="obit">Paste a page I copied (Ancestry, Find a Grave, obituary)</button>
    <button id="add">A person</button>
  </div></details>
  <button id="m_how">How to</button>
  <button id="m_about">About</button>
  <span class="navsep"></span>
  <button id="n_back" title="Back to the person you were looking at before (Alt+Left)">&#9664; Back</button>
  <button id="n_fwd" title="Forward (Alt+Right)">Forward &#9654;</button>
  <button id="n_home" title="Back to the adoptee at the root of the tree (Alt+Home)">&#8962; Home</button>
  <span class="mute" id="saved"></span>
</nav>
<div id="prog" hidden><div id="prog_bar"><div id="prog_fill"></div></div><span id="prog_txt"></span></div>
<div id="live" class="goal" hidden><strong>The collecting window is open.</strong> Keep copying the pages it asks for. This tree updates by itself a few seconds after each page: the closest matches are placed, shared ancestors appear, and the DNA counts on each family grow.</div>
<h1>Birth family tree: <span id="case"></span></h1>
<details class="goal"><summary><strong>The goal, and a caution</strong> (click to read)</summary>Build a family tree file (a GEDCOM), upload it to Ancestry as a private tree, put the adoptee at its root, and link the adoptee's DNA test to them. Ancestry then shows which DNA matches fit which branch. Each time the matches confirm a branch, add to this tree, export again, and repeat one generation closer to the birth parents. <strong>This tree is a working guess, not proof.</strong> Confirm every person with records before relying on it.</details>
<p class="key"><span style="background:var(--me)">Adoptee</span><span style="background:var(--hit)">DNA match</span><span style="background:var(--unk)">Unknown, still to find</span><span style="background:var(--card)">Named from a tree or typed by you</span></p>
<dialog id="how" style="width:min(820px, 94vw)">
  <h2>How to</h2>
  <details open><summary>Get the tree into Ancestry</summary><ol>
    <li>File, then Export GEDCOM file. It saves to your Downloads folder.</li>
    <li>On Ancestry choose Trees, then Create &amp; Manage Trees, then Upload a GEDCOM file.</li>
    <li>Keep the new tree private and unsearchable.</li>
    <li>Open the adoptee's DNA settings and link the DNA test to the adoptee in the new tree.</li>
    <li>Give Ancestry a day or two. Common-ancestor hints and record hints then appear against this tree.</li></ol></details>
  <details><summary>Update the tree after new findings</summary><p>Make your changes here, export again, and upload the new file as a new tree. Link the DNA test to the adoptee in the newest tree and delete the older one. Ancestry cannot merge a second upload into an existing tree.</p></details>
  <details><summary>Change, join or remove a person</summary><p>Click the person. You can rename them, set sex and dates, add parents, a spouse or a child, or remove them. If two boxes are the same person, open one and pick the other under "This is the same person as".</p></details>
  <details><summary>Add a family from an obituary</summary><p>Click the person, use the links to find their obituary, and copy its text. Then Add, Family from an obituary or Find a Grave, choose who it is about, paste, and press Read it. Check each line before adding. A daughter listed under another surname is entered under the family name, married to that surname.</p></details>
  <details><summary>Add a family from Find a Grave</summary><p>Open the memorial page, press Ctrl+A then Ctrl+C, and paste it the same way. It reads the Family Members panel: parents, spouses, siblings and children with their years. It only lists relatives who have their own memorial, so living people are usually missing.</p></details>
  <details><summary>Search old newspapers</summary><p>Click the person and press Search old newspapers. It searches the Library of Congress collection, which ends in 1963 and is thin after the 1920s. Finding nothing does not mean there was no obituary.</p></details>
  <details><summary>Decide which branch to follow next</summary><p>At each couple, list their children and who each child married. Then check whether the spouse's family, including the spouse's mother's maiden name, shows up among the adoptee's DNA matches (menu choice 5 in the tool, "Does a surname show up in the matches?"). A spouse whose family has matches is the branch to follow. Repeat one generation down.</p></details>
  <details><summary>Collect more DNA matches</summary><p>Run "1 Capture a case" again with the same case name. It keeps what you have and adds new match pages and shared matches.</p></details>
  <details><summary>Keep case data private</summary><p>This page and the case folder hold real names. Keep them on your own computer and delete the folder when the case is closed.</p></details>
  <div class="bar"><button type="button" id="how_close">Close</button></div>
</dialog>
<dialog id="about" style="width:min(820px, 94vw)">
  <h2>About</h2>
  <p>Birth Family Tool, for Search Angels researchers. It never visits Ancestry or any site you copy from: you browse and copy, and it reads what you paste. Suggestions and problems: support@searchangels.org</p>
  <h3>Sources and credits</h3>
  <ul>
    <li><strong>Historic newspapers to 1963:</strong> Library of Congress, <em>Chronicling America: Historic American Newspapers</em>, loc.gov/collections/chronicling-america. Searched only when you press Search old newspapers.</li>
    <li><strong>1950 census:</strong> National Archives and Records Administration, 1950census.archives.gov. Public domain; the name index was read by software and corrected by volunteers. Looked up only when you press Look in the 1950 census.</li>
    <li><strong>Social Security applications and claims, 1936–2007:</strong> Wray, Anthony; U.S. Federal Security Agency, Social Security Board; and U.S. Social Security Administration. "Numerical Identification Files (NUMIDENT), 1936–2007." ICPSR [distributor], doi.org/10.3886/E207202V3. Creative Commons Attribution 4.0. Reduced to names, dates, places and parents' names; Social Security numbers removed. Please also cite Hollingsworth, Karbownik, Thomasson and Wray (2024), <em>American Economic Review</em> 114 (7): 2201–38.</li>
    <li><strong>Texas marriage and divorce indexes, 1966–2024:</strong> Texas Department of State Health Services, Vital Statistics Section.</li>
    <li><strong>Deceased veterans (BIRLS death file):</strong> U.S. Department of Veterans Affairs, obtained and published by Reclaim The Records, birls.org. Public domain.</li>
    <li><strong>Veterans cemetery burials:</strong> U.S. Department of Veterans Affairs, National Cemetery Administration, Gravesite Locator. CC0 public domain.</li>
    <li><strong>Read from text you paste, and not stored:</strong> Ancestry match pages, Find a Grave memorial pages, and obituaries.</li>
  </ul>
  <p class="mute">The Social Security, Texas and veterans record files are prepared for this tool but are not yet searchable from this page.</p>
  <div class="bar"><button type="button" id="about_close">Close</button></div>
</dialog>
<div class="scroll" id="out"></div>
<dialog id="dlg"><form method="dialog">
  <label>Name <input id="f_name"></label>
  <label>Sex <select id="f_sex"><option value="">Not known</option><option value="M">Male</option><option value="F">Female</option></select></label>
  <label>Born (year or date) <input id="f_b"></label>
  <label>Died (year or date) <input id="f_d"></label>
  <label>Note <input id="f_note"></label>
  <label>This is the same person as <select id="f_same"></select></label>
  <p class="mute" id="f_links"></p>
  <p><button type="button" id="f_census">Look in the 1950 census</button>
     <span class="mute">Free, from the National Archives. Shows the household they were living in. This sends this one person's name to 1950census.archives.gov.</span></p>
  <p><button type="button" id="f_news">Search old newspapers for this person</button>
     <span class="mute">Free, from the Library of Congress, up to 1963. This sends this one person's name and years to loc.gov.</span></p>
  <div class="bar">
    <button class="go" value="save">Save</button>
    <button value="parents">Add parents</button><button value="spouse">Add spouse</button><button value="child">Add child</button>
    <button value="remove">Remove</button><button value="cancel">Cancel</button>
  </div>
</form></dialog>
<div id="mini"></div>
<dialog id="exp" style="width:min(760px, 94vw)">
  <h2>GEDCOM file exported</h2>
  <p>The file was saved as <strong id="exp_name"></strong>. It holds the tree exactly as it looks now.</p>
  <h3>Upload it to Ancestry as a new tree</h3>
  <ol>
    <li>Sign in to Ancestry on the account that has access to the adoptee's DNA.</li>
    <li>Open the <strong>Trees</strong> menu and choose <strong>Create &amp; Manage Trees</strong>.</li>
    <li>Choose <strong>Upload a GEDCOM file</strong>. Do not open an existing tree: a GEDCOM always makes a new one.</li>
    <li>Press <strong>Choose file</strong> and pick the file from your Downloads folder.</li>
    <li>Give the tree a name, for example the case name and today's date.</li>
    <li>Untick <strong>Allow others to see my tree</strong> so it is private, then accept the submission agreement and press <strong>Upload</strong>.</li>
    <li>In the new tree's <strong>Tree Settings</strong>, under Privacy, also tick the box that keeps it out of search results.</li>
  </ol>
  <h3>Link the adoptee's DNA to it</h3>
  <ol>
    <li>Open the <strong>DNA</strong> menu and choose the adoptee's test, then <strong>Settings</strong>.</li>
    <li>Find <strong>Family tree linking</strong> and choose <strong>Link to tree</strong>.</li>
    <li>Pick the tree you just uploaded and select the adoptee as the person.</li>
    <li>Allow a day or two for Ancestry's common-ancestor hints to appear against the new tree.</li>
  </ol>
  <p class="mute">Ancestry changes its wording from time to time. If a step looks different, pick the closest match. Uploading again later makes another new tree: link the DNA to the newest one and delete the older one.</p>
  <div class="bar"><button type="button" id="exp_close">Close</button></div>
</dialog>
<dialog id="rec" style="width:min(900px, 95vw)">
  <p><strong id="rec_head"></strong></p>
  <div id="rec_rows"></div>
  <div class="bar"><button type="button" id="rec_close">Close</button></div>
</dialog>
<dialog id="busy" style="width:min(460px,90vw);box-shadow:0 8px 30px rgba(0,0,0,.35)">
  <p style="margin:0 0 6px"><strong id="busy_title">Working</strong></p>
  <p id="busy_txt" style="margin:0 0 8px"></p>
  <p class="mute" style="margin:0 0 8px">Please wait. The tree will update, and then you can carry on.</p>
  <button type="button" id="busy_stop" hidden>Stop the census import</button>
</dialog>
<dialog id="cen" style="width:min(860px, 95vw)">
  <p><strong id="cen_head"></strong></p>
  <div class="bar"><label style="margin:0">State, two letters, if you know it <input id="cen_state" size="3" maxlength="2" style="width:60px"></label><button type="button" id="cen_go">Search again</button></div>
  <div id="cen_rows"></div>
  <p class="mute">Source: 1950 Census, National Archives and Records Administration. The names were read from handwriting by software and are often misspelled; lines marked "checked" were corrected by volunteers. The census gives the head of the household first, then usually the wife, then the children. Open the page to see ages and relationships before trusting a household.</p>
  <div class="bar"><button type="button" id="cen_close">Close</button></div>
</dialog>
<dialog id="news" style="width:min(820px, 94vw)">
  <p><strong id="news_head"></strong></p>
  <div id="news_rows"></div>
  <p class="mute">Source: Chronicling America, Library of Congress. The text is machine-read from old print, so expect misspellings and misses. Open a page, find the obituary, copy its text, then use Add family from an obituary or Find a Grave.</p>
  <div class="bar"><button type="button" id="news_close">Close</button></div>
</dialog>
<dialog id="ob">
  <div id="ob_head">Paste a page <span class="mute">drag this bar to move the box</span><button type="button" id="ob_x" title="Close">&#10005;</button></div>
  <div id="ob_next"></div>
  <p><strong>Paste a page you copied.</strong> Open the page yourself, press Ctrl+A then Ctrl+C, and paste it below with Ctrl+V. The tool works out what it is and adds what it finds to this case: an Ancestry DNA match list, a match's page or Shared Matches list, a person's page from someone's Ancestry family tree, a Find a Grave memorial, or an obituary. It only reads what you paste.</p>
  <label>For an obituary or memorial page: who is it about? <select id="ob_who"></select></label>
  <label>Pasted text <textarea id="ob_text" rows="7" style="width:100%;box-sizing:border-box;font:inherit"></textarea></label>
  <div class="bar"><button class="go" id="ob_read" type="button">Read it</button><span class="mute">Check every line before adding. The tool guesses, and obituaries are worded many ways.</span></div>
  <div id="ob_rows"></div>
  <div class="bar"><button class="go" id="ob_add" type="button" disabled>Add the ticked people to the tree</button><button id="ob_close" type="button">Close</button></div>
  <details id="ob_rep"><summary>The page was not read, or was read wrongly? Send us a report</summary>
    <p class="mute">A report lets us teach the tool this kind of page and send out an update. It is saved as a file on your computer; nothing is sent until you email it.</p>
    <label>Where did you copy this from? (the website and which page) <input id="rep_src" style="width:100%;box-sizing:border-box"></label>
    <label>What went wrong? <input id="rep_what" style="width:100%;box-sizing:border-box"></label>
    <label><input type="checkbox" id="rep_mask" checked style="width:auto"> Blank out names before saving (recommended for Ancestry pages; the layout is what we need)</label>
    <div class="bar"><button type="button" id="rep_go">Save the report</button><span id="rep_done"></span></div>
  </details>
</dialog>
<script>
const DATA = /*DATA*/;
const KEY = "birthfamilytree:" + DATA.case + ":" + DATA.stamp;
// LIVE: the page was opened by the tool (menu choice 3), which can save the tree in the case folder, compare
// surnames with the DNA matches and look people up in the record files. Opened as a plain file it still works
// for viewing, editing and exporting, with edits kept in the browser.
const LIVE = typeof location !== "undefined" && /^http/.test(location.protocol), APIKEY = LIVE ? new URLSearchParams(location.search).get("k") : "";
const api = (path, body) => fetch(path, {method: body ? "POST" : "GET", headers: {"X-Key": APIKEY}, body: body ? JSON.stringify(body) : undefined}).then(r => r.json());
let S = JSON.parse(JSON.stringify(DATA)), cur = null, focus = null, trail = [], DNA = {};
if (!LIVE) try { const t = localStorage.getItem(KEY); if (t) S = JSON.parse(t); } catch (e) {}
const $ = id => document.getElementById(id);
const person = id => S.people.find(p => p.id === id);
const esc = s => String(s).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const famsOf = id => S.fams.filter(f => f.a === id || f.b === id);
const hasParents = id => S.fams.some(f => f.kids.includes(id));
function fresh(name) { const p = {id: Math.max(0, ...S.people.map(p => p.id)) + 1, name: name || "Unknown", sex: "", b: "", d: "", note: ""}; S.people.push(p); return p; }
// The surname to compare with DNA matches: a wife's or mother's likely maiden name when known, else the last word.
function surOf(p) {
  const m = /maiden name probably (\S+)/.exec(p.note || ""), w = (p.name || "").replace(/["().,]/g, "").trim().split(/\s+/).filter(x => !/^Unknown$/i.test(x));
  return (m ? m[1] : w.length > 1 || (w.length === 1 && /^Unknown/i.test(p.name)) ? w[w.length - 1] : w.length === 1 ? w[0] : "").toUpperCase();
}
function box(id) {
  const p = person(id); if (!p) return "";
  const isMatch = /^DNA match/.test(p.note), cls = p.root ? "me" : /^Unknown/.test(p.name) && !surOf(p) ? "unk" : isMatch ? "hit" : "";
  const dates = (p.b || p.d) ? (p.b || "?") + " - " + (p.d || "") : "", d = !p.root && !isMatch && DNA[surOf(p)];
  return `<div class="p ${cls}${id === focus ? " focus" : ""}" tabindex="0" data-id="${id}">${esc(p.name || "Unknown")}<small>${esc(dates)}</small><small>${esc(p.note)}</small>`
    + (d && d.n ? `<small class="dna">${esc(surOf(p))}: ${d.n} DNA match${d.n > 1 ? "es" : ""}, ${d.cm} cM</small>` : "") + "</div>";
}
// A person with their parents above them, their parents above those, and so on: the two sides always stay joined.
// In the whole-tree view each person on the line also has their brothers and sisters beside them, each with
// who they married and their own children and grandchildren below. That is where imported children and the
// DNA matches sit, so an import shows on the diagram.
// ponytail: children of an ancestor's OTHER marriages are not drawn here; they show when you click that ancestor.
function desc(id, depth, seen) {
  const mine = famsOf(id), kids = depth > 0 && !seen.includes(id) ? mine.flatMap(f => f.kids) : [];
  return `<div class="desc"><div class="couple">${box(id)}${mine.map(f => partnerOf(f, id)).filter(Boolean).map(box).join("")}</div>`
    + (kids.length ? `<div class="sibs below">${kids.map(k => desc(k, depth - 1, seen.concat(id))).join("")}</div>` : "") + "</div>";
}
function anc(id, depth, seen, left) {
  const f = !seen.includes(id) && S.fams.find(f => f.kids.includes(id));
  const up = !f ? "" : depth > 0 ? `<div class="par">${[f.a, f.b].filter(Boolean).map((x, i) => anc(x, depth - 1, seen.concat(id), i === 0)).join("")}</div>`
                                 : `<div class="more">more above: click to move up</div>`;
  const sibs = wholeTree && f ? f.kids.filter(k => k !== id).map(k => desc(k, 6, [id])).join("") : "";
  // a husband's brothers and sisters go to his left and a wife's to her right, so the couple stay side by side
  return `<div class="anc">${up}${sibs ? `<div class="sibs">${left ? sibs : ""}<div class="stem">${box(id)}</div>${left ? "" : sibs}</div>` : `<div class="stem">${box(id)}</div>`}</div>`;
}
const partnerOf = (f, id) => f.a === id ? f.b : f.a;
function connected() {   // everyone reachable from the adoptee through any family link
  const root = (S.people.find(p => p.root) || {}).id, seen = new Set(root ? [root] : []), todo = [...seen];
  while (todo.length) { const id = todo.pop();
    for (const f of S.fams) if (f.a === id || f.b === id || f.kids.includes(id))
      for (const x of [f.a, f.b, ...f.kids]) if (x && !seen.has(x)) { seen.add(x); todo.push(x); } }
  return seen;
}
// Moving around the tree, with Back, Forward and Home like a browser. Home is always the adoptee.
let navList = [], navAt = -1;
function go(id, viaButtons) {
  if (id === focus || !person(id)) return;
  if (!viaButtons) { if (navAt < 0) { navList = [focus]; navAt = 0; } navList = navList.slice(0, navAt + 1).concat(id).slice(-60); navAt = navList.length - 1; }
  trail = trail.filter(x => x !== id && x !== focus).concat(focus).slice(-6); focus = id; draw();
  if (typeof scrollTo === "function") scrollTo(0, 0);
}
function navStep(by) {
  let to = navAt + by;
  while (to >= 0 && to < navList.length && !person(navList[to])) to += by;   // skip anyone since removed
  if (to < 0 || to >= navList.length) return;
  navAt = to; go(navList[to], true);
}
function navButtons() {
  const root = S.people.find(p => p.root);
  $("n_back").disabled = navAt <= 0; $("n_fwd").disabled = navAt < 0 || navAt >= navList.length - 1; $("n_home").disabled = !root || root.id === focus;
}
const parentsOfId = id => { const f = S.fams.find(f => f.kids.includes(id)); return f ? [f.a, f.b] : [null, null]; };
const named = p => p && !/^Unknown$/.test(p.name);
const yrOf = s => +(String(s || "").match(/\d{4}/) || [0])[0];
// What to work on next, worked out from the tree and the DNA matches each time the page is drawn. The strongest
// lead comes first: a family that married into the adoptee's line and also shows up among the DNA matches.
function nextSteps() {
  const root = S.people.find(p => p.root), out = [], seen = new Set();
  if (!root) return out;
  const add = (score, id, text) => { if (!seen.has(id + text)) { seen.add(id + text); out.push({score, id, text}); } };
  const dna = p => (p && !p.root && !/^DNA match/.test(p.note) && DNA[surOf(p)]) || null;
  const line = [];                                   // the adoptee's direct ancestors, nearest first
  (function up(id, depth) { for (const x of parentsOfId(id)) if (x && depth < 6) { line.push(x); up(x, depth + 1); } })(root.id, 0);
  if (!yrOf(root.b) || /^Adoptee$/.test(root.name)) add(999, root.id, "Start here: click the adoptee's box and enter their name and birth year. The birth year is what lets the tool put shared ancestors at the right generation.");
  const [p1, p2] = parentsOfId(root.id).map(person);
  if (!named(p1) || !named(p2)) add(1, root.id, "One or both birth parents are still unknown. Work down each side from the closest DNA match until the two sides meet.");
  for (const id of line) {
    const p = person(id); if (!named(p)) continue;
    const d = dna(p), hasPar = hasParents(id);
    if (d && d.n && !hasPar) add(100 + d.cm, id, `${surOf(p)} has ${d.n} DNA match${d.n > 1 ? "es" : ""} (${d.cm} cM) but ${p.name}'s parents are not in the tree. Add them: look ${p.name} up in the records, or find an obituary.`);
    for (const f of famsOf(id)) {
      const others = f.kids.filter(k => !line.includes(k) && k !== root.id);
      if (!others.length) add(40, id, `No brothers or sisters are entered for the children of ${p.name}. List them and who each married: a spouse's family that has DNA matches shows which child to follow.`);
      for (const k of f.kids) for (const g of famsOf(k)) {
        const sp = person(partnerOf(g, k)), sd = dna(sp);
        if (sd && sd.n && !line.includes(k) && k !== root.id) add(60 + sd.cm, k, `${person(k).name} married into the ${surOf(sp)} family, which has ${sd.n} DNA match${sd.n > 1 ? "es" : ""} (${sd.cm} cM). That points to this child's branch: follow it down a generation.`);
      }
    }
  }
  if (LEADS && LEADS.unplaced && LEADS.unplaced.length) add(95, root.id, yrOf(root.b)
    ? `${LEADS.unplaced.length} shared ancestor${LEADS.unplaced.length > 1 ? "s were" : " was"} found but not put on the tree: ${LEADS.unplaced.slice(0, 4).map(c => c.name + " b. " + c.b + " (" + (c.why || "their matches are not assigned to a parent's side") + ")").join("; ")}. They are not added until they can be placed properly. Copy more of those matches' pages, or work out the line above them.`
    : "Enter the adoptee's birth year (click the adoptee). The tool needs it to place shared ancestors at the right generation.");
  const guessed = S.people.filter(p => /which line of that generation is not known yet/.test(p.note || "")).sort((a, b) => (b.b || "").localeCompare(a.b || ""));
  if (guessed.length) add(85, guessed[0].id, `${guessed.length} shared ancestor${guessed.length > 1 ? "s are" : " is"} on the right side and generation, but the exact line is not proven. Start with the most recent, ${guessed[0].name} b. ${guessed[0].b}: look them up in the records and list their children, to find which child leads down to the adoptee.`);
  if (LEADS) {
    const pg = LEADS.pages;
    if (pg.trees.length) add(80, root.id, `Copy more from Ancestry: ${pg.trees.length + pg.more_trees} DNA match${pg.trees.length + pg.more_trees > 1 ? "es have" : " has"} a public tree that has not been copied yet (${pg.trees.join(", ")}${pg.more_trees ? ", and more" : ""}). On Ancestry open DNA, Matches, click each name and copy the page.`);
    if (pg.shared.length) add(50, root.id, `Copy the Shared Matches list for: ${pg.shared.join(", ")}. It shows which matches belong to the same family group.`);
    if (pg.thin_sides.length) add(70, root.id, `No DNA match of 200 cM or more is assigned to ${pg.thin_sides.join(" or ")}. On Ancestry, check the match list's "By parent" view and the Unassigned group for closer matches on that side, and copy any you find.`);
  }
  for (const f of S.fams) {   // families known to be bigger than the tree shows: every missing child may lead to more DNA matches
    const miss = (+f.count || 0) - f.kids.length, head = person(f.a) || person(f.b);
    if (miss > 0 && head) add(75 + miss + (line.includes(head.id) ? 20 : 0), head.id, `${[f.a, f.b].map(person).filter(named).map(p => p.name).join(" and ") || "This couple"} had ${f.count} children and only ${f.kids.length} ${f.kids.length === 1 ? "is" : "are"} on the tree. Find the other ${miss}: their page lists who to search for and where.`);
  }
  const loose = S.people.filter(p => /^DNA match/.test(p.note) && !connected().has(p.id));
  if (loose.length) add(30, loose[0].id, `${loose.length} DNA match${loose.length > 1 ? "es are" : " is"} in the tree but not yet joined to the adoptee's line. Open one and use "This is the same person as" on the shared ancestors once you know them.`);
  if (line.some(id => named(person(id)))) add(5, root.id, S.exported ? "Changed the tree since your last export? Export the GEDCOM again and upload it to Ancestry as a new tree." : "When a generation is filled in, export the GEDCOM (File menu), upload it to Ancestry, and link the adoptee's DNA to it.");
  return out.sort((a, b) => b.score - a.score).slice(0, 7);
}
const lineOf = () => { const root = S.people.find(p => p.root), line = new Set(); if (root) (function up(id, d) { for (const x of parentsOfId(id)) if (x && d < 9) { line.add(x); up(x, d + 1); } })(root.id, 0); return line; };
// The pendulum. A couple on the adoptee's line has a child whose spouse's family also shows up in the DNA matches:
// that child and spouse are the next generation down. They take the place of the Unknown pair that stood there,
// keeping everything already hung above and below it, so the line moves one generation closer to the birth parents.
function followBranch(k) {
  const pf = S.fams.find(f => f.kids.includes(k)), line = lineOf(), root = S.people.find(p => p.root);
  const L = pf && pf.kids.find(x => x !== k && line.has(x)), old = person(L);
  if (!old) return $("saved").textContent = "This couple is not on the adoptee's line, so there is nothing to move here.";
  if (!/^Unknown$/.test(old.name)) return $("saved").textContent = `The adoptee's line already runs through ${old.name} here. Rename or remove that person first if this is wrong.`;
  const lf = famsOf(L).find(f => f.kids.some(x => line.has(x) || x === root.id)), kf = famsOf(k)[0];
  const sp = kf && partnerOf(kf, k), osp = lf && partnerOf(lf, L);
  const swap = (from, to) => S.fams.forEach(f => { if (f.a === from) f.a = to; if (f.b === from) f.b = to; f.kids = f.kids.map(x => x === from ? to : x); });
  pf.kids = pf.kids.filter(x => x !== L);                        // the Unknown child gives way to the named one
  swap(L, k);
  if (sp && osp && sp !== osp && /^Unknown$/.test(person(osp).name)) {   // and the Unknown spouse to the named spouse
    if (hasParents(sp)) S.fams.forEach(f => f.kids = f.kids.filter(x => x !== osp)); // the named spouse's own parents win
    swap(osp, sp);
    S.people = S.people.filter(p => p.id !== osp);
  }
  if (kf && lf && kf !== lf) { lf.kids = [...new Set(lf.kids.concat(kf.kids))]; S.fams = S.fams.filter(f => f !== kf); }
  S.people = S.people.filter(p => p.id !== L);
  person(k).note = (person(k).note ? person(k).note + ". " : "") + "On the adoptee's line: chosen because the spouse's family has DNA matches";
  focus = k; draw();
}
function collectMore() {   // open the copy-and-paste window for this case, then pick up what it found
  api("/api/capture", {name: LEADS ? LEADS.name : DATA.case}).then(() => { $("saved").textContent = "The collecting window is open. When you close it, press Refresh findings."; });
}
// A small map of the adoptee's ancestors, fixed in the lower left: where you are, and where the DNA matches point.
function miniMap() {
  const root = S.people.find(p => p.root); if (!root) return "";
  const W = 230, H = 150, G = 5, dots = [], lines = [];
  (function place(id, g, i) {
    if (!id || g > G) return;
    const p = person(id); if (!p) return;
    const x = (i + 0.5) / Math.pow(2, g) * W, y = H - 16 - g * ((H - 28) / G), d = !p.root && !/^DNA match/.test(p.note) && DNA[surOf(p)];
    const [a, b] = parentsOfId(id);
    [a, b].forEach((q, k) => { if (q && g < G) lines.push(`<line x1="${x}" y1="${y}" x2="${(2 * i + k + 0.5) / Math.pow(2, g + 1) * W}" y2="${H - 16 - (g + 1) * ((H - 28) / G)}"/>`); });
    dots.push(`<circle class="${p.root ? "me" : !named(p) ? "unk" : d && d.n ? "dna" : "known"}${id === focus ? " here" : ""}" cx="${x}" cy="${y}" r="${id === focus ? 7 : 5}" data-go="${id}"><title>${esc(p.name)}${d && d.n ? " - " + d.n + " DNA matches" : ""}</title></circle>`);
    place(a, g + 1, 2 * i); place(b, g + 1, 2 * i + 1);
  })(root.id, 0, 0);
  const here = person(focus), off = here && !dots.some(d => d.includes(`data-go="${focus}"`));
  return `<div class="minihead">Where you are${off ? ": " + esc(here.name) + " (a side branch)" : ""}</div><svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${lines.join("")}${dots.join("")}</svg>`
    + '<div class="minikey"><i class="dna"></i>DNA matches <i class="known"></i>named <i class="unk"></i>unknown</div>';
}
function draw() {
  S.fams = S.fams.filter(f => (f.a || f.b) && (f.kids.length || (f.a && f.b)));
  const root = S.people.find(p => p.root);
  if (!person(focus)) focus = person(S.focus) ? S.focus : (root || S.people[0] || {}).id;   // reopen where the researcher left off
  S.focus = focus;
  const me = person(focus);
  if (!me) { $("out").innerHTML = '<p class="mute">The tree is empty. Use Add, A person.</p>'; return save(); }
  const crumbs = trail.filter(person).map(id => `<a href="#" data-go="${id}">${esc(person(id).name)}</a>`).concat(`<strong>${esc(me.name)}</strong>`).join(" &rsaquo; ");
  const fams = famsOf(focus), pf = S.fams.find(f => f.kids.includes(focus)), onLine = lineOf();
  const sibs = pf ? pf.kids.filter(k => k !== focus) : [];
  const marriages = fams.map(f => { const o = partnerOf(f, focus), of = o && S.fams.find(g => g.kids.includes(o));
    return `<div class="fam"><div class="couple">${box(focus)}${o ? box(o) : '<div class="p unk">spouse not known</div>'}</div>`
      + (of ? `<p class="mute">Parents of ${esc(person(o).name)}:</p><div class="couple">${[of.a, of.b].filter(Boolean).map(box).join("")}</div>` : "")
      + (f.kids.length ? `<p class="mute">Children, each with who they married:</p><div class="kids">` + f.kids.map(k => {
          const mates = famsOf(k).map(g => partnerOf(g, k)).filter(Boolean), hit = mates.map(person).map(m => m && !m.root && DNA[surOf(m)]).find(d => d && d.n);
          const can = onLine.has(focus) && !onLine.has(k) && !person(k).root && f.kids.some(x => onLine.has(x) && /^Unknown$/.test(person(x).name));
          return `<div><div class="couple">${box(k)}${mates.map(box).join("")}</div>` + (can ? `<button data-follow="${k}"${hit ? ' class="go"' : ""}>${hit ? "DNA points here: " : ""}move the adoptee's line under ${esc(person(k).name)}</button>` : "") + "</div>"; }).join("") + "</div>" : '<p class="mute">No children entered.</p>')
      + `<p>${f.kids.length} on the tree. How many children did this couple have in all, if you know? <input type="number" min="0" max="30" style="width:70px" data-count="${S.fams.indexOf(f)}" value="${f.count || ""}"></p>` + missingKids(f) + "</div>"; }).join("");
  // surnames worth checking against the DNA matches around this person: their own, each spouse's, each mother's
  const near = [focus, ...(pf ? [pf.a, pf.b] : []), ...fams.flatMap(f => { const o = partnerOf(f, focus), of = o && S.fams.find(g => g.kids.includes(o)); return [o, ...(of ? [of.a, of.b] : [])]; })];
  // the adoptee's own surname is their adoptive one, so it says nothing about the birth family
  const surs = [...new Set(near.filter(Boolean).map(person).filter(p => !p.root).map(surOf).filter(Boolean))];
  const dna = !LIVE ? '<p class="mute">Open this page from the tool (choice 3 in "2 Work the case") to see which DNA matches fit each family and to look people up in the record files.</p>'
    : "<table>" + surs.map(s => { const d = DNA[s]; return `<tr><td><strong>${esc(s)}</strong></td><td>${!d ? "checking..." : d.n ? `${d.n} match${d.n > 1 ? "es" : ""}, ${d.cm} cM, mostly ${esc(d.side || "side not assigned")}` : "no matches"}</td><td class="mute">${d && d.n ? esc(d.who.join("; ")) : ""}</td></tr>`; }).join("") + "</table>";
  const linked = connected(), loose = S.people.filter(p => !linked.has(p.id));
  const steps = nextSteps();
  $("mini").innerHTML = miniMap();
  const stepsHtml = (steps.length ? `<div class="next"><strong>What to do next</strong><ol>${steps.map(s =>
      `<li>${esc(s.text)} ${s.id !== focus ? `<a href="#" data-go="${s.id}">Go to ${esc(person(s.id).name)}</a>` : ""}</li>`).join("")}</ol>`
      + (APP ? '<div class="bar"><button data-act="capture">Collect more DNA matches from Ancestry</button><button data-act="refresh">Refresh findings</button></div>' : "") + "</div>" : "");
  $("out").innerHTML = `<p class="trail">${root && root.id !== focus ? `<a href="#" data-go="${root.id}">Back to the adoptee</a> &nbsp;|&nbsp; ` : ""}${crumbs}</p>`
    + `<h2>${wholeTree && root ? "The whole tree" : esc(me.name) + " and the generations above"}</h2>`
    + `<div class="bar viewbar">${APP ? '<button class="go" data-act="next">What should I copy next?</button><span class="navsep"></span>' : ""}<button data-view="whole"${wholeTree ? ' class="go"' : ""}>Whole tree</button><button data-view="line"${wholeTree ? "" : ' class="go"'}>Just this person's line</button>`
    + `<span class="navsep"></span><button data-view="out" title="Zoom out">&minus;</button><button data-view="in" title="Zoom in">+</button><button data-view="fit">Fit all on screen</button>`
    + `<span class="mute">Drag the tree to move it. Scroll the mouse wheel to zoom. Click a person to go to them; click the highlighted person to change them.</span></div>`
    + `<div class="view" id="view"><div class="pedigree" id="ped">${anc(wholeTree && root ? root.id : focus, wholeTree ? 12 : 3, [])}</div></div>`
    + stepsHtml
    + `<div class="bar"><button data-act="edit">Change ${esc(me.name)}</button>${named(me) && !me.root ? `<button data-act="older">Look up ${esc(me.name)} in older records, step by step</button>` : ""}${LIVE ? `<button data-act="records">Look up ${esc(me.name)} in the record files</button>` : ""}</div>`
    + `<h2>Marriages and children</h2>${marriages || '<p class="mute">No marriage entered. Click the highlighted person to add a spouse or child, or use the Add menu.</p>'}`
    + (sibs.length ? `<h2>Brothers and sisters</h2><div class="kids">${sibs.map(k => `<div class="couple">${box(k)}${famsOf(k).map(g => partnerOf(g, k)).filter(Boolean).map(box).join("")}</div>`).join("")}</div>` : "")
    + `<h2>DNA matches around this family</h2><p class="mute">A spouse's family with DNA matches is the branch to follow. A family with none is probably married in.</p>${dna}`
    + (loose.length ? `<h2>Not connected to the adoptee yet</h2><p class="mute">Open one and use "This is the same person as" to join it to the tree.</p><div class="kids">${loose.map(p => box(p.id)).join("")}</div>` : "");
  $("case").textContent = DATA.case;
  if (S.people.length !== fitCount) { fitCount = S.people.length; fitNext = true; }   // people were added or removed: show all of it again
  placeView();
  remember();
  navButtons();
  guide(false);
  // entering or correcting the adoptee's birth year changes where shared ancestors belong: place them again
  const by = root ? yrOf(root.b) : 0;
  if (lastRootBirth === null) lastRootBirth = by;
  else if (by !== lastRootBirth) { lastRootBirth = by; if (LIVE) setTimeout(() => api("/api/leads").then(foldIn, () => {}), 0); }
  save();
  if (LIVE) { loadDNA(); if (LEADS) { clearTimeout(fillTimer); fillTimer = setTimeout(autoFill, 700); } }
}
// Undo and redo. Every change to the tree, whether made by the researcher or added by the tool, is a step that can be
// taken back and put back. Moving the view around is not a step. Up to 200 steps are kept while the tree is open.
let lastRootBirth = null;
let hist = [], hpos = -1, restoring = false, auto = false;   // auto: the tool, not the researcher, is changing the tree
const treeOnly = () => { const {focus: skip, ...rest} = S; return JSON.stringify(rest); };
function remember() {
  if (restoring) return;
  if (hpos >= 0 && !auto && !S.touched && treeOnly() !== hist[hpos]) S.touched = true;   // from here on the tool only adds; it never rebuilds over the researcher's work
  const now = treeOnly();
  if (hpos >= 0 && now === hist[hpos]) return;
  hist = hist.slice(0, hpos + 1).concat(now).slice(-200); hpos = hist.length - 1;
  undoState();
}
function undoState() {
  $("m_undo").disabled = hpos <= 0; $("m_redo").disabled = hpos >= hist.length - 1;
}
function step(by) {
  const to = hpos + by;
  if (to < 0 || to >= hist.length) return;
  hpos = to; restoring = true;
  S = JSON.parse(hist[to]); S.focus = focus;
  fillSeen = S.people.length;          // don't let the automatic lookup immediately redo what was just undone
  draw(); restoring = false; undoState();
  $("saved").textContent = by < 0 ? "Undone. Edit, Redo puts it back." : "Redone.";
}
$("n_back").onclick = () => navStep(-1);
$("n_fwd").onclick = () => navStep(1);
$("n_home").onclick = () => { const root = S.people.find(p => p.root); if (root) go(root.id); };
$("m_undo").onclick = () => step(-1);
$("m_redo").onclick = () => step(1);
document.addEventListener("keydown", e => {
  if (e.altKey && e.key === "ArrowLeft") { e.preventDefault(); return navStep(-1); }
  if (e.altKey && e.key === "ArrowRight") { e.preventDefault(); return navStep(1); }
  if (e.altKey && e.key === "Home") { e.preventDefault(); return $("n_home").onclick(); }
  if (!(e.ctrlKey || e.metaKey) || /^(INPUT|TEXTAREA|SELECT)$/.test((e.target && e.target.tagName) || "")) return;   // typing keeps its own undo
  const k = (e.key || "").toLowerCase();
  if (k === "z" && !e.shiftKey) { e.preventDefault(); step(-1); }
  if (k === "y" || (k === "z" && e.shiftKey)) { e.preventDefault(); step(1); }
});
(function () {   // drag the paste box by its title bar
  const box = $("ob"), bar = $("ob_head"); let sx = 0, sy = 0, ox = 0, oy = 0;
  bar.addEventListener("pointerdown", e => {
    if (e.target.id === "ob_x") return;
    const r = box.getBoundingClientRect(); sx = e.clientX; sy = e.clientY; ox = r.left; oy = r.top; bar.setPointerCapture(e.pointerId);
    bar.onpointermove = m => { box.style.left = Math.max(0, Math.min(innerWidth - 80, ox + m.clientX - sx)) + "px"; box.style.top = Math.max(40, Math.min(innerHeight - 40, oy + m.clientY - sy)) + "px"; box.style.right = "auto"; };
    bar.onpointerup = () => { bar.onpointermove = null; };
  });
  $("ob_x").onclick = () => box.close();
  // pasting is enough: the page is read as soon as it lands, with no button to press
  $("ob_text").addEventListener("paste", () => setTimeout(() => $("ob_read").onclick(), 60));
  box.addEventListener("change", e => { if (e.target.id === "st_site") { stepSite = +e.target.value; guide(false); } });
  box.addEventListener("click", e => {
    if (e.target.id === "st_end") { e.preventDefault(); stepFor = null; return guide(false); }
    if (e.target.id !== "st_copy") return;
    const name = person(stepFor).name.replace(/["().,]/g, "").trim(), ok = () => $("st_done").textContent = `"${name}" is on your clipboard.`;
    const old = () => { const t = document.createElement("textarea"); t.value = name; document.body.appendChild(t); t.select(); document.execCommand("copy"); t.remove(); ok(); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(name).then(ok, old); else old();
  });
  box.addEventListener("click", e => { const s = e.target.closest("[data-skipgrave]"); if (s) { e.preventDefault(); const p = person(+s.dataset.skipgrave); if (p) p.grave = "skipped"; draw(); } });
})();
let saveTimer = 0;
function save() {
  S.rules = RULES;
  if (LIVE) { clearTimeout(saveTimer); saveTimer = setTimeout(() => api("/api/tree", S).then(() => $("saved").textContent = "Saved in the case folder.",
                                                                                           () => $("saved").textContent = "NOT saved: the tool's window was closed. Export the GEDCOM now."), 400); return; }
  try { localStorage.setItem(KEY, JSON.stringify(S)); $("saved").textContent = "Built " + DATA.stamp + ". Edits are remembered in this browser only."; }
  catch (e) { $("saved").textContent = "This browser can't remember edits, so export the GEDCOM before closing."; }
}
function loadDNA() {   // ask the tool which surnames in the tree show up among the adoptee's matches; only new ones are asked
  const want = [...new Set(S.people.map(surOf).filter(s => s.length > 2 && !(s in DNA)))];
  if (!want.length) return;
  want.forEach(s => DNA[s] = null);
  api("/api/dna", want).then(r => { Object.assign(DNA, r); draw(); }, () => {});
}
function edit(id) {
  cur = person(id);
  $("f_name").value = cur.name; $("f_sex").value = cur.sex; $("f_b").value = cur.b; $("f_d").value = cur.d; $("f_note").value = cur.note;
  $("f_same").innerHTML = '<option value="">(nobody else)</option>' + S.people.filter(p => p !== cur)
    .sort((x, y) => x.name.localeCompare(y.name)).map(p => `<option value="${p.id}">${esc(p.name)}${p.b ? " b. " + esc(p.b) : ""} (#${p.id})</option>`).join("");
  // Ready-made searches the searcher opens in their own browser. These are plain links: the tool never visits them.
  $("f_links").innerHTML = !searchLinks(cur) ? "" : "Look for this person (opens in your browser): " + searchLinks(cur)
    + ". Open the right page, copy its text, then use Add family from an obituary or Find a Grave.";
  $("dlg").showModal();
}
function searchLinks(p) {
  const w = p.name.replace(/^Unknown\s*/, "").replace(/["().,]/g, "").trim().split(/\s+/).filter(Boolean), q = encodeURIComponent, yr = s => (String(s).match(/\d{4}/) || [""])[0];
  return !w.length ? "" : [
    ["Find a Grave", "https://www.findagrave.com/memorial/search?firstname=" + q(w.length > 1 ? w[0] : "") + "&lastname=" + q(w[w.length - 1])
                     + (yr(p.b) ? "&birthyear=" + yr(p.b) + "&birthyearfilter=5" : yr(p.d) ? "&deathyear=" + yr(p.d) + "&deathyearfilter=5" : "")],   // one year, give or take 5: years taken from other people's trees are often a little off, and two exact years found nothing
    ["obituary search", "https://www.google.com/search?q=" + q('"' + w.join(" ") + '" obituary')],
    ["Legacy.com", "https://www.google.com/search?q=" + q('site:legacy.com "' + w.join(" ") + '"')],
  ].map(([t, u]) => `<a href="${esc(u)}" target="_blank" rel="noopener">${t}</a>`).join(" · ");
}
// A couple known to have had more children than the tree shows: say how many are missing and where to look. An
// obituary or memorial for a parent, or for any child already known, usually names the rest and who they married.
function missingKids(f) {
  const n = +f.count || 0, miss = n - f.kids.length;
  if (miss <= 0) return "";
  const who = [f.a, f.b, ...f.kids].map(person).filter(p => named(p) && !p.root);
  return `<div class="next"><strong>${miss} of their ${n} children ${miss > 1 ? "are" : "is"} not on the tree yet.</strong> Look up the people already known. A parent's obituary lists the children and their spouses; a child's lists brothers and sisters. Copy what you find and use Add, Family from an obituary or Find a Grave. Each name added is then looked up in the records automatically.<ul>`
    + who.map(p => `<li>${esc(p.name)}${p.d ? " (died " + esc(p.d) + ")" : ""}: ${searchLinks(p)}</li>`).join("") + "</ul></div>";
}
// The tree diagram is a window onto a larger sheet: drag to move it, wheel to zoom, "Fit all" to see every
// generation at once. The position and zoom are kept while the tree is redrawn, so a change never jumps the view.
let fitCount = 0, wholeTree = true, zoom = 1, panX = 0, panY = 0, fitNext = true, dragged = false, drag = null;
function placeView() {
  const ped = $("ped"), view = $("view");
  if (!ped || !ped.style || !view) return;
  if (fitNext && view.clientWidth && ped.scrollWidth) {   // show everything, centred, never larger than life
    zoom = Math.max(.12, Math.min(1, (view.clientWidth - 24) / ped.scrollWidth, (view.clientHeight - 24) / ped.scrollHeight));
    panX = (view.clientWidth - ped.scrollWidth * zoom) / 2; panY = (view.clientHeight - ped.scrollHeight * zoom) / 2; fitNext = false;
  }
  ped.style.transform = `translate(${panX}px, ${panY}px) scale(${zoom})`;
}
function zoomBy(factor, cx, cy) {
  const view = $("view"), x = cx === undefined ? view.clientWidth / 2 : cx, y = cy === undefined ? view.clientHeight / 2 : cy, z = Math.max(.12, Math.min(2.5, zoom * factor));
  panX = x - (x - panX) * z / zoom; panY = y - (y - panY) * z / zoom; zoom = z; placeView();   // the point under the cursor stays put
}
$("out").addEventListener("pointerdown", e => { if (e.target.closest("#view") && e.button === 0) drag = {x: e.clientX, y: e.clientY, px: panX, py: panY}; dragged = false; });
$("out").addEventListener("pointermove", e => {
  if (!drag) return;
  if (!dragged && Math.abs(e.clientX - drag.x) + Math.abs(e.clientY - drag.y) < 5) return;   // a small wobble is still a click
  dragged = true; panX = drag.px + e.clientX - drag.x; panY = drag.py + e.clientY - drag.y; placeView();
});
["pointerup", "pointerleave"].forEach(t => $("out").addEventListener(t, () => { drag = null; }));
$("out").addEventListener("wheel", e => { const v = e.target.closest("#view"); if (!v) return; e.preventDefault(); const r = v.getBoundingClientRect(); zoomBy(e.deltaY < 0 ? 1.12 : 1 / 1.12, e.clientX - r.left, e.clientY - r.top); }, {passive: false});
function pick(e) {
  if (dragged) { dragged = false; return; }          // the tree was being moved, not clicked
  const v = e.target.closest("[data-view]");
  if (v) { const k = v.dataset.view;
    if (k === "in") return zoomBy(1.25); if (k === "out") return zoomBy(.8);
    if (k === "fit") { fitNext = true; return placeView(); }
    wholeTree = k === "whole"; fitNext = true; return draw(); }
  const a = e.target.closest("[data-go]"), act = e.target.closest("[data-act]"), b = e.target.closest(".p[data-id]"), up = e.target.closest(".more");
  if (a) { e.preventDefault(); return go(+a.dataset.go); }
  const fol = e.target.closest("[data-follow]");
  if (fol) return confirm("Move the adoptee's line under this person and their spouse? Do this when the spouse's family has DNA matches to the adoptee.") && followBranch(+fol.dataset.follow);
  if (act && act.dataset.act === "next") { stepFor = null; return guide(true); }
  if (act && act.dataset.act === "older") return startLookup(focus);
  if (act && act.dataset.act === "capture") return collectMore();
  if (act && act.dataset.act === "refresh") return api("/api/leads").then(foldIn, () => {});
  if (act) return act.dataset.act === "edit" ? edit(focus) : records(focus);
  if (up) return go(+up.parentNode.querySelector(".p[data-id]").dataset.id);
  if (b) return +b.dataset.id === focus ? edit(focus) : go(+b.dataset.id);
}
$("out").addEventListener("click", pick);
$("mini").addEventListener("click", pick);
$("out").addEventListener("change", e => { const i = e.target.dataset && e.target.dataset.count; if (i !== undefined && S.fams[+i]) { S.fams[+i].count = +e.target.value || 0; draw(); } });
$("out").addEventListener("keydown", e => { if (e.key === "Enter") pick(e); });

// ---------- record files: parents, brothers and sisters, children and marriages from public records ----------
let recRows = [];
function records(id) {
  const p = person(id), w = p.name.replace(/["().,]/g, "").trim().split(/\s+/), pf = S.fams.find(f => f.kids.includes(id));
  const par = pf ? [pf.a, pf.b].map(person).filter(x => x && !/^Unknown$/.test(x.name)) : [], words = s => s.replace(/["().,]/g, "").trim().split(/\s+/);
  const q = new URLSearchParams({first: w[0] || "", last: w[w.length - 1] || "", maiden: (/maiden name probably (\S+)/.exec(p.note || "") || ["", ""])[1], b: p.b || "", d: p.d || ""});
  if (par[0]) { q.set("ff", words(par[0].name)[0]); q.set("fl", words(par[0].name).pop()); }
  if (par[1]) q.set("ml", surOf(par[1]));
  $("rec_head").textContent = "Looking up " + p.name + "..."; $("rec_rows").innerHTML = ""; $("rec").showModal();
  api("/api/records?" + q).then(r => {
    recRows = []; const add = (label, items) => { recRows.push({label, items}); return recRows.length - 1; };
    const line = (text, i) => `<tr><td>${text}</td><td><button type="button" data-rec="${i}">Add</button></td></tr>`, yrs = x => (x.b || "?") + "-" + (x.d || "");
    const sec = (title, rows, note) => `<h3>${title}</h3>` + (rows.length ? `<table>${rows.join("")}</table>` : `<p class="mute">${note || "Nothing found."}</p>`);
    const SS = "From Social Security records", parentsOf = x => [x.father && {role: "parent", name: x.father, src: SS}, x.mother && {role: "parent", name: x.mother, src: SS}].filter(Boolean);
    $("rec_head").textContent = "Records for " + p.name;
    $("rec_rows").innerHTML = (r.have.length ? "" : '<p><strong>No record files were found on this computer.</strong> They are large and are installed separately.</p>')
      + sec("This person, with the parents they named", r.self.map(x => line(`<strong>${esc(x.name)}</strong> ${esc(yrs(x))}${x.also.length ? ", also " + esc(x.also.join(", ")) : ""}${x.pob ? ", born " + esc(x.pob) : ""}`
          + `<br>father ${esc(x.father || "not given")}, mother ${esc(x.mother || "not given")}`, add("self", parentsOf(x).concat({self: true, b: x.b, d: x.d})))), "Not in the Social Security file. It holds people who died by 2007, and about half its records name no parents.")
      + sec("Children: people who named this person as a parent", r.children.map(x => line(`<strong>${esc(x.name)}</strong> ${esc(yrs(x))}${x.also.length ? ", also " + esc(x.also.join(", ")) : ""}${x.pob ? ", born " + esc(x.pob) : ""}<br>other parent ${esc(x.other || "not given")}`,
          add("child", [{role: "child", name: x.name, b: x.b, d: x.d, src: SS}].concat(x.other ? [{role: "spouse", name: x.other, src: SS + ", named as a parent by a child"}] : [])))), "None found. Living children are not in this file.")
      + sec("Brothers and sisters: people who named the same parents", r.siblings.map(x => line(`<strong>${esc(x.name)}</strong> ${esc(yrs(x))}${x.also.length ? ", also " + esc(x.also.join(", ")) : ""}${x.pob ? ", born " + esc(x.pob) : ""}`,
          add("sibling", [{role: "sibling", name: x.name, b: x.b, d: x.d, src: SS}]))), "Needs both parents named in the tree, or none found.")
      + sec("Texas marriages and divorces", r.marriages.map(x => line(`${esc(x.kind)} ${esc(x.year)}, ${esc(x.county)} County: <strong>${esc(x.self)}</strong> (${esc(x.age || "?")}) and <strong>${esc(x.spouse)}</strong> (${esc(x.spouse_age || "?")})`,
          add("spouse", [{role: "spouse", name: x.spouse, src: "From the Texas " + x.kind + " index, " + x.year}]))), "None found. Texas only, 1966 to 2024.")
      + sec("Veterans cemetery burials tied to this person", r.burials.map(x => line(`<strong>${esc(x.name)}</strong> ${esc(yrs(x))}, ${esc(x.relationship)} of ${esc(x.veteran)}, ${esc(x.cemetery)} ${esc(x.state)}`,
          add("burial", /wife|husband/i.test(x.relationship) ? [{role: "spouse", name: x.name, b: x.b, d: x.d, src: "From veterans cemetery records"}]
                      : /son|daughter|child/i.test(x.relationship) ? [{role: "child", name: x.name, b: x.b, d: x.d, src: "From veterans cemetery records"}] : [{self: true, b: x.b, d: x.d}]))), "None found.")
      + sec("Deceased veterans with this name", r.veterans.map(x => line(`<strong>${esc(x.name)}</strong> ${esc(x.bfull || "?")} to ${esc(x.dfull || "?")}, ${esc(x.branch)}`, add("veteran", [{self: true, b: x.b, d: x.d}]))), "None found.")
      + '<p class="mute">A matching name is not proof it is the same person. Check the dates and places before adding. Sources are credited under About.</p>';
  }, () => { $("rec_head").textContent = "The record files could not be reached. Is the tool's window still open?"; });
}
$("rec_rows").addEventListener("click", e => {
  const b = e.target.closest("[data-rec]"); if (!b) return;
  const items = recRows[+b.dataset.rec].items, me = person(focus), own = items.find(i => i.self);
  if (own) { if (own.b && !me.b) me.b = own.b; if (own.d && !me.d) me.d = own.d; }
  addObit(focus, items.filter(i => !i.self)); b.textContent = "Added"; b.disabled = true; draw();
});
$("rec_close").onclick = () => $("rec").close();
$("dlg").addEventListener("close", () => {
  const act = $("dlg").returnValue, p = cur;
  if (!p || act === "cancel" || !act) return;
  if (act !== "remove") Object.assign(p, {name: $("f_name").value.trim() || "Unknown", sex: $("f_sex").value, b: $("f_b").value.trim(), d: $("f_d").value.trim(), note: $("f_note").value.trim()});
  const same = +$("f_same").value;
  if (act === "save" && same) {   // p and `same` are one person: keep `same`, move p's place in the tree onto it
    const q = person(same);
    for (const k of ["sex", "b", "d", "note"]) if (!q[k]) q[k] = p[k];
    if (/^Unknown/.test(q.name)) q.name = p.name;
    q.root = q.root || p.root;
    const hadParents = hasParents(q.id);
    S.fams.forEach(f => { if (f.a === p.id) f.a = q.id; if (f.b === p.id) f.b = q.id;
      f.kids = f.kids.map(k => k === p.id ? (hadParents ? null : q.id) : k).filter(Boolean); });
    S.people = S.people.filter(x => x !== p);
  }
  if (act === "parents" && !hasParents(p.id)) S.fams.push({a: fresh().id, b: fresh().id, kids: [p.id]});
  if (act === "spouse") S.fams.push({a: p.id, b: fresh().id, kids: []});
  if (act === "child") { let f = famsOf(p.id)[0]; if (!f) { f = {a: p.id, b: fresh().id, kids: []}; S.fams.push(f); } f.kids.push(fresh().id); }
  if (act === "remove" && confirmRemove(p)) {
    if (/^(Shared ancestor|Suggested by)/.test(p.note)) (S.dismissed = S.dismissed || []).push(p.name + "|" + String(p.b).slice(0, 4));   // don't bring it back on the next import
    S.fams.forEach(f => { if (f.a === p.id) f.a = null; if (f.b === p.id) f.b = null; f.kids = f.kids.filter(k => k !== p.id); });
    S.people = S.people.filter(x => x !== p);
  }
  draw();
});
function confirmRemove(p) { return !p.root || confirm("This is the adoptee. Remove anyway?"); }
$("add").onclick = () => { const p = fresh(); draw(); edit(p.id); };

// ---------- old newspapers: the one place this page contacts a website, and only when the button is pressed ----------
function newsQuery(p) {
  const w = p.name.replace(/^Unknown\s*/, "").replace(/["().]/g, "").trim().split(/\s+/).filter(Boolean), yr = s => +(String(s).match(/\d{4}/) || [0])[0];
  if (w.length < 2) return null;
  const d = yr(p.d), b = yr(p.b);
  // around the death year when known; otherwise adult life from the birth year; the collection ends in 1963
  const from = d ? d - 1 : b ? b + 18 : 1850, to = Math.min(d ? d + 1 : b ? b + 95 : 1963, 1963);
  if (from > 1963) return null;
  const phrase = w[0] + " " + w[w.length - 1];
  return {phrase, from, to, url: "https://www.loc.gov/collections/chronicling-america/?dl=page&ops=PHRASE&searchType=advanced&qs=" + encodeURIComponent(phrase)
          + "&start_date=" + from + "-01-01&end_date=" + to + "-12-31"};
}
function newsRows(json, phrase) {
  const last = phrase.split(" ").pop().toLowerCase();
  return (json.results || []).map(r => {
    const text = [].concat(r.description || []).join(" ").replace(/\s+/g, " "), at = text.toLowerCase().indexOf(last);
    return {date: r.date || "", paper: [].concat(r.partof_title || r.title || [""])[0], state: [].concat(r.location_state || [""])[0],
            url: r.url || r.id || "", snip: at < 0 ? "" : text.slice(Math.max(0, at - 140), at + 200)};
  });
}
$("f_news").onclick = async () => {
  const p = cur, q = p && newsQuery(p);
  $("dlg").close("cancel"); $("news_rows").innerHTML = ""; $("news").showModal();
  if (!q) { $("news_head").textContent = "Needs a first and last name, and a death or birth year before 1963, to search."; return; }
  $("news_head").textContent = `Searching for "${q.phrase}", ${q.from} to ${q.to}...`;
  try {
    const j = await (await fetch(q.url + "&fo=json&c=25")).json(), rows = newsRows(j, q.phrase);
    $("news_head").textContent = `"${q.phrase}", ${q.from} to ${q.to}: ${(j.pagination || {}).of || rows.length} newspaper pages` + (rows.length ? ", first " + rows.length + " shown" : "");
    $("news_rows").innerHTML = rows.map(r => `<p><a href="${esc(r.url)}" target="_blank" rel="noopener">${esc(r.date)} · ${esc(r.paper)}</a> <span class="mute">${esc(r.state)}</span><br>${esc(r.snip)}</p>`).join("")
      || `<p class="mute">Nothing found. Most newspapers are not in this collection, so this does not mean there was no obituary. <a href="${esc(q.url)}" target="_blank" rel="noopener">Open the search on loc.gov</a> to change the years or spelling.</p>`;
  } catch (e) {
    $("news_head").textContent = "The search could not be reached.";
    $("news_rows").innerHTML = `<p><a href="${esc(q.url)}" target="_blank" rel="noopener">Open the same search on loc.gov</a></p>`;
  }
};
$("news_close").onclick = () => $("news").close();
// ---------- the 1950 census: the household a person was living in ----------
let cenFor = null, cenList = [];
function censusSearch() {
  const p = cenFor; if (!p) return;
  const w = p.name.replace(/["().,]/g, "").trim().split(/\s+/).filter(x => !/^Unknown$/i.test(x));
  if (w.length < 2) { $("cen_head").textContent = "Needs a first and last name to search."; $("cen_rows").innerHTML = ""; return; }
  $("cen_head").textContent = `Looking for ${w[0]} ${w[w.length - 1]} in the 1950 census...`; $("cen_rows").innerHTML = "";
  api("/api/census?" + new URLSearchParams({first: w[0], last: w[w.length - 1], state: $("cen_state").value.trim()})).then(r => {
    cenList = r.results || [];
    $("cen_head").textContent = r.error || `${w[0]} ${w[w.length - 1]} in the 1950 census${r.state ? ", " + r.state : ""}: ${cenList.length} sheet${cenList.length === 1 ? "" : "s"} with that name`;
    $("cen_rows").innerHTML = cenList.map((s, i) => {
      // the household: the line with this person, then the lines straight after it
      const fam = s.lines.slice(s.at, s.at + 12).filter((n, k, all) => k === 0 || (s.checked && n.row != null && all[k - 1].row != null && n.row === all[k - 1].row + 1));
      s.fam = fam;
      return `<div class="fam"><strong>${esc(s.county)} County, ${esc(s.state)}</strong> <span class="mute">district ${esc(s.ed)}${s.checked ? ", checked by a volunteer" : ", machine-read only"}</span>`
        + `<div>${fam.map((n, k) => (k ? "" : "<strong>") + esc(n.name) + (k ? "" : "</strong>")).join(", ")}</div>`
        + `<div class="bar"><a href="${esc(s.view)}" target="_blank" rel="noopener">Open this page on the census site</a>`
        + (fam.length > 1 ? `<button type="button" data-cen="${i}">Review this household to add it</button>` : '<span class="mute">Only this line was read. Open the page to see the rest of the household.</span>') + "</div></div>"; }).join("")
      || '<p class="mute">Nothing found under that spelling. Try adding or removing the state, or open the 1950 census site and search by hand: the index misreads many names.</p>';
  }, () => { $("cen_head").textContent = "The lookup could not be made. Is the tool still running?"; });
}
$("f_census").onclick = () => { cenFor = cur; $("dlg").close("cancel"); $("cen").showModal(); censusSearch(); };
$("cen_go").onclick = censusSearch;
$("busy_stop").onclick = () => { censusStop = true; $("busy_txt").textContent = "Stopping after this lookup..."; };
$("busy").addEventListener("cancel", e => e.preventDefault());   // Esc must not dismiss it part-way through an import
$("cen_close").onclick = () => $("cen").close();
$("cen_rows").addEventListener("click", e => {
  const b = e.target.closest("[data-cen]"); if (!b || !cenFor) return;
  const s = cenList[+b.dataset.cen], fam = cenFor.name.trim().split(/\s+/).pop(), src = `From the 1950 census, ${s.county} County, ${s.state}, check`;
  // names may be written surname first, or with the surname left off for the rest of the family
  const tidy = raw => censusName(raw, fam);
  obSelf = null;
  obItems = s.fam.slice(1).map((n, k) => ({role: k === 0 ? "spouse" : "child", name: tidy(n.name), spouse: "", dead: false, src}));
  $("cen").close(); $("obit").onclick(); $("ob_who").value = String(cenFor.id); focus === cenFor.id || go(cenFor.id);
  showRows(`<p>From the 1950 census: the lines written after <strong>${esc(cenFor.name)}</strong> in ${esc(s.county)} County, ${esc(s.state)}. The index gives names only, so the roles below are a first guess (next line as the spouse, the rest as children). <strong>Open the census page and correct each role</strong> before adding.</p>`);
});
$("m_how").onclick = () => $("how").showModal(); $("how_close").onclick = () => $("how").close();
$("m_about").onclick = () => $("about").showModal(); $("about_close").onclick = () => $("about").close();
document.addEventListener("click", e => document.querySelectorAll(".menu[open]").forEach(m => { if (!m.contains(e.target) || e.target.tagName === "BUTTON") m.open = false; }));

// ---------- obituaries: the searcher pastes the text, this proposes relatives, the searcher confirms ----------
const STATES = "Alabama|Alaska|Arizona|Arkansas|California|Colorado|Connecticut|Delaware|Florida|Georgia|Hawaii|Idaho|Illinois|Indiana|Iowa|Kansas|Kentucky|Louisiana|Maine|Maryland|Massachusetts|Michigan|Minnesota|Mississippi|Missouri|Montana|Nebraska|Nevada|New Hampshire|New Jersey|New Mexico|New York|North Carolina|North Dakota|Ohio|Oklahoma|Oregon|Pennsylvania|Rhode Island|South Carolina|South Dakota|Tennessee|Texas|Utah|Vermont|Virginia|Washington|West Virginia|Wisconsin|Wyoming|[A-Z]{2}|[A-Z][a-z]{2,5}\\.";
const PLACE = new RegExp(",?\\s+(?:all |both )?of\\s+(?:the\\s+)?[A-Z][\\w.'-]*(?:\\s+[A-Z][\\w.'-]*){0,2}(?:,\\s*(?:" + STATES + "))?", "g");
const ROLE = {wife: "spouse", husband: "spouse", spouse: "spouse", son: "child", daughter: "child", children: "child",
              brother: "sibling", sister: "sibling", siblings: "sibling", parents: "parent", mother: "parent", father: "parent"};
function obitPerson(raw, kw, fam) {
  let s = raw.replace(/\b(?:Mr|Dr|Rev|Ms|Miss)\.?\s+/g, "").replace(/^[\s,;:]+|[\s,;:.]+$/g, ""), spouse = "", m;
  if (!/^[A-Z("]/.test(s) || /\d|Church|Funeral|Cemetery|Service|Hospital|Home\b|Chapel|Memorial/.test(s) || s.split(" ").length > 6) return null;
  if ((m = s.match(/^Mrs\.?\s+(.+)$/))) return {name: "Unknown", spouse: m[1]};   // a woman listed under her husband's name
  if ((m = s.match(/^(.*?)\s*\(([^)]+)\)\s*(.*)$/))) { s = (m[1] + " " + m[3]).trim(); spouse = m[2].trim(); }
  if (!s) return null;
  if (!/ /.test(s) && fam) s += " " + fam;                    // "sons, John and Bob": family surname understood
  const w = s.split(" "), last = w[w.length - 1], first = w.slice(0, -1).join(" ");
  if (/daughter|children/.test(kw) && fam && last.toLowerCase() !== fam.toLowerCase()) {
    // a daughter with another surname: taken as her married name, so she goes in under the family name
    return {name: first + " " + fam, spouse: (spouse || "Unknown") + " " + last};
  }
  if (spouse && !/ /.test(spouse) && /daughter|sister/.test(kw)) spouse += " " + last;
  return {name: s, spouse};
}
function parseObit(text, fam) {
  const t = text.replace(/\s+/g, " "), items = [];
  const cue = /survived by|survivors (?:include|are)|surviving (?:are|is)|(?:he|she) leaves|preceded in death by|predeceased by|preceded (?:him|her) in death/gi;
  const cuts = [...t.matchAll(cue)].map(m => ({at: m.index, end: m.index + m[0].length, dead: /preced|predeceas/i.test(m[0])}));
  const kwRe = /\b(sons?-in-law|daughters?-in-law|brothers?-in-law|sisters?-in-law|great[- ]grand\w+|grand\w+|step\w+|nieces?|nephews?|cousins?|friends|wife|husband|spouse|sons?|daughters?|children|brothers?|sisters?|siblings|parents|mother|father)\b/gi;
  cuts.forEach((c, i) => {
    // "Susan Ortiz and her husband Mike" becomes "Susan Ortiz (Mike)" before the relationship words are looked for
    const sec = t.slice(c.end, i + 1 < cuts.length ? cuts[i + 1].at : t.length)
      .replace(/\s+and\s+(?:his |her )?(?:wife|husband),?\s+([A-Z][\w'-]+)/g, " ($1)"), kws = [...sec.matchAll(kwRe)];
    kws.forEach((k, j) => {
      const kw = k[1].toLowerCase().replace(/s$/, "").replace(/^(children|parents|siblings)$/, "$1"), role = ROLE[k[1].toLowerCase()] || ROLE[kw];
      if (!role) return;
      let chunk = sec.slice(k.index + k[0].length, j + 1 < kws.length ? kws[j + 1].index : sec.length);
      chunk = chunk.split(/(?<!\b(?:Mrs|Mr|Dr|Jr|Sr|Rev|St|Ms|[A-Z]))\.\s+(?=[A-Z])/)[0];        // stop at the end of the sentence
      chunk = chunk.replace(PLACE, "");
      for (const raw of chunk.split(/;|,(?!\s*(?:Jr|Sr|II|III)\b)|\band\b|&/)) {
        const p = obitPerson(raw, k[1].toLowerCase(), role === "child" ? fam : role === "sibling" && /brother/.test(kw) ? "" : "");
        if (p) items.push({role, dead: c.dead, ...p});
      }
    });
  });
  let m = t.match(/\b(?:son|daughter|child) of (?:the late )?([A-Z][^.;]*?)(?=[.;]|,? (?:in|on|at|and was|was|he|she)\b)/)
       || t.match(/\bborn\b[^.;]*?\bto ((?:the late )?[A-Z][\w.'-]*(?: [A-Z(][\w.'()-]*)* and [A-Z][^.;,]*)/);
  if (m) {   // "son of Edwin and Mary (Quinn) Parker": in this spot the bracket is her maiden name
    let [a, b] = m[1].replace(/^the late /i, "").replace(PLACE, "").split(/\s+and\s+/), maiden = (b || "").match(/\(([^)]+)\)/);
    const sur = (b || a).replace(/\s*\([^)]*\)/, "").trim().split(" ").pop();
    if (a && !/ /.test(a.trim())) a = a.trim() + " " + sur;
    if (b) b = maiden ? b.replace(/\s*\([^)]*\).*/, "").trim() + " " + maiden[1] : b.trim();
    [a, b].filter(Boolean).forEach(n => items.push({role: "parent", dead: false, name: n.trim(), spouse: ""}));
  }
  m = t.match(/\bmarried (?:to )?(?:the former )?([A-Z][\w.'-]+(?: [A-Z][\w.'-]+){0,3})/);
  const sp = items.find(i => i.role === "spouse");
  if (m && !sp) items.push({role: "spouse", dead: false, name: m[1], spouse: ""});
  else if (m && !/ /.test(sp.name) && m[1].startsWith(sp.name + " ")) sp.name = m[1];   // "his wife, Ida" + "married Ida Lane"
  return items;
}
// Two names fit the same person when every word of the shorter is in the longer: "Fred Parker" and
// "Fred Calvin Parker", or "Marsh" and "Unknown Marsh". Only used among one person's own close family.
function same(a, b) {
  // a nickname in quotes or brackets is left out, and an initial stands for a name starting with it, so
  // 'Edwin B. "Ned" PARKER' and 'Edwin Blake Parker' are one person and are not put on the tree twice
  const words = s => s.toLowerCase().replace(/"[^"]*"|“[^”]*”|\([^)]*\)/g, " ").replace(/[.,]/g, "").split(/\s+/).filter(w => w && w !== "unknown");
  const [x, y] = [words(a), words(b)].sort((p, q) => p.length - q.length);
  const hit = (w, v) => w === v || (w.length === 1 && v[0] === w) || (v.length === 1 && w[0] === v);
  return x.length > 0 && x.every(w => y.some(v => hit(w, v)));
}
const fuller = (p, name) => { if (name.replace(/^Unknown\s*/, "").length > p.name.replace(/^Unknown\s*/, "").length) p.name = name; };
function addObit(pid, items) {   // hang the confirmed people on the person whose obituary it is
  const note = it => (it.src || "From obituary") + (it.dead ? ", died before them" : "") + (it.maiden ? ", maiden name probably " + it.maiden : "");
  const stamp = (p, it) => { p.note = p.note || note(it); if (it.b && !p.b) p.b = it.b; if (it.d && !p.d) p.d = it.d; };
  const named = (ids, name) => { const hit = ids.map(person).find(p => p && same(p.name, name)); if (hit) { fuller(hit, name); return hit; }
                                 const p = fresh(name); p.note = ""; return p; };
  const partner = f => f.a === pid ? f.b : f.a;
  const marry = (p, spouseName) => { if (!spouseName) return;
    const had = famsOf(p.id).map(f => person(f.a === p.id ? f.b : f.a)).find(o => o && same(o.name, spouseName));
    if (had) return fuller(had, spouseName);
    const s = fresh(spouseName); s.note = "From obituary, by marriage"; S.fams.push({a: p.id, b: s.id, kids: []}); };
  const parentFam = () => { let f = S.fams.find(f => f.kids.includes(pid));
    if (!f) { f = {a: fresh().id, b: fresh().id, kids: [pid]}; S.fams.push(f); } return f; };
  for (const it of items) {
    if (it.role === "spouse") {
      let f = famsOf(pid).find(f => { const o = person(partner(f)); return !o || /^Unknown/.test(o.name) || same(o.name, it.name); });
      if (!f) { f = {a: pid, b: null, kids: []}; S.fams.push(f); }
      let o = person(partner(f));
      if (!o) { o = fresh(); if (f.a === pid) f.b = o.id; else f.a = o.id; }
      if (/^Unknown/.test(o.name)) o.name = it.name; else fuller(o, it.name);
      stamp(o, it);
    } else if (it.role === "child") {
      let f = (it.partner && famsOf(pid).find(f => partnerOf(f, pid) === it.partner)) || famsOf(pid)[0];   // the right marriage, when known
      if (!f) { f = {a: pid, b: fresh().id, kids: []}; S.fams.push(f); }
      const kids = famsOf(pid).flatMap(f => f.kids);
      // "Mrs. Harold Smith" has no first name: recognise her on a second paste by her husband
      const wifeOf = it.spouse && kids.map(person).find(k => k && famsOf(k.id).some(f => { const o = person(f.a === k.id ? f.b : f.a); return o && same(o.name, it.spouse); }));
      const kid = /^Unknown$/.test(it.name) && wifeOf ? wifeOf : named(kids, it.name);
      if (!famsOf(pid).some(f => f.kids.includes(kid.id))) f.kids.push(kid.id);
      stamp(kid, it); marry(kid, it.spouse);
    } else if (it.role === "sibling") {
      const f = parentFam(), sib = named(f.kids, it.name);
      if (!f.kids.includes(sib.id)) f.kids.push(sib.id);
      stamp(sib, it); marry(sib, it.spouse);
    } else if (it.role === "parent") {
      const f = parentFam(), cur = [f.a, f.b].map(person);
      const hit = cur.find(p => p && !/^Unknown$/.test(p.name) && same(p.name, it.name)), slot = cur.find(p => p && /^Unknown$/.test(p.name));
      if (hit) { fuller(hit, it.name); stamp(hit, it); }
      else if (slot) { slot.name = it.name; slot.note = ""; stamp(slot, it); }
    }
  }
}
// A Find a Grave memorial page copied with Ctrl+A, Ctrl+C. Its "Family Members" panel lists each relative as a
// name line followed by a years line ("1897–1985"), under the headings Parents, Spouses, Siblings, Children.
// Built from one real copied page (Nov 2026); "Half Siblings" and "Children" headings are assumed, not yet seen.
function parseGrave(text, fam) {
  const L = text.split(/\r?\n/).map(s => s.trim()), at = L.indexOf("Family Members");
  if (at < 0) return null;
  const heads = {Parents: "parent", Parent: "parent", Spouse: "spouse", Spouses: "spouse", Siblings: "sibling", Sibling: "sibling",
                 "Half Siblings": "sibling", Children: "child", Child: "child"};
  const years = /^(\d{4}|unknown)?\s*[–—-]\s*(\d{4}|unknown)?(?:\s*\(m\.\s*(\d{4})\))?$/i, items = [];
  let role = null;
  for (let i = at + 1; i < L.length; i++) {
    if (!L[i]) continue;
    if (heads[L[i]]) { role = heads[L[i]]; continue; }
    if (/^(Flowers|See more|Sponsored|Inscription|Gravesite Details)\b/.test(L[i])) break;
    const j = L.findIndex((x, k) => k > i && x), y = j > 0 && L[j].match(years);
    if (!role || !y) continue;
    const w = L[i].replace(/\s+/g, " ").split(" "), k = fam ? w.map(x => x.toLowerCase()).lastIndexOf(fam.toLowerCase()) : -1;
    items.push({role, name: w.join(" "), spouse: "", dead: false, src: "From Find a Grave", b: /^\d/.test(y[1] || "") ? y[1] : "", d: /^\d/.test(y[2] || "") ? y[2] : "",
                // a wife or mother shown as "June Carter Cash": the word before the family surname is likely her maiden name
                maiden: (role === "spouse" || role === "parent") && k >= 2 && !/^"/.test(w[k - 1]) ? w[k - 1] : ""});
    i = j;
  }
  const after = h => { const n = L.indexOf(h); return n < 0 ? "" : (L.slice(n + 1).find(x => x) || "").replace(/\s*\(aged[^)]*\)/, ""); };
  const m = L.indexOf("Maintained by: Find a Grave") > 0 ? L.indexOf("Maintained by: Find a Grave") : L.findIndex(x => /^Maintained by:/.test(x));
  return {name: m > 0 ? L[m - 1] : "", b: after("Birth"), d: after("Death"), items};
}
// A person's profile page from someone's Ancestry family tree (learned from a real copied page): the name with
// "Birth ..." and "Death ..." lines under it, and near the end "Relationships" with Parents, Siblings, and
// "Spouse and children", each person a name line and a years line. Marriage records listed under Sources
// ("... in entry for A and B") say who a child married.
// ponytail: under "Spouse and children" a second wife cannot be told from a child by layout; anyone born less
// than 14 years after the person is taken as another spouse. Half-siblings in the timeline are not read.
const yearOf = s => (String(s || "").match(/\d{4}/) || [""])[0];
const nearName = (a, b) => { if (same(a, b)) return true;      // also "Frances Lola Dalton" and "Frances Lola Dalten"
  const x = a.trim().split(/\s+/), y = b.trim().split(/\s+/), lx = x.pop().toLowerCase(), ly = y.pop().toLowerCase();
  return x.length > 0 && y.length > 0 && lx.slice(0, 4) === ly.slice(0, 4) && Math.abs(lx.length - ly.length) <= 1 && same(x.join(" "), y.join(" ")); };
function parseProfile(text) {
  const L = text.split(/\r?\n/).map(s => s.trim()), at = L.indexOf("Relationships"), bi = L.findIndex(x => /^Birth \S/.test(x));
  if (at < 0 || !L.some(x => /^(Spouses? and children|Parents|Biological parents)$/.test(x))) return null;
  const years = /^(\d{4})?\s*[–—-]\s*(\d{4})?$/, heads = {Parents: "parent", "Biological parents": "parent", "Additional parent relationships": "skip", Siblings: "sibling", "Spouse and children": "family", "Spouses and children": "family"};
  const self = {name: bi > 0 ? L[bi - 1] : "", b: bi > 0 ? L[bi].replace(/^Birth\s*/, "").split("•")[0].trim() : "", d: "", profile: true, items: []};
  const di = L.findIndex(x => /^Death \S/.test(x)); if (di > 0 && di < at) self.d = L[di].replace(/^Death\s*/, "").split("•")[0].trim();
  const place = bi > 0 && L[bi].includes("•") ? L[bi].split("•")[1].trim() : ""; if (place) self.place = place;
  let role = null, first = true;
  for (let i = at + 1; i < L.length; i++) {
    if (!L[i]) continue;
    if (heads[L[i]]) { role = heads[L[i]]; first = true; continue; }
    if (/^(Support Center|Notifications|©)/.test(L[i])) break;
    const j = L.findIndex((x, k) => k > i && x), y = j > 0 && L[j].match(years);
    if (!role || role === "skip" || !y) continue;
    const b = y[1] || "", me = +yearOf(self.b), r = role !== "family" ? role : first || (me && b && +b < me + 14) ? "spouse" : "child";
    first = false;
    const name = L[i].replace(/\s+/g, " "); i = j;
    // trees often hold the same child twice ("Lena Dalten Stone 1913" and "Lena Dalten 1914"): keep the first
    if (self.items.some(it => it.role === r && same(it.name, name) && (!it.b || !b || Math.abs(it.b - b) <= 1))) continue;
    const it = {role: r, name, spouse: "", dead: false, b, d: y[2] || "", src: "From a family tree on Ancestry, check"};
    // a parent must be old enough: this is how a wrongly built tree shows itself
    if (r === "parent" && me && b && me - b < 13) { it.off = true; it.doubt = `born ${b}, only ${me - b} years before ${self.name.split(" ")[0]}: this tree has the wrong parent here`; }
    self.items.push(it);
  }
  for (const m of text.matchAll(/in entry for (.+?) and (.+?), "/g))      // a marriage record: one of the two is a child here
    for (const [a, b] of [[m[1], m[2]], [m[2], m[1]]]) { const kid = self.items.find(it => it.role === "child" && !it.spouse && nearName(it.name, a)); if (kid) kid.spouse = b.trim(); }
  return self.items.length ? self : null;
}
// The same people are often spelled a little differently in someone else's tree ("Adeline Francis Reeve" for
// "Addy Frances Reeve"). A relative already on the tree whose birth and death years both agree, and whose first
// or last name agrees, is the same person: the line takes the tree's spelling so they are not added twice.
function useTreeNames(who, items) {
  const near = [...famsOf(who.id).flatMap(f => [partnerOf(f, who.id), ...f.kids]), ...S.fams.filter(f => f.kids.includes(who.id)).flatMap(f => [f.a, f.b, ...f.kids])].map(person).filter(p => p && p.id !== who.id && !/^Unknown/.test(p.name));
  const w = s => s.toLowerCase().trim().split(/\s+/);
  for (const it of items) { if (!it.b || !it.d) continue;
    const hit = near.find(p => yearOf(p.b) && Math.abs(yearOf(p.b) - it.b) <= 1 && yearOf(p.d) === it.d && (w(p.name)[0] === w(it.name)[0] || w(p.name).pop() === w(it.name).pop()));
    if (hit && !same(hit.name, it.name)) it.name = hit.name; }
}
let obItems = [], obSelf = null;
$("obit").onclick = () => {
  $("ob_who").innerHTML = S.people.filter(p => !/^Unknown/.test(p.name)).sort((x, y) => x.name.localeCompare(y.name))
    .map(p => `<option value="${p.id}">${esc(p.name)}${p.b ? " b. " + esc(p.b) : ""}</option>`).join("");
  $("ob_rows").innerHTML = ""; $("ob").show();   // not modal: the tree behind stays usable and shows each change
};
$("ob_close").onclick = () => $("ob").close();
$("rep_go").onclick = () => {
  if (!LIVE) return $("rep_done").textContent = "Open the tree from the tool to save a report.";
  api("/api/report", {source: $("rep_src").value, what: $("rep_what").value, text: $("ob_text").value, mask: $("rep_mask").checked}).then(r => {
    $("rep_done").innerHTML = `Saved as <strong>${esc(r.path)}</strong>. Please attach that file to an email to <a href="mailto:support@searchangels.org?subject=${encodeURIComponent("Birth Family Tool report")}">support@searchangels.org</a>.`;
  }, () => $("rep_done").textContent = "Could not save the report: the tool has stopped.");
};
$("ob_read").onclick = () => {
  const text = $("ob_text").value;
  // An Ancestry page (match list, a match's page, Shared Matches) goes to the case itself, then the tree picks up what it adds.
  if (LIVE && ((/All matches/.test(text) && /By cluster/.test(text)) || /shared DNA/.test(text) || /^Shared Matches$/m.test(text))) {
    $("ob_rows").innerHTML = "<p>Reading the Ancestry page...</p>"; obItems = [];
    lock("Reading the page you pasted", "Adding what it holds to the case.");
    return api("/api/paste", {text}).then(r => {
      if (!r.ok) unlock();
      $("ob_rows").innerHTML = `<p class="${r.ok ? "" : "mute"}"><strong>${r.ok ? "Added to this case." : "Not read."}</strong> ${esc(r.msg)}</p>`
        + (r.ok ? '<p class="mute">Paste the next page, or close this box. The tree updates with anything new it found.</p>' : "<p>If this is the right page, open the report section below and send it to us.</p>");
      if (!r.ok) $("ob_rep").open = true;
      if (r.ok) { $("ob_text").value = ""; fillRounds = 0; return api("/api/leads").then(foldIn, unlock); }
    }, () => { unlock(); $("ob_rows").innerHTML = "<p>The tool has stopped, so the page could not be read.</p>"; });
  }
  let who = person(+$("ob_who").value);
  const prof = parseProfile($("ob_text").value);
  if (prof) {   // the page says who it is about: find that person on the tree by both years, then by name
    const by = S.people.filter(p => !p.root && yearOf(p.b) && yearOf(p.b) === yearOf(prof.b) && yearOf(p.d) === yearOf(prof.d) && nearName(p.name, prof.name));
    const hit = by.length === 1 ? by[0] : who && nearName(who.name, prof.name) ? who : S.people.filter(p => nearName(p.name, prof.name)).length === 1 ? S.people.find(p => nearName(p.name, prof.name)) : null;
    if (!hit) { obItems = []; obSelf = null; $("ob_rows").innerHTML = `<p><strong>Read as the Ancestry tree page of ${esc(prof.name)}</strong> (${esc(prof.b || "?")} to ${esc(prof.d || "?")}), but nobody on this tree fits that name and those years, so nothing is added. If they are on the tree under another spelling, choose them under "who is it about" and press Read it again.</p>`; return; }
    who = hit; $("ob_who").value = String(hit.id); obSelf = prof; obItems = prof.items; useTreeNames(who, obItems);
    // parents that contradict the ones already on the tree are shown but not ticked: this tree and that one disagree
    const has = parentsOfId(who.id).map(person).filter(q => q && named(q)), odd = obItems.filter(it => it.role === "parent" && !has.some(q => nearName(q.name, it.name)));
    if (has.length + odd.length > 2) odd.forEach(it => { it.off = true; it.doubt = it.doubt || `your tree already has ${has.map(q => q.name).join(" and ")} as a parent, and this tree names someone else`; });
    return showRows(`<p>Read as the Ancestry tree page of <strong>${esc(prof.name)}</strong> (${esc(prof.b || "?")}${prof.place ? ", " + esc(prof.place) : ""}; died ${esc(prof.d || "?")}), taken as <strong>${esc(who.name)}</strong> on this tree. Someone else built that tree, so check each line.</p>`);
  }
  const fam = who ? who.name.trim().split(" ").pop() : "";
  obSelf = parseGrave($("ob_text").value, fam);
  const house = !obSelf && who ? parseHousehold($("ob_text").value, who) : null;
  if (house) { obItems = house.items; if (house.selfBorn && !who.b) who.b = house.selfBorn;
    return showRows(`<p>Read as a ${house.year ? house.year + " " : ""}census household of ${house.rows} people. <strong>${esc(who.name)}</strong> is listed as ${esc(house.as)}, and the others are placed from that. Ages give the birth years. Check each line before adding.</p>`); }
  obItems = obSelf ? obSelf.items : parseObit($("ob_text").value, fam);
  const head = obSelf ? `<p>Read as a Find a Grave page for <strong>${esc(obSelf.name)}</strong> (born ${esc(obSelf.b || "?")}, died ${esc(obSelf.d || "?")}). ` +
    "Make sure that is the person chosen above. Find a Grave only lists relatives who have their own memorial, so living people are usually missing.</p>" : "";
  showRows(head);
};
function showRows(head) {   // the proposed people, one per line, each to be checked, corrected or unticked before adding
  $("ob_add").disabled = !obItems.length;
  $("ob_rows").innerHTML = head + (obItems.length ? `<p class="todo"><strong>Now check these ${obItems.length} lines.</strong> Untick anyone who is wrong, correct a name if needed, then press the green <strong>Add the ticked people to the tree</strong> button under the list. Nothing is added until you press it.</p>` : "") + (obItems.length ? "<table>" + obItems.map((it, i) =>
    `<tr><td><input type="checkbox" id="ob_on${i}"${it.off ? "" : " checked"}></td><td><select id="ob_role${i}">${["spouse", "child", "sibling", "parent"].map(r =>
      `<option${r === it.role ? " selected" : ""}>${r}</option>`).join("")}</select></td><td><input id="ob_name${i}" value="${esc(it.name)}"></td>` +
    `<td>married to <input id="ob_sp${i}" value="${esc(it.spouse || "")}"></td><td class="mute">${it.dead ? "died before" : ""} ${esc([it.b, it.d].some(Boolean) ? (it.b || "?") + "-" + (it.d || "") : "")}${it.maiden ? " maiden name probably " + esc(it.maiden) : ""}${it.doubt ? " <strong>Not ticked: " + esc(it.doubt) + ".</strong>" : ""}</td></tr>`).join("") + "</table>"
    : '<p><strong>Nothing was recognised in that page.</strong> Check that you copied the whole page (Ctrl+A, then Ctrl+C). If you did, open the report section below and send it to us, so the tool can learn this kind of page. You can still add people by hand.</p>');
  if (!obItems.length) $("ob_rep").open = true;
};
$("ob_add").onclick = () => {
  const picked = obItems.map((it, i) => $("ob_on" + i).checked ? {...it, role: $("ob_role" + i).value, name: $("ob_name" + i).value.trim(),
                                                                  spouse: $("ob_sp" + i).value.trim()} : null).filter(it => it && it.name);
  const who = person(+$("ob_who").value);
  if (obSelf && who) { if (!who.b) who.b = obSelf.b; if (!who.d) who.d = obSelf.d; who.grave = "done"; }   // the page's own dates, only where the tree has none
  // "one of seven children" tells us how big the family was, even when the page names only some of them
  const num = t => { const i = "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen".split(" ").indexOf(String(t).toLowerCase()); return i >= 0 ? i : +t || 0; };
  const sib = /\bone of (\w+) children\b/i.exec($("ob_text").value), kid = /\b(?:survived by|leaves|parents? of|father of|mother of)\s+(?:(?:his|her|their)\s+)?(\w+) children\b/i.exec($("ob_text").value);
  var obCounts = {sib: sib ? num(sib[1]) : 0, kid: kid ? num(kid[1]) : 0};
  const obCount = S.people.length;
  addObit(+$("ob_who").value, picked);
  if (who && obCounts.sib) { const pf = S.fams.find(f => f.kids.includes(who.id)) || (S.fams.push({a: fresh().id, b: fresh().id, kids: [who.id]}), S.fams[S.fams.length - 1]); pf.count = Math.max(+pf.count || 0, obCounts.sib); }
  if (who && obCounts.kid && famsOf(who.id)[0]) famsOf(who.id)[0].count = Math.max(+famsOf(who.id)[0].count || 0, obCounts.kid);
  const n0 = obCount; draw();
  const added = S.people.length - n0;
  $("ob_text").value = ""; obItems = []; obSelf = null; $("ob_add").disabled = true;
  $("ob_rows").innerHTML = `<p class="done"><strong>Done.</strong> ${added > 0 ? added + (added === 1 ? " person was" : " people were") + " added to the tree" : "Nobody new was added (they were already on the tree)"}${who ? " around " + esc(who.name) : ""}. Edit, Undo takes it back. The next page to copy is at the top of this box; or close the box to look at the tree.</p>`;
  guide(false);
};
$("reset").onclick = () => { if (confirm("Throw away every edit made on this page and go back to what the tool built?")) { S = JSON.parse(JSON.stringify(DATA)); draw(); } };
function gedcom() {
  const L = ["0 HEAD", "1 SOUR birthfamily", "1 GEDC", "2 VERS 5.5.1", "2 FORM LINEAGE-LINKED", "1 CHAR UTF-8"];
  const one = s => String(s).replace(/\s+/g, " ").trim();
  for (const p of [...S.people].sort((x, y) => (y.root ? 1 : 0) - (x.root ? 1 : 0))) {   // the adoptee first: Ancestry takes the first person as the tree's home person
    const w = one(p.name).split(" ");
    // a single word is taken as a surname: that is what the walk has when only a spouse's family is known
    const name = /^Unknown$/i.test(one(p.name)) || !w[0] ? "Unknown //" : (w.slice(0, -1).join(" ") + " /" + w[w.length - 1] + "/").trim();
    L.push(`0 @I${p.id}@ INDI`, "1 NAME " + name);
    if (p.sex) L.push("1 SEX " + p.sex);
    if (p.b) L.push("1 BIRT", "2 DATE " + one(p.b));
    if (p.d) L.push("1 DEAT", "2 DATE " + one(p.d));
    if (p.note) L.push("1 NOTE " + one(p.note));
    S.fams.forEach((f, i) => { if (f.kids.includes(p.id)) L.push(`1 FAMC @F${i + 1}@`); });
    S.fams.forEach((f, i) => { if (f.a === p.id || f.b === p.id) L.push(`1 FAMS @F${i + 1}@`); });
  }
  S.fams.forEach((f, i) => {
    let h = f.a, w = f.b;   // a known woman goes in the wife slot, a known man in the husband slot
    if ((h && person(h).sex === "F") || (w && person(w).sex === "M")) [h, w] = [w, h];
    L.push(`0 @F${i + 1}@ FAM`);
    if (h) L.push(`1 HUSB @I${h}@`);
    if (w) L.push(`1 WIFE @I${w}@`);
    f.kids.forEach(k => L.push(`1 CHIL @I${k}@`));
  });
  return L.concat("0 TRLR").join("\n") + "\n";
}
$("dl").onclick = () => {
  S.exported = true;
  if (LIVE) return api("/api/export", {name: DATA.case + "-birth-family-tree.ged", text: gedcom()}).then(r => { $("exp_name").textContent = r.path; $("exp").showModal(); save(); },
                                                                                                     () => $("saved").textContent = "Could not export: the tool has stopped.");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([gedcom()], {type: "text/plain"}));
  a.download = DATA.case + "-birth-family-tree.ged";
  a.click();
  $("exp_name").textContent = a.download; $("exp").showModal();
};
$("exp_close").onclick = () => $("exp").close();
// Changes are saved as they are made. Closing writes the tree once more, then hands control back to the tool's window.
const APP = LIVE && new URLSearchParams(location.search).get("app") === "1";   // opened from the home page: no command window
if (LIVE) setInterval(() => api("/api/ping").catch(() => {}), 20000);          // tells the tool the page is still open
function closeCase(next) {
  if (!LIVE) { $("saved").textContent = "Opened as a plain file: export the GEDCOM to keep your work, then close this tab."; return; }
  clearTimeout(saveTimer);
  if (APP) return api("/api/tree", S).then(() => location = "/home?k=" + APIKEY,
                                           () => $("saved").textContent = "Could not save: the tool has stopped. Export the GEDCOM to keep your work.");
  api("/api/tree", S).then(() => api("/api/close", {next})).then(() => {
    document.body.innerHTML = "<h1>Saved and closed</h1><p>This case is saved in its folder. You can close this tab.</p><p>"
      + (next ? "Go to the tool's window to open another case or start a new one." : "The tool's window is back at this case's menu.") + "</p>";
  }, () => { $("saved").textContent = "Could not reach the tool's window. Export the GEDCOM to keep your work."; });
}
$("m_close").onclick = () => closeCase("");
$("m_other").onclick = () => closeCase("again");
// Opened by the tool: start from the copy saved in the case folder, if there is one.
// Then fold in what the latest import found: anyone who appears in two or more matches' trees is added as a
// candidate shared ancestor (once; removing one keeps it removed), so every new batch of copied pages grows the tree.
//
// Nothing is left floating. A shared ancestor goes ON the tree, on the parent's side its matches are assigned to:
//   1. in the shared-couple box above that side's closest match, when that match's own tree contains them; else
//   2. at the generation their birth year implies (about 28 years a generation), joined down to that parent through
//      Unknown boxes. The side and the generation are evidence; WHICH line of that generation is not known yet,
//      and the note on the person says so. One whose side or generation can't be worked out is not added at all.
let LEADS = null;
function ensureParents(id) { let f = S.fams.find(f => f.kids.includes(id)); if (!f) { f = {a: fresh().id, b: fresh().id, kids: [id]}; S.fams.push(f); } return f; }
function depthAbove(from, target, d) {   // generations from one person up to an ancestor, or -1
  if (from === target) return d || 0;
  if ((d || 0) > 7) return -1;
  for (const x of parentsOfId(from)) if (x) { const r = depthAbove(x, target, (d || 0) + 1); if (r >= 0) return r; }
  return -1;
}
function freeSlot(id, depth) {           // an unclaimed Unknown box exactly `depth` generations above `id`, made if need be
  if (depth === 0) { const x = person(id); return x && /^Unknown$/.test(x.name) && !x.ph ? x : null; }
  const f = ensureParents(id);
  for (const par of [f.a, f.b]) if (par) { const r = freeSlot(par, depth - 1); if (r) return r; }
  return null;
}
function placeAncestor(c) {
  const root = S.people.find(x => x.root), side = /^Parent 2/.test(c.side) ? 1 : /^Parent 1/.test(c.side) ? 0 : -1;
  const yr = s => +(String(s || "").match(/\d{4}/) || [0])[0], rb = root && yr(root.b), cb = yr(c.b);
  if (!root || side < 0 || !rb || !cb) return false;
  const gens = Math.max(2, Math.min(6, Math.round((rb - cb) / 28))), start = [ensureParents(root.id).a, ensureParents(root.id).b][side];
  let slot = null, how = "";
  const sameTrees = p => p && p.inTrees && p.inTrees.length === c.matches.length && c.matches.every(m => p.inTrees.includes(m)) && Math.abs(yr(p.b) - cb) <= 15;
  const closest = LEADS.closest && LEADS.closest["Parent " + (side + 1)];
  const couple = S.people.find(x => x.ph === `?Parent ${side + 1} ancestor A` && /^Unknown$/.test(x.name));
  // ponytail: no couple is inferred. Appearing in the same matches' trees does not make two people husband and wife
  // (every shared ancestor of those matches does), so each one keeps an Unknown spouse until a record names them.
  if (couple && closest && c.matches.includes(closest) && Math.abs(depthAbove(root.id, couple.id) - gens) <= 1) {
    slot = couple; how = `placed as the ancestor shared with ${closest}, whose tree lists them`;
  } else if (start) {
    for (const g of [gens, gens + 1]) if (!slot) { slot = freeSlot(start, g - 1); how = `placed about ${g} generations above the adoptee by birth year; which line of that generation is not known yet`; }
  }
  if (!slot) { c.why = "no free place at that generation on that side"; return false; }
  Object.assign(slot, {name: c.name, b: c.b, inTrees: c.matches, note: c.note + ". " + how[0].toUpperCase() + how.slice(1) + "."});
  const f = S.fams.find(f => f.a === slot.id || f.b === slot.id), mate = f && person(f.a === slot.id ? f.b : f.a);
  if (mate && /^Unknown$/.test(mate.name) && !mate.ph) mate.ph = "spouse of a placed ancestor";   // only a spouse proven by the same trees may sit here
  return true;
}
// A couple Ancestry names as shared ancestors, with each one's exact relationship to the adoptee. They are seated
// together at exactly that generation on their side: in the box above that side's closest match when that match is
// one of theirs, otherwise in a free pair of boxes at that level. The side and generation are certain; which of the
// lines at that generation is still for the researcher to settle, and the note says so.
function placeCouple(h) {
  const root = S.people.find(x => x.root), side = /^Parent 2/.test(h.side) ? 1 : /^Parent 1/.test(h.side) ? 0 : -1;
  if (!root || side < 0 || !h.people.length) return;
  const linked = connected(), onTree = c => S.people.find(p => linked.has(p.id) && same(p.name, c.name) && (!p.b || !c.b || String(p.b).slice(0, 4) === c.b));
  if (h.people.every(onTree)) return h.people.forEach(c => { const p = onTree(c); if (!p.b) p.b = c.b; if (!p.d) p.d = c.d; if (!p.sex) p.sex = c.sex; });
  if (h.people.some(c => (S.dismissed || []).includes(c.name + "|" + c.b))) return;
  const gen = h.people[0].gen, start = [ensureParents(root.id).a, ensureParents(root.id).b][side];
  const A = S.people.find(x => x.ph === `?Parent ${side + 1} ancestor A` && /^Unknown$/.test(x.name));
  const closest = LEADS.closest && LEADS.closest["Parent " + (side + 1)];
  let slot = A && h.matches.includes(closest) && depthAbove(root.id, A.id) === gen ? A : start && gen >= 2 ? freeSlot(start, gen - 1) : null;
  if (!slot) return;
  const f = S.fams.find(f => f.a === slot.id || f.b === slot.id), mate = f && person(f.a === slot.id ? f.b : f.a);
  const fill = (p, c) => Object.assign(p, {name: c.name, b: c.b, d: c.d, sex: c.sex, ph: "", inTrees: h.matches,
    note: h.note + ". Ancestry gives this person as the adoptee's " + c.rel + (slot === A ? "." : "; which line of that generation is not known yet.")});
  fill(slot, h.people[0]);
  if (h.people[1] && mate && /^Unknown$/.test(mate.name)) fill(mate, h.people[1]);
}
// The closest match's branch must hang on the tree too. If the walk named that parent's line, the branch is tied to
// whoever sits at the right generation and still has no parents, preferring someone with the match's own surname.
function attachBranches() {
  const root = S.people.find(x => x.root); if (!root) return;
  for (const side of [0, 1]) {
    const A = S.people.find(x => x.ph === `?Parent ${side + 1} ancestor A`), f = A && S.fams.find(f => f.a === A.id || f.b === A.id);
    const depth = (DATA.couple_depth || {})["Parent " + (side + 1)], start = parentsOfId(root.id)[side];
    if (!f || !depth || !start || connected().has(A.id)) continue;
    let level = [start];
    for (let i = 1; i < depth; i++) level = level.flatMap(id => parentsOfId(id).filter(Boolean));
    const open = level.map(person).filter(p => p && !hasParents(p.id)), kin = new Set();
    (function down(g) { for (const k of g.kids) { const p = person(k); if (p && /^DNA match/.test(p.note)) kin.add(surOf(p)); famsOf(k).forEach(down); } })(f);
    const pick = open.find(p => kin.has(surOf(p))) || (open.length === 1 ? open[0] : null);
    if (pick) f.kids.push(pick.id);
  }
}
function foldIn(r) {
  LEADS = r; S.dismissed = S.dismissed || []; LEADS.unplaced = [];
  lock("Building the tree"); status(`Building the tree from ${r.pages.matches} DNA matches`, .35);
  attachBranches();
  (r.couples || []).forEach(placeCouple);   // Ancestry's own named ancestors first: their generation is exact
  const linked = connected();
  for (const c of r.ancestors) {
    const key = c.name + "|" + c.b, have = S.people.find(p => same(p.name, c.name) && (!p.b || !c.b || String(p.b).slice(0, 4) === c.b));
    if (have && linked.has(have.id)) { if (!have.b) have.b = c.b; if (!have.note) have.note = c.note; continue; }
    if (S.dismissed.includes(key)) continue;
    if (have) S.people = S.people.filter(p => p !== have);     // an earlier version left it floating: place it properly now
    if (!placeAncestor(c)) LEADS.unplaced.push(c);
  }
  (r.pedigrees || []).forEach(applyPedigree);
  auto = true; draw(); auto = false;
}
// A DNA match's own family tree, with its parent links, laid onto ours. Wherever someone in their tree is already
// on our tree (same name, and the same birth year when both have one), the people around them are filled in:
// UP, their parents and grandparents; and DOWN, the child of each couple and who that child married, generation
// by generation toward the match, who is then hung under their own parents. This is how the families between a
// shared ancestor and each DNA match appear. It comes from a public tree, so every person is marked to be checked.
function applyPedigree(g) {
  const P = g.people, src = `From the public tree of DNA match ${g.match} (${g.cm} cM), check`, at = {};
  const linked = connected(), yr4 = s => String(s || "").slice(0, 4);
  const fits = (p, c) => !/^Unknown/.test(p.name) && c.name.split(" ").length > 1 && same(p.name, c.name) && (!p.b || !c.b || yr4(p.b) === c.b);
  P.forEach((c, i) => { if (c.name) { const hit = S.people.find(p => linked.has(p.id) && !p.root && fits(p, c)); if (hit) at[i] = hit.id; } });
  if (!Object.keys(at).length) return;
  const item = (i, role, extra) => ({role, name: P[i].name, b: P[i].b, d: P[i].d, src, ...extra});
  const find = (ids, c) => { const p = ids.map(person).find(p => p && fits(p, c)); return p && p.id; };
  for (let i = 0; i < P.length; i++) {      // up: a person's parents come later in the list, so one pass reaches them all
    if (!at[i] || g.father[i] === undefined) continue;
    // Both parents are already on the tree as a couple, and this person is on it too but with no parents: this
    // tree says whose child they are. They take the place of that couple's "Unknown" child if there is one, so
    // the couple is joined to the named person and no second, empty line is left standing in for them.
    const pa = at[g.father[i]], ma = at[g.mother[i]], fam0 = pa && ma && S.fams.find(x => (x.a === pa && x.b === ma) || (x.a === ma && x.b === pa));
    if (fam0 && !hasParents(at[i]) && at[i] !== pa && at[i] !== ma) {
      const u = fam0.kids.findIndex(k => /^Unknown$/.test((person(k) || {}).name || ""));
      if (u >= 0) fam0.kids.splice(u, 1);
      fam0.kids.push(at[i]);
    }
    for (const [k, sex] of [[g.father[i], "M"], [g.mother[i], "F"]]) {
      if (k === undefined || at[k] || !P[k].name) continue;
      addObit(at[i], [item(k, "parent")]);
      at[k] = find(parentsOfId(at[i]).filter(Boolean), P[k]);
      if (at[k] && !person(at[k]).sex) person(at[k]).sex = sex;
    }
  }
  for (let i = P.length - 1; i >= 0; i--) { // down: a couple's child, then that child's spouse, toward the match
    const f = g.father[i], m = g.mother[i];
    if (f === undefined || m === undefined) continue;   // only a full preview says who the mother is
    for (const [have, need] of [[f, m], [m, f]]) {   // one parent on the tree, the other named: they are a couple
      if (at[have] && !at[need] && P[need].name) { addObit(at[have], [item(need, "spouse")]); at[need] = find(famsOf(at[have]).map(x => partnerOf(x, at[have])).filter(Boolean), P[need]); }
    }
    if (!at[f] || !at[m] || at[i]) continue;
    const child = i === 0 ? {role: "child", partner: at[m], name: g.match, b: P[0].b, d: P[0].d, src: `DNA match, ${g.cm} cM, placed from their own public tree`} : P[i].name ? item(i, "child", {partner: at[m]}) : null;
    if (!child) continue;                   // a living person, shown as Private: the line stops here
    addObit(at[f], [child]);
    const fam = S.fams.find(x => (x.a === at[f] && x.b === at[m]) || (x.a === at[m] && x.b === at[f]));
    at[i] = fam && find(fam.kids, {name: child.name, b: child.b});
  }
}
// While pages are being copied, the tree keeps up. Every few seconds the page asks whether anything new has arrived
// in the case. If so it takes in the new findings; and until the researcher has changed the tree by hand, it also
// re-seats the closest match on each side from the fuller match list, so the tree visibly builds as copying goes on.
let lastStamp = null, collecting = false;
function refreshFromCase() {
  const base = S.touched ? Promise.resolve(null) : api("/api/base");
  return base.then(b => { if (b && b.people) { b.focus = S.focus; b.dismissed = S.dismissed; S = b; fillRounds = 0; auto = true; draw(); auto = false; } return api("/api/leads"); })
             .then(r => { if (r && r.ancestors) foldIn(r); }, () => {});
}
function pulse() {
  api("/api/pulse").then(p => {
    if (p.busy !== collecting) { collecting = p.busy; $("live").hidden = !collecting; }
    if (lastStamp !== null && p.stamp !== lastStamp) refreshFromCase();
    lastStamp = p.stamp;
  }, () => {});
}
// Add what the records settle beyond doubt (see auto_fill in the tool): parents, children, and for each marriage the
// wife's maiden name, birth year and parents. It runs again whenever people are added, so a newly entered marriage
// is looked up straight away and each answer can lead to the next generation. Capped so it cannot run away.
let fillRounds = 0, fillSeen = 0, fillTimer = 0;
// The bar under the menus says what the tool is doing, so the researcher can see progress and knows when it is their turn.
function status(text, part) {
  $("prog").hidden = false; $("prog_txt").textContent = text; $("prog_fill").style.width = Math.round(part * 100) + "%";
  if ($("busy").open) $("busy_txt").textContent = text;
}
// While the tool is importing, the screen is held: nothing can be clicked until the import has finished and the
// tree has been redrawn with the new people. A stuck import lets go by itself after two minutes.
let lockTimer = 0;
function lock(title, text) {
  const b = $("busy"); $("busy_title").textContent = title; if (text) $("busy_txt").textContent = text;
  $("busy_stop").hidden = !/census/i.test(title);
  if (!b.open) b.showModal();
  clearTimeout(lockTimer); lockTimer = setTimeout(unlock, 120000);
}
function unlock() { clearTimeout(lockTimer); const b = $("busy"); if (b.open) b.close(); }
// The next thing to copy, most useful first: match pages with trees, then Shared Matches lists, then memorial pages
// for the ancestors on the tree. It is shown at the top of the paste box, which opens by itself when the tool is done.
// Any DNA match's family tree that holds this person, linked or unlinked, has a page for them listing their
// parents, spouse and every child. The person's note says which match's tree named them.
function inTreeOf(p) {
  const m = /shared with ([^(]+?) \(/.exec(p.note || "") || /DNA match ([^(]+?) \(/.exec(p.note || "");
  return `<strong>Best first, on Ancestry:</strong> ${m ? "open DNA, Matches, click <strong>" + esc(m[1].trim()) + "</strong>, and on the Trees tab look in their Linked Tree:" : "open any DNA match's family tree that has them (a linked or an unlinked tree):"} find <strong>${esc(p.name)}</strong>${p.b ? " (" + esc(yearOf(p.b)) + (p.d ? " to " + esc(yearOf(p.d)) : "") + "; their tree may spell the name differently)" : ""} and click them, then click <strong>View in tree</strong> (on an unlinked tree the button is <strong>Profile</strong>). Copy the page that opens; an advert at the top does not matter. It lists the parents, the spouse and every child, and the tool works out who it is about. Or find `;
}
function nextCopy() {
  const pg = LEADS && LEADS.pages, wait = " Wait until the page has finished loading, press Ctrl+A then Ctrl+C, click in the box below and press Ctrl+V.";
  // First the couples Ancestry names as shared ancestors: their children and who those children married are what
  // lead to more DNA matches, and the public records often have nothing on them.
  const couple = S.people.find(p => /^Common ancestor named/.test(p.note || "") && !p.grave && named(p));
  if (couple) { const f = famsOf(couple.id)[0], mate = f && person(partnerOf(f, couple.id)), n = f ? f.kids.length : 0;
    return `<strong>Next: the children of ${esc(couple.name)}${mate && named(mate) ? " and " + esc(mate.name) : ""}.</strong> ${n <= 1 ? "Only the line to the adoptee is on the tree so far, and the record files have no other children for them." : n + " children are on the tree."} `
      + inTreeOf(couple) + `${esc(couple.name)}'s memorial or obituary: ${searchLinks(couple)}. Copy the page, choose ${esc(couple.name)} under "who is it about", and paste it below. It usually lists the children and who they married. `
      + `<a href="#" data-skipgrave="${couple.id}">Skip this person</a>`; }
  if (pg && pg.trees.length) return `<strong>Next on Ancestry (${pg.trees.length + pg.more_trees} to go):</strong> open DNA, Matches and click <strong>${esc(pg.trees[0])}</strong>. On their page click the <strong>Trees</strong> tab (Ancestry may open on a different tab), so their family tree is showing.` + wait;
  if (pg && pg.shared.length) return `<strong>Next on Ancestry (${pg.shared.length} to go):</strong> open <strong>${esc(pg.shared[0])}</strong>, click the Shared Matches tab, and wait until you can see "Items per page" at the bottom.` + wait;
  const who = S.people.find(p => /^(Common ancestor named|Shared ancestor|Suggested by)/.test(p.note || "") && !p.grave && named(p));
  if (who) return `<strong>Next: ${esc(who.name)}</strong>${who.b ? " (born " + esc(who.b) + ")" : ""}. ` + inTreeOf(who) + `a memorial: ${searchLinks(who)}. Copy the memorial page, choose ${esc(who.name)} under "who is it about", and paste it below. `
    + `<a href="#" data-skipgrave="${who.id}">Skip this person</a>`;
  return "<strong>Nothing more is waiting to be copied.</strong> Paste any page you find, or export the GEDCOM from the File menu.";
}
// A census record page copied from Ancestry, FamilySearch or a library site: a household listed as name, age and
// relationship to the head, either on one line or on lines of their own. Everyone is placed from how the chosen
// person is listed (head, wife, son, daughter).
// ponytail: written from how census households are laid out in general, NOT yet from a real copied page. It only
// proposes lines for the researcher to tick; the first real page that reads wrongly should come in as a report.
function parseHousehold(text, who) {
  const REL = /^(head|self|wife|husband|son|daughter|mother|father|brother|sister|step-?son|step-?daughter|grandson|granddaughter|mother-in-law|father-in-law|son-in-law|daughter-in-law|nephew|niece|boarder|lodger|servant|roomer)$/i;
  const isName = s => /^[A-Z][A-Za-z.'-]*(?: [A-Z][A-Za-z.'-]*){0,4}$/.test(s) && !REL.test(s) && !/^(Name|Age|Relationship|Household|Members|Gender|Male|Female|Birth|Year|Census|View|Record|Role)$/i.test(s);
  const L = text.split(/\r?\n/).map(s => s.trim()), rows = [];
  L.forEach((line, i) => {
    const cells = line.split(/\t+|\s{2,}/).map(c => c.trim()).filter(Boolean);
    const rel = cells.find(c => REL.test(c));
    if (!rel) return;
    let name = cells.find(isName), age = cells.find(c => /^\d{1,3}$/.test(c) && +c <= 110);
    for (let k = i - 1; k >= Math.max(0, i - 4) && !name; k--) {   // the name is on a line above; an age may sit between them
      if (isName(L[k])) name = L[k];
      else if (age === undefined && /^\d{1,3}$/.test(L[k]) && +L[k] <= 110) age = L[k];
      else if (REL.test(L[k])) break;
    }
    if (age === undefined && /^\d{1,3}$/.test(L[i + 1] || "") && +L[i + 1] <= 110) age = L[i + 1];
    if (name && !rows.some(r => r.name === name)) rows.push({name, age: age === undefined ? null : +age, rel: rel.toLowerCase().replace("self", "head")});
  });
  if (rows.length < 2 || !rows.some(r => r.rel === "head")) return null;
  const y = /\b(1[789]\d0|19[0-5]0)\b[^\n]{0,40}census|census[^\n]{0,40}\b(1[789]\d0|19[0-5]0)\b/i.exec(text), year = y ? +(y[1] || y[2]) : 0;
  const first = who.name.trim().split(/\s+/)[0].toLowerCase(), me = rows.find(r => r.name.toLowerCase().split(" ")[0] === first);
  if (!me) return null;
  const born = r => year && r.age !== null ? String(year - r.age) : "", src = `From the ${year || ""} census (a copied household), check`.replace("  ", " ");
  const child = /son|daughter/.test(me.rel) && !/in-law|grand|step/.test(me.rel), roleOf = r =>
    me.rel === "head" ? (/^(wife|husband)$/.test(r.rel) ? "spouse" : /^(son|daughter)$/.test(r.rel) ? "child" : /^(mother|father)$/.test(r.rel) ? "parent" : /^(brother|sister)$/.test(r.rel) ? "sibling" : "")
    : /^(wife|husband)$/.test(me.rel) ? (r.rel === "head" ? "spouse" : /^(son|daughter)$/.test(r.rel) ? "child" : "")
    : child ? (/^(head|wife|husband)$/.test(r.rel) ? "parent" : /^(son|daughter)$/.test(r.rel) ? "sibling" : "") : "";
  const items = rows.filter(r => r !== me && roleOf(r)).map(r => ({role: roleOf(r), name: r.name, spouse: "", dead: false, b: born(r), d: "", src}));
  return items.length ? {items, year, rows: rows.length, as: me.rel === "head" ? "the head of the household" : "the head's " + me.rel, selfBorn: born(me)} : null;
}
// Looking someone up in older records, one step at a time. The tool picks the site, puts the name on the clipboard
// and opens the search; the researcher finds the record, copies the page and pastes it back here.
const LOOKUPS = [
  ["A census, 1790 to 1940, on FamilySearch (free account)", "https://www.familysearch.org/en/search/", "Type the name, add a birth year or place if you know it, and choose a census from the results. Open the record so the whole household is listed."],
  ["A census, 1790 to 1950, on Ancestry", "https://www.ancestry.com/search/categories/usfedcen/", "Search the name, open the census record, and make sure the household members are showing."],
  ["Find a Grave", "", "The search allows the birth year to be 5 years either way. Nothing listed? Remove the year on their search form, or add a state. Open the memorial that matches (copy the memorial itself, not the list of results) and scroll to the bottom once so the family section loads."]];
let stepFor = null, stepSite = 0;
function lookupSteps() {
  const p = person(stepFor); if (!p) return "";
  const [, url, hint] = LOOKUPS[stepSite], name = p.name.replace(/["().,]/g, "").trim(), link = stepSite === 2 ? (searchLinks(p).match(/href="([^"]+)"/) || [0, ""])[1].replace(/&amp;/g, "&") : url;
  return `<strong>Look up ${esc(p.name)}${p.b ? ", born " + esc(p.b) : ""} in older records</strong>`
    + `<ol style="margin:6px 0 0;padding-left:20px"><li>Where to look: <select id="st_site">${LOOKUPS.map((s, i) => `<option value="${i}"${i === stepSite ? " selected" : ""}>${esc(s[0])}</option>`).join("")}</select></li>`
    + `<li><button type="button" id="st_copy">Copy the name</button> <span id="st_done" class="mute">puts "${esc(name)}" on your clipboard, ready to paste into the site's search box.</span></li>`
    + `<li><a href="${esc(link)}" target="_blank" rel="noopener"><strong>Open the site</strong></a>, paste the name with Ctrl+V and search. ${esc(hint)}</li>`
    + `<li>When the right record is on screen and has finished loading, press <strong>Ctrl+A</strong> then <strong>Ctrl+C</strong>.</li>`
    + `<li>Come back here, click in the box below and press <strong>Ctrl+V</strong>. The tree updates and you can undo it with Ctrl+Z.</li></ol>`
    + `<p style="margin:6px 0 0"><a href="#" id="st_end">Finished with ${esc(p.name)}</a></p>`;
}
function startLookup(id) {
  stepFor = id; guide(false);
  if (!$("ob").open) $("obit").onclick();
  $("ob_who").value = String(id); $("ob_rows").innerHTML = "";
}
function guide(open) {
  $("ob_next").innerHTML = person(stepFor) ? lookupSteps() : LIVE && LEADS ? nextCopy() : "";
  if (open && APP && !$("ob").open) { $("ob_who").innerHTML = S.people.filter(p => !/^Unknown/.test(p.name)).sort((x, y) => x.name.localeCompare(y.name))
      .map(p => `<option value="${p.id}">${esc(p.name)}${p.b ? " b. " + esc(p.b) : ""}</option>`).join(""); $("ob").show(); }
}
// names in the census index may be written surname first, or with the surname left off for the rest of a family
function censusName(raw, fam) {
  const w = raw.replace(/[^A-Za-z.' -]/g, " ").trim().split(/\s+/).filter(Boolean), i = w.findIndex(x => x.toLowerCase() === fam.toLowerCase());
  if (i >= 0) w.splice(i, 1);
  return w.concat(fam).map(x => x[0].toUpperCase() + x.slice(1).toLowerCase()).join(" ");
}
// Once the tree is built and the record files have been applied, every couple who would have been keeping house in
// 1950 is looked up in the 1950 census, one at a time with the pause the Archives' site asks for. A household is
// added by itself only when it cannot be anyone else: a volunteer-checked sheet on which the husband's line is
// followed by a line with his wife's first name, and only one such sheet exists. The other lines of that household
// are added as their children. Anything less certain is left for the researcher's own "Look in the 1950 census".
let censusBusy = false, censusStop = false;
async function censusPass() {
  const wordsOf = s => s.toLowerCase().replace(/[^a-z ]/g, " ").split(/\s+/).filter(Boolean), full = p => named(p) && wordsOf(p.name).length >= 2;
  const home1950 = p => (!yrOf(p.b) || yrOf(p.b) <= 1932) && (!yrOf(p.d) || yrOf(p.d) >= 1950);
  const todo = S.fams.filter(f => !f.c50 && full(person(f.a)) && named(person(f.b)) && home1950(person(f.a)) && home1950(person(f.b))).slice(0, 40);
  if (!LIVE || censusBusy || !todo.length) return 0;
  censusBusy = true; censusStop = false;
  let added = 0;
  for (let i = 0; i < todo.length && !censusStop; i++) {
    const f = todo[i], a = person(f.a), b = person(f.b); if (!a || !b) continue;
    const [head, wife] = a.sex === "F" ? [b, a] : [a, b], hw = wordsOf(head.name), sp = wordsOf(wife.name)[0];
    // A husband's and wife's names alone fit strangers all over the country, so a household is only accepted when
    // it also holds a child this couple is already known to have. With no named child there is nothing to check
    // against, and the couple is left for the researcher to look up by hand.
    const known = f.kids.map(person).filter(k => k && !k.root && named(k)).map(k => wordsOf(k.name)[0]);
    if (!known.length) { f.c50 = "none"; continue; }
    const text = `Checking the 1950 census: couple ${i + 1} of ${todo.length} (${head.name})`;
    lock("Importing from the 1950 census", text); status(text + (added ? ". " + added + " people added so far" : ""), .9 + .1 * i / todo.length);
    try {
      const r = await api("/api/census?" + new URLSearchParams({first: hw[0], last: hw[hw.length - 1], state: ""}));
      const homes = (r.results || []).filter(s => s.checked).map(s => ({s, fam: s.lines.slice(s.at, s.at + 14).filter((n, k, all) => k === 0 || (n.row != null && all[k - 1].row != null && n.row === all[k - 1].row + 1))}))
                                     .filter(h => h.fam.slice(1).some(n => wordsOf(n.name).includes(sp)) && h.fam.slice(1).some(n => known.some(k => wordsOf(n.name).includes(k))));
      f.c50 = homes.length === 1 ? "found" : "none";
      if (homes.length === 1) {
        const {s, fam} = homes[0], sur = head.name.trim().split(/\s+/).pop(), src = `From the 1950 census, ${s.county} County, ${s.state} (a household with ${head.name} and ${wife.name}), check`;
        const kids = fam.slice(1).filter(n => !wordsOf(n.name).includes(sp)).map(n => ({role: "child", partner: wife.id, name: censusName(n.name, sur), src}));
        const before = S.people.length; addObit(head.id, kids); added += S.people.length - before;
        auto = true; draw(); auto = false;
      }
    } catch (e) { break; }
  }
  censusBusy = false;
  if (added) $("saved").textContent = `${added} people added from the 1950 census. Edit, Undo takes it back.`;
  return added;
}
// The tool has done what it can by itself: hand over to the researcher.
function settle() {
  censusPass().then(added => { if (added) { clearTimeout(fillTimer); fillTimer = setTimeout(autoFill, 50); return; }   // new people start another round of record lookups, which ends back here
    unlock(); status("Up to date. Your turn: see the box on the right for what to copy next.", 1); guide(true); });
}
function autoFill() {
  if (censusBusy) return;   // the census import is still adding people; it starts the record lookups again when it ends
  if (!LIVE || fillRounds >= 8 || S.people.length === fillSeen) { if (LIVE && LEADS) settle(); return; }
  fillSeen = S.people.length; fillRounds++;
  lock("Importing Social Security records");
  status(`Importing Social Security records for ${S.people.filter(named).length} people (pass ${fillRounds})`, Math.min(.95, .55 + .08 * fillRounds));
  api("/api/autofill", S).then(acts => {
    let changed = false;
    for (const a of acts) {
      const p = person(a.pid); if (!p) continue;
      const own = a.items.find(i => i.self), before = S.people.length + (p.b || "") + (p.note || "");
      if (own) {
        if (own.b && !p.b) p.b = own.b;
        if (own.d && !p.d) p.d = own.d;
        if (own.maiden && !/maiden name probably/.test(p.note || "")) p.note = (p.note ? p.note + ", " : "") + "maiden name probably " + own.maiden;
      }
      addObit(a.pid, a.items.filter(i => !i.self && !(S.dismissed || []).includes(i.name + "|" + (i.b || ""))));
      changed = changed || before !== S.people.length + (p.b || "") + (p.note || "");
    }
    if (changed) { auto = true; draw(); auto = false; $("saved").textContent = "Added from Social Security records. Edit, Undo takes it back."; }
    else settle();
  }, () => { unlock(); status("Up to date.", 1); });
}
const RULES = 2;   // the same number as TREE_RULES in the tool
// A tree saved by an older version is rebuilt from the collected pages when it is opened (the tool has already
// set the old one aside as a backup). The adoptee's own details are carried over.
async function fixOld(t) {
  lock("Updating this tree", "This tree was built by an older version of the tool. Building it again from the pages you collected.");
  const b = await api("/api/base"), was = (t.people || []).find(p => p.root), now = b.people.find(p => p.root);
  if (was && now) for (const k of ["name", "sex", "b", "d"]) if (was[k]) now[k] = was[k];
  S = b; fillRounds = 0; hist.length = 0; hpos = -1;
  $("saved").textContent = "This tree was rebuilt by the newer version. The old one is kept in the case folder as a backup.";
}
if (LIVE) api("/api/tree").then(t => (t && t.outdated ? fixOld(t).catch(() => { S = t; }) : Promise.resolve(t && t.people && t.people.length && (S = t))).then(() => { auto = true; draw(); auto = false; return api("/api/leads"); }), draw).then(r => r && r.ancestors && foldIn(r), () => {}); else draw();
if (LIVE) { pulse(); setInterval(pulse, 4000); }
</script></body></html>
"""


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


FROZEN = getattr(sys, "frozen", False)   # True in the installed app (packaged with Python inside)


def program_dir():
    return os.path.dirname(sys.executable if FROZEN else os.path.abspath(__file__))


def data_dir():
    """Where cases live. The installed app keeps them outside its program folder, so installing an update or
    repairing the app never touches them, and off OneDrive, so match names are not synced to the cloud."""
    if not FROZEN:
        return program_dir()
    d = os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "Birth Family Tool")
    os.makedirs(d, exist_ok=True)
    return d


def case_folder(name):
    """Cases live in one 'cases' folder, so researchers never type a path."""
    return os.path.join(data_dir(), "cases", "_".join(name.split()))


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


VERSION = "0.3.0"
REPO = "searchangels/searchangel-tools"            # where updates come from
REPO_FILE = "BirthFamilyTool/birthfamily.py"


def newer_version():
    """(version, new file text) when GitHub has a newer copy of this file, else None. Download only: nothing
    about any case is sent. A private repo needs a read-only token saved next to this file as
    update_token.txt; without one (or offline) the check quietly finds nothing."""
    import urllib.request
    here = os.path.abspath(__file__)
    tok = os.path.join(os.path.dirname(here), "update_token.txt")
    ver = lambda s: tuple(int(x) for x in s.split("."))
    if FROZEN:
        # The installed app updates by running the newer installer from the latest GitHub release.
        # ponytail: not yet exercised against a real published release; verify on the first one.
        try:
            import json
            with urllib.request.urlopen(f"https://api.github.com/repos/{REPO}/releases/latest", timeout=6) as r:
                rel = json.load(r)
            tag = rel.get("tag_name", "").lstrip("v")
            setup = next((a["browser_download_url"] for a in rel.get("assets", [])
                          if re.fullmatch(r"BirthFamilyTool-Setup-[\d.]+\.exe", a.get("name", ""))), None)
            return (tag, setup) if setup and ver(tag) > ver(VERSION) else None
        except Exception:
            return None
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
            return None
        compile(new, "birthfamily.py", "exec")               # never install a download that isn't valid code
        return m.group(1), new
    except Exception:
        return None


def install_update(new):
    import shutil
    if FROZEN:   # `new` is the installer's address: download it and let it replace this version; cases are not touched
        import subprocess, tempfile, urllib.request
        if not new.startswith(f"https://github.com/{REPO}/releases/download/"):
            return
        path = os.path.join(tempfile.mkdtemp(), "BirthFamilyTool-Setup.exe")
        with urllib.request.urlopen(new, timeout=120) as r, open(path, "wb") as f:
            shutil.copyfileobj(r, f)
        subprocess.Popen([path, "/SILENT", "/RELAUNCH=1"])
        return
    here = os.path.abspath(__file__)
    shutil.copyfile(here, here + ".previous")                # one step back, in case the new one misbehaves
    with open(here, "w", encoding="utf-8", newline="\n") as f:
        f.write(new)


def check_update():
    """Command-line launchers: offer a newer version and wait for the answer."""
    found = newer_version()
    if not found:
        return
    ans = ask(f"\nVersion {found[0]} of this tool is available (you have {VERSION}).\n"
              "Press Enter to update now, or type n to keep this version: ")
    if ans.lower().startswith("n"):
        return
    install_update(found[1])
    print(f"Updated to version {found[0]}. Close this window and start the tool again.")
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
    print(f"Birth Family Tool, version {VERSION}")
    while True:   # one case at a time: open an existing one or start a new one, work it, close it, pick again
        root = os.path.dirname(case_folder("x"))
        have = sorted(d for d in os.listdir(root) if os.path.exists(os.path.join(root, d, "matches.csv"))) if os.path.isdir(root) else []
        print("\nCases on this computer: " + (", ".join(have) if have else "none yet"))
        name = ask("Type a case name to open it, or a new name to start a new case (just Enter to quit): ")
        if not name:
            return
        folder = case_folder(name)
        if not os.path.exists(os.path.join(folder, "matches.csv")):
            if not open_files:
                print(f'Nothing has been captured for "{name}" yet.')
                return
            print(f'\n"{name}" is a new case. A small window is opening to collect its DNA matches from Ancestry.\n'
                  "Keep it next to your browser and follow what it says. Close it when you are finished or want a break.\n")
            cmd_guide(folder, 90.0)
            if not os.path.exists(os.path.join(folder, "matches.csv")):
                print("Nothing was captured, so the case was not started.")
                continue
            show_report(folder)
        if case_menu(folder, open_files) == "quit":
            return


def case_menu(folder, open_files=True):
    """The numbered menu for one open case. Returns 'again' to pick another case, or 'quit'."""
    while True:
        print(f"\n------------- case: {os.path.basename(folder)} -------------\n"
              "  1  Save the adoptee's birth year and place\n"
              "  2  Show the report\n"
              "  3  Build the tree: see it, edit it, download the file for Ancestry\n"
              "  4  Walk down the generations\n"
              "  5  Does a surname show up in the matches?\n"
              "  6  Does a cM amount fit a relationship?\n"
              "  7  Collect more DNA matches for this case\n"
              "  8  Close this case and open or start another\n"
              "  9  Quit\n"
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
                page = os.path.join(folder, "tree.html")
                if os.path.exists(page) and not ask("A tree page already exists. Press Enter to open it with your edits, "
                                                    "or type new to rebuild it from the latest findings: ").lower().startswith("n"):
                    pass
                else:
                    who = ask("Adoptee's name as it should appear in the tree: ")
                    page = cmd_skeleton(folder, who) if who else None
                    if page and os.path.exists(os.path.join(folder, "tree_saved.json")):
                        os.replace(os.path.join(folder, "tree_saved.json"), os.path.join(folder, "tree_saved.previous.json"))
                if page and open_files and show_tree(folder) == "again":
                    return "again"   # the researcher chose "Open or start another case" on the tree page
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
            elif ch == "7" and open_files:
                cmd_guide(folder, 90.0)
                show_report(folder)
            elif ch == "8":
                return "again"
            elif ch == "9":
                return "quit"
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
        early = g.on_copy("You and Ann Lee\nShared Matches\nBo Wu\n250 cM\n")           # copied before the list finished loading
        assert early[0] is False and "not finished loading" in early[1], early
        assert g.on_copy("You and Ann Lee\nShared Matches\nBo Wu\n250 cM\n20\nItems per page\n1\n")[0]
        # a test shared with the searcher: header names the test taker, not "You"
        assert g.on_copy("Pat Doe and Zed Nobody\nHow are Pat Doe and Zed Nobody related?\nSecret Name\n")[0] is False
        lay = open(os.path.join(d, "case", "layout_report.txt"), encoding="utf-8").read()
        assert "Secret" not in lay and "Zed" not in lay and "How are xxx xxx and" in lay, lay
        assert g.on_copy("Pat Doe and Bo Wu\nHow are Pat Doe and Bo Wu related?\n")[0] is False   # first tab refused
        ok, m = g.on_copy("Pat Doe and Bo Wu\nShared Matches\nAnn Lee\n300 cM\n20\nItems per page\n1\n")
        assert ok and "1 shared" in m, m
        g.finish_step()                                                                       # leave Dee Fox
        assert g.stage == "done", g.stage
        rep = open(os.path.join(d, "case", "clusters.txt"), encoding="utf-8").read()
        assert "Cluster 1" in rep and "Ann Lee" in rep and "Bo Wu" in rep, rep
        assert g.rows()[0]["tree"] == "Ann_Lee.txt"
        with open(os.path.join(cd, "tree_parent_1_s_side.txt"), "w", encoding="utf-8") as f:   # as the walk leaves it
            f.write("Ned_Parker = Edwin_Parker Mary_Quinn\n# born Ned_Parker 1920\nDot_Parker = Ned_Parker Ida_Lane\n"
                    "@parents|Parent 1's side|Dot_Parker|Al_Marsh\n")
        import json
        html = open(cmd_skeleton(cd, "Test Adoptee"), encoding="utf-8").read()   # case has Ann Lee, '2nd cousin', Parent 1's side
        td = json.loads(html.split("const DATA = ", 1)[1].split(";\nconst KEY", 1)[0])
        who = {p["id"]: p for p in td["people"]}
        kid = next(p for p in td["people"] if p.get("root"))
        fam = next(f for f in td["fams"] if kid["id"] in f["kids"])
        assert (who[fam["a"]]["name"], who[fam["b"]]["name"]) == ("Dot Parker", "Al Marsh"), fam     # walk's birth parents
        dot = next(f for f in td["fams"] if fam["a"] in f["kids"])
        assert who[dot["a"]]["name"] == "Ned Parker" and who[dot["a"]]["b"] == "1920" and who[dot["a"]]["sex"] == "", dot
        assert next(p for p in td["people"] if p["name"] == "Edwin Parker")["sex"] == "M"             # typed top couple only
        assert "300 cM" in next(p for p in td["people"] if p["name"] == "Ann Lee")["note"]
        # the page server: refuses a request without the key, saves a tree, reports DNA surnames, survives no record files
        import threading, urllib.request, urllib.error
        srv, address = tree_server(cd)
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        base, key = address.split("/?k=")
        call = lambda path, body=None, k=key: json.loads(urllib.request.urlopen(urllib.request.Request(
            base + path, data=json.dumps(body).encode() if body is not None else None, headers={"X-Key": k}), timeout=10).read())
        assert b"Birth family tree" in urllib.request.urlopen(address, timeout=10).read()
        try:
            call("/api/tree", k="wrong")
            raise AssertionError("a request with the wrong key was answered")
        except urllib.error.HTTPError as e:
            assert e.code == 403, e.code
        assert call("/api/tree") == {} and call("/api/tree", td)["ok"] and call("/api/tree")["people"] == td["people"]
        dn = call("/api/dna", ["Lee", "Zz"])
        assert dn["LEE"]["n"] == 1 and dn["LEE"]["cm"] == 300 and dn["ZZ"]["n"] == 0, dn
        assert record_lookup({"first": "Ann", "last": "Lee"}, files={}) == {"have": [], "self": [], "children": [], "siblings": [],
                                                                           "marriages": [], "veterans": [], "burials": []}
        srv.shutdown()
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
    if FROZEN:
        if sys.stdout is None or sys.stderr is None:   # the installed app has no console: keep messages in a log file
            sys.stdout = sys.stderr = open(os.path.join(data_dir(), "app.log"), "a", encoding="utf-8", buffering=1)
        a = a or ["app"]                               # double-clicking the app opens its window
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
    elif a[:1] == ["app"]:
        try:
            cmd_app()
        except Exception:   # there is no window to show an error in, so leave it where it can be found and sent
            import traceback
            with open(os.path.join(data_dir(), "app_error.txt"), "w", encoding="utf-8") as f:
                f.write(f"Birth Family Tool {VERSION}\n" + traceback.format_exc())
            raise
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
