# Search Angels tools

Tools for Search Angels volunteers who help adoptees find their birth families.

## Birth Family Tool (test version)

### [Download the latest version](https://github.com/searchangels/searchangel-tools/releases/latest)

Click the zip under **Assets**. Free, and you don't need an account.

The tool helps a searcher work an adoptee's AncestryDNA matches: collect the matches, sort them into the two parents' sides and into family groups, find the most recent shared ancestors, and step down one generation at a time toward the birth parents.

Ancestry doesn't allow bots, AI agents or other automated tools to collect data from its site. This tool never touches Ancestry. You browse as usual, press Ctrl+A then Ctrl+C on each page, and the tool reads your clipboard. Case data stays on your computer.

### Install (Windows)

1. Install Python: open the Microsoft Store, search for **Python**, and install the newest one from the Python Software Foundation.
2. Download the zip from the link above and unzip it somewhere easy to find, such as Documents.
3. Open **READ ME FIRST.html** in the unzipped folder and follow it.

### Use

- Double-click **1 Capture a case** to collect a case. A small window tells you which page to copy next.
- Double-click **2 Work the case** for the report, the starting tree file, the generation walk and the checks.

Full instructions are in `BirthFamilyTool/READ ME FIRST.html`.

### Warnings you may see

- Your browser may say the zip "isn't commonly downloaded". Choose Keep.
- The first time you double-click a launcher, Windows may show "Windows protected your PC". Click **More info**, then **Run anyway**.

These appear because the tool is new and isn't signed by a software publisher.

### Updates

The tool checks here for a newer version each time it starts and asks before changing anything. Your cases aren't touched by an update.

### Suggestions and bugs

Email support@searchangels.org with "Birth Family Tool" in the subject line. Say which step you were on and what the window said.

### Privacy

Each case you capture holds real names of an adoptee's DNA matches. Keep case folders on your own computer and delete them when the case is finished. Never post one here or attach one to a bug report.

### For maintainers

- The whole tool is `BirthFamilyTool/birthfamily.py`. Run `python birthfamily.py selftest` before publishing a change.
- To publish an update, raise `VERSION` in that file and push to `main`. Everyone is offered it the next time they start the tool. Then publish a new release with a fresh zip so new downloads get it too.
- Never commit a `cases` folder or pasted pages. The `.gitignore` blocks them.
