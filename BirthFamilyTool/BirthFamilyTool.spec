# -*- mode: python ; coding: utf-8 -*-
# Birth Family Tool: PyInstaller build spec (one folder, no console window).
# Build from this folder:  python -m PyInstaller BirthFamilyTool.spec --noconfirm
import re

version = re.search(r'^VERSION = "([\d.]+)"', open("birthfamily.py", encoding="utf-8").read(), re.M).group(1)
parts = (version.split(".") + ["0", "0", "0"])[:4]
with open("version_info.txt", "w", encoding="utf-8") as f:   # what Windows shows under Properties > Details
    f.write(f"""VSVersionInfo(
  ffi=FixedFileInfo(filevers=({', '.join(parts)}), prodvers=({', '.join(parts)})),
  kids=[StringFileInfo([StringTable('040904B0', [
      StringStruct('CompanyName', 'Search Angels'),
      StringStruct('FileDescription', 'Birth Family Tool'),
      StringStruct('FileVersion', '{version}'),
      StringStruct('ProductName', 'Birth Family Tool'),
      StringStruct('ProductVersion', '{version}')])]),
    VarFileInfo([VarStruct('Translation', [1033, 1200])])])
""")

a = Analysis(["birthfamily.py"], pathex=[], binaries=[],
             datas=[("searchangels.ico", "."), ("searchangels_mark.png", "."), ("READ ME FIRST.html", "."), ("../SOURCES.md", ".")],
             hiddenimports=[], excludes=["numpy", "PIL", "matplotlib", "pandas", "scipy", "pytest"], noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="BirthFamilyTool", console=False,
          icon="searchangels.ico", version="version_info.txt")
coll = COLLECT(exe, a.binaries, a.datas, name="BirthFamilyTool")
