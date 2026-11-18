"""Build the signed Windows installer for the Birth Family Tool.

  python build_installer.py            build, sign, verify
  python build_installer.py --unsigned build only (Windows Smart App Control will block the result)

Steps: self-test, package with PyInstaller, sign the program, build the installer with Inno Setup (which signs
its own setup program and uninstaller through the same tool), sign the installer, verify every signature.
Signing uses the Azure Trusted Signing account already set up on this computer (`az login` + the `sign` tool).
"""
import glob, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ACCOUNT = os.environ.get("TRUSTED_SIGNING_ACCOUNT", "herocompanionsign")
PROFILE = os.environ.get("TRUSTED_SIGNING_PROFILE", "hero-companion-public")
ENDPOINT = os.environ.get("TRUSTED_SIGNING_ENDPOINT", "https://eus.codesigning.azure.net/")
ISCC = os.path.join(os.environ.get("LOCALAPPDATA", ""), "Programs", "Inno Setup 6", "ISCC.exe")


def run(cmd, **kw):
    print(">", cmd if isinstance(cmd, str) else " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=HERE, **kw)


def sign_args(sign):
    return [sign, "code", "trusted-signing", "--trusted-signing-account", ACCOUNT, "--trusted-signing-certificate-profile", PROFILE,
            "--trusted-signing-endpoint", ENDPOINT, "--azure-credential-type", "azure-cli"]


def signed(path):
    out = subprocess.run(["powershell", "-NoProfile", "-Command",
                          f"$s = Get-AuthenticodeSignature -LiteralPath '{path}'; \"$($s.Status)|$($s.SignerCertificate.Subject)\""],
                         capture_output=True, text=True).stdout.strip()
    print(f"  {os.path.basename(path)}: {out.split(',')[0]}")
    return out.startswith("Valid|")


def main():
    sign = None if "--unsigned" in sys.argv else shutil.which("sign")
    if "--unsigned" not in sys.argv and not sign:
        sys.exit("The `sign` tool was not found. Install it with: dotnet tool install --global sign")
    version = re.search(r'^VERSION = "([\d.]+)"', open(os.path.join(HERE, "birthfamily.py"), encoding="utf-8").read(), re.M).group(1)
    run([sys.executable, "birthfamily.py", "selftest"], stdout=subprocess.DEVNULL)
    run([sys.executable, "-m", "PyInstaller", "BirthFamilyTool.spec", "--noconfirm", "--log-level", "ERROR"])
    exe = os.path.join(HERE, "dist", "BirthFamilyTool", "BirthFamilyTool.exe")
    setup = os.path.join(HERE, "dist", f"BirthFamilyTool-Setup-{version}.exe")
    iscc = [ISCC, "/Q", f"/DAppVersion={version}"]
    if sign:
        run(sign_args(sign) + [exe], stdout=subprocess.DEVNULL)
        # Inno Setup runs this for its setup program and uninstaller; $f is the file, $q a quote mark
        iscc += ["/DSigned=1", "/Sazure=" + " ".join(f"$q{a}$q" for a in sign_args(sign)) + " $q$f$q"]
    run(iscc + [os.path.join("installer", "BirthFamilyTool.iss")])
    if sign:
        run(sign_args(sign) + [setup], stdout=subprocess.DEVNULL)
        print("signatures:")
        if not all([signed(exe), signed(setup)]):
            sys.exit("A signature did not verify.")
    print(f"\nBuilt {setup} ({os.path.getsize(setup) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
