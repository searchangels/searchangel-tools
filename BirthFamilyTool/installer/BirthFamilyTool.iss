; Birth Family Tool: Inno Setup installer script.
; Build (after PyInstaller has produced dist\BirthFamilyTool):
;   ISCC.exe /DAppVersion=0.3.0 installer\BirthFamilyTool.iss
;
; Installs for the current user with no administrator prompt, adds an entry to
; "Add or remove programs" (Uninstall, and Repair by running this installer again),
; a Start menu entry and an optional desktop icon. Installing over an older version
; replaces it. Cases live in %LOCALAPPDATA%\Birth Family Tool and are never touched
; by an install or update; uninstalling asks whether to delete them.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#define AppName "Birth Family Tool"
#define AppPublisher "Search Angels"
#define AppURL "https://github.com/searchangels/searchangel-tools"

[Setup]
AppId={{5E0C7C1B-6B0B-4F4E-9D1E-2B6C2A7F3A51}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
AppPublisherURL={#AppURL}
AppSupportURL={#AppURL}
AppUpdatesURL={#AppURL}/releases
DefaultDirName={userpf}\BirthFamilyTool
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=..\dist
OutputBaseFilename=BirthFamilyTool-Setup-{#AppVersion}
SetupIconFile=..\searchangels.ico
UninstallDisplayIcon={app}\BirthFamilyTool.exe
UninstallDisplayName={#AppName}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
VersionInfoVersion={#AppVersion}
CloseApplications=force
RestartApplications=no
DirExistsWarning=no
; Windows Smart App Control blocks any unsigned program, including the setup program this installer unpacks
; into the temporary folder and the uninstaller it leaves behind. With a sign tool named "azure" given on the
; command line (see build_installer.py) Inno Setup signs those too.
#ifdef Signed
SignTool=azure
SignedUninstaller=yes
#endif
; the "Ready to install" page lists what is installed and the disk space it needs
DisableReadyPage=no
InfoBeforeFile=before_install.txt

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop icon"; GroupDescription: "Additional shortcuts:"

[InstallDelete]
; an update must leave exactly the new version's files: clear the old program files first (cases are elsewhere)
Type: filesandordirs; Name: "{app}\_internal"

[Files]
Source: "..\dist\BirthFamilyTool\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{userprograms}\{#AppName}"; Filename: "{app}\BirthFamilyTool.exe"
Name: "{userdesktop}\{#AppName}"; Filename: "{app}\BirthFamilyTool.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\BirthFamilyTool.exe"; Description: "Open {#AppName} now"; Flags: nowait postinstall skipifsilent
Filename: "{app}\BirthFamilyTool.exe"; Flags: nowait; Check: WantRelaunch

[UninstallDelete]
Type: filesandordirs; Name: "{app}"

[Code]
procedure CloseApp();
var
  R: Integer;
begin
  Exec(ExpandConstant('{sys}\taskkill.exe'), '/F /IM BirthFamilyTool.exe', '', SW_HIDE, ewWaitUntilTerminated, R);
end;

function PrepareToInstall(var NeedsRestart: Boolean): String;
begin
  CloseApp();   // nothing is lost: the tree is saved after every change
  Result := '';
end;

function InitializeUninstall(): Boolean;
begin
  CloseApp();
  Result := True;
end;

function WantRelaunch(): Boolean;
begin
  Result := ExpandConstant('{param:RELAUNCH|0}') = '1';
end;

// Uninstall keeps the cases unless the researcher says to delete them.
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  Data: String;
begin
  if CurUninstallStep = usPostUninstall then
  begin
    Data := ExpandConstant('{localappdata}\Birth Family Tool');
    if DirExists(Data) and not UninstallSilent() then
      if MsgBox('Birth Family Tool has been removed.' + #13#10 + #13#10 +
                'Do you also want to delete the cases and trees you created with it?' + #13#10 + Data + #13#10 + #13#10 +
                'Choose No to keep them. They will be there if you install the tool again.',
                mbConfirmation, MB_YESNO or MB_DEFBUTTON2) = IDYES then
        DelTree(Data, True, True, True);
  end;
end;
