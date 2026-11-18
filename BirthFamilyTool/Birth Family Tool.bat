@echo off
rem Starts the Birth Family Tool with no command window: it opens in your web browser.
cd /d "%~dp0"
pythonw -c "import tkinter" 2>nul
if errorlevel 1 goto nopython
start "" pythonw birthfamily.py app
exit /b
:nopython
echo.
echo Python is not installed on this computer yet.
echo Please follow Step 1 in "READ ME FIRST.html", then run this again.
echo.
pause
