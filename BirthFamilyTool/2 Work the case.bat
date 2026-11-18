@echo off
title Birth Family Tool
cd /d "%~dp0"
set PYTHONIOENCODING=utf-8
python -c "import tkinter" 2>nul
if errorlevel 1 goto nopython
python birthfamily.py menu
goto done
:nopython
echo.
echo Python is not installed on this computer yet.
echo Please follow Step 1 in "READ ME FIRST.html", then run this again.
:done
echo.
pause
