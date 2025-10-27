@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 4736)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 16752)

del /F cleanup-ansys-windows-laptop-16752.bat
