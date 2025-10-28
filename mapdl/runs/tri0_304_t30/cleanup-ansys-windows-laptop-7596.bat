@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 15832)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 7596)

del /F cleanup-ansys-windows-laptop-7596.bat
