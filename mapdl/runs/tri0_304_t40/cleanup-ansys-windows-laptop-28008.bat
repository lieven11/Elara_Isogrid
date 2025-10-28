@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 22864)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 28008)

del /F cleanup-ansys-windows-laptop-28008.bat
