@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 27664)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 27608)

del /F cleanup-ansys-windows-laptop-27608.bat
