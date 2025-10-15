@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 27596)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 19564)

del /F cleanup-ansys-windows-laptop-19564.bat
