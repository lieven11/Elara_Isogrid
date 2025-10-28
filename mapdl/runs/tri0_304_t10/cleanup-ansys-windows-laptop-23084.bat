@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 8584)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 23084)

del /F cleanup-ansys-windows-laptop-23084.bat
