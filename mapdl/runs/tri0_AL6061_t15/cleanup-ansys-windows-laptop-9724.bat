@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 16732)
if /i "%LOCALHOST%"=="windows-laptop" (taskkill /f /pid 9724)

del /F cleanup-ansys-windows-laptop-9724.bat
