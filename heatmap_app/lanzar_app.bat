@echo off
cd /d "%~dp0"
"..\R-4.5.3\bin\Rscript.exe" run_app.R
set EXITCODE=%ERRORLEVEL%

if %EXITCODE% NEQ 0 (
    echo.
    echo ==== The app exited unexpectedly ^(error code %EXITCODE%^). ====
    echo ==== Press any key to close this window and review the message above. ====
    pause >nul
)