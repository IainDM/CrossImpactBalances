@echo off
rem ============================================================
rem  Cross-Impact Balances desktop app launcher.
rem  Starts the app using all CPU cores, then the app opens
rem  your browser. Close this window (or click Quit in the
rem  page) to stop it.
rem ============================================================
set "JULIA_NUM_THREADS=auto"
"%~dp0bin\CrossImpactBalances.exe" %*
