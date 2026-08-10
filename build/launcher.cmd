@echo off
rem ============================================================
rem  Cross-Impact Balances desktop app launcher.
rem  Starts the app using all CPU cores, then the app opens
rem  your browser. Close this window (or click Quit in the
rem  page) to stop it.
rem ============================================================
rem  "auto,1" is all cores for the analyses plus one interactive thread, which
rem  is what keeps the page answering (and Quit working) while a long analysis
rem  is running. Recent Julia provides that thread by default; asking for it
rem  explicitly means the app does not depend on which Julia built it.
set "JULIA_NUM_THREADS=auto,1"
"%~dp0bin\CrossImpactBalances.exe" %*
