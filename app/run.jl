# Launch the Cross-Impact Balances desktop app from source.
#
#   julia --project=app -t auto,1 app/run.jl
#
# Starts a localhost server and opens the default browser. `auto,1` is all CPU
# cores for the exhaustive search and basin analysis, plus one interactive
# thread — that thread is what keeps the page answering, and Quit working,
# while a long analysis has every other core busy.
using CIBApp
exit(CIBApp.julia_main())
