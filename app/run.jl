# Launch the Cross-Impact Balances desktop app from source.
#
#   julia --project=app -t auto app/run.jl
#
# Starts a localhost server and opens the default browser. Use `-t auto`
# so the exhaustive search and basin analysis use all CPU cores.
using CIBApp
exit(CIBApp.julia_main())
