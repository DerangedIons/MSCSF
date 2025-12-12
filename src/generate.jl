# Example: julia --project=. src/generate.jl 300 1

BCL, ISO = parse.(Int, ARGS)

using MSCSF

sims = MSCSF.all_simulations()

sim = sims[BCL][ISO]

MSCSF.run(sim)
