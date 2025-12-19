# Example: julia --project=. src/generate.jl 300 1

BCL, ISO = parse.(Int, ARGS)

using CellSims

sims = CellSims.all_simulations()

sim = only(filter(sims) do sim
    sim.ISO == ISO && sim.BCL == BCL
end)

CellSims.run(sim)
