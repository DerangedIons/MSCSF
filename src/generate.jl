# Example: julia --project=. src/generate.jl 300 1

BCL = parse(Int, ARGS[1])
ISO = parse(Int, ARGS[2])
N = parse(Int, get(ARGS, 3, "1"))

using CellSims

sims = CellSims.all_simulations()

sim = only(filter(sims) do sim
    sim.runner.ISO == ISO && sim.runner.BCL == BCL
end)

for _ in 1:N
    CellSims.run(sim)
end
