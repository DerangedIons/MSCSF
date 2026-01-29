# Usage:
#   Model3D:   julia --project=. src/generate.jl BCL ISO [N]
#   CaClamp3D: julia --project=. src/generate.jl --clamp Cai CaSR RyR_Po [N]
#
# Examples:
#   julia --project=. src/generate.jl 300 1 10
#   julia --project=. src/generate.jl --clamp 0.1 1000 1.0 5

using CellSims

if first(ARGS) == "--clamp"
    Cai = parse(Float64, ARGS[2])
    CaSR = parse(Float64, ARGS[3])
    RyR_Po = parse(Float64, ARGS[4])
    N = parse(Int, get(ARGS, 5, "1"))

    sims = CellSims.all_clamp_simulations()

    sim = only(filter(sims) do sim
        sim.runner.Cai == Cai &&
        sim.runner.CaSR == CaSR &&
        sim.runner.RyR_Po == RyR_Po
    end)
else
    BCL = parse(Int, ARGS[1])
    ISO = parse(Int, ARGS[2])
    N = parse(Int, get(ARGS, 3, "1"))

    sims = CellSims.all_simulations()

    sim = only(filter(sims) do sim
        sim.runner.ISO == ISO && sim.runner.BCL == BCL
    end)
end

for _ in 1:N
    CellSims.run(sim)
end
