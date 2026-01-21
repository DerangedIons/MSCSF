# Usage:
#   Model3D:   julia --project=. src/generate_all.jl N
#   CaClamp3D: julia --project=. src/generate_all.jl --clamp N
#
# Examples:
#   julia --project=. src/generate_all.jl 10
#   julia --project=. src/generate_all.jl --clamp 5

using CellSims

if first(ARGS) == "--clamp"
    N = parse(Int, ARGS[2])
    sims = CellSims.all_clamp_simulations()

    @info "Running $(length(sims)) Ca clamp simulation settings, $N times each ($(length(sims) * N) total runs)"

    for i in 1:N
        @info "████████████████████████████████████████████████ Run $i/$N..."
        for sim in sims
            r = sim.runner
            @info "Running setting: Cai=$(r.Cai), CaSR=$(r.CaSR), ISO=$(r.ISO), RyR_Po=$(r.RyR_Po)"
            CellSims.run(sim)
        end
    end
else
    N = parse(Int, ARGS[1])
    sims = CellSims.all_simulations()

    @info "Running $(length(sims)) simulation settings, $N times each ($(length(sims) * N) total runs)"

    for i in 1:N
        @info "████████████████████████████████████████████████ Run $i/$N..."
        for sim in sims
            ISO, BCL = sim.runner.ISO, sim.runner.BCL
            @info "Running setting: BCL=$BCL, ISO=$ISO"
            CellSims.run(sim)
        end
    end
end

@info "All simulations complete!"
