# Example: julia --project=. src/generate_all.jl 10
# Runs all simulation settings N times each

N = parse(Int, ARGS[1])

using CellSims

sims = CellSims.all_simulations()

@info "Running $(length(sims)) simulation settings, $N times each ($(length(sims) * N) total runs)"

for i in 1:N
    @info "████████████████████████████████████████████████ Run $i/$N..."
    for sim in sims
        ISO, BCL = sim.runner.ISO, sim.runner.BCL
        @info "Running setting: BCL=$BCL, ISO=$ISO ($N times)"
        CellSims.run(sim)
    end
end

@info "All simulations complete!"
