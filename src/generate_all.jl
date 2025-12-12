# Example: julia --project=. src/generate_all.jl 10
# Runs all simulation settings N times each

N = parse(Int, ARGS[1])

using MSCSF

sims = MSCSF.all_simulations()

@info "Running $(length(sims)) simulation settings, $N times each ($(length(sims) * N) total runs)"

for (settings, sim) in sims
    ISO, BCL = settings.ISO, settings.BCL
    @info "Running setting: BCL=$BCL, ISO=$ISO ($N times)"

    for i in 1:N
        @info"████████████████████████████████████████████████ Run $i/$N..."
        MSCSF.run(sim)
    end
end

@info "All simulations complete!"
