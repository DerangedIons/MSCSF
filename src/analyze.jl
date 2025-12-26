# This assumes the `generate_all.jl` script has already been run to generate simulation data

using CellSims, ProgressMeter, OrderedCollections, Plots

sims = all_simulations()

out = OrderedDict()

for (i, sim) in enumerate(sims)
    @info "Analyzing simulation $(i)/$(length(sims)): BCL=$(sim.runner.BCL), ISO=$(sim.runner.ISO)..."
    n = length(results(sim))

    @info "    Number of runs: $n"

    if n == 0
        @warn "    No runs!  Skipping analysis."
        continue
    end


    out[sim.Reference] = CellSims.analyze(sim)
end

@info "Analysis complete!  $(length(out)) / $(length(sims)) simulations analyzed."


# init_Ca_cyto = OrderedDict(
#     k => v.df.Ca_cyto[400] for (k, v) in results
# )

# v = collect(values(init_Ca_cyto))

# barplot(collect(keys(init_Ca_cyto)), collect(values(init_Ca_cyto)), xlim=extrema(v))


# results = Dict(
#     sim.Reference => CellSims.analyze(sim) for sim in sims
# )


# init = Dict(
#     CellSims.get_df(sim).
# )

# results = map(sims) do sim
#     @info sim
#     df = CellSims.get_df(sim)
#     stats = CellSims.stats(df)
#     (; df, stats)
# end
