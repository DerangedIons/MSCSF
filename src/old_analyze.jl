# This assumes the `generate_all.jl` script has already been run to generate simulation data

using CellSims, ProgressMeter, OrderedCollections, Plots, DataFrames

sims = all_simulations()

#-----------------------------------------------------------------------------# Get analysis dict
# Suppress warnings for when generative PCA fitting fails
CellSims.QUIET = true

dict = CellSims.analyze(sims)

#-----------------------------------------------------------------------------# Make summary DataFrame
df = DataFrame(
    sim = collect(keys(dict)),
    bcl = [v[:bcl] for v in values(dict)],
    iso = [v[:iso] for v in values(dict)],
    n = [v[:n] for v in values(dict)],
    n_sr = [v[:n_sr] for v in values(dict)],
    pscr = [v[:pscr] for v in values(dict)],
    prepace_final_ca_nsr = [v[:prepace_full][end, :Ca_NSR] for v in values(dict)],
)

#-----------------------------------------------------------------------------# Plots
dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "summary"))

savefig(scatter(df.bcl, df.pscr, group=df.iso, xlab="BCL", ylab="P(SCR)", legendtitle="ISO"), "$dir/pscr_vs_bcl.png")

savefig(scatter(df.bcl, df.prepace_final_ca_nsr, group=df.iso, xlab="BCL", ylab="Final Ca_NSR (prepace)", legendtitle="ISO"), "$dir/ca_nsr_vs_bcl.png")

df
