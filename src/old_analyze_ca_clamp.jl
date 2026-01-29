# This assumes the `generate_all.jl` script has already been run to generate simulation data

using CellSims, ProgressMeter, OrderedCollections, Plots, DataFrames

sims = all_clamp_simulations()

#-----------------------------------------------------------------------------# Get analysis dict
# Suppress warnings for when generative PCA fitting fails
CellSims.QUIET = true

dict = CellSims.analyze(sims)

#-----------------------------------------------------------------------------# Make summary DataFrame
df = DataFrame(
    sim = collect(keys(dict)),
    cai = [v[:cai] for v in values(dict)],
    casr = [v[:casr] for v in values(dict)],
    ryr_po = [v[:ryr_po] for v in values(dict)],
    iso = [v[:iso] for v in values(dict)],
    n = [v[:n] for v in values(dict)],
    n_sr = [v[:n_sr] for v in values(dict)],
    pscr = [v[:pscr] for v in values(dict)],
)

#-----------------------------------------------------------------------------# Plots
dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "ca_clamp_summary"))

savefig(scatter(df.casr, df.pscr, group=df.ryr_po, xlab="CaSR (µM)", ylab="P(SCR)", legendtitle="RyR Po"), "$dir/pscr_vs_casr.png")

savefig(scatter(df.ryr_po, df.pscr, group=df.casr, xlab="RyR Po", ylab="P(SCR)", legendtitle="CaSR (µM)"), "$dir/pscr_vs_ryr_po.png")

savefig(scatter(df.casr, df.n_sr, group=df.ryr_po, xlab="CaSR (µM)", ylab="N(SCR)", legendtitle="RyR Po"), "$dir/n_sr_vs_casr.png")

df
