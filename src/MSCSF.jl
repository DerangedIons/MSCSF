module MSCSF

using CSV, DataFrames, Distributions, Statistics, StatsPlots

dir::String = ""

function __init__()
    global dir = mkpath(abspath(joinpath(@__DIR__, "..", "data")))
end

output_dir(ref = "SR") = joinpath(dir, "Outputs_3Dcell_$ref")

#-----------------------------------------------------------------------------# sr_df
# Step 4: Get all runs into single DataFrame
function get_df(ref = "SR")
    out = output_dir(ref)
    dirs = filter(readdir(out, join=true)) do path
        startswith(basename(path), "Results_run_")
    end
    files = joinpath.(dirs, Ref("CRU.txt"))
    dfs = [CSV.read(file, DataFrame, header=false) for file in files]
    rename!.(dfs, Ref(cru_cols))
    df = vcat(dfs...; cols=:union, source="run")
    select!(df, "run", All())
end

# Field names for CRU.txt output files
cru_cols = [
    "t",               # ms
    "Vm",              # mV
    "Ca_ds",           # [Ca²⁺]ds (µM)
    "Ca_ss",           # [Ca²⁺]ss (µM)
    "Ca_cyto",         # [Ca²⁺]cyto (µM)
    "Ca_JSR",          # [Ca²⁺]JSR (µM)
    "Ca_NSR",          # [Ca²⁺]NSR (µM)
    "Jrel",            # release flux (µM/ms)
    "RyR_OA",          # fraction RyR open activated
    "RyR_OI",          # fraction RyR open inactivated
    "RyR_CA",          # fraction RyR closed activated
    "RyR_CI",          # fraction RyR closed inactivated
    "Monomer_state",   # csqn monomer state
    "N_active_CRUs",   # number of active CRUs
    "N_active_frac",   # N_active / N_tot
    "JCaL",            # LTCC flux (µM/ms)
    "LTCC_open_frac",  # fraction open LTCC
    "JSERCA",          # uptake flux (µM/ms)
    "Jleak",           # leak flux (µM/ms)
    "JNCX_cyto",       # NaCa exchanger, cyto
    "JNCX_ss",         # NaCa exchanger, subspace
    "JCaP_cyto",       # PMCA, cyto
    "JCaP_ss",         # PMCA, subspace
    "JCab_cyto",       # background Ca²⁺ current, cyto
    "JCab_ss"          # background Ca²⁺ current, subspace
]

#-----------------------------------------------------------------------------# stats
function stats(df, col = "Ca_cyto")
    # Determine SR threshold
    x = df[!, col]
    cutoff = 0.1 * (maximum(x) - minimum(x)) + minimum(x)

    @info "Getting stats for $col with cutoff $cutoff."
    inits = combine(groupby(df, "run"), first)

    # Drop the initial 300ms
    df = filter(row -> row.t > 300, df)

    # Get stats
    g = groupby(df, "run")
    out = combine(g,
        col => (x -> findfirst(>(cutoff), x)) => :i,  # Start time of SR
        col => (x -> findlast(>(cutoff), x)) => :j,  # End time of SR
        col => findmax => :max,  # (max value, index)
    )
    out[!, :init] = inits[!, col]
    out[!, :duration] = out.j .- out.i
    out
end

end
