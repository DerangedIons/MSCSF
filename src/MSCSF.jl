module MSCSF

using CSV, DataFrames, Distributions, Statistics, StatsPlots, Scratch, DefaultApplication,
    MultivariateStats, LinearAlgebra

export Model3D, Model3DSimulations, sims_low, sims_high, get_df,
    reference, results, open_reference, open_results, stats

#-----------------------------------------------------------------------------# init
DIR::String = ""  # where data goes, e.g. $DIR/$Reference/$Results_Reference/
STATE_AND_GEOMETRY_FILES = ""  # Where things get saved to/read from when using `Read_state` and `Write_state`

function __init__()
    global DIR = @get_scratch!("data")
    global STATE_AND_GEOMETRY_FILES = joinpath(DIR, "state_and_geometry_files")

    # Directory needs to exist or model_single_3D code crashes
    write(joinpath(DIR, "PATH.txt"), STATE_AND_GEOMETRY_FILES)
    mkpath(joinpath(STATE_AND_GEOMETRY_FILES, "State_files", "Single_cell"))
end

#-----------------------------------------------------------------------------# Model3D
# Specification of a single run of the 3D model
@kwdef mutable struct Model3D
    bin::String = joinpath(@__DIR__, "..", "CODE", "model_single_3D")
    Model::String = "minimal"
    ISO::Int = 0
    Jup_scale::Float64 = 1.0
    BCL::Int = 1000
    tau_ss_type::String = "medium_fast"
    Spatial_output_interval_data::Int = 0
    Spatial_output_interval_vtk::Int = 0
    Reference::String = "temp"
    Results_Reference::String = "temp"
    Beats::Int = 1
    Read_state::String = "Off"
    Write_state::String = "Off"
    Total_time::Int = Beats * BCL
    Sim_cell_size::String = "full"
end

function Base.Cmd(o::Model3D)
    out = [o.bin]
    ref_dir = joinpath(DIR, "Outputs_3Dcell_$(o.Reference)")
    res_dir = joinpath(ref_dir, "Results_$(o.Results_Reference)")
    mkpath(res_dir)

    for name in setdiff(fieldnames(Model3D), (:bin,))
        value = getfield(o, name)
        push!(out, string(name), string(value))
    end
    Cmd(out)
end

Base.run(o::Model3D) = cd(() -> run(Cmd(o)), DIR)

reference(o::Model3D) = mkpath(joinpath(DIR, "Outputs_3Dcell_$(o.Reference)"))
results(o::Model3D) = mkpath(joinpath(DIR, "Outputs_3Dcell_$(o.Reference)", "Results_$(o.Results_Reference)"))

# Opens the "Reference" or "Results" directory in your file browser
open_reference(o) = DefaultApplication.open(reference(o))
open_results(o) = DefaultApplication.open(results(o))

#-----------------------------------------------------------------------------# Model3DSimulations
# Specification of multiple runs for generating SR
struct Model3DSimulations
    Reference::String
    prepace::Model3D
    prepace_full::Model3D
    runner::Model3D

    function Model3DSimulations(Reference::String, prepace, prepace_full, runner)
        prepace.Reference = Reference
        prepace.Results_Reference = "prepace"
        prepace.Write_state = "ave"

        prepace_full.Reference = Reference
        prepace_full.Results_Reference = "prepace_full"
        prepace_full.Read_state = "ave"
        prepace_full.Write_state = "On"

        runner.Reference = Reference
        runner.Results_Reference = "run_" * lpad(1, 4, '0')
        runner.Read_state = "On"
        new(Reference, prepace, prepace_full, runner)
    end
end

@recipe function f(o::Model3DSimulations, col="RyR_OA"; runs=nothing)
    title --> "$(o.Reference) Simulations"
    label --> ""
    xlabel --> "Time (ms)"
    ylabel --> col
    linewidth --> 1
    for (i, df) in enumerate(groupby(get_df(o), "run"))
        if isnothing(runs) || i in runs
            @series df.t, df[!, col]
        end
    end
end


reference(o::Model3DSimulations) = reference(o.runner)
function results(o::Model3DSimulations)
    filter(readdir(reference(o))) do dir
        startswith(dir, "Results_run_")
    end
end

Base.length(o::Model3DSimulations) = length(results(o))

function Base.run(o::Model3DSimulations)
    n = length(o) + 1
    o.runner.Results_Reference = "run_" * lpad(n, 4, '0')
    if !isfile(joinpath(reference(o), "Results_prepace", "CRU.txt"))
        @info "Running prepace for $(o.Reference)..."
        run(o.prepace)
    end
    if !isfile(joinpath(reference(o), "Results_prepace_full", "CRU.txt"))
        @info "Running full prepace for $(o.Reference)..."
        run(o.prepace_full)
    end
    @info "Running simulation run $n for $(o.Reference)..."
    run(o.runner)
end

#-----------------------------------------------------------------------------# sims_low
# Total_time is passed to the `runner`
function sims_low(; Total_time=1300)
    common = (ISO = 0, Jup_scale = 1, tau_ss_type = "medium_fast", BCL = 1000)
    prepace = Model3D(; Beats=40, Sim_cell_size="testing", common...)
    prepace_full = Model3D(; Beats=4, Sim_cell_size="full", common...)
    runner = Model3D(; Beats=1, Sim_cell_size="full", Total_time, common...)
    Model3DSimulations("sr_low", prepace, prepace_full, runner)
end

#-----------------------------------------------------------------------------# sims_high
function sims_high(; Total_time=1300)
    common = (ISO = 1, Jup_scale = 1, tau_ss_type = "medium_fast", BCL = 350)
    prepace = Model3D(; Beats=40, Sim_cell_size="testing", common...)
    prepace_full = Model3D(; Beats=4, Sim_cell_size="full", common...)
    runner = Model3D(; Beats=1, Sim_cell_size="full", Total_time, common...)
    Model3DSimulations("sr_high", prepace, prepace_full, runner)
end

#-----------------------------------------------------------------------------# get_df
function get_df(Reference::String, select=1:12)
    result_dirs = filter(readdir(joinpath(DIR, "Outputs_3Dcell_$Reference"))) do dir
        startswith(dir, "Results_run_")
    end
    files = joinpath.(DIR, "Outputs_3Dcell_$Reference", result_dirs, Ref("CRU.txt"))
    dfs = [CSV.read(file, DataFrame; select, header=false) for file in files]
    rename!.(dfs, Ref(cru_cols[select]))
    df = vcat(dfs...; cols=:union, source="run")
    select!(df, "run", All())
end

get_df(o::Model3D) = get_df(o.Reference)
get_df(o::Model3DSimulations) = get_df(o.Reference)

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
# Applied to a SubDataFrame (single run)
function _stats(df, threshold = 0.05, skip = 100)
    df
    # Determine SR threshold for APD90, start/stop of SR
    cutoff = threshold * (maximum(df.Ca_cyto) - minimum(df.Ca_cyto)) + minimum(df.Ca_cyto)

    # Drop APD90
    i = findfirst(>(skip), df.t)
    i = findnext(<(cutoff), df.Ca_cyto, i)  # Find index of first value below cutoff
    df2 = @view df[i:end, :]

    # Get stats
    ti_i = findfirst(>(cutoff), df2.Ca_cyto)
    ti = isnothing(ti_i) ? -1 : df2.t[ti_i]
    tf_i = findlast(>(cutoff), df2.Ca_cyto)
    tf = isnothing(tf_i) ? -1 : df2.t[tf_i]
    λ = isnothing(ti_i) ? -1 : tf - ti
    peak = maximum(df2.Ca_cyto)
    plat = median(df2.Ca_cyto)
    tp = df2.t[findfirst(==(peak), df2.Ca_cyto)]
    (; ti, tf, λ, tp, peak, plat)
end

stats(df) = combine(groupby(df, "run")) do sdf
    _stats(sdf)
end

stats(o::Model3DSimulations) = stats(get_df(o))

end
