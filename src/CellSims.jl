module CellSims

using CSV, DataFrames, Distributions, Statistics, StatsPlots, Scratch, DefaultApplication,
    MultivariateStats, LinearAlgebra, OrderedCollections

export Model3D, Model3DSimulations,
    reference, results, open_reference, open_results, all_simulations, get_df, preview, analyze

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

# for displaying plots on server
function preview(x)
    file = joinpath(@__DIR__, "..", "data", "preview.png")
    savefig(x, file)
    file
end

#-------------------------------------------------all----------------------------# Model3D
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
    State_Reference_read::String = "Off"
    State_Reference_write::String = "Off"
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

function Base.run(o::Model3D)
    # Optimal threads: total cores ÷ N simultaneous simulations
    n_threads = max(1, div(Sys.CPU_THREADS, length(all_simulations())))
    cd(DIR) do
        withenv("OMP_NUM_THREADS" => string(n_threads)) do
            run(Cmd(o))
        end
    end
end

reference(o::Model3D) = mkpath(joinpath(DIR, "Outputs_3Dcell_$(o.Reference)"))
results(o::Model3D) = mkpath(joinpath(DIR, "Outputs_3Dcell_$(o.Reference)", "Results_$(o.Results_Reference)"))

# Opens the "Reference" or "Results" directory in your file browser
open_reference(o::Model3D) = DefaultApplication.open(reference(o))
open_results(o::Model3D) = DefaultApplication.open(results(o))

#-----------------------------------------------------------------------------# Model3DSimulations
# Specification for generating SR datasets
struct Model3DSimulations
    Reference::String
    prepace::Model3D
    prepace_full::Model3D
    runner::Model3D

    function Model3DSimulations(Reference::String, prepace, prepace_full, runner)
        prepace.Reference = Reference
        prepace.Results_Reference = "prepace"
        prepace.Write_state = "ave"
        prepace.State_Reference_write = Reference

        prepace_full.Reference = Reference
        prepace_full.Results_Reference = "prepace_full"
        prepace_full.Read_state = "ave"
        prepace_full.Write_state = "On"
        prepace_full.State_Reference_read = Reference
        prepace_full.State_Reference_write = Reference

        runner.Reference = Reference
        runner.Results_Reference = "run_" * lpad(1, 4, '0')
        runner.Read_state = "On"
        runner.State_Reference_read = Reference
        new(Reference, prepace, prepace_full, runner)
    end
end

Base.show(io::IO, o::Model3DSimulations) = print(io, "Model3DSimulations: $(repr(o.Reference))")

open_reference(o::Model3DSimulations) = DefaultApplication.open(reference(o))

# Default Model is "minimal".  Should we be using ord model?
# Note: Ca input to generative PCA model is last row of prepace_full *Ca_NSR*
function Model3DSimulations(; ISO=1, BCL=1000, Total_time=2000)
    common = (Jup_scale=2, tau_ss_type="medium_fast", ISO, BCL)
    prepace =      Model3D(; Beats=50, Sim_cell_size="testing", common...)
    prepace_full = Model3D(; Beats=4, Sim_cell_size="full", common...)
    runner =       Model3D(; Beats=1, Sim_cell_size="full", Total_time, common...)
    Model3DSimulations("sr_$(BCL)_$(ISO)", prepace, prepace_full, runner)
end


all_simulations() = vec([Model3DSimulations(; ISO, BCL) for (ISO, BCL) in Iterators.product(0:1, 300:200:1500)])

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

#-----------------------------------------------------------------------------# loading DataFrames
function load_cru_file(path::String)
    isfile(path) || error("File does not exist: $path")
    df = CSV.read(path, DataFrame; header=false)
    rename!(df, cru_cols)
end
function load_cru_files(paths::Vector{String})
    dfs = [load_cru_file(file) for file in paths]
    vcat(dfs...; cols=:union, source="run")
end

function get_prepace_df(o::Model3DSimulations)
    file = joinpath(reference(o), "Results_prepace", "CRU.txt")
    load_cru_file(file)
end

function get_prepace_full_df(o::Model3DSimulations)
    file = joinpath(reference(o), "Results_prepace_full", "CRU.txt")
    load_cru_file(file)
end

function get_df(o::Model3DSimulations)
    dirs = filter(readdir(reference(o))) do dir
        startswith(dir, "Results_run_")
    end
    files = joinpath.(reference(o), dirs, Ref("CRU.txt"))
    load_cru_files(files)
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
function _get_apd90(df, threshold = 0.05)
    a, b = extrema(df.Vm)
    cutoff = threshold * (b - a) + a

    i = findfirst(==(b), df.Vm)  # Find index of max Vm
    ti = findnext(<(cutoff), df.Vm, i)  # Find index when Vm falls below cutoff after max
    isnothing(ti) ? -1 : df.t[ti]
end
get_apd90(df::DataFrame) = rename!(combine(_get_apd90, groupby(df, "run")), :x1 => :apd90)


# Applied to a SubDataFrame (single run)
function _get_colman_stats(df, threshold = 0.05)
    i = findfirst(>(_get_apd90(df, threshold)), df.t)
    df2 = @view df[i:end, :]

    a, b = extrema(df.RyR_OA)
    cutoff = threshold * (b - a) + a

    ti_i = findfirst(>(cutoff), df2.RyR_OA)
    ti = isnothing(ti_i) ? -1 : df2.t[ti_i]
    tf_i = findlast(>(cutoff), df2.RyR_OA)
    tf = isnothing(tf_i) ? -1 : df2.t[tf_i]
    λ = isnothing(ti_i) ? -1 : tf - ti
    peak = maximum(df2.RyR_OA)
    plat = median(df2.RyR_OA)
    tp = df2.t[findfirst(==(peak), df2.RyR_OA)]
    (; ti, tf, λ, tp, peak, plat)
end

get_colman_stats(df::DataFrame) = combine(groupby(df, "run")) do sdf
    _get_colman_stats(sdf)
end

#-----------------------------------------------------------------------------# analyze
function analyze(sim::Model3DSimulations)
    out = OrderedDict()

    df = get_df(sim)
    prepace = get_prepace_df(sim)
    prepace_full = get_prepace_full_df(sim)

    out[:df] = df
    out[:prepace] = prepace
    out[:prepace_full] = prepace_full
    out[:colman] = get_colman_stats(df)
    out[:apd90] = get_apd90(df)
    out[:n] = length(unique(df.run))
    out[:n_sr] = count(row -> row.λ != -1, eachrow(out[:colman]))

    out[:pscr] = out[:n_sr] / out[:n]

    out[:plot] = plot(df.t, df.RyR_OA, g=df.run, xlab="t", ylab="RyR_OA", title="$(sim.Reference)")

    runs_with_sr = filter(!=(-1), out[:colman].ti)
    out[:df_sr_only] = filter(row -> row.run in runs_with_sr, df)

    return out
end

function analyze(sims::Vector{Model3DSimulations})
    OrderedDict(
        sim.Reference => analyze(sim) for sim in sims
    )
end

#-----------------------------------------------------------------------------# Generative PCA
struct ApproxLog
    eps::Float64
end
(o::ApproxLog)(x) = log(x + o.eps)
inverse(o::ApproxLog) = x -> exp(x) - o.eps

inverse(::typeof(identity)) = identity


function fit_generative_pca(df_sr::AbstractDataFrame, f = ApproxLog(0.1); k=25, pratio=.999)
    isempty(df_sr) && error("df_sr is empty; cannot fit PCA model")
    f2 = inverse(f)

    w = reduce(hcat, f.(sub.RyR_OA) for sub in groupby(df_sr, "run"))

    model = fit(PCA, w; maxoutdim=k, pratio);

    z = MultivariateStats.transform(model, w)

    μ = vec(mean(z, dims=2))
    Σ = cov(z; dims=2)
    dist = MvNormal(μ, Symmetric(Σ))

    minval = minimum(df_sr.RyR_OA)

    rand_waveform() = max.(minval, reconstruct(model, rand(dist)))
end

end
