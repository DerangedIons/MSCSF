module CellSims

using CSV, DataFrames, Distributions, Statistics, StatsPlots, Scratch, DefaultApplication,
    MultivariateStats, LinearAlgebra, OrderedCollections

export Model3D, Model3DSimulations, CaClamp3D, CaClamp3DSimulations,
    reference, results, open_reference, open_results, all_simulations, all_clamp_simulations, get_df, preview,
    generate, generate_all, generate_all_clamp, analyze, summarize, summarize_clamp, @trycatch,
    plot_all_clamp_simulations

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

QUIET::Bool = false

macro trycatch(ex)
    quote
        try
            $(esc(ex))
        catch e
            if !CellSims.QUIET
                @warn "Error occurred in expression" exception=(e, catch_backtrace()) location=$(__source__)
            end
        end
    end
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

#-----------------------------------------------------------------------------# CaClamp3D
# Specification of a single run of the 3D calcium clamp model
@kwdef mutable struct CaClamp3D
    bin::String = joinpath(@__DIR__, "..", "CODE", "model_Ca_clamp_3D")
    Model::String = "minimal"
    Jup_scale::Float64 = 1.0
    Jrel_scale::Float64 = 1.0   # Scale factor for release flux magnitude
    RyR_Po::Float64 = 1.0       # Scale factor for RyR Ca-dependent open rate (sensitivity)
    tau_ss_type::String = "medium_fast"
    Spatial_output_interval_data::Int = 0
    Spatial_output_interval_vtk::Int = 0
    Reference::String = "temp"
    Results_Reference::String = "temp"
    Total_time::Int = 1500
    Sim_cell_size::String = "full"
    Cai::Float64 = 0.1      # Initial cytosolic Ca (µM) - clamp value
    CaSR::Float64 = 1000.0  # Initial SR Ca (µM) - clamp value
end

function Base.Cmd(o::CaClamp3D)
    out = [o.bin]
    ref_dir = joinpath(DIR, "Outputs_Ca_clamp_3Dcell_$(o.Reference)")
    res_dir = joinpath(ref_dir, "Results_$(o.Results_Reference)")
    mkpath(res_dir)

    for name in setdiff(fieldnames(CaClamp3D), (:bin,))
        value = getfield(o, name)
        push!(out, string(name), string(value))
    end
    Cmd(out)
end

function Base.run(o::CaClamp3D)
    n_threads = max(1, div(Sys.CPU_THREADS, length(all_simulations())))
    cd(DIR) do
        withenv("OMP_NUM_THREADS" => string(n_threads)) do
            run(Cmd(o))
        end
    end
end

reference(o::CaClamp3D) = mkpath(joinpath(DIR, "Outputs_Ca_clamp_3Dcell_$(o.Reference)"))
results(o::CaClamp3D) = mkpath(joinpath(DIR, "Outputs_Ca_clamp_3Dcell_$(o.Reference)", "Results_$(o.Results_Reference)"))

open_reference(o::CaClamp3D) = DefaultApplication.open(reference(o))
open_results(o::CaClamp3D) = DefaultApplication.open(results(o))

#-----------------------------------------------------------------------------# CaClamp3DSimulations
# Specification for generating calcium clamp datasets (no prepacing needed due to clamping)
struct CaClamp3DSimulations
    Reference::String
    runner::CaClamp3D

    function CaClamp3DSimulations(Reference::String, runner::CaClamp3D)
        runner.Reference = Reference
        runner.Results_Reference = "run_" * lpad(1, 4, '0')
        new(Reference, runner)
    end
end

Base.show(io::IO, o::CaClamp3DSimulations) = print(io, "CaClamp3DSimulations: $(repr(o.Reference))")

function CaClamp3DSimulations(; CaSR=1000.0, RyR_Po=1.0, Total_time=1500)
    runner = CaClamp3D(;
        CaSR,
        RyR_Po,
        Total_time,
        Jup_scale=2.0,
        tau_ss_type="medium_fast",
        Sim_cell_size="full"
    )
    CaClamp3DSimulations("ca_clamp_CaSR$(CaSR)_Po$(RyR_Po)", runner)
end

function all_clamp_simulations(; CaSR = [400:100:1600; 1800; 2000], RyR_Po = [1])
    vec([
        CaClamp3DSimulations(; CaSR, RyR_Po) for (CaSR, RyR_Po) in Iterators.product(CaSR, RyR_Po)
    ])
end

reference(o::CaClamp3DSimulations) = reference(o.runner)

function results(o::CaClamp3DSimulations)
    filter(readdir(reference(o))) do dir
        startswith(dir, "Results_run_")
    end
end

Base.length(o::CaClamp3DSimulations) = length(results(o))

function Base.run(o::CaClamp3DSimulations)
    n = length(o) + 1
    o.runner.Results_Reference = "run_" * lpad(n, 4, '0')
    @info "Running Ca clamp simulation run $n for $(o.Reference)..."
    run(o.runner)
end

open_reference(o::CaClamp3DSimulations) = DefaultApplication.open(reference(o))

@recipe function f(o::CaClamp3DSimulations, col="RyR_OA"; runs=nothing)
    title --> "$(o.Reference) Ca Clamp Simulations"
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

function get_df(o::CaClamp3DSimulations)
    dirs = filter(readdir(reference(o))) do dir
        startswith(dir, "Results_run_")
    end
    files = joinpath.(reference(o), dirs, Ref("CRU.txt"))
    load_cru_files(files)
end

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

#-----------------------------------------------------------------------------# generate
# Run a simulation N times
function generate(sim::Union{Model3DSimulations, CaClamp3DSimulations}, N::Int=1)
    for _ in 1:N
        run(sim)
    end
end

# Run multiple simulations N times each
function generate(sims::Vector{Model3DSimulations}, N::Int=1)
    @info "Running $(length(sims)) simulation settings, $N times each ($(length(sims) * N) total runs)"
    for i in 1:N
        @info "████████████████████████████████████████████████ Run $i/$N..."
        for sim in sims
            ISO, BCL = sim.runner.ISO, sim.runner.BCL
            @info "Running setting: BCL=$BCL, ISO=$ISO"
            run(sim)
        end
    end
    @info "All simulations complete!"
end

function generate(sims::Vector{CaClamp3DSimulations}, N::Int=1)
    @info "Running $(length(sims)) Ca clamp simulation settings, $N times each ($(length(sims) * N) total runs)"
    for i in 1:N
        @info "████████████████████████████████████████████████ Run $i/$N..."
        for sim in sims
            r = sim.runner
            @info "Running setting: CaSR=$(r.CaSR), RyR_Po=$(r.RyR_Po)"
            run(sim)
        end
    end
    @info "All simulations complete!"
end

# Convenience: generate all Model3D simulations
generate_all(N::Int=1) = generate(all_simulations(), N)

# Convenience: generate all CaClamp3D simulations
generate_all_clamp(N::Int=1) = generate(all_clamp_simulations(), N)

# Find and run a specific Model3D simulation
function generate(N::Int, BCL::Int, ISO::Int)
    sims = all_simulations()
    sim = only(filter(sims) do sim
        sim.runner.ISO == ISO && sim.runner.BCL == BCL
    end)
    generate(sim, N)
end

# Find and run a specific CaClamp3D simulation
function generate(N::Int; CaSR::Float64, RyR_Po::Float64)
    sims = all_clamp_simulations()
    sim = only(filter(sims) do sim
        sim.runner.CaSR == CaSR &&
        sim.runner.RyR_Po == RyR_Po
    end)
    generate(sim, N)
end

#-----------------------------------------------------------------------------# loading DataFrames
function load_cru_file(path::String)
    isfile(path) || error("File does not exist: $path")
    df = CSV.read(path, DataFrame; header=false)
    rename!(df, cru_cols)
end
function load_cru_files(paths::Vector{String})
    dfs = DataFrame[]
    for file in paths
        filesize(file) > 0 || continue
        push!(dfs, load_cru_file(file))
    end
    isempty(dfs) && return DataFrame()
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

    dir = mkpath(joinpath(@__DIR__, "..", "data", "results", sim.Reference))

    df = get_df(sim)
    prepace = get_prepace_df(sim)
    prepace_full = get_prepace_full_df(sim)

    out[:bcl] = sim.runner.BCL
    out[:iso] = sim.runner.ISO
    out[:df] = df
    out[:prepace] = prepace
    out[:prepace_full] = prepace_full
    out[:colman] = get_colman_stats(df)
    out[:apd90] = get_apd90(df)
    out[:n] = length(unique(df.run))
    out[:n_sr] = count(row -> row.λ != -1, eachrow(out[:colman]))

    out[:pscr] = out[:n_sr] / out[:n]

    out[:plot] = plot(df.t, df.RyR_OA, g=df.run, xlab="t", ylab="RyR_OA", lab="", title="$(sim.Reference)")
    savefig(out[:plot], "$dir/plot.png")

    runs_with_sr = filter(row -> row.ti != -1, out[:colman]).run
    out[:df_sr_only] = filter(row -> row.run in runs_with_sr, df)

    @trycatch out[:gen_pca] = GenerativePCA(out[:df_sr_only], identity)
    @trycatch out[:gen_pca_log] = GenerativePCA(out[:df_sr_only], ApproxLog(0.1))

    if haskey(out, :gen_pca)
        waves = [out[:gen_pca]() for _ in 1:10]
        p = plot(waves, xlab="Time (ms)", ylab="RyR_OA", title="Generative PCA Samples for $(sim.Reference)")
        savefig(p, "$dir/gen_pca_samples.png")
    end
    if haskey(out, :gen_pca_log)
        waves = [out[:gen_pca_log]() for _ in 1:10]
        p = plot(waves, xlab="Time (ms)", ylab="RyR_OA", title="Generative Log PCA Samples for $(sim.Reference)")
        savefig(p, "$dir/gen_pca_log_samples.png")
    end

    return out
end

function analyze(sims::Vector{Model3DSimulations})
    OrderedDict(
        sim.Reference => (@info("Analyzing: $(sim.Reference)"); analyze(sim)) for sim in sims
    )
end

function analyze(sim::CaClamp3DSimulations)
    out = OrderedDict()

    dir = mkpath(joinpath(@__DIR__, "..", "data", "results", sim.Reference))

    df = get_df(sim)

    out[:casr] = sim.runner.CaSR
    out[:ryr_po] = sim.runner.RyR_Po
    out[:df] = df
    out[:colman] = get_colman_stats(df)
    out[:n] = length(unique(df.run))
    out[:n_sr] = count(row -> row.λ != -1, eachrow(out[:colman]))

    out[:pscr] = out[:n_sr] / out[:n]

    out[:plot] = plot(df.t, df.RyR_OA, g=df.run, xlab="t", ylab="RyR_OA", lab="", title="$(sim.Reference)")
    savefig(out[:plot], "$dir/plot.png")

    runs_with_sr = filter(row -> row.ti != -1, out[:colman]).run
    out[:df_sr_only] = filter(row -> row.run in runs_with_sr, df)

    @trycatch out[:gen_pca] = GenerativePCA(out[:df_sr_only], identity)
    @trycatch out[:gen_pca_log] = GenerativePCA(out[:df_sr_only], ApproxLog(0.1))

    if haskey(out, :gen_pca)
        waves = [out[:gen_pca]() for _ in 1:10]
        p = plot(waves, xlab="Time (ms)", ylab="RyR_OA", title="Generative PCA Samples for $(sim.Reference)")
        savefig(p, "$dir/gen_pca_samples.png")
    end
    if haskey(out, :gen_pca_log)
        waves = [out[:gen_pca_log]() for _ in 1:10]
        p = plot(waves, xlab="Time (ms)", ylab="RyR_OA", title="Generative Log PCA Samples for $(sim.Reference)")
        savefig(p, "$dir/gen_pca_log_samples.png")
    end

    return out
end

function analyze(sims::Vector{CaClamp3DSimulations})
    OrderedDict(
        sim.Reference => (@info("Analyzing: $(sim.Reference)"); analyze(sim)) for sim in sims
    )
end

#-----------------------------------------------------------------------------# summarize
function summarize(sims::Vector{Model3DSimulations}; quiet=true)
    old_quiet = QUIET
    QUIET = quiet

    dict = analyze(sims)

    # Make summary DataFrame
    df = DataFrame(
        sim = collect(keys(dict)),
        bcl = [v[:bcl] for v in values(dict)],
        iso = [v[:iso] for v in values(dict)],
        n = [v[:n] for v in values(dict)],
        n_sr = [v[:n_sr] for v in values(dict)],
        pscr = [v[:pscr] for v in values(dict)],
        prepace_final_ca_nsr = [v[:prepace_full][end, :Ca_NSR] for v in values(dict)],
    )

    # Plots
    dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "summary"))
    savefig(scatter(df.bcl, df.pscr, group=df.iso, xlab="BCL", ylab="P(SCR)", legendtitle="ISO"), "$dir/pscr_vs_bcl.png")
    savefig(scatter(df.bcl, df.prepace_final_ca_nsr, group=df.iso, xlab="BCL", ylab="Final Ca_NSR (prepace)", legendtitle="ISO"), "$dir/ca_nsr_vs_bcl.png")

    QUIET = old_quiet
    return (; dict, df)
end
summarize(; kw...) = summarize(all_simulations(); kw...)

function summarize(sims::Vector{CaClamp3DSimulations}; quiet=true)
    old_quiet = QUIET
    QUIET = quiet

    dict = analyze(sims)

    # Make summary DataFrame
    df = DataFrame(
        sim = collect(keys(dict)),
        casr = [v[:casr] for v in values(dict)],
        ryr_po = [v[:ryr_po] for v in values(dict)],
        n = [v[:n] for v in values(dict)],
        n_sr = [v[:n_sr] for v in values(dict)],
        pscr = [v[:pscr] for v in values(dict)],
    )

    # Plots
    dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "ca_clamp_summary"))
    savefig(scatter(df.casr, df.pscr, group=df.ryr_po, xlab="CaSR (µM)", ylab="P(SCR)", legendtitle="RyR Po"), "$dir/pscr_vs_casr.png")
    savefig(scatter(df.ryr_po, df.pscr, group=df.casr, xlab="RyR Po", ylab="P(SCR)", legendtitle="CaSR (µM)"), "$dir/pscr_vs_ryr_po.png")
    savefig(scatter(df.casr, df.n_sr, group=df.ryr_po, xlab="CaSR (µM)", ylab="N(SCR)", legendtitle="RyR Po"), "$dir/n_sr_vs_casr.png")

    QUIET = old_quiet
    return (; dict, df)
end
summarize_clamp(; kw...) = summarize(all_clamp_simulations(); kw...)

#-----------------------------------------------------------------------------# plot_all_clamp_simulations
"""
    plot_all_clamp_simulations(sims=all_clamp_simulations(); col="RyR_OA", size=(1600, 1200))

Create a grid plot with subplots for each calcium clamp simulation.
All subplots share the same y-axis limits for easy comparison.

# Arguments
- `sims`: Vector of CaClamp3DSimulations (default: all_clamp_simulations())
- `col`: Column to plot (default: "RyR_OA")
- `size`: Figure size tuple (default: (1600, 1200))
"""
function plot_all_clamp_simulations(sims::Vector{CaClamp3DSimulations}=all_clamp_simulations();
                                     col::String="RyR_OA", size=(1600, 1200))
    # Load all data first to determine global y-axis limits
    data = [(sim, get_df(sim)) for sim in sims]
    filter!(x -> !isempty(x[2]), data)

    isempty(data) && error("No simulation data found")

    # Get global y-axis limits across all simulations
    ymin = minimum(minimum(df[!, col]) for (_, df) in data)
    ymax = maximum(maximum(df[!, col]) for (_, df) in data)

    n = length(data)
    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)

    plots = Plots.Plot[]
    for (i, (sim, df)) in enumerate(data)
        p = plot(; title="CaSR=$(sim.runner.CaSR)", titlefontsize=10,
                   xlabel=(i > n - ncols) ? "Time (ms)" : "",
                   ylabel=(mod1(i, ncols) == 1) ? col : "",
                   ylims=(ymin, ymax), legend=false)
        for sub in groupby(df, "run")
            plot!(p, sub.t, sub[!, col], linewidth=1.5, alpha=0.7)
        end
        push!(plots, p)
    end

    p = plot(plots...; layout=(nrows, ncols), size=size,
             plot_title="Calcium Clamp Simulations", plot_titlefontsize=14)

    # Save to data/results folder
    dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "ca_clamp_summary"))
    path = joinpath(dir, "all_clamp_simulations_$(col).png")
    savefig(p, path)
    @info "Saved plot to $path"

    return p
end

#-----------------------------------------------------------------------------# Generative PCA
struct ApproxLog
    eps::Float64
end
(o::ApproxLog)(x) = log(x + o.eps)
inverse(o::ApproxLog) = x -> exp(x) - o.eps

inverse(::typeof(identity)) = identity


@kwdef struct GenerativePCA{F, D}
    fun::F = ApproxLog(0.1)
    model::PCA{Float64}
    dist::D
    min::Float64
end
function GenerativePCA(df_sr::AbstractDataFrame, f = ApproxLog(0.1); k=25, pratio=.999)
    isempty(df_sr) && error("df_sr is empty; cannot fit PCA model")

    gdf = groupby(df_sr, "run")

    # If simulations are currently running, the last run might have a different number of rows
    # We'll drop simulations that don't have the same n as the first run.
    n = nrow(first(gdf))
    gdf2 = filter(sub -> nrow(sub) == n, gdf)

    w = reduce(hcat, f.(sub.RyR_OA) for sub in gdf2)
    model = fit(PCA, w; maxoutdim=k, pratio);
    z = MultivariateStats.transform(model, w)
    μ = vec(mean(z, dims=2))
    Σ = Symmetric(cov(z; dims=2))
    dist = MvNormal(μ, Σ)
    minval = minimum(df_sr.RyR_OA)
    GenerativePCA(f, model, dist, minval)
end

function (o::GenerativePCA)()
    f_inv = inverse(o.fun)
    w = reconstruct(o.model, rand(o.dist))
    max.(o.min, f_inv.(w))
end

#-----------------------------------------------------------------------------# generate_vtk_with_sr
function generate_vtk(BCL=300, ISO=1, Total_time=2000)
    sim = Model3DSimulations(; ISO, BCL)
    sim.runner.Spatial_output_interval_data = 1
    sim.runner.Spatial_output_interval_vtk = 1
    sim.runner.Total_time = Total_time
    run(sim)
    return joinpath(reference(sim), "Spatial_" * results(sim)[end])
end

end
