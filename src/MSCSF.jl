module MSCSF

using CSV, DataFrames

dir::String = ""

function __init__()
    global dir = mkpath(abspath(joinpath(@__DIR__, "..", "data")))

    # copy binaries to data/
    for file in readdir(joinpath(@__DIR__, "..", "CODE"); join=true)
        startswith(basename(file), "model_") && cp(file, joinpath(dir, basename(file)); force=true)
    end

    # Create PATH.txt and the state_and_geometry_files directory
    write(joinpath(dir, "PATH.txt"), mkpath(joinpath(dir, "state_and_geometry_files")))
    mkpath(joinpath(dir, "state_and_geometry_files", "State_files", "Single_cell"))
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

function _cmd(kw)
    out = [joinpath(dir, "model_single_3D")]
    for (k, v) in kw
        push!(out, string(k), string(v))
    end
    return Cmd(out)
end

# Run model and return kw arguments
_run(; kw...) = cd(dir) do
    @time run(_cmd(kw))
    NamedTuple(kw)
end

#-----------------------------------------------------------------------------# sr_prepace
# Step 1 (~2.5 minutes for 50 beats at 250ms BCL)
function sr_prepace(; Beats=50, BCL=350, kw...)
    kw = _run(;
        Model = "minimal",
        Beats = Beats,
        BCL = BCL,
        ISO = 1,
        Jup_scale = 2,
        tau_ss_type = "medium_fast",
        Sim_cell_size = "testing",
        Total_time = Beats * BCL,
        Write_state = "ave",
        Reference = "SR",
        Results_Reference = "prepace",
        Spatial_output_interval_data = 0,
        Spatial_output_interval_vtk = 0,
        kw...
    )
    joinpath(dir, "Outputs_3Dcell_$(kw.Reference)", "Results_$(kw.Results_Reference)")
end


#-----------------------------------------------------------------------------# sr_prepace_full
# Step 2 (~8 minutes for 4 beats at 350ms BCL)
function sr_prepace_full(; Beats=4, BCL=350, kw...)
    kw = _run(;
        Model = "minimal",
        Beats = Beats,
        BCL = BCL,
        ISO = 1,
        Jup_scale = 2,
        tau_ss_type = "medium_fast",
        Sim_cell_size = "full",
        Total_time = Beats * BCL,
        Read_state = "ave",
        Write_state = "On",
        Reference = "SR",
        Results_Reference = "prepace_full",
        Spatial_output_interval_data = 0,
        Spatial_output_interval_vtk = 0,
        kw...
    )
    joinpath(dir, "Outputs_3Dcell_$(kw.Reference)", "Results_$(kw.Results_Reference)")
end


#-----------------------------------------------------------------------------# sr_run
# Step 3 (~5 minutes/per 1000ms)
function sr_run(i::Integer; Beats=1, BCL=350, kw...)
    kw = _run(;
        Model = "minimal",
        Beats = Beats,
        BCL = BCL,
        ISO = 1,
        Jup_scale = 2,
        tau_ss_type = "medium_fast",
        Sim_cell_size = "full",
        Read_state = "On",
        Write_state = "On",
        Total_time = 1000,
        Reference = "SR",
        Results_Reference = "Run_$i",
        Spatial_output_interval_data = 0,
        Spatial_output_interval_vtk = 0,
        kw...
    )
    joinpath(dir, "Outputs_3Dcell_$(kw.Reference)", "Results_$(kw.Results_Reference)")
end


#-----------------------------------------------------------------------------# sr_df
# Step 4: Get all runs into single DataFrame
function sr_df(Reference = "SR")
    dirs = filter(readdir(joinpath(dir, "Outputs_3Dcell_$Reference"), join=true)) do path
        startswith(basename(path), "Results_Run_")
    end
    files = joinpath.(dirs, Ref("CRU.txt"))
    dfs = [CSV.read(file, DataFrame, header=false) for file in files]
    rename!.(dfs, Ref(cru_cols))
    df = vcat(dfs...; cols=:union, source="run")
    select!(df, "run", All())
end


end
