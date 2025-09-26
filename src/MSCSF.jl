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
# Step 3 ()
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
        startswith(path, "Results_Run_")
    end
    files = joinpath.(dirs, Ref("CRU.txt"))
    dfs = [CSV.read(file, DataFrame, header=false) for file in files]
    rename!.(dfs, Ref(cru_cols))
    df = vcat(dfs...; cols=:union, source="run")
    select!(df, "run", All())
end









# # Get all the runs into a single DataFrame
# function single_3d_runs_df()
#     dirs = filter(readdir(project_path("Outputs_3Dcell_SR"), join=true)) do path
#         endswith(path, r"Results_Run_\d+") && !occursin("Spatial", path)
#     end
#     files = joinpath.(dirs, Ref("CRU.txt"))
#     dfs = [CSV.read(file, DataFrame, header=false) for file in files]
#     rename!.(dfs, Ref(cru_cols))
#     df = vcat(dfs...; cols=:union, source="run")
#     select!(df, "run", All())
# end

# # Get all the runs into a single DataFrame
# function single_3d_runs_df()
#     dirs = filter(readdir(project_path("Outputs_3Dcell_SR"), join=true)) do path
#         endswith(path, r"Results_Run_\d+") && !occursin("Spatial", path)
#     end
#     files = joinpath.(dirs, Ref("CRU.txt"))
#     dfs = [CSV.read(file, DataFrame, header=false) for file in files]
#     rename!.(dfs, Ref(cru_cols))
#     df = vcat(dfs...; cols=:union, source="run")
#     select!(df, "run", All())
# end

# @enum MODEL ca_clamp_0d ca_clamp_3d single_0d single_3d single_native tissue_0d tissue_native

# output_dir(x::MODEL) =
#     x == ca_clamp_0d ? "Outputs_Ca_clamp_0D" :
#     x == ca_clamp_3d ? "Outputs_Ca_clamp_3D" :
#     x == single_0d ? "Outputs_0Dcell" :
#     x == single_3d ? "Outputs_3Dcell" :
#     x == single_native ? "Outputs_single_native" :
#     x == tissue_native ? "Outputs_tissue_native" :
#     x == tissue_0d ? "Outputs_0Dtissue" :
#     error("Unknown model: $x")



# #-----------------------------------------------------------------------------# Project
# function run_model(model::MODEL, project::String; kw...)
#     cd(mkpath(dir)) do
#         # copy binary to `$dir/`
#         bin = joinpath(@__DIR__, "..", "CODE", string("model_", model))
#         isfile(bin) || error("Model binary does not exist: $bin")
#         cp(bin, joinpath(dir, string(model)); force=true)
#         mkpath(joinpath(dir, project, "State_files", "Single_cell"))
#         # Add PATH.txt file to `$dir/`
#         write("PATH.txt", abspath(project))
#         args = String["Spatial_output_interval_data", "0", "Spatial_output_interval_vtk", "0"]  # Don't save spatial output files by default
#         for (k, v) in kw
#             push!(args, string(k), string(v))
#         end
#         cmd = Cmd([joinpath(dir, string(model)), args...])
#         @info "Running command: $cmd"
#         run(cmd)
#     end
#     out = output_dir(model)
#     kwd = Dict(kw)
#     if haskey(kwd, :Reference)
#         out *= "_" * string(kwd[:Reference])
#     end
#     out = joinpath(out, "Results")
#     if haskey(kwd, :Results_Reference)
#         out *= "_" * string(kwd[:Results_Reference])
#     end
#     return joinpath(dir, out)
# end




# #-----------------------------------------------------------------------------# run_cmd
# function make_cmd(bin; kw...)
#     out = [string(bin)]
#     for (k, v) in kw
#         push!(out, string(k), string(v))
#     end
#     return Cmd(out)
# end

# run_cmd(bin; kw...) = run(make_cmd(bin; kw...))


# #-----------------------------------------------------------------------------# single_3d
# run_single_3d(project; kw...) = run_model(single_3d, project; kw...)

# # Takes ~2.5 minutes
# function run_single_3d_prepace(project; Model="minimal", Beats=50, BCL=350, ISO=1, Jup_scale=2,
#     tau_ss_type="medium_fast", Sim_cell_size="testing", Total_time = BCL * Beats, Write_state="ave",
#     Reference="SR", Results_Reference="prepace", kw...)
#     @time run_single_3d(project; Model, Beats, BCL, ISO, Jup_scale, tau_ss_type, Sim_cell_size, Total_time, Write_state, Reference, Results_Reference, kw...)
# end

# # Must call run_single_3d_prepace() first
# # Take 8 minutes  14:15:23 to ????
# function run_single_3d_prepace_full(project; Model="minimal", Beats=4, BCL=350, ISO=1, Jup_scale=2,
#     tau_ss_type = "medium_fast", Sim_cell_size="full", Total_time = BCL * Beats, Write_state="On",
#     Reference="SR", Results_Reference="prepace_full", kw...)
#     @time run_single_3d(project; Model, Beats, BCL, ISO, Jup_scale, tau_ss_type, Sim_cell_size, Total_time, Write_state, Reference, Results_Reference, kw...)
# end

# function run_single_3d_sr(project; Beats=1, BCL=350)
#     files = readdir(joinpath(dir, project, "Outputs_3dcell", "Results"))
# end

# # Takes ~3.5 minutes/run
# function single_3d_run(i::Integer)
#     # ./model_single_3D Model $model Beats 1 BCL ${BCL} Read_state On ISO 1 Jup_scale 2 tau_ss_type medium_fast Sim_cell_size full Spatial_output_interval_data 0 Reference Spontaneous_release Results_Reference Run_${run}
#     single_3d(;
#         Model = "minimal",
#         Beats = 1,
#         BCL = 350,
#         # Read_state = "On",
#         ISO = 1,
#         Jup_scale = 2,
#         tau_ss_type = "medium_fast",
#         Sim_cell_size = "full",
#         Read_state = "ave",
#         Write_state = "On",
#         Total_time = 1000,
#         Reference = "SR",
#         Results_Reference = "Run_$i",
#         # Other settings:
#         Spatial_output_interval_data = 0,  # Don't save `Spatial_Results` .data files
#         Spatial_output_interval_vtk = 0   # Don't save `Spatial_Results` .vtk files
#     )
# end

# # Get all the runs into a single DataFrame
# function single_3d_runs_df()
#     dirs = filter(readdir(project_path("Outputs_3Dcell_SR"), join=true)) do path
#         endswith(path, r"Results_Run_\d+") && !occursin("Spatial", path)
#     end
#     files = joinpath.(dirs, Ref("CRU.txt"))
#     dfs = [CSV.read(file, DataFrame, header=false) for file in files]
#     rename!.(dfs, Ref(cru_cols))
#     df = vcat(dfs...; cols=:union, source="run")
#     select!(df, "run", All())
# end

end
