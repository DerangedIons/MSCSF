using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using CSV, DataFrames

DIR = mkpath(abspath(joinpath(@__DIR__, "..", "data")))
write(joinpath(DIR, "PATH.txt"), mkpath(joinpath(DIR, "state_and_geometry_files")))
mkpath(joinpath(DIR, "state_and_geometry_files", "State_files", "Single_cell"))

cp(joinpath(@__DIR__, "..", "CODE", "model_single_3D"), joinpath(DIR, "model_single_3D"); force=true)

# Number of runs to perform; default is 10 if not provided as command line argument
N = length(ARGS) == 1 ? parse(Int, ARGS[1]) : 10

#-----------------------------------------------------------------------------# common_args
common_args = (
    Model = "minimal",
    ISO = 1,
    Jup_scale = 2,
    tau_ss_type = "medium_fast",
    BCL = 350,
    Reference = "SR",
    Spatial_output_interval_data = 0,
    Spatial_output_interval_vtk = 0,
)

#-----------------------------------------------------------------------------# prepace
prepace_args = (
    Beats = 50,
    Sim_cell_size = "testing",
    Write_state = "ave",
    Results_Reference = "prepace",
    Total_time = 50 * 350,
)

#-----------------------------------------------------------------------------# prepace_full
prepace_full_args = (
    Beats = 4,
    Sim_cell_size = "full",
    Read_state = "ave",
    Write_state = "On",
    Results_Reference = "prepace_full",
    Total_time = 4 * 350,
)

#-----------------------------------------------------------------------------# run_args
run_args = (
    Beats = 1,
    Sim_cell_size = "full",
    Read_state = "On",
    total_time = 1000
)


#-----------------------------------------------------------------------------# functions
function make_cmd(args)
    out = ["./model_single_3D"]
    for (k, v) in pairs(merge(common_args, args))
        push!(out, string(k), string(v))
    end
    return Cmd(out)
end

function run_model(args)
    cd(joinpath(@__DIR__, "..", "data")) do
        @time run(make_cmd(args))
    end
end

#-----------------------------------------------------------------------------# runs
if isdir(joinpath(DIR, "Outputs_3Dcell_SR", "Results_prepace"))
    @info "`prepace` already done, skipping..."
else
    run_model(prepace_args)  # ~2.5 minutes for 50 beats at 250ms BCL
end

if isdir(joinpath(DIR, "Outputs_3Dcell_SR", "Results_prepace_full"))
    @info "`prepace_full` already done, skipping..."
else
    run_model(prepace_full_args)  # ~7.5 minutes for 4 beats at 350ms BCL
end

for i in 1:N
    n = length(filter(x -> startswith(x, r"Results_run"), readdir(joinpath(@__DIR__, "..", "data", "Outputs_3Dcell_SR")))) + 1
    @info "Run: $n"
    run_model(merge(run_args, (; Results_Reference = "run_$(lpad(n, 3, '0'))")))
end
