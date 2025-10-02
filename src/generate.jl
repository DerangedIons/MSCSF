using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using CSV, DataFrames

DIR = mkpath(abspath(joinpath(@__DIR__, "..", "data")))
write(joinpath(DIR, "PATH.txt"), mkpath(joinpath(DIR, "state_and_geometry_files")))
mkpath(joinpath(DIR, "state_and_geometry_files", "State_files", "Single_cell"))

cp(joinpath(@__DIR__, "..", "CODE", "model_single_3D"), joinpath(DIR, "model_single_3D"); force=true)

# Number of runs to perform; default is 10 if not provided as command line argument
N = length(ARGS) == 1 ? parse(Int, ARGS[1]) : 10

#-----------------------------------------------------------------------------# notes
# Low adrenergic / Low SR: ISO 0, BCL 1000, Jup_scale 1.0, tau_ss_type medium_low
# High adrenergic / High SR: ISO 1, BCL 350, Jup_scale 2.0, tau_ss_type medium_fast,


#-----------------------------------------------------------------------------# low
low = (
    common_args = (
        Model = "minimal",
        ISO = 0,
        Jup_scale = 1,
        tau_ss_type = "medium_low",
        BCL = 1000,
        Reference = "sr_low",
        Spatial_output_interval_data = 0,
        Spatial_output_interval_vtk = 0,
    ),
    prepace_args = (
        Beats = 40,
        Sim_cell_size = "testing",
        Write_state = "ave",
        Results_Reference = "prepace",
        Total_time = 40 * 1000,
    ),
    prepace_full_args = (
        Beats = 4,
        Sim_cell_size = "full",
        Read_state = "ave",
        Write_state = "On",
        Results_Reference = "prepace_full",
        Total_time = 4 * 1000,
    ),
    run_args = (
        Beats = 1,
        Sim_cell_size = "full",
        Read_state = "On",
        total_time = 1100
    )
)

#-----------------------------------------------------------------------------# high
high = (
    common_args = (
        Model = "minimal",
        ISO = 1,
        Jup_scale = 1,
        tau_ss_type = "medium_fast",
        BCL = 350,
        Reference = "sr_high",
        Spatial_output_interval_data = 0,
        Spatial_output_interval_vtk = 0,
    ),
    prepace_args = (
        Beats = 40,
        Sim_cell_size = "testing",
        Write_state = "ave",
        Results_Reference = "prepace",
        Total_time = 40 * 350,
    ),
    prepace_full_args = (
        Beats = 4,
        Sim_cell_size = "full",
        Read_state = "ave",
        Write_state = "On",
        Results_Reference = "prepace_full",
        Total_time = 4 * 350,
    ),
    run_args = (
        Beats = 1,
        Sim_cell_size = "full",
        Read_state = "On",
        total_time = 1100
    )
)

#-----------------------------------------------------------------------------# init
# low
run_model(merge(low.common_args, low.prepace_args))
run_model(merge(low.common_args, low.prepace_full_args))
# high
run_model(merge(high.common_args, high.prepace_args))
run_model(merge(high.common_args, high.prepace_full_args))

#-----------------------------------------------------------------------------# runs
for _ in 1:N
    # Make sure we don't overwrite existing runs
    n = length(filter(x -> startswith(x, r"Results_run"), readdir(joinpath(@__DIR__, "..", "data", "Outputs_3Dcell_sr_low")))) + 1
    args = merge(low.common_args, low.run_args, (; Results_Reference = "run_$(lpad(n, 3, '0'))"))
    @info "Running low SR Run $n" args
    run_model(args)

    n = length(filter(x -> startswith(x, r"Results_run"), readdir(joinpath(@__DIR__, "..", "data", "Outputs_3Dcell_sr_high")))) + 1
    args = merge(high.common_args, high.run_args, (; Results_Reference = "run_$(lpad(n, 3, '0'))"))
    @info "Running high SR Run $n" args
    run_model(args)
end


#-----------------------------------------------------------------------------# functions
function make_cmd(args)
    out = ["./model_single_3D"]
    for (k, v) in pairs(merge(common_args, args))
        push!(out, string(k), string(v))
    end
    return Cmd(out)
end

function run_model(args)
    results_dir = joinpath(DIR, "Outputs_3Dcell_$(args.Reference)")
    if isdir(results_dir)
        @info "`$(args.Reference)` results directory already exists, skipping..."
    else
        cd(joinpath(@__DIR__, "..", "data")) do
            @time run(make_cmd(args))
        end
    end
end
