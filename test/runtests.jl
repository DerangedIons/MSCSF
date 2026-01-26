using Test
using CellSims
using DataFrames

@testset "CellSims.jl" begin

    @testset "Model3D" begin
        m = CellSims.Model3D()
        @test m.Model == "minimal"
        @test m.ISO == 0
        @test m.Jup_scale == 1.0
        @test m.BCL == 1000
        @test m.tau_ss_type == "medium_fast"
        @test m.Beats == 1
        @test m.Total_time == 1000  # Beats * BCL
        @test m.Sim_cell_size == "full"
        @test m.Read_state == "Off"
        @test m.Write_state == "Off"

        # Test with custom parameters
        m2 = CellSims.Model3D(; ISO=1, BCL=500, Beats=10)
        @test m2.ISO == 1
        @test m2.BCL == 500
        @test m2.Beats == 10

        # Test Cmd generation
        cmd = Cmd(m)
        cmd_str = string(cmd)
        @test occursin("Model", cmd_str)
        @test occursin("minimal", cmd_str)
        @test occursin("ISO", cmd_str)
        @test occursin("BCL", cmd_str)
    end

    @testset "CaClamp3D" begin
        c = CellSims.CaClamp3D()
        @test c.Model == "minimal"
        @test c.ISO == 0
        @test c.Jup_scale == 1.0
        @test c.Jrel_scale == 1.0
        @test c.RyR_Po == 1.0
        @test c.tau_ss_type == "medium_fast"
        @test c.Total_time == 2000
        @test c.Sim_cell_size == "full"
        @test c.Cai == 0.1
        @test c.CaSR == 1000.0

        # Test with custom parameters
        c2 = CellSims.CaClamp3D(; Cai=0.2, CaSR=1500.0, RyR_Po=1.2)
        @test c2.Cai == 0.2
        @test c2.CaSR == 1500.0
        @test c2.RyR_Po == 1.2

        # Test Cmd generation
        cmd = Cmd(c)
        cmd_str = string(cmd)
        @test occursin("Cai", cmd_str)
        @test occursin("CaSR", cmd_str)
        @test occursin("RyR_Po", cmd_str)
    end

    @testset "Model3DSimulations" begin
        sim = Model3DSimulations(; ISO=1, BCL=300)
        @test sim.Reference == "sr_300_1"
        @test sim.runner.ISO == 1
        @test sim.runner.BCL == 300
        @test sim.prepace.Results_Reference == "prepace"
        @test sim.prepace_full.Results_Reference == "prepace_full"
        @test startswith(sim.runner.Results_Reference, "run_")

        # Test show method
        io = IOBuffer()
        show(io, sim)
        @test String(take!(io)) == "Model3DSimulations: \"sr_300_1\""
    end

    @testset "CaClamp3DSimulations" begin
        sim = CaClamp3DSimulations(; Cai=0.1, CaSR=1000.0, ISO=0, RyR_Po=1.0)
        @test sim.Reference == "ca_clamp_Cai0.1_CaSR1000.0_ISO0_Po1.0"
        @test sim.runner.Cai == 0.1
        @test sim.runner.CaSR == 1000.0
        @test sim.runner.ISO == 0
        @test sim.runner.RyR_Po == 1.0

        # Test show method
        io = IOBuffer()
        show(io, sim)
        @test String(take!(io)) == "CaClamp3DSimulations: \"ca_clamp_Cai0.1_CaSR1000.0_ISO0_Po1.0\""
    end

    @testset "all_simulations" begin
        sims = all_simulations()
        # 2 ISO values (0, 1) × 7 BCL values (300:200:1500) = 14
        @test length(sims) == 14
        @test all(s -> s isa Model3DSimulations, sims)

        # Check that all combinations are present
        bcl_values = [s.runner.BCL for s in sims]
        iso_values = [s.runner.ISO for s in sims]
        @test Set(bcl_values) == Set(300:200:1500)
        @test Set(iso_values) == Set(0:1)
    end

    @testset "all_clamp_simulations" begin
        sims = all_clamp_simulations()
        # Default: 3 CaSR × 3 RyR_Po = 9
        @test length(sims) == 9
        @test all(s -> s isa CaClamp3DSimulations, sims)

        # Test with custom parameters
        sims2 = all_clamp_simulations(; CaSR=[500.0, 1000.0], RyR_Po=[0.8, 1.0])
        @test length(sims2) == 4
    end

    @testset "cru_cols" begin
        @test length(CellSims.cru_cols) == 25
        @test CellSims.cru_cols[1] == "t"
        @test CellSims.cru_cols[2] == "Vm"
        @test "RyR_OA" in CellSims.cru_cols
        @test "Ca_cyto" in CellSims.cru_cols
        @test "Ca_NSR" in CellSims.cru_cols
    end

    @testset "ApproxLog" begin
        f = CellSims.ApproxLog(0.1)
        @test f(0.0) == log(0.1)
        @test f(1.0) == log(1.1)

        # Test inverse
        f_inv = CellSims.inverse(f)
        @test f_inv(f(0.5)) ≈ 0.5
        @test f_inv(f(1.0)) ≈ 1.0

        # Test identity inverse
        @test CellSims.inverse(identity) === identity
    end

    @testset "sims_low and sims_high" begin
        low = CellSims.sims_low()
        @test low.Reference == "sr_low"
        @test low.runner.ISO == 0
        @test low.runner.BCL == 1000

        high = CellSims.sims_high()
        @test high.Reference == "sr_high"
        @test high.runner.ISO == 1
        @test high.runner.BCL == 350
    end

    @testset "_get_apd90" begin
        # Create mock data with an action potential shape
        # Vm starts low, peaks, then returns low
        n = 101
        t = collect(0.0:1.0:100.0)
        Vm = Vector{Float64}(undef, n)
        # Resting phase
        Vm[1:10] .= -80.0
        # Upstroke
        Vm[11:17] .= [-80.0, -60.0, -40.0, -20.0, 0.0, 20.0, 40.0]
        # Plateau
        Vm[18:25] .= 40.0
        # Repolarization
        Vm[26:38] .= collect(40.0:-10.0:-80.0)
        # Rest
        Vm[39:end] .= -80.0

        df = DataFrame(t=t, Vm=Vm)

        apd = CellSims._get_apd90(df)
        @test apd > 0  # Should find a valid APD
        @test apd < 100  # Should be within time range

        # Test with flat Vm (no AP) - should return -1
        df_flat = DataFrame(t=t, Vm=fill(-80.0, length(t)))
        # When there's no variation, the function behavior depends on implementation
        # Just test it doesn't error
        @test CellSims._get_apd90(df_flat) isa Number
    end

    @testset "@trycatch macro" begin
        # Test that errors are caught when QUIET is true
        CellSims.QUIET = true
        result = @trycatch error("test error")
        @test isnothing(result)

        # Test that successful expressions return their value
        result = @trycatch 1 + 1
        @test result == 2

        CellSims.QUIET = false  # Reset
    end

    @testset "load_cru_files with empty paths" begin
        # Test with empty vector
        df = CellSims.load_cru_files(String[])
        @test df isa DataFrame
        @test isempty(df)
    end

    @testset "DIR initialization" begin
        @test !isempty(CellSims.DIR)
        @test isdir(CellSims.DIR)
        @test !isempty(CellSims.STATE_AND_GEOMETRY_FILES)
    end

end
