# CaSR=1000 Reconstruction Comparison
#
# Compares Colman SRF and PCA reconstruction methods on 20 example SCR
# waveforms: 10 from the Ca-clamp model (CaSR=1000 µM) and 10 from the 3D
# cell model (ISO=1, BCL=300 ms).  Produces one plot per reconstruction
# method, each showing original (black) vs. reconstruction (red/blue).

using CellSims, DataFrames, StatsPlots, MultivariateStats, Statistics

#-----------------------------------------------------------------------------# Colman SRF reconstruction
function colman_srf(t, stats)
    k_from_tau(τ) = 0.16980607 * τ + 0.00254852

    if stats.λ < 300
        t_to_peak = stats.tp - stats.ti
        decay_time = stats.tf - stats.tp
        k1 = k_from_tau(t_to_peak)
        k2 = k_from_tau(decay_time)
        thalf_1 = stats.ti + t_to_peak / 2
        thalf_2 = stats.tp + decay_time / 2
        peak = stats.peak
        map(t) do ti
            F1 = 1 / (1 + exp(-(ti - thalf_1) / k1))
            F2 = 1 / (1 + exp((ti - thalf_2) / k2))
            peak * F1 * F2
        end
    else
        tau_spike_rise = 50.0
        tau_spike_decay = 35.0
        tau_plat = 50.0
        k1_spike = k_from_tau(tau_spike_rise)
        k2_spike = k_from_tau(tau_spike_decay)
        thalf_spike_1 = stats.tp - tau_spike_rise / 2
        thalf_spike_2 = stats.tp + tau_spike_decay / 2
        plat_amp = 31.09 * (0.01 * stats.λ)^(-7.39) + 0.034
        spike_amp = stats.peak - plat_amp
        k_plat = k_from_tau(tau_plat)
        thalf_p1 = stats.ti + tau_plat / 2
        thalf_p2 = stats.ti + stats.λ - tau_plat / 2
        map(t) do ti
            F1 = 1 / (1 + exp(-(ti - thalf_spike_1) / k1_spike))
            F2 = 1 / (1 + exp((ti - thalf_spike_2) / k2_spike))
            F3 = 1 / (1 + exp(-(ti - thalf_p1) / k_plat))
            F4 = 1 / (1 + exp((ti - thalf_p2) / k_plat))
            spike_amp * F1 * F2 + plat_amp * F3 * F4
        end
    end
end

#-----------------------------------------------------------------------------# PCA hold-one-out reconstruction (sqrt transform)
function pca_reconstruct(gdf2, run_id)
    run_ids = [first(sub.run) for sub in gdf2]
    held_out_idx = findfirst(==(run_id), run_ids)
    train_idx = setdiff(1:length(gdf2), held_out_idx)
    f = sqrt
    f_inv = abs2
    w_train = reduce(hcat, f.(gdf2[i].RyR_OA) for i in train_idx)
    model = fit(PCA, w_train; maxoutdim=25, pratio=0.999)
    waveform = f.(collect(gdf2[held_out_idx].RyR_OA))
    z = MultivariateStats.transform(model, waveform)
    pca_y = f_inv.(MultivariateStats.reconstruct(model, z))
    (; pca_y, model)
end

#-----------------------------------------------------------------------------# Example selection helpers
function smooth(y, w=20)
    [mean(y[max(1,i-w):min(end,i+w)]) for i in eachindex(y)]
end

function is_good_example(run_df, stats; quiet_threshold=0.1, peak_threshold=0.3)
    outside = filter(row -> row.t < stats.ti || row.t > stats.tf, run_df)
    if !isempty(outside) && maximum(outside.RyR_OA) >= quiet_threshold * stats.peak
        return false
    end
    active = filter(row -> stats.ti <= row.t <= stats.tf, run_df)
    y = smooth(active.RyR_OA)
    isempty(y) && return false
    cutoff = peak_threshold * stats.peak
    n_peaks = sum(y[i] > y[i-1] && y[i] > y[i+1] && y[i] > cutoff for i in 2:length(y)-1)
    n_peaks <= 1
end

#-----------------------------------------------------------------------------# Data loading
function load_scr_data(df)
    colman = CellSims.get_colman_stats(df)
    scr_runs = filter(row -> row.λ != -1, colman)
    isempty(scr_runs) && return nothing
    df_sr = filter(row -> row.run in scr_runs.run, df)
    gdf = groupby(df_sr, "run")
    n = nrow(first(gdf))
    gdf2 = filter(sub -> nrow(sub) == n, gdf)
    length(gdf2) < 3 && return nothing
    (; df, scr_runs, gdf2)
end

function make_example(d, stats, label)
    run_df = filter(row -> row.run == stats.run, d.df)
    pca_y, model = pca_reconstruct(d.gdf2, stats.run)
    (; t=run_df.t, y=run_df.RyR_OA, colman=colman_srf(run_df.t, stats),
     pca=pca_y, model, stats, label)
end

# Find up to N good examples, sorted by peak (highest first)
function find_n_examples(d, label_fn; N=10, quiet=0.1, peak=0.3)
    examples = NamedTuple[]
    used_runs = Set{Int}()
    sorted = sort(d.scr_runs, :peak, rev=true)
    for qt in [quiet, Inf]
        for row in eachrow(sorted)
            row.run in used_runs && continue
            run_df = filter(r -> r.run == row.run, d.df)
            if is_good_example(run_df, row; quiet_threshold=qt, peak_threshold=peak)
                push!(examples, make_example(d, row, label_fn(row)))
                push!(used_runs, row.run)
                length(examples) >= N && return examples
            end
        end
    end
    length(examples) < N && @warn "Only found $(length(examples)) of $N requested examples"
    examples
end

#-----------------------------------------------------------------------------# Collect 10 clamp + 10 3D cell examples

# Ca-clamp CaSR=1000
d_clamp = load_scr_data(get_df(only(filter(s -> s.runner.CaSR == 1000, all_clamp_simulations()))))
clamp_examples = find_n_examples(d_clamp, row -> "Clamp, λ=$(round(Int, row.λ))"; N=10)

# 3D cell model ISO=1 BCL=300
d_cell = load_scr_data(get_df(Model3DSimulations(; ISO=1, BCL=300)))
cell_examples = find_n_examples(d_cell, row -> "3D, λ=$(round(Int, row.λ))"; N=10)

all_examples = vcat(clamp_examples, cell_examples)
@info "Collected $(length(all_examples)) examples: $(length(clamp_examples)) clamp + $(length(cell_examples)) 3D cell"

#-----------------------------------------------------------------------------# Generate and save individual plots

dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "comparison_plots"))

for (i, e) in enumerate(all_examples)
    margin = 0.2 * e.stats.λ
    xl = (e.stats.ti - margin, e.stats.tf + margin)
    yl = (min(minimum(e.y), minimum(e.colman), minimum(e.pca)),
          max(maximum(e.y), maximum(e.colman), maximum(e.pca)) * 1.05)
    colman_title = e.stats.λ > 300 ? "Colman SRF + Plateau" : "Colman SRF"
    kw = (; lw=2, label="", ylabel="RyR OA", ylims=yl, xlims=xl)

    p1 = plot(e.t, e.y;       color=:black, title="Original: $(e.label)", xformatter=x -> "", kw...)
    p2 = plot(e.t, e.colman;  color=:red,   title=colman_title, xformatter=x -> "", kw...)
    p3 = plot(e.t, e.pca;     color=:blue,  title="PCA ($(outdim(e.model)) comp.)", xlabel="Time (ms)", kw...)

    p = plot(p1, p2, p3; layout=(3, 1), size=(700, 600))
    fname = "waveform_$(lpad(i, 2, '0')).png"
    savefig(p, joinpath(dir, fname))
end

@info "Saved $(length(all_examples)) individual waveform plots to $dir"
