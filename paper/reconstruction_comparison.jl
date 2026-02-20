using CellSims, DataFrames, StatsPlots, MultivariateStats, Statistics

#-----------------------------------------------------------------------------# Colman SRF reconstruction
# Product-of-sigmoids model from Colman (2019) Spontaneous Release Functions
#   NRyRo(t) = NRyRo_peak * F1(t) * F2(t)    [spike]
#            + NRyRo_plat * F3(t) * F4(t)     [plateau, if λ > 300]
# where Fn(t) = 1/(1 + exp(±(t - thalf)/k))   [sigmoids]
# k = 0.16980607 * τ + 0.00254852             [slope from time constant]
function colman_srf(t, stats)
    k_from_tau(τ) = 0.16980607 * τ + 0.00254852

    if stats.λ < 300
        # Short event: spike spans full duration
        t_to_peak = stats.tp - stats.ti
        decay_time = stats.tf - stats.tp
        k1 = k_from_tau(t_to_peak)
        k2 = k_from_tau(decay_time)
        thalf_1 = stats.ti + t_to_peak / 2   # midpoint of rise
        thalf_2 = stats.tp + decay_time / 2   # midpoint of decay
        peak = stats.peak

        map(t) do ti
            F1 = 1 / (1 + exp(-(ti - thalf_1) / k1))
            F2 = 1 / (1 + exp((ti - thalf_2) / k2))
            peak * F1 * F2
        end
    else
        # Long event: plateau spanning full duration + embedded spike (fixed 50ms/35ms)
        tau_spike_rise = 50.0
        tau_spike_decay = 35.0
        tau_plat = 50.0

        # Spike: fixed timing, centered on measured peak time
        k1_spike = k_from_tau(tau_spike_rise)
        k2_spike = k_from_tau(tau_spike_decay)
        thalf_spike_1 = stats.tp - tau_spike_rise / 2
        thalf_spike_2 = stats.tp + tau_spike_decay / 2

        # Plateau amplitude from Colman power law (not measured median)
        plat_amp = 31.09 * (0.01 * stats.λ)^(-7.39) + 0.034
        spike_amp = stats.peak - plat_amp

        # Plateau: spans from ti to ti + λ with tau_plat rise/decay
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
    f_inv = abs2  # x² (inverse of √x); abs2 guarantees non-negative

    w_train = reduce(hcat, f.(gdf2[i].RyR_OA) for i in train_idx)
    model = fit(PCA, w_train; maxoutdim=25, pratio=0.999)

    waveform = f.(collect(gdf2[held_out_idx].RyR_OA))
    z = MultivariateStats.transform(model, waveform)
    pca_y = f_inv.(MultivariateStats.reconstruct(model, z))
    (; pca_y, model)
end

#-----------------------------------------------------------------------------# Check if a waveform has a single isolated, unimodal SCR event
function smooth(y, w=20)
    [mean(y[max(1,i-w):min(end,i+w)]) for i in eachindex(y)]
end

function is_good_example(run_df, stats; quiet_threshold=0.1, peak_threshold=0.3)
    # Check that waveform outside [ti, tf] is quiet relative to peak
    outside = filter(row -> row.t < stats.ti || row.t > stats.tf, run_df)
    if !isempty(outside) && maximum(outside.RyR_OA) >= quiet_threshold * stats.peak
        return false
    end
    # Check unimodal: smooth then count peaks above threshold
    active = filter(row -> stats.ti <= row.t <= stats.tf, run_df)
    y = smooth(active.RyR_OA)
    isempty(y) && return false
    cutoff = peak_threshold * stats.peak
    n_peaks = sum(y[i] > y[i-1] && y[i] > y[i+1] && y[i] > cutoff for i in 2:length(y)-1)
    n_peaks <= 1
end

#-----------------------------------------------------------------------------# Helper: load data and prepare PCA groups for a given CaSR
function load_condition(casr)
    sims = all_clamp_simulations()
    sim = only(filter(s -> s.runner.CaSR == casr, sims))
    df = get_df(sim)
    colman = CellSims.get_colman_stats(df)
    scr_runs = filter(row -> row.λ != -1, colman)
    df_sr = filter(row -> row.run in scr_runs.run, df)
    gdf = groupby(df_sr, "run")
    n = nrow(first(gdf))
    gdf2 = filter(sub -> nrow(sub) == n, gdf)
    (; df, scr_runs, gdf2)
end

function find_good_example(candidates, df; quiet=0.1, peak=0.3)
    for row in eachrow(sort(candidates, :peak, rev=true))
        run_df = filter(r -> r.run == row.run, df)
        is_good_example(run_df, row; quiet_threshold=quiet, peak_threshold=peak) && return row
    end
    error("No good example waveform found")
end

#-----------------------------------------------------------------------------# Short event (λ ≤ 300) from CaSR=1000
d_s = load_condition(1000)
short_stats = find_good_example(filter(row -> row.λ <= 300, d_s.scr_runs), d_s.df; peak=0.5)
short_id = short_stats.run
short_df = filter(row -> row.run == short_id, d_s.df)

t_s = short_df.t
y_s = short_df.RyR_OA
colman_s = colman_srf(t_s, short_stats)
pca_s, model_s = pca_reconstruct(d_s.gdf2, short_id)

#-----------------------------------------------------------------------------# Long event (λ > 300) from CaSR=1000
d_l = load_condition(1000)
long_stats = find_good_example(filter(row -> row.λ > 300, d_l.scr_runs), d_l.df; peak=0.5)
long_id = long_stats.run
long_df = filter(row -> row.run == long_id, d_l.df)

t_l = long_df.t
y_l = long_df.RyR_OA
colman_l = colman_srf(t_l, long_stats)
pca_l, model_l = pca_reconstruct(d_l.gdf2, long_id)

#-----------------------------------------------------------------------------# Plot (3×2 grid)
yl_s = (min(minimum(y_s), minimum(colman_s), minimum(pca_s)),
        max(maximum(y_s), maximum(colman_s), maximum(pca_s)) * 1.05)
yl_l = (min(minimum(y_l), minimum(colman_l), minimum(pca_l)),
        max(maximum(y_l), maximum(colman_l), maximum(pca_l)) * 1.05)

kw = (; lw=2, label="", xformatter=x -> "")

p1 = plot(t_s, y_s; color=:black, ylabel="RyR OA", title="Original (λ=$(round(Int, short_stats.λ)) ms)", ylims=yl_s, kw...)
p2 = plot(t_l, y_l; color=:black, title="Original (λ=$(round(Int, long_stats.λ)) ms)", ylims=yl_l, kw...)

p3 = plot(t_s, colman_s; color=:red, ylabel="RyR OA", title="Colman SRF", ylims=yl_s, kw...)
p4 = plot(t_l, colman_l; color=:red, title="Colman SRF + Plateau", ylims=yl_l, kw...)

p5 = plot(t_s, pca_s; color=:blue, ylabel="RyR OA", title="PCA ($(outdim(model_s)) comp.)", ylims=yl_s, lw=2, label="", xlabel="Time (ms)")
p6 = plot(t_l, pca_l; color=:blue, title="PCA ($(outdim(model_l)) comp.)", ylims=yl_l, lw=2, label="", xlabel="Time (ms)")

p = plot(p1, p2, p3, p4, p5, p6; layout=(3, 2), size=(1000, 750),
    plot_title="Waveform Reconstruction Comparison (CaSR = 1000 µM)")

dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "ca_clamp_summary"))
savefig(p, joinpath(dir, "reconstruction_comparison.png"))
@info "Saved reconstruction_comparison.png"

p
