# Reconstruction Comparison
#
# Compares two methods for reconstructing RyR open-probability waveforms:
#
#   1. Colman SRF — A product-of-sigmoids parametric model (Colman 2019).
#      Short events (λ ≤ 300 ms) use a single spike (two sigmoids for rise/decay).
#      Long events (λ > 300 ms) add a plateau component (four sigmoids total).
#      Parameters are derived from summary statistics (ti, tp, tf, peak, λ).
#
#   2. PCA — Hold-one-out principal component analysis on sqrt-transformed
#      waveforms.  A PCA basis is fit on all SCR runs except the target, then
#      the held-out waveform is projected and reconstructed.
#
# Examples are drawn from two model types:
#   - Ca-clamp simulations (CaSR = 900, 1000, 1100, 1200 µM; short + long events)
#   - 3D cell model simulations (ISO=1, varying BCL)
#
# Output: two separate plots (one per method), each showing original vs.
# reconstruction across all examples.

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

function find_good_example(candidates, df; quiet=0.1, peak=0.3)
    sorted = sort(candidates, :peak, rev=true)
    # Try strict (quiet + unimodal), then fall back to unimodal-only
    for row in eachrow(sorted)
        run_df = filter(r -> r.run == row.run, df)
        is_good_example(run_df, row; quiet_threshold=quiet, peak_threshold=peak) && return row
    end
    for row in eachrow(sorted)
        run_df = filter(r -> r.run == row.run, df)
        is_good_example(run_df, row; quiet_threshold=Inf, peak_threshold=peak) && return row
    end
    nothing
end

#-----------------------------------------------------------------------------# Data loading helpers
function load_scr_data(df)
    colman = CellSims.get_colman_stats(df)
    scr_runs = filter(row -> row.λ != -1, colman)
    isempty(scr_runs) && return nothing
    df_sr = filter(row -> row.run in scr_runs.run, df)
    gdf = groupby(df_sr, "run")
    n = nrow(first(gdf))
    gdf2 = filter(sub -> nrow(sub) == n, gdf)
    length(gdf2) < 3 && return nothing  # need enough runs for PCA
    (; df, scr_runs, gdf2)
end

function load_condition(casr)
    sims = all_clamp_simulations()
    sim = only(filter(s -> s.runner.CaSR == casr, sims))
    load_scr_data(get_df(sim))
end

function load_condition_3d(sim)
    load_scr_data(get_df(sim))
end

#-----------------------------------------------------------------------------# Example extraction

# Make a single example named tuple from a run
function make_example(d, stats, label)
    run_df = filter(row -> row.run == stats.run, d.df)
    pca_y, model = pca_reconstruct(d.gdf2, stats.run)
    (; t=run_df.t, y=run_df.RyR_OA, colman=colman_srf(run_df.t, stats),
     pca=pca_y, model, stats, label)
end

# Extract short + long examples from a Ca-clamp condition
function extract_clamp_examples(casr; quiet=0.1, peak=0.3)
    d = load_condition(casr)
    isnothing(d) && return NamedTuple[]
    examples = NamedTuple[]

    short_cands = filter(row -> row.λ <= 300, d.scr_runs)
    if !isempty(short_cands)
        ss = find_good_example(short_cands, d.df; quiet, peak)
        !isnothing(ss) && push!(examples, make_example(d, ss, "Clamp $(casr), λ=$(round(Int, ss.λ))"))
    end

    long_cands = filter(row -> row.λ > 300, d.scr_runs)
    if !isempty(long_cands)
        ls = find_good_example(long_cands, d.df; quiet, peak)
        !isnothing(ls) && push!(examples, make_example(d, ls, "Clamp $(casr), λ=$(round(Int, ls.λ))"))
    end

    examples
end

# Extract up to N examples from a 3D cell model simulation
function extract_3d_examples(sim; N=1, quiet=0.1, peak=0.3)
    d = load_condition_3d(sim)
    isnothing(d) && return NamedTuple[]

    examples = NamedTuple[]
    used_runs = Set{Int}()
    sorted = sort(d.scr_runs, :peak, rev=true)
    iso, bcl = sim.runner.ISO, sim.runner.BCL

    # Strict pass, then relaxed pass
    for qt in [quiet, Inf]
        for row in eachrow(sorted)
            row.run in used_runs && continue
            run_df = filter(r -> r.run == row.run, d.df)
            if is_good_example(run_df, row; quiet_threshold=qt, peak_threshold=peak)
                push!(examples, make_example(d, row, "3D ISO=$(iso) BCL=$(bcl), λ=$(round(Int, row.λ))"))
                push!(used_runs, row.run)
                length(examples) >= N && return examples
            end
        end
    end

    length(examples) < N && @warn "Only found $(length(examples)) examples for ISO=$iso BCL=$bcl"
    examples
end

#-----------------------------------------------------------------------------# Collect all examples

# Ca-clamp: 4 CaSR values × short + long
clamp_examples = reduce(vcat, extract_clamp_examples(casr) for casr in [900, 1000, 1100, 1200])

# 3D cell model: ISO=1 with 4 BCL values
sims_3d = filter(s -> s.runner.ISO == 1 && s.runner.BCL in [300, 500, 700, 1100], all_simulations())
cell_examples = reduce(vcat, extract_3d_examples(sim) for sim in sims_3d)

all_examples = vcat(clamp_examples, cell_examples)
@info "Collected $(length(all_examples)) examples: $(length(clamp_examples)) clamp + $(length(cell_examples)) 3D cell"

#-----------------------------------------------------------------------------# Plotting

function ylims_for(e)
    lo = min(minimum(e.y), minimum(e.colman), minimum(e.pca))
    hi = max(maximum(e.y), maximum(e.colman), maximum(e.pca)) * 1.05
    (lo, hi)
end

function make_reconstruction_plot(examples, method::Symbol; ncols=4)
    N = length(examples)
    nrow_pairs = ceil(Int, N / ncols)
    nrows = 2 * nrow_pairs

    recon_color = method == :colman ? :red : :blue
    method_name = method == :colman ? "Colman SRF" : "PCA"

    subplots = []

    for pair in 1:nrow_pairs
        # Original row
        for col in 1:ncols
            idx = (pair - 1) * ncols + col
            if idx > N
                push!(subplots, plot(; framestyle=:none, label=""))
                continue
            end
            e = examples[idx]
            yl = ylims_for(e)
            is_first = col == 1
            push!(subplots, plot(e.t, e.y;
                color=:black, lw=2, label="",
                title=e.label, titlefontsize=8,
                ylabel=is_first ? "RyR OA" : "",
                xformatter=x -> "", ylims=yl,
            ))
        end
        # Reconstruction row
        for col in 1:ncols
            idx = (pair - 1) * ncols + col
            if idx > N
                push!(subplots, plot(; framestyle=:none, label=""))
                continue
            end
            e = examples[idx]
            yl = ylims_for(e)
            recon_data = method == :colman ? e.colman : e.pca
            is_first = col == 1
            is_last_pair = pair == nrow_pairs
            subtitle = method == :colman ?
                (e.stats.λ > 300 ? "SRF + Plateau" : "SRF") :
                "PCA ($(outdim(e.model)) comp.)"
            push!(subplots, plot(e.t, recon_data;
                color=recon_color, lw=2, label="",
                title=subtitle, titlefontsize=8,
                ylabel=is_first ? "RyR OA" : "",
                xlabel=is_last_pair ? "Time (ms)" : "",
                xformatter=is_last_pair ? :auto : (x -> ""),
                ylims=yl,
            ))
        end
    end

    plot(subplots...; layout=(nrows, ncols),
        size=(ncols * 400, nrows * 200),
        plot_title="$method_name Reconstruction Comparison")
end

#-----------------------------------------------------------------------------# Generate and save plots
dir = mkpath(joinpath(@__DIR__, "..", "data", "results", "ca_clamp_summary"))

p_colman = make_reconstruction_plot(all_examples, :colman)
savefig(p_colman, joinpath(dir, "reconstruction_colman.png"))
@info "Saved reconstruction_colman.png"

p_pca = make_reconstruction_plot(all_examples, :pca)
savefig(p_pca, joinpath(dir, "reconstruction_pca.png"))
@info "Saved reconstruction_pca.png"
