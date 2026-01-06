using CellSims, Distributions, HypothesisTests, StatsPlots, DataFrames, MultivariateStats, LinearAlgebra


# PREREQUISITE:
# - You need to have compiled the model_single_3D code (good task for Claude Code).

#--------------------------------------------------------------------------# Create Simulation Runs
# See also: sims_low() and the lower-level Model3DSimulations struct
# sims_high akes about 6.5 minutes per run
sim = CellSims.sims_high()

# How many runs to simulate
n = 25

# !!! This will take a while depending on n !!!
@elapsed for i in 1:n
    run(sim)
end

plot(sim, runs=1:10)  # quick look at waveforms.  Remove `runs=...` to plot all

# Opens Reference directory in case you want to check what directories/files were created
open_reference(h)

#-----------------------------------------------------------------------------# Analyze Results
# DataFrame with columns: run, t, Vm, ..., Ca_cyto, ..., RyR_OA, RyR_OI
df = get_df(sim)

# Statistics DataFrame with columns: run, ti, tf, λ, tp, peak, plat
# If no SR: ti, tf, λ = -1
df_stats = stats(df)
stats_dict = Dict(row.run => (; ti=row.ti, tf=row.tf, λ=row.λ, tp=row.tp, peak=row.peak, plat=row.plat) for row in eachrow(df_stats))
df_stats_sr = filter(row -> row.ti != -1, df_stats)

# Estimate of Prob(SR)
pscr = nrow(df_stats_sr) / nrow(df_stats)

# Get df with only runs that had SR
sr_runs = df_stats.run[df_stats.ti .!= -1]
df_sr = filter(x -> x.run in sr_runs, df)

# Get only the SR part from df_sr (drop timesteps before ti, after tf)
df_sr_filtered = filter(df_sr) do row
    (; ti, tf) = stats_dict[row.run]
    ti == -1 ? false : ti ≤ row.t ≤ tf
end

plot(df_stats.ti, df_stats.tf, seriestype=:scatter, xlabel="ti", ylabel="tf", title="SR start and end times")


#-----------------------------------------------------------------------------# Fit Distributions
# For fitting distributions for Colman's parameterization.  Not necessary if we go the generative PCA route.
# # Fit Distribution to e.g. ti
# x = df_stats_sr.ti

# dist = fit(Gamma, x)

# # Evaluate the fit (plot)
# plt = histogram(x, nbin=25, normalize=true, lab="")
# plot!(plt, x -> pdf(dist, x), link=:x, label="Gamma fit")
# plot(plt, qqplot(dist, x, title="Q-Q Plot"))

# # Evaluate fit: Kolmogorov-Smirnov Test
# @info pvalue(ExactOneSampleKSTest(x, dist), tail = :right)

# # If you decide this distribution is a good fit, then randomly generate values via:
# rand(dist)

#-----------------------------------------------------------------------------# Generative PCA
# Get log(waveforms) as matrix
ϵ = 0.1  # small offset to avoid log(0)
_log(x) = log(x + ϵ)
_log_inverse(x) = exp(x) - ϵ



w = reduce(hcat, sub.RyR_OA for sub in groupby(df_sr, "run"))
w = reduce(hcat, _log.(sub.RyR_OA) for sub in groupby(df_sr, "run"))

# Max number of principal components to keep
k = 25

# Fit PCA model
model = fit(PCA, w; maxoutdim=k, pratio=0.999);
@info "N components: $(outdim(model))"

z = MultivariateStats.transform(model, w)

μ = vec(mean(z, dims=2))
Σ = cov(z; dims=2)

dist2 = MvNormal(μ, Symmetric(Σ))

rand_waveform(pca::PCA=model, dist::MvNormal = dist2) = max.(df_sr.RyR_OA[1], reconstruct(pca, rand(dist)))
rand_waveform(pca::PCA=model, dist::MvNormal = dist2) = max.(df_sr.RyR_OA[1], _log_inverse.(reconstruct(pca, rand(dist))))

# Test out plotting random waveforms
n_sr_runs = length(unique(df_sr.run))

waveforms = reduce(hcat, rand_waveform() for _ in 1:n_sr_runs)
plot(waveforms, lab="")


# Compare generated waveforms to actual
plot(
    plot(h, title="Actual"),
    plot(waveforms, lab="", title="Generated"),
    link = :all, linewidth=0.5
)

# Closer look at generated waveforms
plot(
    plot(sim, runs=1:12, layout=12, title="Sim"),
    plot(waveforms[:, 1:12], lab="", title="Gen", layout=12),
    xlab="", ylab="", link=:all, linewidth=0.5, ticks=false, xlim=(400, Inf)
)
