using MSCSF, Distributions, HypothesisTests, StatsPlots


# PREREQUISITE:
# - You need to have compiled the model_single_3D code (good task for Claude Code).

#--------------------------------------------------------------------------# Create Simulation Runs
 # See also: sims_low() and the lower-level Model3DSimulations struct
h = sims_high()

# How many runs to simulate
n = 100

# !!! This will take a good long while !!!
for i in 1:n
    run(h)
end

plot(h, runs=[1,2,3])  # quick look at waveforms.  Remove `runs=...` to plot all

# Opens Reference directory in case you want to check what directories/files were created
open_reference(h)

#-----------------------------------------------------------------------------# Analyze Results
# DataFrame with columns: run, t, Vm, ..., Ca_cyto, ..., RyR_OA, RyR_OI
df = get_df(h)

# Statistics DataFrame with columns: run, ti, tf, λ, tp, peak, plat
# If no SR: ti, tf, λ = -1
df_stats = stats(df)
df_stats_sr = filter(row -> row.ti != -1, df_stats)

# Estimate of Prob(SR)
pscr = nrow(df_stats_sr) / nrow(df_stats)

# Get df with non-SR runs removed
sr_runs = df_stats.run[df_stats.ti .!= -1]
df_sr = filter(x -> x.run in sr_runs, df)


#-----------------------------------------------------------------------------# Fit Distributions
# Fit Distribution to e.g. ti
x = df_stats_sr.ti

dist = fit(Gamma, x)

# Evaluate the fit (plot)
plt = histogram(x, nbin=25, normalize=true, lab="")
plot!(plt, x -> pdf(dist, x), link=:x, label="Gamma fit")
plot(plt, qqplot(dist, x, title="Q-Q Plot"))

# Evaluate fit: Kolmogorov-Smirnov Test
@info pvalue(ExactOneSampleKSTest(x, dist), tail = :right)

# If you decide this distribution is a good fit, then randomly generate values via:
rand(dist)









# EXPERIMENTAL:

#-------------------------------------------# Idea: PCA + Multivariate Normal to generate waveforms
using MultivariateStats, LinearAlgebra

# Get log(waveforms) as matrix
ϵ = 0.1  # small offset to avoid log(0)

w = reduce(hcat, log.(sub.RyR_OA .+ ϵ) for sub in groupby(df_sr, "run"))

# Max number of principal components to keep
k = 25

# Fit PCA model
model = fit(PCA, w; maxoutdim=k);
@info "N components: $(outdim(model))"

z = MultivariateStats.transform(model, w)

μ = vec(mean(z, dims=2))
Σ = cov(z; dims=2)

dist2 = MvNormal(μ, Symmetric(Σ))

rand_waveform(pca::PCA=model, dist::MvNormal = dist2) = max.(df_sr.RyR_OA[1], exp.(reconstruct(pca, rand(dist))) .- ϵ)

# Test out plotting random waveforms
waveforms = reduce(hcat, rand_waveform() for _ in 1:20)
plot(waveforms, lab="")


# Compare generated waveforms to actual
plot(
    plot(h, title="Actual"),
    plot(waveforms, lab="", title="Generated"),
    link = :all, linewidth=0.5
)
