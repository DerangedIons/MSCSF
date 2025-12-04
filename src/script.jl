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

#-------------------------------------------------------------------# New Waveform Parameterization
