# Generative PCA for Synthetic Calcium Release Waveforms in CellSims

## 1. Introduction

CellSims is a Julia framework for running and analyzing stochastic 3D cardiac cell
simulations. A core capability is its **generative PCA** model, which learns the
statistical structure of simulated spontaneous calcium release (SCR) waveforms and
samples new, synthetic waveforms that are statistically consistent with the training
data but were never directly produced by the simulator.

The target variable is **RyR_OA** (the fraction of ryanodine receptors in the
open-activated state), a time series recorded at each calcium release unit (CRU) during
a simulation beat. RyR_OA traces capture the onset, amplitude, and duration of
spontaneous calcium release events---key determinants of arrhythmogenic risk.

## 2. Simulation Pipeline

Each experiment is managed by a `Model3DSimulations` object that orchestrates three
stages:

1. **Prepace** -- A reduced-cell-size run (e.g. 50 beats) that equilibrates the model
   and saves an averaged intracellular state.
2. **Prepace Full** -- A full-cell-size continuation (e.g. 4 beats) that reads the
   averaged state and writes a final state file.
3. **Runner** -- Multiple independent full-cell-size runs (each 1 beat) that read the
   final state and record time series to `CRU.txt` files. Each run is an independent
   realization of the stochastic calcium dynamics.

The `CRU.txt` output contains 25 columns per time step (membrane voltage, calcium
concentrations in five compartments, ion fluxes, RyR and LTCC channel states, etc.).
After loading all runner outputs into a single DataFrame, each run is identified by a
`run` column.

## 3. Preprocessing

Before fitting the generative model, the raw data is filtered:

1. **Completeness filter.** Runs whose CRU output has a different number of time
   points than the first run (e.g. from an interrupted simulation) are discarded.
2. **SCR detection.** For each run, an APD90-based cutoff is computed on the membrane
   voltage to identify the post-repolarization window. Within that window, an
   amplitude-based threshold (5% of the RyR_OA range) determines whether an SCR event
   occurred. Runs without a detected SCR ($t_i = -1$) are excluded from the generative
   model.

The remaining "SR-only" runs form the training set of $N$ RyR_OA waveforms, each of
length $T$ (the number of time points per run).

## 4. Generative PCA Model

### 4.1 Overview

Standard PCA is a tool for dimensionality reduction: it finds an orthogonal basis that
maximizes explained variance and can reconstruct the original data from a reduced set of
coefficients. Generative PCA goes further by **modeling the distribution of those
coefficients** and **sampling from that distribution** to produce novel waveforms.

The pipeline is:

```
Training waveforms ──> (optional transform) ──> PCA fit ──> latent scores
                                                                  │
                                                          fit 𝒩(μ, Σ)
                                                                  │
New waveform <── (inverse transform + clamp) <── reconstruct <── sample z
```

### 4.2 Optional Nonlinear Transform

Because RyR_OA values are non-negative and may span several orders of magnitude, the
framework supports an optional variance-stabilizing transform $f$ applied element-wise
before PCA. Two choices are provided:

**Identity transform:**

$$f(x) = x, \qquad f^{-1}(x) = x$$

**Approximate log transform** with offset $\varepsilon > 0$ (default $\varepsilon = 0.1$):

$$f(x) = \log(x + \varepsilon), \qquad f^{-1}(y) = e^{y} - \varepsilon$$

The offset $\varepsilon$ prevents $\log(0)$. Working in log-space stabilizes variance
across the dynamic range of RyR_OA and can make the distribution of PC scores closer to
Gaussian.

### 4.3 PCA Decomposition

Let $\mathbf{w}_i \in \mathbb{R}^{T}$ denote the $i$-th (possibly transformed)
training waveform, for $i = 1, \dots, N$. Arrange these as columns of the data matrix

$$\mathbf{W} = \begin{bmatrix} \mathbf{w}_1 & \cdots & \mathbf{w}_N \end{bmatrix} \in \mathbb{R}^{T \times N}.$$

PCA is fit via `MultivariateStats.fit(PCA, W)` with two controls:

- `maxoutdim` $= k$ (default 25): upper bound on retained components.
- `pratio` $= 0.999$: retain the fewest components that explain at least 99.9% of the
  total variance.

The fit produces:

- The sample mean $\bar{\mathbf{w}} = \frac{1}{N} \sum_{i=1}^{N} \mathbf{w}_i \in \mathbb{R}^{T}$.
- The loading matrix $\mathbf{P} \in \mathbb{R}^{T \times d}$, whose columns are the
  top $d \leq k$ eigenvectors of the sample covariance, ordered by decreasing
  eigenvalue.

Each training waveform is projected into the $d$-dimensional latent space:

$$\mathbf{z}_i = \mathbf{P}^\top (\mathbf{w}_i - \bar{\mathbf{w}}), \qquad i = 1, \dots, N.$$

### 4.4 Latent Distribution

The score vectors $\{\mathbf{z}_1, \dots, \mathbf{z}_N\} \subset \mathbb{R}^{d}$ are
modeled as i.i.d. draws from a multivariate normal distribution:

$$\mathbf{z} \sim \mathcal{N}(\boldsymbol{\mu},\, \boldsymbol{\Sigma})$$

where the parameters are estimated from the scores:

$$\boldsymbol{\mu} = \frac{1}{N} \sum_{i=1}^{N} \mathbf{z}_i \in \mathbb{R}^{d}$$

$$\boldsymbol{\Sigma} = \frac{1}{N-1} \sum_{i=1}^{N} (\mathbf{z}_i - \boldsymbol{\mu})(\mathbf{z}_i - \boldsymbol{\mu})^\top \in \mathbb{R}^{d \times d}$$

The covariance $\boldsymbol{\Sigma}$ is symmetrized explicitly
($\boldsymbol{\Sigma} \leftarrow \tfrac{1}{2}(\boldsymbol{\Sigma} + \boldsymbol{\Sigma}^\top)$)
to guard against floating-point asymmetry.

This captures both the marginal variance along each principal axis and the covariance
structure between components---information that would be lost by treating components as
independent (i.e. using a diagonal $\boldsymbol{\Sigma}$).

### 4.5 Sampling New Waveforms

To generate a synthetic waveform:

1. **Sample** a latent vector from the fitted distribution:

   $$\mathbf{z}_{\text{new}} \sim \mathcal{N}(\boldsymbol{\mu},\, \boldsymbol{\Sigma})$$

2. **Reconstruct** in (transformed) waveform space using the PCA loading matrix:

   $$\tilde{\mathbf{w}} = \mathbf{P}\,\mathbf{z}_{\text{new}} + \bar{\mathbf{w}}$$

3. **Inverse-transform** element-wise:

   $$\hat{w}_t = f^{-1}(\tilde{w}_t), \qquad t = 1, \dots, T$$

4. **Clamp** to enforce physical bounds:

   $$w_t^{*} = \max\!\bigl(w_{\min},\; \hat{w}_t\bigr)$$

   where $w_{\min} = \min_{i,t}\, \text{RyR\_OA}_{i,t}$ is the minimum value observed
   in the (untransformed) training data.

The result $\mathbf{w}^{*} \in \mathbb{R}^{T}$ is a synthetic RyR_OA time series that
can be plotted or fed into downstream models.

### 4.6 Relationship to Probabilistic PCA

The generative PCA model described here is closely related to **Probabilistic PCA**
(PPCA; Tipping & Bishop, 1999). PPCA assumes a latent variable model

$$\mathbf{w} = \mathbf{P}\,\mathbf{z} + \bar{\mathbf{w}} + \boldsymbol{\epsilon}, \qquad \mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}), \quad \boldsymbol{\epsilon} \sim \mathcal{N}(\mathbf{0}, \sigma^2 \mathbf{I})$$

and recovers standard PCA in the $\sigma^2 \to 0$ limit. The CellSims approach differs
in two ways: (i) it fits a **full covariance** $\boldsymbol{\Sigma}$ over latent scores
rather than assuming $\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I})$, and (ii)
it omits the observation noise term $\boldsymbol{\epsilon}$, reconstructing
deterministically from the sampled scores.

## 5. Implementation

The model is encapsulated in two types:

```julia
struct ApproxLog
    eps::Float64
end
(o::ApproxLog)(x) = log(x + o.eps)
inverse(o::ApproxLog) = x -> exp(x) - o.eps

struct GenerativePCA{F, D}
    fun::F              # transform f: ApproxLog or identity
    model::PCA{Float64} # fitted PCA model (stores P, w̄)
    dist::D             # MvNormal(μ, Σ) over latent scores
    min::Float64        # w_min from training data
end
```

Construction fits the full pipeline from a filtered DataFrame:

```julia
gen = GenerativePCA(df_sr_only, ApproxLog(0.1))       # log-space model
gen_id = GenerativePCA(df_sr_only, identity)           # linear model
```

The constructor internally performs:

```julia
function GenerativePCA(df_sr, f; k=25, pratio=0.999)
    # Group by run, filter incomplete runs
    gdf = groupby(df_sr, "run")
    n = nrow(first(gdf))
    gdf2 = filter(sub -> nrow(sub) == n, gdf)

    # Build W matrix and fit PCA
    w = reduce(hcat, f.(sub.RyR_OA) for sub in gdf2)   # T × N
    model = fit(PCA, w; maxoutdim=k, pratio)            # P, w̄

    # Fit latent distribution
    z = MultivariateStats.transform(model, w)            # d × N
    μ = vec(mean(z, dims=2))                             # d × 1
    Σ = Symmetric(cov(z; dims=2))                        # d × d
    dist = MvNormal(μ, Σ)

    GenerativePCA(f, model, dist, minimum(df_sr.RyR_OA))
end
```

Sampling (the struct is callable):

```julia
function (o::GenerativePCA)()
    f_inv = inverse(o.fun)
    w = reconstruct(o.model, rand(o.dist))   # P * z_new + w̄
    max.(o.min, f_inv.(w))                   # inverse transform + clamp
end

waveform = gen()   # returns Vector{Float64} of length T
```

## 6. Integration with the Analysis Pipeline

The `analyze` function automates the full workflow for a set of simulations:

1. Load all runner CRU outputs into a DataFrame.
2. Compute per-run SCR statistics (onset time $t_i$, offset time $t_f$, duration
   $\lambda = t_f - t_i$, peak, plateau).
3. Filter to SR-only runs.
4. Fit two generative models---one with $f = \text{identity}$, one with
   $f = \text{ApproxLog}(0.1)$.
5. Draw 10 samples from each model and save diagnostic plots.
6. Return all results in an `OrderedDict` for downstream use.

The probability of SCR is computed as

$$\hat{p}_{\text{SCR}} = \frac{N_{\text{SR}}}{N_{\text{total}}}$$

alongside the generative models, providing both a scalar summary statistic and a full
waveform-level generative capability for each experimental condition (BCL, ISO level,
SERCA scaling, etc.).

## 7. Assumptions and Limitations

- **Gaussian latent space.** The $\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})$
  assumption is exact when the true waveform distribution is Gaussian in the (possibly
  transformed) feature space. For strongly non-Gaussian distributions, the log
  transform helps but does not guarantee normality.
- **Linear reconstruction.** PCA captures only linear structure in $f$-transformed
  space. Nonlinear variability (e.g. bimodal waveform shapes) may be poorly represented
  unless $f$ approximately linearizes it.
- **Fixed dimensionality.** The number of retained components $d$ is chosen by variance
  explained ($\geq 99.9\%$), not by generative quality. In principle, cross-validated
  log-likelihood could yield a better choice of $d$.
- **Minimum clamp.** Clamping at $w_{\min}$ prevents negative values but can introduce
  a point mass at the boundary if the Gaussian tail extends significantly below
  $f(w_{\min})$ in latent space.
- **Sample size.** The covariance estimate $\boldsymbol{\Sigma}$ requires $N > d$.
  When few runs exhibit SCR, the estimated covariance may be noisy or singular.

## 8. Dependencies

| Package            | Role                                              |
|--------------------|---------------------------------------------------|
| MultivariateStats  | PCA fitting ($\mathbf{P}$, $\bar{\mathbf{w}}$), projection, reconstruction |
| Distributions      | $\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\Sigma})$ sampling |
| DataFrames / CSV   | Data loading and manipulation                     |
| Statistics         | $\boldsymbol{\mu}$, $\boldsymbol{\Sigma}$ estimation |
| LinearAlgebra      | `Symmetric` wrapper for $\boldsymbol{\Sigma}$     |
| StatsPlots         | Visualization of waveforms and diagnostics        |
