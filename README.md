# README for Julia code:

## Usage

- REQUIREMENT: You've run `CODE/Makefile` to generate the `CODE/model_single_3D` executable.
  - Makefiles are finnicky...this is the kind of thing Claude Code is good at figuring out.
- Change directory to root of repo, then:

```julia
using CellSims

# All settings used to generate different Ca_NSR inputs
sims = all_simulations()
```

## Generating Runs

```
# Generating simulations for a single setting of BCL/ISO:
julia --project=. src/generate.jl BCL ISO [N]

# Example: 10 runs with BCL=300 and ISO=1
julia --project=. src/generate.jl 300 1 10

# Generating simulations across multiple settings
julia --project=. src/generate_all.jl N

# Example: 5 runs of every setting in all_simulations()
julia --project=. src/generate_all.jl 5
```

## Analyzing Simulation Results

```julia
include("src/analyze.jl")
```

- This script will populate `data/results/sr_$BCL_$ISO` with plots of the results.
- Summary plots in `data/results/summary`
- A DataFrame `df` is also created with all the grouped results.
