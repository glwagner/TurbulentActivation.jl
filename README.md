# TurbulentActivation.jl

Reproducing [Anderson et al. (2023)](https://doi.org/10.1029/2022GL102635) — droplet
activation enhanced by turbulent supersaturation fluctuations in the Michigan Tech Pi
Chamber — with *online* Lagrangian droplets in [Breeze.jl](https://github.com/NumericalEarth/Breeze.jl).

The infrastructure (six-wall bulk fluxes, κ-Köhler droplet physics, Lagrangian droplet
dynamics) lives in Breeze on the branch `glw/anderson-chamber`. This repository holds the
reproduction itself: the chamber configuration, the aerosol, the experiment drivers, the
reference audit, and the analysis.

## Layout

- `src/chamber.jl` — `PiChamber` configuration and `pi_chamber_model` (walls, dynamics, advection)
- `src/aerosol.jl` — Anderson's NaCl aerosol and droplet seeding at equilibrium
- `scripts/` — runnable drivers (`run_chamber.jl`, `run_droplets.jl`) and a Slurm helper
- `reference/audit.md` — paper vs. archive vs. live-code discrepancies (kept current)
- `analysis/` — post-processing

## Status

Field parity with the reference LES reached for the one-way experiment, with a physical wall
model and the reference's own condensation sink; no reproduction claim against the paper's
numbers yet.

The chamber (six-wall bulk fluxes, WENO(5) implicit LES with bounds-preserving WENO on the
moisture species, adaptive time step at CFL 0.7) convects, 10⁴ κ-Köhler droplets are advected
and grown online, and every droplet carries Anderson's counterfactual replicas for 19 target
mean supersaturations. The same droplet physics and the same replica construction are also
replayed offline through Anderson's SAM fields (`analysis/reference_trajectories.jl`), which
is the reference the online chamber is judged against.

Two quantities were calibrated from the reference fields themselves. The condensation rate in
the SAM snapshots correlates at 0.96 with the local supersaturation excess and gives a phase
relaxation time of 72 s (`analysis/reference_condensation.jl`), which the chamber's warm-only
one-moment host reproduces with `relaxation_time = 28` (Breeze's rate carries the latent-heat
factor Γ ≈ 2.5). The wall model is a log law with the Monin–Obukhov stability correction on the
floor and ceiling and neutral side walls (`log_law_coefficient`), and a 1.25 mm roughness
length gives, at 64 × 64 × 32 with a 15 min spin-up and three 60 s windows:

| quantity | Breeze chamber | reference SAM fields |
|---|---:|---:|
| interior spread of 𝒮 | 1.24 % | 1.23 % |
| spread of 𝒮 over the box | 1.96 % | 1.94 % |
| Lagrangian spread of 𝒮 along droplets | 1.89–1.93 % | 1.85–1.96 % |
| correlation time τₛ | 3.0–3.1 s | 2.9–3.1 s |
| σ(T), σ(qᵛ) | 0.99 K, 0.79 g kg⁻¹ | 0.85 K, 0.72 g kg⁻¹ |
| mean 𝒮, interior | +0.19 % | +0.71 % |
| mean T | 288.6 K | 287.6 K |
| activated fraction after 60 s at a −1 % target (fluctuating) | 0.36–0.41 | 0.355 |
| at 0 % (fluctuating / instantaneous) | 0.82–0.84 / 0.59–0.62 | 0.85 / 0.51 |
| at −2 % (fluctuating) | 0.085–0.09 | 0.11 |

![Default chamber against the reference replay](figures/default_chamber_vs_reference_activation.png)

The roughness sweep brackets it (1 mm undershoots, 1.5 mm overshoots; `figures/wall_model_activation.png`).
The residuals are a 1 K warm bias of the mean temperature and a mean supersaturation of
+0.2 % against +0.7 %, both a flux-partition question between floor, ceiling, and side walls.
A constant coefficient of 2 × 10⁻² reproduces the fluctuations equally well but leaves the
mean near zero, and the neutral log-law value 6 × 10⁻³ gives fluctuations three times too
weak. The 128 × 128 × 64 chamber reproduces the 64 × 64 × 32 statistics at the same settings
(`figures/ladder_activation.png`). WENO(9) gives an interior spread about 12 % larger than
WENO(5) at the same roughness.

Open items: the warm bias and the mean; the correlation-time definition (Anderson reports
7.5 s where the integral estimator gives 3 s in his own fields); second-order particle
advection; checkpointing; the Breeze example.

## Installing

The package depends on Breeze's `glw/anderson-chamber` branch, declared in `[sources]`, so
`Pkg.instantiate()` fetches it. To develop against a local Breeze checkout instead:

```julia
using Pkg; Pkg.develop(path="/path/to/Breeze.jl")
```

## Running

```julia
julia --project=. scripts/run_windows.jl 64 64 32 25 3 10000 parity
```
spins the chamber up for 25 minutes and then runs three 60 s windows of Anderson's experiment
with 10⁴ droplets, writing `parity_activation.jld2` (activation curves, Lagrangian samples,
chamber statistics) and `parity_profiles.jld2` (horizontal-mean profiles). Environment
variables: `WALL_C` (wall transfer coefficient, default the `PiChamber` value), `DT` (time
step, 0.02 s), `HOST=bulk` with `HOST_TAU` (warm-only one-moment host microphysics with the
given relaxation time in seconds). `scripts/gpu.sbatch` wraps any script for Slurm.

Analysis: `analysis/activation_curve.jl` (curves and bands of one run),
`analysis/supersaturation_statistics.jl` (Lagrangian PDF, autocorrelation, correlation time),
`analysis/compare_curves.jl online.jld2 reference.jld2 out.png` (online against the replay),
`analysis/coefficient_sweep.jl reference.jld2 out.png C=file.jld2 ...`,
`analysis/reference_statistics.jl` (Eulerian statistics of the mirrored SAM fields), and
`analysis/reference_trajectories.jl N windows first_step out.jld2` (the replay of Anderson's
experiment through the SAM fields with Breeze's droplet physics).

Tests: `julia --project=. test/runtests.jl` (`Pkg.test()` currently fails inside Pkg with the
git-sourced Breeze dependency).
