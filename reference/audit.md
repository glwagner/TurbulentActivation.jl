# Reference audit: Anderson et al. (2023)

Live code `lfierce2/LagrangianDroplets` read at commit `5ca80674d173481bc1c51b96d2f2d5c10a728ee9`
(2026-09-03). The archived Zenodo release v1.0.0 (10.5281/zenodo.7259494, commit `488336d`,
mirrored in `~/anderson_reference/zenodo/`) was diffed against it on 2026-09-04: `parcels.py`,
`microphysics.py`, and `run.py` differ only in file headers and one directory-creation line, so
the archive and the live code are physically identical. In particular both use κ = 1.0
(`main_particle-traces.py`, "0.65 better?"), a 130 nm dry diameter, N = 10⁶ m⁻³, and the 19
target means from −4 % to +5 %; the archived figure script `main_make-figs.py` plots 15 targets
from −4 % to +3 %, which is presumably the abscissa of the published figure.

| Item | Paper / draft | Live code | Consequence |
|---|---|---|---|
| Hygroscopicity | κ = 1.2 | `kappa = 1.0  # 0.65 better?` (`main_particle-traces.py`) | keep both as named cases |
| Dry diameter, density, N | 130 nm | 130 nm, 2160 kg m⁻³, 10⁶ m⁻³ | match |
| Mean-supersaturation sweep | 19 values, −4 % to +5 % | `np.linspace(-0.04, 0.05, 19)` | use exactly this abscissa |
| Counterfactual construction | shift a common fluctuation field | `SS + avg_SS − mean(SS)` per trajectory | the mean is per trajectory, not domain-wide |
| Growth equation | Maxwell–Mason with κ-Köhler 𝒮ᵉ | `dr/dt = (G/r)(S − S_eq)`, `G = 1/(G_a + G_b)`, kinetic corrections `dv(T, r, P, accom)` and `ka(T, r, ρ)`, accommodation 0.3, thermal accommodation 0.96; state in ln D; particle temperature ignored | ported verbatim into `Breeze.LagrangianMicrophysics` |
| Integrator | — | SciPy BDF, `max_step` 10⁻² or 10⁻³ s, on 0.5 s trajectory samples | oracle for the implicit D² step |
| Initial wet size | equilibrium at the initial S | bisection of `Seq(r) − SS` between `r_dry` and `r_crit` | `equilibrium_diameter` in Breeze |
| Supersaturation from fields | S = qᵛ/qᵛ⁺(T, p) − 1 | Magnus `es` with fixed P = 1000 hPa in `parcels.py`; P₀ = 101325 Pa in `run.py` | Breeze uses its own thermodynamics at the reference pressure; quantify the difference |
| Trajectories | passive tracers | adaptive RK45 on 4-D linear interpolation of U, V, W, T, QV every 0.5 s | Breeze advects online, forward Euler per step (second-order integrator planned) |
| Settling | "first-order settling model" | none in `parcels.py`; wall handling via `RBC` | locate in the archive or `process.py`; do not assume |
| Dependencies | — | pyrcel (constants, thermo), numba.pycc, netCDF4 | freeze in the R0 environment |
| Wall model | SAM H-M case (Yang et al. 2022) | not in this repository | **Resolved from Yang et al. (2022, JAMES, §2):** "The floor and ceiling surfaces are water saturated, and the sensible and latent heat fluxes are calculated based on the Monin-Obukhov similarity theory. Side-wall sensible and latent heat fluxes are calculated using the same parameterization. The wall surface wetness is set to be 0.78 (0.78 times its saturated value), such that the domain-averaged supersaturation in steady state without cloud droplets is about 2.5%". Roughness lengths are not given there (Thomas et al. 2019, paywalled for us). The fluxes are grid-spacing dependent: "A smaller grid spacing leads to a larger gradient in temperature and water vapor mixing ratio between surface and the adjacent air, leading to larger calculated sensible and latent heat fluxes." Breeze uses constant bulk coefficients; C = 6e-3 gives a chamber that is too dry and too quiet (⟨𝒮⟩ −1.1 %, σ 0.6 %), C = 2e-2 reproduces the reference (⟨𝒮⟩ +1.4 %, σ 1.7 %), see `figures/online_C2e-2_vs_reference_activation.png` |
| Correlation time τₛ | "Using s along the ensemble of Lagrangian parcels, τₛ ≈ 7.5 s" (preprint), no definition given | not computed anywhere in the public code (no autocorrelation in `process.py`, `parcels.py`, `main_*.py`) | our integral of the Lagrangian autocorrelation to its first zero crossing gives 2.9–3.1 s in the same fields (`analysis/supersaturation_statistics.jl`); the autocorrelation first crosses zero near 9 s and reaches 1/e near 4 s, so 7.5 s is probably a different estimator — treat τₛ as a shape comparison, not a number to match |

## From the preprint (arXiv:2210.15766v1, mirrored as `~/anderson_reference/anderson2023_arxiv.txt`)

| Item | Preprint | Code | Consequence |
|---|---|---|---|
| Condensation (mass accommodation) coefficient | `ac = 1`; thermal accommodation `aT = 0.96` | `accom = 0.3` (`main_particle-traces.py` → `run.main`) | named sensitivity case; Breeze default follows the code (0.3) |
| Hygroscopicity | κ = 1 (v1 preprint) | κ = 1.0 | the draft's κ = 1.2 must come from the published GRL version; check |
| LES mean supersaturation | with cloud (bin microphysics, injected NaCl): s̄ = −2.7 % in the domain (Fig. 1 case); without droplets ≈ +2.5 % "consistent with chamber observations" | — | gate 2 target for the cloud-free chamber: ⟨𝒮⟩ ≈ +2.5 %; the cloudy reference has ⟨𝒮⟩ ≈ −2.7 % |
| Steady state | "after about 5 min"; 1 h total; fields every 0.5 s | files from step 90000 (30 min) | our spin-up must be checked, not assumed |
| Lagrangian statistics | τₛ ≈ 7.5 s along the parcel ensemble; τ_evap = 9.9 s at s̄ = −1.2 % and 120 s at −0.2 % | — | targets for the L1 cloudy gate |
| Trajectories | passive tracers, SciPy IVP solver with max step 0.1 s on linear 4-D interpolation | `solve_ivp` default RK45 | second-order integrator (A1c) desirable |
| Settling | Stokes decay of number per parcel available, but NOT applied when computing the activated fraction | none in `parcels.py` | ignore settling for the activation curves |
| Supersaturation | `s = e/e_sat − 1`, `e = p q/(0.622 + q)`, Magnus `e_sat`, p in kPa | fixed P = 1000 hPa in `parcels.py` | Breeze uses its own thermodynamics at the reference pressure |
| Sweep and curves | activated fraction after 60 s vs s̄; uniform = step at s_crit; fluctuating = `D ≥ D_crit`; instantaneous = `s ≥ s_crit` (τ_evap = 0) | 19 targets −4 %…+5 %; figure script 15 targets −4 %…+3 % | matches the online replicas |
