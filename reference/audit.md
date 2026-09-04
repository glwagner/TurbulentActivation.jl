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
| Wall model | SAM H-M case (Yang et al. 2022) | not in this repository | audit Yang et al. for roughness lengths and stability form |
