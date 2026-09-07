"""Grow droplets along prescribed supersaturation traces with Anderson's own code.

Reads the CSV traces written by `analysis/export_traces.jl` (column 1 is time, every other
column is one droplet's supersaturation), integrates each with `microphysics.one_particle`
from `lfierce2/LagrangianDroplets` exactly as `run.one_particle` does (initial size from
`equilibrate_h2o`, `solve_ivp` with BDF and a 10 ms maximum step, `ignore_Tp=True`), and
writes the diameters back as CSV. Activation is judged against `pyrcel.thermo.kohler_crit`.

Usage: python3 analysis/anderson_growth.py <tracedir> [dry_diameter] [kappa] [T] [P]
"""
import os, sys, glob
import numpy as np

CODE = os.path.expanduser("~/anderson_reference/zenodo/lfierce2-LagrangianDroplets-488336d")
sys.path.insert(0, os.environ.get("NUMBA_SHIM", ""))
sys.path.insert(0, CODE)
import microphysics                     # noqa: E402
from pyrcel.thermo import kohler_crit   # noqa: E402
from scipy.integrate import solve_ivp   # noqa: E402

tracedir = sys.argv[1]
Dd = float(sys.argv[2]) if len(sys.argv) > 2 else 130e-9
kappa = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
T = float(sys.argv[4]) if len(sys.argv) > 4 else 287.6
P = float(sys.argv[5]) if len(sys.argv) > 5 else 1e5
accom, rho_aero, N_conc = 0.3, 2160.0, 1e6

r_crit, S_crit = kohler_crit(T, Dd / 2.0, kappa)
print("critical diameter %.4e m, critical supersaturation %.5f %%" % (2 * r_crit, 100 * S_crit))

for path in sorted(glob.glob(os.path.join(tracedir, "traces_*.csv"))):
    A = np.loadtxt(path, delimiter=",")
    t, traces = A[:, 0], A[:, 1:]
    out = np.zeros_like(traces)
    activated = 0
    for n in range(traces.shape[1]):
        S = traces[:, n]
        S0 = min(S[0], 0.0)            # as in run.one_particle: equilibrate at a subsaturated value
        D0, _, _ = microphysics.equilibrate_h2o(np.array([Dd]), np.array([kappa]), np.array([N_conc]), T, S0, P)
        S_of_t = lambda τ: np.interp(τ, t, S)
        rhs = lambda τ, y: microphysics.one_particle(y, τ, 0.0, S_of_t(τ), T, P, Dd, kappa, N_conc,
                                                     accom, rho_aero, ignore_Tp=True)
        soln = solve_ivp(rhs, [t[0], t[-1]], np.array([D0[0], T]), method="BDF",
                         max_step=1e-2, t_eval=t)
        out[:, n] = soln["y"][0]
        activated += out[-1, n] >= 2 * r_crit
    name = os.path.basename(path).replace("traces_", "diameters_")
    np.savetxt(os.path.join(tracedir, name), np.column_stack([t, out]), delimiter=",")
    print("%s: %d droplets, activated fraction %.3f" % (name, traces.shape[1], activated / traces.shape[1]))
print("ANDERSON GROWTH OK")
