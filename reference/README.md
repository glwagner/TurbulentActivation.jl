# Reference data (R0)

The released SAM fields, the DNS case, and the Zenodo code archive are mirrored outside
the repository in `~/anderson_reference/` (not committed):

- `LES/`  3601 SAM `bin3D` snapshots `PiChamber_huji_19K_trj_32_<step>.bin3D`, steps 90000 to
  180000 every 25 (0.5 s at Δt = 0.02 s; the second half-hour of the run), 9.4 MB each,
  plus one converted NetCDF snapshot and the `bin3D2nc.f` converter, from
  https://portal.nersc.gov/project/m1657/LagrangianDroplets/LES/
- `DNS/case3.mat` (1.3 GB) from the `DNS/` directory of the same portal
- `zenodo/` the archived code release 10.5281/zenodo.7259494

Checksums (`sha256sum`) are recorded in `checksums.txt` once the mirror is complete.
