# The effective phase relaxation time of the reference fields: SAM's `QLDot` (the condensation
# rate, g kg⁻¹ s⁻¹) regressed on the local supersaturation excess qᵛ⁺ 𝒮 in the interior,
# QLDot = qᵛ⁺ 𝒮 / τ. The result calibrates the relaxation time of the chamber's host
# microphysics (`chamber_microphysics(; relaxation_time)`, `HOST_TAU` in `run_windows.jl`).
#
# Usage: julia --project analysis/reference_condensation.jl [stride]
using Statistics, Printf
include(joinpath(homedir(), "TurbulentActivation.jl", "reference", "read_bin3d.jl"))
using .SAMBin3D
dir = joinpath(homedir(), "anderson_reference", "LES", "LES")
esat(T) = 610.94 * exp(17.625 * (T - 273.15) / (T - 30.11))
qsat(T, p) = 0.622 * esat(T) / (p - esat(T))
p₀ = 1e5
num = 0.0; den = 0.0; n = 0
allQ = Float64[]; allX = Float64[]
stride = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 5000
for step in 90000:stride:180000
    s = read_bin3d(joinpath(dir, @sprintf("PiChamber_huji_19K_trj_32_%010d.bin3D", step)))
    T = s.TABS; q = s.QV ./ 1000; qn = s.QN ./ 1000; Qdot = s.QLDot
    if n == 0
        println("QLDot: mean ", mean(Qdot), " min ", minimum(Qdot), " max ", maximum(Qdot), "  QP mean ", mean(s.QP), " max ", maximum(s.QP), "  QN mean ", mean(qn))
        println("fields: ", s.names)
    end
    S = @. p₀ * q / (0.622 + q) / esat(T) - 1
    X = @. qsat(T, p₀) * S           # kg/kg supersaturation excess
    ix = (3:62, 3:62, 3:30)
    x = vec(view(X, ix...)); y = vec(view(Qdot, ix...))
    append!(allX, x); append!(allQ, y)
    global num += sum(x .* y); global den += sum(x .^ 2); global n += 1
end
# Try several unit hypotheses for QLDot (kg/kg/s, g/kg/s, g/kg/day, kg/kg/day)
slope = num / den
println("snapshots: ", n, "  correlation(QLDot, qsat S) = ", cor(allX, allQ))
for (name, f) in (("kg/kg/s", 1.0), ("g/kg/s", 1e-3), ("g/kg/day", 1e-3 / 86400), ("kg/kg/day", 1 / 86400), ("g/kg/hour", 1e-3 / 3600))
    @printf("  if QLDot is in %-10s → 1/τ = %.4g s⁻¹, τ = %.3g s\n", name, slope * f, 1 / (slope * f))
end
# Sign check and conditional means
pos = allX .> 0
@printf("mean QLDot where S>0: %.4g ; where S<0: %.4g ; fraction S>0: %.3f\n", mean(allQ[pos]), mean(allQ[.!pos]), mean(pos))
# Domain-mean condensation vs the LWC turnover: QN mean / QLDot mean = residence time
