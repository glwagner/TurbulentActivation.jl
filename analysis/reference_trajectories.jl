# Anderson's experiment replayed through the reference SAM fields, in Julia and with Breeze's
# droplet physics: passive tracers are advected through the 0.5 s snapshots with linear
# interpolation in space and time (as in `parcels.py`), the supersaturation is sampled along
# each trajectory (Anderson's Magnus-based definition), and the droplet and its replicas are
# grown with the same implicit κ-Köhler step the online experiment uses. Output: the
# activated fraction after 60 s for the fluctuating, uniform, and instantaneous replicas at
# every target mean, the Lagrangian PDF and autocorrelation time of the supersaturation, and
# sample trajectories — the reference against which the online chamber is judged.
#
# Usage: julia --project analysis/reference_trajectories.jl [N] [n_windows] [first_step] [output.jld2]
using Breeze, JLD2, Statistics, Random, Printf
using Breeze.LagrangianMicrophysics: implicit_growth_step
include(joinpath(@__DIR__, "..", "reference", "read_bin3d.jl"))
using .SAMBin3D

N = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 10_000
n_windows = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 3
first_step = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 90000
output = length(ARGS) ≥ 4 ? ARGS[4] : joinpath(@__DIR__, "..", "runs", "reference_trajectories.jld2")

dir = joinpath(homedir(), "anderson_reference", "LES", "LES")
snapshot(step) = read_bin3d(joinpath(dir, @sprintf("PiChamber_huji_19K_trj_32_%010d.bin3D", step)))

#####
##### Anderson's aerosol and thermodynamics
#####

FT = Float64
constants = ThermodynamicConstants(FT)
Dᵈ, κ, T₀ = FT(130e-9), FT(1), FT(289)
Dᶜ = critical_diameter(Dᵈ, κ, T₀, constants)
𝒮ᶜ = critical_supersaturation(Dᵈ, κ, T₀, constants)
targets = collect(range(-0.04, 0.05, length=19))
p₀ = FT(100000)                                       # Pa; Anderson uses a fixed 1000 hPa
parameters = DropletDynamics(FT; substeps=25)         # 0.02 s substeps within each 0.5 s sample
esat(T) = 610.94 * exp(17.625 * (T - 273.15) / (T - 30.11))
supersaturation(T, q) = p₀ * q / (0.622 + q) / esat(T) - 1
@info "Reference replay" N n_windows first_step Dᶜ 𝒮ᶜ

#####
##### Trilinear interpolation on the collocated snapshot grid (cell centres, Δ = 3.125 cm)
#####

const Δ = 0.03125
const Lx, Ly, Lz = 2.0, 2.0, 1.0
const nx, ny, nz = 64, 64, 32

@inline function interpolate(field, x, y, z)
    # centres at Δ/2 + Δ(i-1); fractional index of the cell centre below
    fi = (x - Δ/2) / Δ + 1; fj = (y - Δ/2) / Δ + 1; fk = (z - Δ/2) / Δ + 1
    i = clamp(floor(Int, fi), 1, nx - 1); j = clamp(floor(Int, fj), 1, ny - 1); k = clamp(floor(Int, fk), 1, nz - 1)
    a = clamp(fi - i, 0, 1); b = clamp(fj - j, 0, 1); c = clamp(fk - k, 0, 1)
    @inbounds begin
        f00 = field[i, j, k] * (1 - a) + field[i+1, j, k] * a
        f10 = field[i, j+1, k] * (1 - a) + field[i+1, j+1, k] * a
        f01 = field[i, j, k+1] * (1 - a) + field[i+1, j, k+1] * a
        f11 = field[i, j+1, k+1] * (1 - a) + field[i+1, j+1, k+1] * a
    end
    return ((f00 * (1 - b) + f10 * b) * (1 - c) + (f01 * (1 - b) + f11 * b) * c)
end

# Reflect a position back into the box
reflect(x, L) = (x = mod(x, 2L); x > L ? 2L - x : x)

# One trajectory step of δt between snapshots A (time 0) and B (time 0.5 s), midpoint rule
function advance!(X, A, B, t, δt)
    θ = t / 0.5
    for n in eachindex(X)
        x, y, z = X[n]
        u = (1 - θ) * interpolate(A.U, x, y, z) + θ * interpolate(B.U, x, y, z)
        v = (1 - θ) * interpolate(A.V, x, y, z) + θ * interpolate(B.V, x, y, z)
        w = (1 - θ) * interpolate(A.W, x, y, z) + θ * interpolate(B.W, x, y, z)
        xm, ym, zm = reflect(x + u * δt/2, Lx), reflect(y + v * δt/2, Ly), reflect(z + w * δt/2, Lz)
        θm = (t + δt/2) / 0.5
        u = (1 - θm) * interpolate(A.U, xm, ym, zm) + θm * interpolate(B.U, xm, ym, zm)
        v = (1 - θm) * interpolate(A.V, xm, ym, zm) + θm * interpolate(B.V, xm, ym, zm)
        w = (1 - θm) * interpolate(A.W, xm, ym, zm) + θm * interpolate(B.W, xm, ym, zm)
        X[n] = (reflect(x + u * δt, Lx), reflect(y + v * δt, Ly), reflect(z + w * δt, Lz))
    end
    return nothing
end

sample_supersaturation(X, A) = [supersaturation(interpolate(A.TABS, x, y, z), interpolate(A.QV, x, y, z) / 1000) for (x, y, z) in X]

#####
##### Windows
#####

M = length(targets)
rng = MersenneTwister(2023)
window_steps = 60 / 0.5                                # 120 snapshots per 60 s window
results = Dict{String, Any}("targets" => targets, "N" => N, "window" => 60.0, "first_step" => first_step)
series = Dict{String, Any}()

for w in 1:n_windows
    step₀ = first_step + (w - 1) * 25 * Int(window_steps)
    X = [(Lx * rand(rng), Ly * rand(rng), Lz * rand(rng)) for _ in 1:N]
    A = snapshot(step₀)
    S = sample_supersaturation(X, A)
    Sₜ = zeros(N, Int(window_steps) + 1); Sₜ[:, 1] .= S
    # Every droplet starts at haze equilibrium with its own initial supersaturation, as in `equilibrate_h2o`;
    # replicas start at equilibrium with their shifted initial value
    D² = [equilibrium_diameter(min(s, 0.9𝒮ᶜ), Dᵈ, κ, T₀, constants)^2 for s in S]
    D²ᶠ = [[equilibrium_diameter(min(targets[m] + S[n] - mean(S), 0.9𝒮ᶜ), Dᵈ, κ, T₀, constants)^2 for m in 1:M] for n in 1:N]
    D²ᵘ = [[equilibrium_diameter(min(targets[m], 0.9𝒮ᶜ), Dᵈ, κ, T₀, constants)^2 for m in 1:M] for n in 1:N]
    instant = [falses(M) for _ in 1:N]
    @info "window $w" step₀ mean_S=mean(S) std_S=std(S)
    for s in 1:Int(window_steps)
        B = snapshot(step₀ + 25s)
        for sub in 1:5                                 # 0.1 s trajectory steps, as in Anderson
            advance!(X, A, B, (sub - 1) * 0.1, 0.1)
        end
        Sⁿ = sample_supersaturation(X, B)
        Sₜ[:, s + 1] .= Sⁿ
        S̄ = mean(Sⁿ)                                  # domain-mean stand-in for the per-trajectory mean (audit)
        Tₙ = [interpolate(B.TABS, x, y, z) for (x, y, z) in X]
        for n in 1:N
            T = Tₙ[n]
            D²[n] = implicit_growth_step(D²[n], Sⁿ[n], T, p₀, Dᵈ, κ, 0.5, parameters, constants)
            for m in 1:M
                sᶠ = targets[m] + (Sⁿ[n] - S̄)
                D²ᶠ[n][m] = implicit_growth_step(D²ᶠ[n][m], sᶠ, T, p₀, Dᵈ, κ, 0.5, parameters, constants)
                D²ᵘ[n][m] = implicit_growth_step(D²ᵘ[n][m], targets[m], T, p₀, Dᵈ, κ, 0.5, parameters, constants)
                instant[n][m] = sᶠ ≥ 𝒮ᶜ
            end
        end
        A = B
    end
    fluctuating = [count(n -> D²ᶠ[n][m] ≥ Dᶜ^2, 1:N) / N for m in 1:M]
    uniform = [count(n -> D²ᵘ[n][m] ≥ Dᶜ^2, 1:N) / N for m in 1:M]
    instantaneous = [count(n -> instant[n][m], 1:N) / N for m in 1:M]
    results["window_$w"] = Dict("step" => step₀, "fluctuating" => fluctuating, "uniform" => uniform, "instantaneous" => instantaneous,
                                "S0" => mean(S), "true_activated" => count(≥(Dᶜ^2), D²) / N)
    series["window_$w"] = Dict("sample_times" => collect(0:0.5:60), "S" => Sₜ')
    @printf("window %d: ⟨S⟩ = %+.3f %%  σ(S) = %.3f %%  true active %.3f\n   target  fluctuating  uniform  instantaneous\n",
            w, 100mean(Sₜ), 100std(Sₜ), results["window_$w"]["true_activated"])
    for m in 1:M
        @printf("   %+.3f   %.3f   %.3f   %.3f\n", targets[m], fluctuating[m], uniform[m], instantaneous[m])
    end
end

jldsave(output; results, series)
@info "saved" output
println("REFERENCE REPLAY OK")
