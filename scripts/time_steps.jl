# Per-step cost of the chamber on the current device: dynamics alone, with plain κ-Köhler droplets,
# and with Anderson's replica droplets; then a short stability check at a larger time step.
#
# Usage: julia --project scripts/time_steps.jl [Nx Ny Nz] [N_droplets] [steps]
using TurbulentActivation, Breeze, Oceananigans, CUDA, Printf, Random, Statistics
using Oceananigans.Units

Nx, Ny, Nz = length(ARGS) ≥ 3 ? parse.(Int, ARGS[1:3]) : (64, 64, 32)
N = length(ARGS) ≥ 4 ? parse(Int, ARGS[4]) : 10_000
steps = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 500
arch = CUDA.functional() ? GPU() : CPU()
sync() = arch isa GPU ? CUDA.synchronize() : nothing

chamber = PiChamber()
grid = pi_chamber_grid(chamber, arch; size=(Nx, Ny, Nz))
aerosol = anderson_aerosol()
targets = anderson_targets()
temperature = (chamber.bottom_temperature + chamber.top_temperature) / 2
rng = MersenneTwister(1)

function time_model(label, model; Δt=0.02)
    initialize_chamber!(model, chamber; temperature, relative_humidity=0.8, rng=MersenneTwister(1))
    for _ in 1:50; time_step!(model, Δt); end            # warm-up and compilation
    sync(); t₀ = time_ns()
    for _ in 1:steps; time_step!(model, Δt); end
    sync(); ms = (time_ns() - t₀) / 1e6 / steps
    @printf("%-32s %7.2f ms/step\n", label, ms)
    return ms
end

t_dyn = time_model("dynamics only", pi_chamber_model(chamber, grid))

droplets = seed_droplets(aerosol, grid, N; temperature, relative_humidity=0.8, rng)
particles = LagrangianParticles(droplets; dynamics=DropletDynamics())
t_drop = time_model("$(N) droplets", pi_chamber_model(chamber, grid; particles))

replicas = seed_replica_droplets(aerosol, grid, N, targets; temperature, supersaturation=-0.2, rng)
particles = LagrangianParticles(replicas; dynamics=ReplicaDynamics(targets))
model = pi_chamber_model(chamber, grid; particles)
t_rep = time_model("$(N) droplets × 19 replicas", model)

microphysics = chamber_microphysics(; relaxation_time=72)
t_cloud = time_model("replicas + one-moment host", pi_chamber_model(chamber, grid; particles, microphysics))

@printf("particles add %.2f ms (plain) and %.2f ms (replicas) per step; host microphysics adds %.2f ms\n",
        t_drop - t_dyn, t_rep - t_dyn, t_cloud - t_rep)

# Stability at larger time steps: two minutes of the replica chamber at each Δt
u, v, w = model.velocities
for Δt in (0.04, 0.06)
    initialize_chamber!(model, chamber; temperature, relative_humidity=0.8, rng=MersenneTwister(1))
    simulation = Simulation(model; Δt, stop_time=2minutes)
    Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
    try
        run!(simulation)
        @printf("Δt = %.2f s: stable for 2 min, max|w| = %.3f, CFL = %.2f\n", Δt, maximum(abs, w), Δt * maximum(abs, w) / minimum_xspacing(grid))
    catch err
        @printf("Δt = %.2f s: FAILED (%s)\n", Δt, sprint(showerror, err)[1:min(end, 80)])
    end
end
println("TIMING OK")
