# The Pi Chamber with online one-way κ-Köhler droplets.
# Usage: julia --project scripts/run_droplets.jl [Nx Ny Nz] [stop_minutes] [N_droplets] [prefix]
using TurbulentActivation, Breeze, Oceananigans, CUDA, Printf, Statistics, Random
using Oceananigans.Units

Nx, Ny, Nz = length(ARGS) ≥ 3 ? parse.(Int, ARGS[1:3]) : (32, 32, 16)
stop_minutes = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 0.5
N = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 1000
prefix = length(ARGS) ≥ 6 ? ARGS[6] : "droplets"
arch = CUDA.functional() ? GPU() : CPU()

# Environment overrides as in run_windows.jl: WALL_C, DT, HOST=bulk with HOST_TAU, PARTICLES=0 for none
chamber = PiChamber(; transfer_coefficient=parse(Float64, get(ENV, "WALL_C", string(PiChamber().transfer_coefficient))))
Δt = parse(Float64, get(ENV, "DT", "0.02"))
host = get(ENV, "HOST", "none")
microphysics = host == "bulk" ? chamber_microphysics(; relaxation_time=parse(Float64, get(ENV, "HOST_TAU", "5"))) : nothing
grid = pi_chamber_grid(chamber, arch; size=(Nx, Ny, Nz))
aerosol = anderson_aerosol()
relative_humidity = 0.8
temperature = (chamber.bottom_temperature + chamber.top_temperature) / 2
droplets = seed_droplets(aerosol, grid, N; temperature, relative_humidity, rng=MersenneTwister(1234))
particles = get(ENV, "PARTICLES", "1") == "1" ? LagrangianParticles(droplets; dynamics=DropletDynamics()) : nothing
model = pi_chamber_model(chamber, grid; particles, microphysics)
initialize_chamber!(model, chamber; temperature, relative_humidity)
@info "Pi Chamber with droplets" chamber aerosol arch size=(Nx, Ny, Nz) stop_minutes N

simulation = Simulation(model; Δt, stop_time=stop_minutes * minutes)
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)

ℋ = RelativeHumidityField(model)
u, v, w = model.velocities

function progress(sim)
    compute!(ℋ)
    s = droplet_statistics(droplets)
    @printf("iter %6d  t = %7.2f s  max|w| = %.3f  ⟨ℋ⟩ = %.3f  max ℋ = %.3f  droplets: ⟨𝒮⟩ = %+.4f  σ(𝒮) = %.4f  ⟨D⟩ = %.2f μm  active = %.3f  wall = %s\n",
            iteration(sim), time(sim), maximum(abs, w), mean(ℋ), maximum(ℋ),
            s.mean_supersaturation, s.std_supersaturation, 1e6 * s.mean_diameter, s.activated_fraction,
            prettytime(sim.run_wall_time))
end
add_callback!(simulation, progress, TimeInterval(10))

profiles = (; T=Average(model.temperature, dims=(1, 2)), ℋ=Average(ℋ, dims=(1, 2)),
              ww=Average(@at((Center, Center, Center), w^2), dims=(1, 2)))
simulation.output_writers[:profiles] = JLD2Writer(model, profiles; filename="$(prefix)_profiles.jld2",
                                                  schedule=TimeInterval(5), overwrite_existing=true)
simulation.output_writers[:particles] = JLD2Writer(model, (; particles=model.particles); filename="$(prefix)_particles.jld2",
                                                   schedule=TimeInterval(0.5), overwrite_existing=true)
run!(simulation)
println(droplet_statistics(droplets))
println("DROPLETS OK")
