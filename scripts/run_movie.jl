# The chamber, saved as slices for a movie: after a spin-up, write a vertical slice through the
# middle of the box (speed, vertical velocity, supersaturation) and the droplets near that slice,
# at movie cadence. `analysis/make_movie.jl` turns the output into an animation.
#
# Usage: julia --project scripts/run_movie.jl [Nx Ny Nz] [spinup_minutes] [movie_seconds] [N_droplets] [prefix]
using TurbulentActivation, Breeze, Oceananigans, CUDA, JLD2, Printf, Statistics, Random
using Oceananigans.Units

Nx, Ny, Nz = length(ARGS) ≥ 3 ? parse.(Int, ARGS[1:3]) : (64, 64, 32)
spinup_minutes = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 15.0
movie_seconds = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 60.0
N = length(ARGS) ≥ 6 ? parse(Int, ARGS[6]) : 5000
prefix = length(ARGS) ≥ 7 ? ARGS[7] : "movie"
frame_interval = parse(Float64, get(ENV, "FRAME_INTERVAL", "0.25"))
arch = CUDA.functional() ? GPU() : CPU()

chamber = PiChamber()
grid = pi_chamber_grid(chamber, arch; size=(Nx, Ny, Nz))
aerosol = anderson_aerosol()
temperature = (chamber.bottom_temperature + chamber.top_temperature) / 2
relative_humidity = 0.8
rng = MersenneTwister(1234)
droplets = seed_droplets(aerosol, grid, N; temperature, relative_humidity, rng)
particles = LagrangianParticles(droplets; dynamics=DropletDynamics())
microphysics = get(ENV, "HOST", "bulk") == "bulk" ? chamber_microphysics(; relaxation_time=parse(Float64, get(ENV, "HOST_TAU", "28"))) : nothing
model = pi_chamber_model(chamber, grid; particles, microphysics)
initialize_chamber!(model, chamber; temperature, relative_humidity, rng)
@info "Chamber movie" chamber arch size=(Nx, Ny, Nz) spinup_minutes movie_seconds N frame_interval

u, v, w = model.velocities
ℋ = RelativeHumidityField(model)
speed = Field(sqrt(u^2 + v^2 + w^2))

function progress(sim)
    compute!(ℋ)
    @printf("iter %6d  t = %7.2f s  Δt = %.3f  max|w| = %.3f  ⟨ℋ⟩ = %.4f  wall = %s\n",
            iteration(sim), time(sim), sim.Δt, maximum(abs, w), mean(ℋ), prettytime(sim.run_wall_time))
end

simulation = Simulation(model; Δt=0.02, stop_time=spinup_minutes * minutes)
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
conjure_time_step_wizard!(simulation; cfl=0.7, max_Δt=0.05)
add_callback!(simulation, progress, TimeInterval(30))
run!(simulation)

# A vertical slice through the middle of the box, and the droplets within one cell of it
j = Ny ÷ 2
slices = (; speed = Field(speed; indices=(:, j, :)),
            w = Field(w; indices=(:, j, :)),
            𝒮 = Field(ℋ - 1; indices=(:, j, :)))
simulation.output_writers[:slices] = JLD2Writer(model, slices; filename="$(prefix)_slices.jld2",
                                                schedule=TimeInterval(frame_interval), overwrite_existing=true,
                                                with_halos=false)
simulation.output_writers[:droplets] = JLD2Writer(model, (; droplets=model.particles); filename="$(prefix)_droplets.jld2",
                                                  schedule=TimeInterval(frame_interval), overwrite_existing=true)
simulation.stop_time = time(simulation) + movie_seconds
run!(simulation)
println("MOVIE OUTPUT OK")
