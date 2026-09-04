# Spin up the Pi Chamber without droplets and write mean profiles and fields.
# Usage: julia --project scripts/run_chamber.jl [Nx Ny Nz] [stop_minutes] [prefix]
using TurbulentActivation, Breeze, Oceananigans, CUDA, Printf, Statistics
using Oceananigans.Units

Nx, Ny, Nz = length(ARGS) ≥ 3 ? parse.(Int, ARGS[1:3]) : (32, 32, 16)
stop_minutes = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 0.5
prefix = length(ARGS) ≥ 5 ? ARGS[5] : "chamber"
arch = CUDA.functional() ? GPU() : CPU()

chamber = PiChamber()
grid = pi_chamber_grid(chamber, arch; size=(Nx, Ny, Nz))
model = pi_chamber_model(chamber, grid)
initialize_chamber!(model, chamber)
@info "Pi Chamber" chamber arch size=(Nx, Ny, Nz) stop_minutes

simulation = Simulation(model; Δt=0.02, stop_time=stop_minutes * minutes)
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)

ℋ = RelativeHumidityField(model)
T = model.temperature
u, v, w = model.velocities

progress(sim) = (compute!(ℋ);
    @printf("iter %6d  t = %7.2f s  max|w| = %.3f m/s  ⟨T⟩ = %.2f K  ⟨ℋ⟩ = %.3f  max ℋ = %.3f  wall = %s\n",
            iteration(sim), time(sim), maximum(abs, w), mean(T), mean(ℋ), maximum(ℋ), prettytime(sim.run_wall_time)))
add_callback!(simulation, progress, TimeInterval(10))

profiles = (; T=Average(T, dims=(1, 2)), ℋ=Average(ℋ, dims=(1, 2)),
              qᵛ=Average(model.microphysical_fields.qᵛ, dims=(1, 2)),
              uu=Average(@at((Center, Center, Center), u^2), dims=(1, 2)),
              vv=Average(@at((Center, Center, Center), v^2), dims=(1, 2)),
              ww=Average(@at((Center, Center, Center), w^2), dims=(1, 2)),
              wT=Average(@at((Center, Center, Center), w * T), dims=(1, 2)))
simulation.output_writers[:profiles] = JLD2Writer(model, profiles; filename="$(prefix)_profiles.jld2",
                                                  schedule=TimeInterval(5), overwrite_existing=true)
simulation.output_writers[:fields] = JLD2Writer(model, (; T, ℋ, u, v, w); filename="$(prefix)_fields.jld2",
                                                schedule=TimeInterval(1minute), overwrite_existing=true)
run!(simulation)
println("CHAMBER OK")
