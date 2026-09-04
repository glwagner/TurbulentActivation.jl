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
chamber = haskey(ENV, "WALL_C") ? PiChamber(; transfer_coefficient=parse(Float64, ENV["WALL_C"])) : PiChamber()
Δt = parse(Float64, get(ENV, "DT", "0.02"))
host = get(ENV, "HOST", "none")
microphysics = host == "bulk" ? chamber_microphysics(; relaxation_time=parse(Float64, get(ENV, "HOST_TAU", "5"))) : nothing
grid = pi_chamber_grid(chamber, arch; size=(Nx, Ny, Nz))
aerosol = anderson_aerosol()
relative_humidity = 0.8
temperature = (chamber.bottom_temperature + chamber.top_temperature) / 2
droplets = seed_droplets(aerosol, grid, N; temperature, relative_humidity, rng=MersenneTwister(1234))
# DEBUG_INTERP=1 wraps the droplet dynamics in a check of the interpolated state after every
# step: an unphysical temperature, pressure, or vapor fraction at any droplet is reported with
# the droplet's position and the run stops there.
struct CheckedDynamics{D}
    inner :: D
end
function (d::CheckedDynamics)(particles, model, Δt)
    droplets = particles.properties
    interpolate_to_droplets!(droplets, model)
    Tmin, Tmax = extrema(droplets.T); pmin, pmax = extrema(droplets.p); qmin, qmax = extrema(droplets.qᵛ)
    if !(250 < Tmin && Tmax < 320 && 9e4 < pmin && pmax < 1.1e5 && 0 ≤ qmin && qmax < 0.05)
        Th, ph, qh = Array(droplets.T), Array(droplets.p), Array(droplets.qᵛ)
        x, y, z = Array(droplets.x), Array(droplets.y), Array(droplets.z)
        bad = findall(n -> !(250 < Th[n] < 320 && 9e4 < ph[n] < 1.1e5 && 0 ≤ qh[n] < 0.05), eachindex(Th))
        @printf("INTERPOLATION CHECK FAILED at iteration %d: %d droplets out of range\n", model.clock.iteration, length(bad))
        for n in bad[1:min(end, 5)]
            @printf("  droplet %d at (%.6f, %.6f, %.6f): T = %g  p = %g  qᵛ = %g\n", n, x[n], y[n], z[n], Th[n], ph[n], qh[n])
        end
        Tf = model.temperature
        @printf("  temperature field extrema (interior) = %s; halo min/max = %s\n", extrema(interior(Tf)), extrema(parent(Tf)))
        error("interpolation check failed")
    end
    return d.inner(particles, model, Δt)
end
dynamics = get(ENV, "DEBUG_INTERP", "0") == "1" ? CheckedDynamics(DropletDynamics()) : DropletDynamics()
particles = get(ENV, "PARTICLES", "1") == "1" ? LagrangianParticles(droplets; dynamics) : nothing
model = pi_chamber_model(chamber, grid; particles, microphysics)
initialize_chamber!(model, chamber; temperature, relative_humidity)
@info "Pi Chamber with droplets" chamber aerosol arch size=(Nx, Ny, Nz) stop_minutes N

simulation = Simulation(model; Δt, stop_time=stop_minutes * minutes)
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
# Adaptive time step: CFL (summed over the three directions, 0.7) and MAX_DT (0.05 s); CFL=0 keeps Δt fixed
cfl = parse(Float64, get(ENV, "CFL", "0.7"))
cfl > 0 && conjure_time_step_wizard!(simulation; cfl, max_Δt=parse(Float64, get(ENV, "MAX_DT", "0.05")))

ℋ = RelativeHumidityField(model)
u, v, w = model.velocities

function progress(sim)
    compute!(ℋ)
    s = droplet_statistics(droplets)
    @printf("iter %6d  t = %7.2f s  Δt = %.3f  max|w| = %.3f  ⟨ℋ⟩ = %.3f  max ℋ = %.3f  droplets: ⟨𝒮⟩ = %+.4f  σ(𝒮) = %.4f  ⟨D⟩ = %.2f μm  active = %.3f  wall = %s\n",
            iteration(sim), time(sim), sim.Δt, maximum(abs, w), mean(ℋ), maximum(ℋ),
            s.mean_supersaturation, s.std_supersaturation, 1e6 * s.mean_diameter, s.activated_fraction,
            prettytime(sim.run_wall_time))
end
add_callback!(simulation, progress, TimeInterval(10))

profiles = (; T=Average(model.temperature, dims=(1, 2)), ℋ=Average(ℋ, dims=(1, 2)),
              ww=Average(@at((Center, Center, Center), w^2), dims=(1, 2)))
simulation.output_writers[:profiles] = JLD2Writer(model, profiles; filename="$(prefix)_profiles.jld2",
                                                  schedule=TimeInterval(5), overwrite_existing=true)
if !isnothing(model.particles)
    simulation.output_writers[:particles] = JLD2Writer(model, (; particles=model.particles); filename="$(prefix)_particles.jld2",
                                                       schedule=TimeInterval(0.5), overwrite_existing=true)
end
run!(simulation)
println(droplet_statistics(droplets))
println("DROPLETS OK")
