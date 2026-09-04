# Anderson's activation experiment, online: spin the chamber up, then run 60 s sampling
# windows. At the start of each window the droplets are redistributed and returned to
# haze equilibrium; at its end the activated fraction of every replica (fluctuating,
# uniform, instantaneous) is recorded for each target mean supersaturation.
#
# Usage: julia --project scripts/run_windows.jl [Nx Ny Nz] [spinup_minutes] [n_windows] [N_droplets] [prefix]
using TurbulentActivation, Breeze, Oceananigans, CUDA, JLD2, Printf, Statistics, Random
using Oceananigans.Units

Nx, Ny, Nz = length(ARGS) ≥ 3 ? parse.(Int, ARGS[1:3]) : (32, 32, 16)
spinup_minutes = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 1.0
n_windows = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 2
N = length(ARGS) ≥ 6 ? parse(Int, ARGS[6]) : 1000
prefix = length(ARGS) ≥ 7 ? ARGS[7] : "windows"
window = 60.0
arch = CUDA.functional() ? GPU() : CPU()

chamber = PiChamber()
grid = pi_chamber_grid(chamber, arch; size=(Nx, Ny, Nz))
aerosol = anderson_aerosol()
targets = anderson_targets()
relative_humidity = 0.8
temperature = (chamber.bottom_temperature + chamber.top_temperature) / 2
rng = MersenneTwister(1234)
droplets = seed_replica_droplets(aerosol, grid, N, targets; temperature, supersaturation=relative_humidity - 1, rng)
dynamics = ReplicaDynamics(targets)
particles = LagrangianParticles(droplets; dynamics)
# HOST=bulk runs the cloudy chamber with the warm-only one-moment host (relaxation time HOST_TAU s)
host = get(ENV, "HOST", "none")
microphysics = host == "bulk" ? chamber_microphysics(; relaxation_time=parse(Float64, get(ENV, "HOST_TAU", "5"))) : nothing
model = pi_chamber_model(chamber, grid; particles, microphysics)
initialize_chamber!(model, chamber; temperature, relative_humidity)
@info "Anderson windows" chamber aerosol dynamics host arch size=(Nx, Ny, Nz) spinup_minutes n_windows N

ℋ = RelativeHumidityField(model)
u, v, w = model.velocities

function progress(sim)
    compute!(ℋ)
    s = droplet_statistics(droplets)
    @printf("iter %6d  t = %7.2f s  max|w| = %.3f  ⟨ℋ⟩ = %.3f  droplets: ⟨𝒮⟩ = %+.4f  σ(𝒮) = %.4f  active = %.3f  wall = %s\n",
            iteration(sim), time(sim), maximum(abs, w), mean(ℋ),
            s.mean_supersaturation, s.std_supersaturation, s.activated_fraction, prettytime(sim.run_wall_time))
end

# Spin-up
simulation = Simulation(model; Δt=0.02, stop_time=spinup_minutes * minutes)
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
add_callback!(simulation, progress, TimeInterval(30))
simulation.output_writers[:profiles] = JLD2Writer(model, (; T=Average(model.temperature, dims=(1, 2)), ℋ=Average(ℋ, dims=(1, 2)),
                                                          ww=Average(@at((Center, Center, Center), w^2), dims=(1, 2)));
                                                  filename="$(prefix)_profiles.jld2", schedule=TimeInterval(5), overwrite_existing=true)
run!(simulation)

# Windows
results = Dict{String, Any}("targets" => collect(targets), "window" => window, "N" => N, "size" => (Nx, Ny, Nz))
series = Dict{String, Any}()
for n in 1:n_windows
    t₀ = time(simulation)
    𝒮₀ = mean(Array(droplets.𝒮))
    reset_window!(droplets, grid, dynamics; supersaturation=𝒮₀, rng)
    @info "Window $n starts at t = $(t₀) s with ⟨𝒮⟩ = $(𝒮₀)"
    times = Float64[]; fluc = Vector{Vector{Float64}}(); unif = Vector{Vector{Float64}}(); inst = Vector{Vector{Float64}}()
    record(sim) = (a = replica_activation(droplets); push!(times, time(sim) - t₀);
                   push!(fluc, collect(a.fluctuating)); push!(unif, collect(a.uniform)); push!(inst, collect(a.instantaneous)))
    add_callback!(simulation, record, TimeInterval(1), name=:record)
    simulation.stop_time = t₀ + window
    run!(simulation)
    delete!(simulation.callbacks, :record)
    a = replica_activation(droplets)
    results["window_$n"] = Dict("start" => t₀, "S0" => 𝒮₀, "fluctuating" => collect(a.fluctuating),
                                "uniform" => collect(a.uniform), "instantaneous" => collect(a.instantaneous))
    series["window_$n"] = Dict("times" => times, "fluctuating" => reduce(hcat, fluc), "uniform" => reduce(hcat, unif), "instantaneous" => reduce(hcat, inst))
    @printf("window %d done: target  fluctuating  uniform  instantaneous\n", n)
    for (m, t) in enumerate(targets)
        @printf("   %+.3f   %.3f   %.3f   %.3f\n", t, a.fluctuating[m], a.uniform[m], a.instantaneous[m])
    end
end
jldsave("$(prefix)_activation.jld2"; results, series)
println("WINDOWS OK")
