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

# WALL_MODEL=loglaw (default) uses the log-law wall model with roughness length ROUGHNESS (1.4e-3 m)
# and the Monin–Obukhov stability correction unless STABILITY=0; WALL_MODEL=constant uses the bulk
# coefficient WALL_C (2e-2) with WALL_C_SIDE on the side walls. DT is the time step (0.02 s).
if get(ENV, "WALL_MODEL", "loglaw") == "loglaw"
    coefficient = log_law_coefficient(; roughness_length=parse(Float64, get(ENV, "ROUGHNESS", "1.4e-3")),
                                        stability=get(ENV, "STABILITY", "1") == "1")
    chamber = PiChamber(; transfer_coefficient=coefficient, side_transfer_coefficient=coefficient)
else
    transfer_coefficient = parse(Float64, get(ENV, "WALL_C", "2e-2"))
    side_transfer_coefficient = parse(Float64, get(ENV, "WALL_C_SIDE", string(transfer_coefficient)))
    chamber = PiChamber(; transfer_coefficient, side_transfer_coefficient)
end
Δt = parse(Float64, get(ENV, "DT", "0.02"))
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
relaxation_time = parse(Float64, get(ENV, "HOST_TAU", "28"))
microphysics = host == "bulk" ? chamber_microphysics(; relaxation_time) :
               host == "twomoment" ? chamber_microphysics(; relaxation_time, scheme=:twomoment,
                                                          aerosol_number=parse(Float64, get(ENV, "AEROSOL_N", "5e6")),
                                                          dry_radius=parse(Float64, get(ENV, "AEROSOL_R", "65e-9")),
                                                          hygroscopicity=parse(Float64, get(ENV, "AEROSOL_KAPPA", "1")),
                                                          geometric_std=parse(Float64, get(ENV, "AEROSOL_SIGMA", "1.5"))) : nothing
model = pi_chamber_model(chamber, grid; particles, microphysics, bounded_moisture=get(ENV, "BOUNDED", "1") == "1")
initialize_chamber!(model, chamber; temperature, relative_humidity)
@info "Anderson windows" chamber aerosol dynamics host arch size=(Nx, Ny, Nz) spinup_minutes n_windows N

ℋ = RelativeHumidityField(model)
u, v, w = model.velocities

function progress(sim)
    compute!(ℋ)
    s = droplet_statistics(droplets)
    @printf("iter %6d  t = %7.2f s  Δt = %.3f  max|w| = %.3f  ⟨ℋ⟩ = %.3f  droplets: ⟨𝒮⟩ = %+.4f  σ(𝒮) = %.4f  active = %.3f  wall = %s\n",
            iteration(sim), time(sim), sim.Δt, maximum(abs, w), mean(ℋ),
            s.mean_supersaturation, s.std_supersaturation, s.activated_fraction, prettytime(sim.run_wall_time))
end

# Spin-up
simulation = Simulation(model; Δt, stop_time=spinup_minutes * minutes)
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
# Adaptive time step: CFL (summed over the three directions, 0.7) and MAX_DT (0.05 s); CFL=0 keeps Δt fixed
cfl = parse(Float64, get(ENV, "CFL", "0.7"))
cfl > 0 && conjure_time_step_wizard!(simulation; cfl, max_Δt=parse(Float64, get(ENV, "MAX_DT", "0.05")))
add_callback!(simulation, progress, TimeInterval(30))
# DIAG=1 prints the extrema of the state every second, to catch a blow-up before it throws
if get(ENV, "DIAG", "0") == "1"
    T = model.temperature
    qᶜˡ = get(model.microphysical_fields, :qᶜˡ, nothing)
    diag(sim) = @printf("t = %7.2f  T ∈ [%.2f, %.2f]  max|u,v,w| = %.3f %.3f %.3f  qᶜˡ ∈ [%.2e, %.2e]  min qᵛ = %.2e\n",
                        time(sim), minimum(T), maximum(T), maximum(abs, u), maximum(abs, v), maximum(abs, w),
                        isnothing(qᶜˡ) ? 0 : minimum(qᶜˡ), isnothing(qᶜˡ) ? 0 : maximum(qᶜˡ), minimum(model.microphysical_fields.qᵛ))
    add_callback!(simulation, diag, TimeInterval(1))
    # droplet extrema: a garbage interpolation or a runaway position shows up here first
    function droplet_diag(sim)
        x, y, z = Array(droplets.x), Array(droplets.y), Array(droplets.z)
        Td, pd, D² = Array(droplets.T), Array(droplets.p), Array(droplets.D²)
        @printf("      droplets: x ∈ [%.3f, %.3f]  y ∈ [%.3f, %.3f]  z ∈ [%.3f, %.3f]  T ∈ [%.2f, %.2f]  p ∈ [%.0f, %.0f]  D² ∈ [%.2e, %.2e]  NaN: %d\n",
                minimum(x), maximum(x), minimum(y), maximum(y), minimum(z), maximum(z), minimum(Td), maximum(Td),
                minimum(pd), maximum(pd), minimum(D²), maximum(D²), count(isnan, x) + count(isnan, Td) + count(isnan, D²))
    end
    add_callback!(simulation, droplet_diag, TimeInterval(1))
end
T = model.temperature
qᵛ = model.microphysical_fields.qᵛ
simulation.output_writers[:profiles] = JLD2Writer(model, (; T=Average(T, dims=(1, 2)), TT=Average(T^2, dims=(1, 2)),
                                                          qᵛ=Average(qᵛ, dims=(1, 2)), qq=Average(qᵛ^2, dims=(1, 2)),
                                                          ℋ=Average(ℋ, dims=(1, 2)), ℋℋ=Average(ℋ^2, dims=(1, 2)),
                                                          ww=Average(@at((Center, Center, Center), w^2), dims=(1, 2)),
                                                          uu=Average(@at((Center, Center, Center), u^2), dims=(1, 2)),
                                                          vv=Average(@at((Center, Center, Center), v^2), dims=(1, 2)));
                                                  filename="$(prefix)_profiles.jld2", schedule=TimeInterval(5), overwrite_existing=true)
run!(simulation)
state = chamber_statistics(model)
print_chamber_statistics(state)

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
    # Lagrangian supersaturation and diameter of the first `n_sampled` droplets every 0.5 s, for the
    # supersaturation PDF, its autocorrelation time, and example trajectories
    n_sampled = min(N, 500)
    sample_times = Float64[]; 𝒮_samples = Vector{Vector{Float64}}(); D_samples = Vector{Vector{Float64}}()
    sample(sim) = (push!(sample_times, time(sim) - t₀); push!(𝒮_samples, Array(droplets.𝒮)[1:n_sampled]);
                   push!(D_samples, sqrt.(Array(droplets.D²)[1:n_sampled])))
    add_callback!(simulation, sample, TimeInterval(0.5), name=:sample)
    simulation.stop_time = t₀ + window
    run!(simulation)
    delete!(simulation.callbacks, :record)
    delete!(simulation.callbacks, :sample)
    a = replica_activation(droplets)
    state = chamber_statistics(model)
    print_chamber_statistics(state)
    results["window_$n"] = Dict("start" => t₀, "S0" => 𝒮₀, "fluctuating" => collect(a.fluctuating),
                                "uniform" => collect(a.uniform), "instantaneous" => collect(a.instantaneous),
                                "state" => state)
    series["window_$n"] = Dict("times" => times, "fluctuating" => reduce(hcat, fluc), "uniform" => reduce(hcat, unif), "instantaneous" => reduce(hcat, inst),
                               "sample_times" => sample_times, "S" => reduce(hcat, 𝒮_samples), "D" => reduce(hcat, D_samples))
    @printf("window %d done: target  fluctuating  uniform  instantaneous\n", n)
    for (m, t) in enumerate(targets)
        @printf("   %+.3f   %.3f   %.3f   %.3f\n", t, a.fluctuating[m], a.uniform[m], a.instantaneous[m])
    end
end
jldsave("$(prefix)_activation.jld2"; results, series)
println("WINDOWS OK")
