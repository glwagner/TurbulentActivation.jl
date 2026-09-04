using TurbulentActivation, Breeze, Oceananigans, CUDA, Random, Test
using Oceananigans.TimeSteppers: step_lagrangian_particles!

arch = CUDA.functional() ? GPU() : CPU()

@testset "TurbulentActivation" begin
    chamber = PiChamber()
    grid = pi_chamber_grid(chamber, arch; size=(8, 8, 8))
    aerosol = anderson_aerosol()
    targets = anderson_targets()
    @test length(targets) == 19
    @test first(targets) ≈ -0.04 && last(targets) ≈ 0.05

    N = 64
    rng = MersenneTwister(1)
    droplets = seed_replica_droplets(aerosol, grid, N, targets; temperature=289, supersaturation=-0.2, rng)
    dynamics = ReplicaDynamics(targets)
    particles = LagrangianParticles(droplets; dynamics)
    model = pi_chamber_model(chamber, grid; particles)
    initialize_chamber!(model, chamber; temperature=289, relative_humidity=0.8, rng)
    @test model.particles === particles

    # Droplets lie inside the chamber and start at equilibrium haze size
    @test all(-1 .≤ Array(droplets.x) .≤ 1) && all(0 .≤ Array(droplets.z) .≤ 1)
    @test all(Array(droplets.D²) .> Array(droplets.Dᵈ) .^ 2)

    # A few steps: the replicas above the critical supersaturation grow beyond the haze branch,
    # those far below stay haze, and every activated fraction is between 0 and 1
    for _ in 1:50
        step_lagrangian_particles!(model, 0.02)
    end
    a = replica_activation(droplets)
    @test all(0 .≤ collect(a.fluctuating) .≤ 1) && all(0 .≤ collect(a.uniform) .≤ 1) && all(0 .≤ collect(a.instantaneous) .≤ 1)
    @test a.uniform[1] == 0                      # −4 %: nothing activates
    @test a.instantaneous[end] == 1              # +5 %: instantaneously above 𝒮ᶜ everywhere
    D²ᵘ = Array(droplets.D²ᵘ)
    @test all(d -> d[end] > d[1], D²ᵘ)          # the +5 % uniform replica is larger than the −4 % one

    # The actual droplets see the ambient (subsaturated) chamber air
    s = droplet_statistics(droplets)
    @test -0.3 < s.mean_supersaturation < 0
    @test s.activated_fraction == 0

    # Resetting a window returns every replica to the same haze size
    reset_window!(droplets, grid, dynamics; supersaturation=-0.2, rng)
    @test all(d -> all(==(d[1]), d), Array(droplets.D²ᶠ))
    @test all(iszero, Array(droplets.instant))
end
