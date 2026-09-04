#####
##### Anderson's counterfactual replicas
#####
##### Anderson et al. (2023) shift a common turbulent supersaturation history by a set of
##### target means 𝒮̄ₘ and grow one droplet per target along the same trajectory, once with
##### the fluctuations ("fluctuating", 𝒮̄ₘ + 𝒮′) and once without ("uniform", 𝒮̄ₘ). Here every
##### Lagrangian droplet carries those replicas as fixed-size tuples of squared diameters, so
##### the whole sweep of target means runs online in one simulation. The fluctuation
##### 𝒮′ = 𝒮 − ⟨𝒮⟩ is taken about the instantaneous mean over all droplets (Anderson
##### removes the per-trajectory mean of the finished 60 s window; see reference/audit.md).
##### A third replica set records instantaneous activation, 𝒮̄ₘ + 𝒮′ ≥ 𝒮ᶜ, as a bit mask.
#####

using Breeze.LagrangianMicrophysics: interpolate_to_droplets!, implicit_growth_step, ambient_supersaturation
using Oceananigans.Architectures: architecture
using Oceananigans.Utils: launch!, KernelParameters
using KernelAbstractions: @kernel, @index

"""
    AndersonDroplet{FT, M}

A [`Droplet`](@ref Breeze.LagrangianMicrophysics.Droplet) that also carries, for each of `M`
target mean supersaturations, a fluctuating replica `D²ᶠ`, a uniform replica `D²ᵘ`, and one
bit of `instant` set while the fluctuating supersaturation exceeds the critical value.
"""
struct AndersonDroplet{FT, M}
    x :: FT
    y :: FT
    z :: FT
    Dᵈ :: FT
    κ :: FT
    D² :: FT
    Dᶜ :: FT
    𝒮ᶜ :: FT
    T :: FT
    qᵛ :: FT
    p :: FT
    𝒮 :: FT
    D²ᶠ :: NTuple{M, FT}
    D²ᵘ :: NTuple{M, FT}
    instant :: UInt32
end

"""
    anderson_targets(FT = Float64)

The 19 target mean supersaturations of the released code, `range(-0.04, 0.05, length=19)`.
"""
anderson_targets(FT = Float64) = ntuple(m -> FT(-0.04 + 0.005 * (m - 1)), 19)

"""
    seed_replica_droplets(aerosol, grid, N, targets; temperature, supersaturation, constants, rng)

`N` droplets of the `aerosol` uniformly distributed in the box of `grid`, every wet diameter
(the droplet and all its replicas) at equilibrium with `supersaturation` at `temperature`.
"""
function seed_replica_droplets(aerosol::Aerosol{FT}, grid, N, targets::NTuple{M, FT};
                               temperature, supersaturation,
                               constants = ThermodynamicConstants(FT),
                               rng::AbstractRNG = default_rng()) where {FT, M}
    arch = grid.architecture
    Dᵈ = aerosol.dry_diameter
    κ = aerosol.hygroscopicity
    T = FT(temperature)
    Dᶜ = critical_diameter(Dᵈ, κ, T, constants)
    𝒮ᶜ = critical_supersaturation(Dᵈ, κ, T, constants)
    Dᵢ = equilibrium_diameter(FT(supersaturation), Dᵈ, κ, T, constants)

    x₁, x₂ = extrema(xnodes(grid, Face()))
    y₁, y₂ = extrema(ynodes(grid, Face()))
    z₁, z₂ = extrema(znodes(grid, Face()))
    uniform(a, b) = on_architecture(arch, FT.(a .+ (b - a) .* rand(rng, N)))
    column(v) = on_architecture(arch, fill(v, N))
    replicas = column(ntuple(m -> Dᵢ^2, M))

    return StructArray{AndersonDroplet{FT, M}}((uniform(x₁, x₂), uniform(y₁, y₂), uniform(z₁, z₂),
                                                column(Dᵈ), column(κ), column(Dᵢ^2), column(Dᶜ), column(𝒮ᶜ),
                                                column(zero(FT)), column(zero(FT)), column(zero(FT)), column(zero(FT)),
                                                replicas, copy(replicas), column(UInt32(0))))
end

"""
    reset_window!(droplets, grid; supersaturation, rng)

Start a new sampling window: redistribute the droplets uniformly and return every wet
diameter to equilibrium with `supersaturation` at the droplets' reference temperature
(the temperature at which they were seeded is not stored, so the critical diameter is kept
and the equilibrium is evaluated with the droplets' last interpolated temperature).
"""
function reset_window!(droplets, grid, dynamics; supersaturation, rng::AbstractRNG = default_rng())
    N = length(droplets)
    FT = eltype(droplets.x)
    x₁, x₂ = extrema(xnodes(grid, Face()))
    y₁, y₂ = extrema(ynodes(grid, Face()))
    z₁, z₂ = extrema(znodes(grid, Face()))
    arch = grid.architecture
    copyto!(droplets.x, on_architecture(arch, FT.(x₁ .+ (x₂ - x₁) .* rand(rng, N))))
    copyto!(droplets.y, on_architecture(arch, FT.(y₁ .+ (y₂ - y₁) .* rand(rng, N))))
    copyto!(droplets.z, on_architecture(arch, FT.(z₁ .+ (z₂ - z₁) .* rand(rng, N))))

    constants = dynamics.droplet.thermodynamic_constants
    Dᵈ = first(Array(droplets.Dᵈ))
    κ = first(Array(droplets.κ))
    T = mean(Array(droplets.T))
    T = iszero(T) ? FT(289) : T
    Dᵢ = equilibrium_diameter(FT(supersaturation), Dᵈ, κ, T, constants)
    M = length(dynamics.targets)
    fill!(droplets.D², Dᵢ^2)
    fill!(droplets.D²ᶠ, ntuple(m -> Dᵢ^2, M))
    fill!(droplets.D²ᵘ, ntuple(m -> Dᵢ^2, M))
    fill!(droplets.instant, UInt32(0))
    return nothing
end

#####
##### Dynamics
#####

"""
    ReplicaDynamics(targets; droplet = DropletDynamics())

The particle `dynamics` that grows an [`AndersonDroplet`](@ref) and all its replicas. `targets`
is the tuple of target mean supersaturations; `droplet` carries the accommodation
coefficients, Newton iterations, substeps, and thermodynamic constants of the growth law.
"""
struct ReplicaDynamics{FT, M, D}
    targets :: NTuple{M, FT}
    droplet :: D
end

ReplicaDynamics(targets::NTuple{M, FT}; droplet = DropletDynamics(FT)) where {FT, M} =
    ReplicaDynamics{FT, M, typeof(droplet)}(targets, droplet)

Base.show(io::IO, d::ReplicaDynamics{FT, M}) where {FT, M} =
    print(io, "ReplicaDynamics{", FT, "}(", M, " targets from ", first(d.targets), " to ", last(d.targets), ")")

function (dynamics::ReplicaDynamics)(particles, model, Δt)
    droplets = particles.properties
    interpolate_to_droplets!(droplets, model)

    grid = model.grid
    arch = architecture(grid)
    launch!(arch, grid, KernelParameters(1:length(droplets)), _replica_supersaturation!, droplets, dynamics.droplet)
    𝒮̄ = mean(droplets.𝒮)
    launch!(arch, grid, KernelParameters(1:length(droplets)), _grow_replicas!, droplets, dynamics, 𝒮̄, Δt)
    return nothing
end

@kernel function _replica_supersaturation!(droplets, parameters)
    n = @index(Global)
    @inbounds droplets.𝒮[n] = ambient_supersaturation(droplets.T[n], droplets.qᵛ[n], droplets.p[n],
                                                      parameters.thermodynamic_constants)
end

@inline function grow(D², 𝒮, T, p, Dᵈ, κ, Δt, parameters, constants)
    δt = Δt / parameters.substeps
    for _ in 1:parameters.substeps
        D² = implicit_growth_step(D², 𝒮, T, p, Dᵈ, κ, δt, parameters, constants)
    end
    return D²
end

@kernel function _grow_replicas!(droplets, dynamics, 𝒮̄, Δt)
    n = @index(Global)
    parameters = dynamics.droplet
    constants = parameters.thermodynamic_constants
    targets = dynamics.targets
    M = length(targets)

    @inbounds begin
        T = droplets.T[n]
        p = droplets.p[n]
        Dᵈ = droplets.Dᵈ[n]
        κ = droplets.κ[n]
        𝒮 = droplets.𝒮[n]
        𝒮ᶜ = droplets.𝒮ᶜ[n]
        D²ᶠ = droplets.D²ᶠ[n]
        D²ᵘ = droplets.D²ᵘ[n]
    end

    𝒮′ = 𝒮 - 𝒮̄

    @inbounds droplets.D²[n] = grow(droplets.D²[n], 𝒮, T, p, Dᵈ, κ, Δt, parameters, constants)

    D²ᶠ⁺ = ntuple(m -> grow(D²ᶠ[m], targets[m] + 𝒮′, T, p, Dᵈ, κ, Δt, parameters, constants), Val(M))
    D²ᵘ⁺ = ntuple(m -> grow(D²ᵘ[m], targets[m], T, p, Dᵈ, κ, Δt, parameters, constants), Val(M))

    instant = UInt32(0)
    for m in 1:M
        instant |= UInt32(targets[m] + 𝒮′ ≥ 𝒮ᶜ) << (m - 1)
    end

    @inbounds begin
        droplets.D²ᶠ[n] = D²ᶠ⁺
        droplets.D²ᵘ[n] = D²ᵘ⁺
        droplets.instant[n] = instant
    end
end

#####
##### Window statistics
#####

"""
    replica_activation(droplets)

The activated fraction for each target mean: `fluctuating` and `uniform` from the replica
diameters against the critical diameter, and `instantaneous` from the activation bits.
"""
function replica_activation(droplets)
    D²ᶠ = Array(droplets.D²ᶠ)
    D²ᵘ = Array(droplets.D²ᵘ)
    Dᶜ² = Array(droplets.Dᶜ) .^ 2
    instant = Array(droplets.instant)
    N = length(Dᶜ²)
    M = length(first(D²ᶠ))
    fluctuating = ntuple(m -> count(n -> D²ᶠ[n][m] ≥ Dᶜ²[n], 1:N) / N, M)
    uniform = ntuple(m -> count(n -> D²ᵘ[n][m] ≥ Dᶜ²[n], 1:N) / N, M)
    instantaneous = ntuple(m -> count(n -> (instant[n] >> (m - 1)) & 1 == 1, 1:N) / N, M)
    return (; fluctuating, uniform, instantaneous)
end
