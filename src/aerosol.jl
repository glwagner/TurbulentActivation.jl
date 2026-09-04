#####
##### Anderson's aerosol and the droplets seeded from it
#####

"""
    Aerosol(FT = Float64; dry_diameter=130e-9, hygroscopicity=1.0, density=2160, number_concentration=1e6)

A monodisperse aerosol population. Anderson et al. (2023) use NaCl with a dry diameter of
130 nm; the paper states κ = 1.2 while the released code uses κ = 1.0 (see `reference/audit.md`).
"""
struct Aerosol{FT}
    dry_diameter :: FT
    hygroscopicity :: FT
    density :: FT
    number_concentration :: FT
end

Aerosol(FT = Float64; dry_diameter = 130e-9, hygroscopicity = 1.0, density = 2160, number_concentration = 1e6) =
    Aerosol{FT}(dry_diameter, hygroscopicity, density, number_concentration)

"""
    anderson_aerosol(FT = Float64; hygroscopicity = 1.0)

Anderson et al. (2023)'s NaCl aerosol, with the hygroscopicity of the released code by default.
"""
anderson_aerosol(FT = Float64; hygroscopicity = 1.0) = Aerosol(FT; hygroscopicity)

Base.show(io::IO, a::Aerosol{FT}) where FT =
    print(io, "Aerosol{", FT, "}(Dᵈ = ", 1e9 * a.dry_diameter, " nm, κ = ", a.hygroscopicity,
          ", N = ", a.number_concentration, " m⁻³)")

"""
    seed_droplets(aerosol, grid, N; temperature, relative_humidity, constants, rng=default_rng())

`N` droplets of the `aerosol`, uniformly distributed in the box of `grid`, with wet diameters in
equilibrium with air at `relative_humidity` and `temperature`, and critical diameters evaluated
at `temperature`. Returns a `StructArray{Droplet}` on the grid's architecture, ready for
`LagrangianParticles(droplets; dynamics=DropletDynamics())`.
"""
function seed_droplets(aerosol::Aerosol{FT}, grid, N;
                       temperature, relative_humidity,
                       constants = ThermodynamicConstants(FT),
                       rng::AbstractRNG = default_rng()) where FT

    arch = grid.architecture
    Dᵈ = aerosol.dry_diameter
    κ = aerosol.hygroscopicity
    T = FT(temperature)
    𝒮 = FT(relative_humidity) - 1
    Dᶜ = critical_diameter(Dᵈ, κ, T, constants)
    Dᵢ = equilibrium_diameter(𝒮, Dᵈ, κ, T, constants)

    x₁, x₂ = extrema(xnodes(grid, Face()))
    y₁, y₂ = extrema(ynodes(grid, Face()))
    z₁, z₂ = extrema(znodes(grid, Face()))
    uniform(a, b) = on_architecture(arch, FT.(a .+ (b - a) .* rand(rng, N)))
    column(v) = on_architecture(arch, fill(FT(v), N))

    return StructArray{Droplet{FT}}((uniform(x₁, x₂), uniform(y₁, y₂), uniform(z₁, z₂),
                                     column(Dᵈ), column(κ), column(Dᵢ^2), column(Dᶜ),
                                     column(0), column(0), column(0), column(0)))
end

"""
    droplet_statistics(droplets)

Mean and standard deviation of the supersaturation seen by the droplets, the mean wet diameter,
and the activated fraction.
"""
function droplet_statistics(droplets)
    𝒮 = Array(droplets.𝒮)
    D = sqrt.(Array(droplets.D²))
    return (; mean_supersaturation = mean(𝒮), std_supersaturation = std(𝒮),
              max_supersaturation = maximum(𝒮), mean_diameter = mean(D),
              activated_fraction = activated_fraction(droplets))
end
