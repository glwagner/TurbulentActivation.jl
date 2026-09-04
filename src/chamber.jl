#####
##### The Pi Chamber: a closed box with a warm, water-saturated floor, a cold, water-saturated
##### ceiling, and side walls at an intermediate temperature and humidity (Anderson et al. 2023,
##### after the H-M case of Yang et al. 2022).
#####

"""
    PiChamber(FT = Float64; kwargs...)

The chamber configuration. Lengths in metres, temperatures in kelvin, pressure in pascals.

- `extent`: the (x, y, z) dimensions of the box (default `(2, 2, 1)`)
- `bottom_temperature`, `top_temperature`, `side_temperature`: wall temperatures (299, 280, 285)
- `side_relative_humidity`: relative humidity of the air in contact with the side walls (0.78);
  the floor and ceiling are water-saturated
- `surface_pressure`: chamber pressure (100000, the pressure of the reference SAM fields)
- `reference_potential_temperature`: potential temperature of the anelastic reference state (290)
- `transfer_coefficient`: the wall law, a number (a constant bulk coefficient) or a Breeze
  `PolynomialCoefficient`. The default is [`log_law_coefficient`](@ref) with a 1.25 mm roughness
  length: the neutral log-law coefficient at the first cell centre (2.5e-2 for 3.1 cm cells)
  with the Monin–Obukhov stability enhancement on the floor and ceiling and neutral side
  walls. Calibrated against the replay of Anderson's experiment through the reference SAM
  fields (`analysis/`): a constant 6e-3 gives fluctuations three times too weak, a constant
  2e-2 reproduces the fluctuations but leaves the mean supersaturation near zero, and with
  bounds-preserving WENO(5) the 1.25 mm log law reproduces the interior and Lagrangian spread,
  the correlation time, and the activation curve, with a 1 K warm bias of the mean temperature
  and a mean supersaturation of +0.2 % against +0.7 %.
- `side_transfer_coefficient`: the coefficient on the four side walls (default: the same). SAM's
  Monin–Obukhov coefficients are enhanced on the unstable floor and ceiling and neutral on the
  side walls, so a smaller side coefficient mimics that partition.
- `advection_order`: order of the WENO advection of the implicit LES (5); momentum and ρθ use
  plain WENO, the moisture species bounds-preserving WENO with bounds (0, 1)
"""
struct PiChamber{FT, C, CS}
    extent :: NTuple{3, FT}
    bottom_temperature :: FT
    top_temperature :: FT
    side_temperature :: FT
    side_relative_humidity :: FT
    surface_pressure :: FT
    reference_potential_temperature :: FT
    transfer_coefficient :: C
    side_transfer_coefficient :: CS
    advection_order :: Int
end

function PiChamber(FT = Float64;
                   extent = (2, 2, 1),
                   bottom_temperature = 299,
                   top_temperature = 280,
                   side_temperature = 285,
                   side_relative_humidity = 0.78,
                   surface_pressure = 100000,
                   reference_potential_temperature = 290,
                   transfer_coefficient = log_law_coefficient(FT),
                   side_transfer_coefficient = transfer_coefficient,
                   advection_order = 5)

    C = transfer_coefficient isa Number ? FT(transfer_coefficient) : transfer_coefficient
    Cₛ = side_transfer_coefficient isa Number ? FT(side_transfer_coefficient) : side_transfer_coefficient
    return PiChamber{FT, typeof(C), typeof(Cₛ)}(FT.(extent), bottom_temperature, top_temperature, side_temperature,
                                                side_relative_humidity, surface_pressure, reference_potential_temperature,
                                                C, Cₛ, advection_order)
end

coefficient_summary(C::Number) = string(C)
coefficient_summary(C) = string(nameof(typeof(C)), "(ℓ = ", C.roughness_length, " m)")

Base.summary(chamber::PiChamber{FT}) where FT =
    string("PiChamber{", FT, "}(", join(chamber.extent, " × "), " m, floor ", chamber.bottom_temperature,
           " K, ceiling ", chamber.top_temperature, " K, walls ", chamber.side_temperature, " K at ",
           100 * chamber.side_relative_humidity, " %, C = ", coefficient_summary(chamber.transfer_coefficient), " / ",
           coefficient_summary(chamber.side_transfer_coefficient), ")")

"""
    log_law_coefficient(FT = Float64; roughness_length = 1.25e-3, scalar_roughness_length = roughness_length / 7.3,
                        stability = true, minimum_wind_speed = 0.05)

A wall law in the form of Breeze's `PolynomialCoefficient`: the neutral coefficient is the
log law `κ² / ln(h / ℓ)²` at the wall distance `h` of the first cell centre (entered through
the 10 m polynomial as a constant, which Breeze transfers to `h` with the same log law), and
on the floor and ceiling the Monin–Obukhov stability correction of `FittedStabilityFunction`
enhances it where the wall layer is unstable, while the side walls stay neutral. With
`ℓ = 1.25 mm` and 3.1 cm cells the neutral value is 2.5e-2; this roughness reproduces the
reference SAM fields' supersaturation statistics and activation curve with WENO(5).
"""
function log_law_coefficient(FT = Float64; roughness_length = 1.25e-3, scalar_roughness_length = roughness_length / 7.3,
                             stability = true, minimum_wind_speed = 0.05)
    ℓ = FT(roughness_length)
    κ = FT(0.4)
    a₀ = 1000 * κ^2 / log(10 / ℓ)^2          # neutral_coefficient_10m = (a₀ + a₁ U + a₂ / U) × 1e-3
    stability_function = stability ? FittedStabilityFunction(FT(scalar_roughness_length)) : nothing
    return PolynomialCoefficient(FT; polynomial=(a₀, zero(FT), zero(FT)), roughness_length=ℓ,
                                 minimum_wind_speed=FT(minimum_wind_speed), stability_function)
end

Base.show(io::IO, chamber::PiChamber) = print(io, summary(chamber))

"""
    pi_chamber_grid(chamber, arch; size)

An all-`Bounded` `RectilinearGrid` filling the chamber, centred horizontally and with the floor at z = 0.
The halo accommodates WENO of the chamber's `advection_order`.
"""
function pi_chamber_grid(chamber::PiChamber{FT}, arch; size) where FT
    Lx, Ly, Lz = chamber.extent
    halo = (chamber.advection_order + 1) ÷ 2
    return RectilinearGrid(arch, FT; size, halo=(halo, halo, halo),
                           x=(-Lx/2, Lx/2), y=(-Ly/2, Ly/2), z=(0, Lz),
                           topology=(Bounded, Bounded, Bounded))
end

"""
    pi_chamber_boundary_conditions(chamber)

Bulk drag on the tangential momentum and bulk heat and vapor fluxes on all six walls, from the
wall temperatures and humidities of the chamber.
"""
function pi_chamber_boundary_conditions(chamber::PiChamber)
    C = chamber.transfer_coefficient
    Cₛ = chamber.side_transfer_coefficient
    T_bottom = chamber.bottom_temperature
    T_top = chamber.top_temperature
    T_side = chamber.side_temperature
    ℋ_side = chamber.side_relative_humidity

    drag(T, C) = BulkDrag(coefficient=C, surface_temperature=T)
    heat(T, C) = BulkSensibleHeatFlux(coefficient=C, surface_temperature=T)
    vapor(T, ℋ, C) = BulkVaporFlux(coefficient=C, surface_temperature=T, surface_relative_humidity=ℋ)

    return (;
        ρu = FieldBoundaryConditions(bottom=drag(T_bottom, C), top=drag(T_top, C), south=drag(T_side, Cₛ), north=drag(T_side, Cₛ)),
        ρv = FieldBoundaryConditions(bottom=drag(T_bottom, C), top=drag(T_top, C), west=drag(T_side, Cₛ), east=drag(T_side, Cₛ)),
        ρw = FieldBoundaryConditions(west=drag(T_side, Cₛ), east=drag(T_side, Cₛ), south=drag(T_side, Cₛ), north=drag(T_side, Cₛ)),
        ρθ = FieldBoundaryConditions(bottom=heat(T_bottom, C), top=heat(T_top, C),
                                     west=heat(T_side, Cₛ), east=heat(T_side, Cₛ), south=heat(T_side, Cₛ), north=heat(T_side, Cₛ)),
        ρqᵛ = FieldBoundaryConditions(bottom=vapor(T_bottom, 1, C), top=vapor(T_top, 1, C),
                                      west=vapor(T_side, ℋ_side, Cₛ), east=vapor(T_side, ℋ_side, Cₛ),
                                      south=vapor(T_side, ℋ_side, Cₛ), north=vapor(T_side, ℋ_side, Cₛ)))
end

"""
    pi_chamber_model(chamber, grid; particles=nothing, microphysics=nothing, kwargs...)

An anelastic `AtmosphereModel` of the chamber on `grid`: WENO implicit LES (bounds-preserving on
the moisture species) with no explicit
closure, the six-wall bulk fluxes of the chamber, and optionally Lagrangian `particles`
(a `LagrangianParticles` of droplets) and a host `microphysics`. Extra keyword arguments are
passed to `AtmosphereModel`.
"""
function pi_chamber_model(chamber::PiChamber{FT}, grid; particles=nothing, microphysics=nothing, bounded_moisture=true, kwargs...) where FT
    constants = ThermodynamicConstants(FT)
    reference_state = ReferenceState(grid, constants;
                                     surface_pressure = chamber.surface_pressure,
                                     potential_temperature = chamber.reference_potential_temperature)
    dynamics = AnelasticDynamics(reference_state)
    boundary_conditions = pi_chamber_boundary_conditions(chamber)
    # WENO of the chamber's order on momentum and on ρθ; bounds-preserving WENO on every moisture
    # species (vapor and the host's prognostic condensate), so that advection undershoot cannot
    # produce negative mass fractions
    order = chamber.advection_order
    moisture_names = filter(name -> startswith(string(name), "ρq"), (:ρqᵛ, prognostic_field_names(microphysics)...))
    bounded = bounded_moisture ? WENO(order=order, bounds=(0, 1)) : WENO(order=order)
    scalar_advection = merge((; ρθ = WENO(order=order)),
                             NamedTuple{moisture_names}(ntuple(_ -> bounded, length(moisture_names))))
    return AtmosphereModel(grid; dynamics, boundary_conditions, momentum_advection=WENO(order=order),
                           scalar_advection, particles, microphysics, thermodynamic_constants=constants, kwargs...)
end

"""
    initialize_chamber!(model, chamber; temperature, relative_humidity=0.8, noise=0.05, rng=default_rng())

Set the chamber air to a uniform `temperature` (default: the mean of the floor and ceiling
temperatures) with a weakly unstable lapse of one kelvin over the chamber height, random
temperature `noise`, and uniform `relative_humidity`.
"""
function initialize_chamber!(model, chamber::PiChamber{FT};
                             temperature = (chamber.bottom_temperature + chamber.top_temperature) / 2,
                             relative_humidity = 0.8, noise = 0.05, rng = default_rng()) where FT
    Lz = chamber.extent[3]
    Tᵢ(x, y, z) = temperature + (1 - z / Lz) + noise * randn(rng)
    set!(model; T=Tᵢ, ℋ=FT(relative_humidity))
    return nothing
end

#####
##### Host microphysics
#####

"""
    chamber_microphysics(FT = Float64; relaxation_time = 28, scheme = :onemoment,
                         aerosol_number = 5e6, dry_radius = 65e-9, hygroscopicity = 1, geometric_std = 1.5)

A warm-only bulk host microphysics for the cloudy chamber. `scheme = :onemoment` is Breeze's
one-moment scheme with every rain process switched off; `scheme = :twomoment` is the
Seifert–Beheng two-moment scheme with Abdul-Razzak–Ghan activation of a single κ-Köhler
aerosol mode (`aerosol_number` per m³, `dry_radius`, `geometric_std`, `hygroscopicity`),
prognostic droplet number, and sedimentation. In both, cloud liquid forms by relaxation
toward saturation on `relaxation_time` (seconds; Breeze multiplies it by the latent-heat
factor Γ ≈ 2.5, so 28 s reproduces the 72 s condensation-rate regression of the reference
SAM fields); the droplet number does not enter the condensation rate in either scheme.
"""
function chamber_microphysics(FT = Float64; relaxation_time = 28, scheme = :onemoment,
                              aerosol_number = 5e6, dry_radius = 65e-9, hygroscopicity = 1,
                              geometric_std = 1.5, molar_mass = 0.058)
    ext = Base.get_extension(Breeze, :BreezeCloudMicrophysicsExt)
    cloud_formation = NonEquilibriumCloudFormation(ConstantRateCondensateFormation(FT(1 / relaxation_time)), nothing)
    negative_moisture_correction = SpeciesBorrowing()
    if scheme == :twomoment
        # Seifert–Beheng two-moment warm microphysics with Abdul-Razzak–Ghan activation of a single
        # κ-Köhler aerosol mode (number per m³, dry radius, geometric standard deviation, hygroscopicity)
        CMAM = CloudMicrophysics.AerosolModel
        mode = CMAM.Mode_κ(FT(dry_radius), FT(geometric_std), FT(aerosol_number), (FT(1),), (FT(1),), (FT(molar_mass),), (FT(hygroscopicity),))
        distribution = CMAM.AerosolDistribution((mode,))
        activation = ext.AerosolActivation(AerosolActivationParameters(FT), distribution, FT(1))
        categories = ext.two_moment_cloud_microphysics_categories(FT; aerosol_activation=activation)
        return ext.TwoMomentCloudMicrophysics(FT; cloud_formation, categories, negative_moisture_correction)
    end
    parameters = Microphysics1MParams(FT; rain_autoconversion=nothing, rain_condensation_evaporation=nothing,
                                      cloud_liquid_rain_accretion=nothing, cloud_ice_formation=nothing, cloud_ice_melt=nothing,
                                      snow_autoconversion=nothing, snow_deposition_sublimation=nothing, snow_melt=nothing,
                                      cloud_liquid_snow_accretion=nothing, cloud_ice_rain_accretion=nothing,
                                      cloud_ice_snow_accretion=nothing, rain_snow_accretion=nothing)
    categories = ext.one_moment_cloud_microphysics_categories(FT; parameters)
    return ext.OneMomentCloudMicrophysics(FT; cloud_formation, categories, negative_moisture_correction)
end

#####
##### Chamber statistics: the Eulerian state of the chamber, to compare with the reference LES
#####

wall_cells(A, ::Val{:bottom}) = view(A, :, :, 1)
wall_cells(A, ::Val{:top}) = view(A, :, :, size(A, 3))
wall_cells(A, ::Val{:west}) = view(A, 1, :, :)
wall_cells(A, ::Val{:east}) = view(A, size(A, 1), :, :)
wall_cells(A, ::Val{:south}) = view(A, :, 1, :)
wall_cells(A, ::Val{:north}) = view(A, :, size(A, 2), :)

"""
    chamber_statistics(model)

Volume statistics of the chamber's Eulerian state: the mean and standard deviation of the
temperature `T`, vapor mass fraction `qᵛ`, and supersaturation `𝒮` over the whole volume and
over the interior (two cells away from every wall), the mean `T`, `qᵛ`, and `𝒮` of the cells
adjacent to each of the six walls, and the root-mean-square velocity components. The
supersaturation is the model's own, `𝒮 = ℋ − 1` with `ℋ` the [`RelativeHumidityField`](@ref).
"""
function chamber_statistics(model)
    T = Array(interior(model.temperature))
    qᵛ = Array(interior(model.microphysical_fields.qᵛ))
    ℋ = RelativeHumidityField(model); compute!(ℋ)
    𝒮 = Array(interior(ℋ)) .- 1
    u, v, w = model.velocities
    Nx, Ny, Nz = size(model.grid)
    inner = (3:Nx-2, 3:Ny-2, 3:Nz-2)
    stats(A) = (mean = mean(A), std = std(A))
    walls = (:bottom, :top, :west, :east, :south, :north)
    near_wall = NamedTuple{walls}(map(w -> (T = mean(wall_cells(T, Val(w))), qᵛ = mean(wall_cells(qᵛ, Val(w))), 𝒮 = mean(wall_cells(𝒮, Val(w)))), walls))
    return (; T = stats(T), qᵛ = stats(qᵛ), 𝒮 = stats(𝒮),
              interior = (T = stats(view(T, inner...)), qᵛ = stats(view(qᵛ, inner...)), 𝒮 = stats(view(𝒮, inner...))),
              near_wall,
              rms_velocity = (u = sqrt(mean(Array(interior(u)).^2)), v = sqrt(mean(Array(interior(v)).^2)), w = sqrt(mean(Array(interior(w)).^2))))
end

"""
    print_chamber_statistics([io], statistics)

Print [`chamber_statistics`](@ref) in the layout of the reference table (`reference/README.md`).
"""
function print_chamber_statistics(io::IO, s)
    @printf(io, "chamber state: ⟨T⟩ = %.3f K  σ(T) = %.3f K  ⟨qᵛ⟩ = %.3f g/kg  σ(qᵛ) = %.3f g/kg  ⟨𝒮⟩ = %+.3f %%  σ(𝒮) = %.3f %%\n",
            s.T.mean, s.T.std, 1000s.qᵛ.mean, 1000s.qᵛ.std, 100s.𝒮.mean, 100s.𝒮.std)
    @printf(io, "  interior:    ⟨T⟩ = %.3f K  σ(T) = %.3f K  ⟨qᵛ⟩ = %.3f g/kg  σ(qᵛ) = %.3f g/kg  ⟨𝒮⟩ = %+.3f %%  σ(𝒮) = %.3f %%\n",
            s.interior.T.mean, s.interior.T.std, 1000s.interior.qᵛ.mean, 1000s.interior.qᵛ.std, 100s.interior.𝒮.mean, 100s.interior.𝒮.std)
    for w in keys(s.near_wall)
        n = s.near_wall[w]
        @printf(io, "  %-7s cell: T = %.3f K  qᵛ = %.3f g/kg  𝒮 = %+.3f %%\n", w, n.T, 1000n.qᵛ, 100n.𝒮)
    end
    @printf(io, "  rms velocity: u = %.3f  v = %.3f  w = %.3f m/s\n", s.rms_velocity.u, s.rms_velocity.v, s.rms_velocity.w)
end
print_chamber_statistics(s) = print_chamber_statistics(stdout, s)
