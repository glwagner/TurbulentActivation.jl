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
- `surface_pressure`: chamber pressure (101325)
- `reference_potential_temperature`: potential temperature of the anelastic reference state (290)
- `transfer_coefficient`: bulk transfer coefficient of the wall laws (6e-3, a neutral log law with
  the first cell centre 1.6 cm from the wall and a 0.1 mm roughness length; to be audited)
- `advection_order`: order of the WENO advection of the implicit LES (9)
"""
struct PiChamber{FT}
    extent :: NTuple{3, FT}
    bottom_temperature :: FT
    top_temperature :: FT
    side_temperature :: FT
    side_relative_humidity :: FT
    surface_pressure :: FT
    reference_potential_temperature :: FT
    transfer_coefficient :: FT
    advection_order :: Int
end

function PiChamber(FT = Float64;
                   extent = (2, 2, 1),
                   bottom_temperature = 299,
                   top_temperature = 280,
                   side_temperature = 285,
                   side_relative_humidity = 0.78,
                   surface_pressure = 101325,
                   reference_potential_temperature = 290,
                   transfer_coefficient = 6e-3,
                   advection_order = 9)

    return PiChamber{FT}(FT.(extent), bottom_temperature, top_temperature, side_temperature,
                         side_relative_humidity, surface_pressure, reference_potential_temperature,
                         transfer_coefficient, advection_order)
end

Base.summary(chamber::PiChamber{FT}) where FT =
    string("PiChamber{", FT, "}(", join(chamber.extent, " × "), " m, floor ", chamber.bottom_temperature,
           " K, ceiling ", chamber.top_temperature, " K, walls ", chamber.side_temperature, " K at ",
           100 * chamber.side_relative_humidity, " %)")

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
    T_bottom = chamber.bottom_temperature
    T_top = chamber.top_temperature
    T_side = chamber.side_temperature
    ℋ_side = chamber.side_relative_humidity

    drag(T) = BulkDrag(coefficient=C, surface_temperature=T)
    heat(T) = BulkSensibleHeatFlux(coefficient=C, surface_temperature=T)
    vapor(T, ℋ) = BulkVaporFlux(coefficient=C, surface_temperature=T, surface_relative_humidity=ℋ)

    return (;
        ρu = FieldBoundaryConditions(bottom=drag(T_bottom), top=drag(T_top), south=drag(T_side), north=drag(T_side)),
        ρv = FieldBoundaryConditions(bottom=drag(T_bottom), top=drag(T_top), west=drag(T_side), east=drag(T_side)),
        ρw = FieldBoundaryConditions(west=drag(T_side), east=drag(T_side), south=drag(T_side), north=drag(T_side)),
        ρθ = FieldBoundaryConditions(bottom=heat(T_bottom), top=heat(T_top),
                                     west=heat(T_side), east=heat(T_side), south=heat(T_side), north=heat(T_side)),
        ρqᵛ = FieldBoundaryConditions(bottom=vapor(T_bottom, 1), top=vapor(T_top, 1),
                                      west=vapor(T_side, ℋ_side), east=vapor(T_side, ℋ_side),
                                      south=vapor(T_side, ℋ_side), north=vapor(T_side, ℋ_side)))
end

"""
    pi_chamber_model(chamber, grid; particles=nothing, microphysics=nothing, kwargs...)

An anelastic `AtmosphereModel` of the chamber on `grid`: WENO implicit LES with no explicit
closure, the six-wall bulk fluxes of the chamber, and optionally Lagrangian `particles`
(a `LagrangianParticles` of droplets) and a host `microphysics`. Extra keyword arguments are
passed to `AtmosphereModel`.
"""
function pi_chamber_model(chamber::PiChamber{FT}, grid; particles=nothing, microphysics=nothing, kwargs...) where FT
    constants = ThermodynamicConstants(FT)
    reference_state = ReferenceState(grid, constants;
                                     surface_pressure = chamber.surface_pressure,
                                     potential_temperature = chamber.reference_potential_temperature)
    dynamics = AnelasticDynamics(reference_state)
    boundary_conditions = pi_chamber_boundary_conditions(chamber)
    advection = WENO(order=chamber.advection_order)
    return AtmosphereModel(grid; dynamics, boundary_conditions, advection, particles, microphysics,
                           thermodynamic_constants=constants, kwargs...)
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
    chamber_microphysics(FT = Float64; relaxation_time = 5)

A warm-only bulk host microphysics for the cloudy chamber: Breeze's one-moment scheme with
prognostic vapor, cloud liquid formed by relaxation toward saturation on `relaxation_time`
(seconds), and every rain process switched off, so that the supersaturation is regulated by
condensation on the cloud but no drizzle forms. The relaxation time stands in for the phase
relaxation time of the droplet population, `1 / (4π Dᵛ N r̄)`; it is a sensitivity parameter
until a supersaturation-driven chamber scheme replaces it.
"""
function chamber_microphysics(FT = Float64; relaxation_time = 5)
    ext = Base.get_extension(Breeze, :BreezeCloudMicrophysicsExt)
    parameters = Microphysics1MParams(FT;
                                      rain_autoconversion = nothing,
                                      rain_condensation_evaporation = nothing,
                                      cloud_liquid_rain_accretion = nothing,
                                      cloud_ice_formation = nothing,
                                      cloud_ice_melt = nothing,
                                      snow_autoconversion = nothing,
                                      snow_deposition_sublimation = nothing,
                                      snow_melt = nothing,
                                      cloud_liquid_snow_accretion = nothing,
                                      cloud_ice_rain_accretion = nothing,
                                      cloud_ice_snow_accretion = nothing,
                                      rain_snow_accretion = nothing)
    categories = ext.one_moment_cloud_microphysics_categories(FT; parameters)
    liquid = ConstantRateCondensateFormation(FT(1 / relaxation_time))
    cloud_formation = NonEquilibriumCloudFormation(liquid, nothing)
    return ext.OneMomentCloudMicrophysics(FT; cloud_formation, categories)
end
