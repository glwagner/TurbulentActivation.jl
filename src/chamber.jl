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
- `transfer_coefficient`: bulk transfer coefficient of the wall laws (2e-2). Calibrated against the
  replay of Anderson's experiment through the reference SAM fields at 64 × 64 × 32: the neutral
  log-law value 6e-3 (first cell centre 1.6 cm from the wall, 0.1 mm roughness) gives a chamber
  whose supersaturation fluctuations are three times too weak, 4e-2 gives fluctuations 50 % too
  strong, and 2e-2 reproduces the reference activation curves to within window scatter
  (`analysis/coefficient_sweep.jl`). The reference itself uses Monin–Obukhov fluxes whose
  magnitude depends on the grid spacing (Yang et al. 2022), so the value belongs to this grid.
- `side_transfer_coefficient`: the coefficient on the four side walls (default: the same). SAM's
  Monin–Obukhov coefficients are enhanced on the unstable floor and ceiling and neutral on the
  side walls, so a smaller side coefficient mimics that partition.
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
    side_transfer_coefficient :: FT
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
                   transfer_coefficient = 2e-2,
                   side_transfer_coefficient = transfer_coefficient,
                   advection_order = 9)

    return PiChamber{FT}(FT.(extent), bottom_temperature, top_temperature, side_temperature,
                         side_relative_humidity, surface_pressure, reference_potential_temperature,
                         transfer_coefficient, side_transfer_coefficient, advection_order)
end

Base.summary(chamber::PiChamber{FT}) where FT =
    string("PiChamber{", FT, "}(", join(chamber.extent, " × "), " m, floor ", chamber.bottom_temperature,
           " K, ceiling ", chamber.top_temperature, " K, walls ", chamber.side_temperature, " K at ",
           100 * chamber.side_relative_humidity, " %, C = ", chamber.transfer_coefficient, " / ",
           chamber.side_transfer_coefficient, ")")

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
    # Advection undershoots leave slightly negative cloud liquid, which the Stokes fall-velocity power
    # law cannot take; borrow it back from vapor at the same level before the auxiliaries are computed.
    negative_moisture_correction = SpeciesBorrowing()
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
