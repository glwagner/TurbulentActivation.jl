module TurbulentActivation

export PiChamber, pi_chamber_grid, pi_chamber_boundary_conditions, pi_chamber_model, initialize_chamber!, chamber_microphysics,
       Aerosol, anderson_aerosol, seed_droplets, droplet_statistics,
       AndersonDroplet, anderson_targets, seed_replica_droplets, reset_window!, ReplicaDynamics, replica_activation

using CloudMicrophysics
using CloudMicrophysics.Parameters: Microphysics1MParams
using Breeze
using Breeze.BoundaryConditions: BulkDrag, BulkSensibleHeatFlux, BulkVaporFlux
using Breeze.Microphysics: ConstantRateCondensateFormation, NonEquilibriumCloudFormation
using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.Grids: xnodes, ynodes, znodes
using Random: AbstractRNG, default_rng
using Statistics: mean, std
using StructArrays: StructArray

include("chamber.jl")
include("aerosol.jl")
include("replicas.jl")

end # module
