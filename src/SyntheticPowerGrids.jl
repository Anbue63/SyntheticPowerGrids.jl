module SyntheticPowerGrids
    using PowerDynamics
    using SyntheticNetworks
    using EmbeddedGraphs
    using OrdinaryDiffEq
    using ForwardDiff
    using Distributions
    using Graphs
    using NetworkDynamics
    using LinearAlgebra
    using ForwardDiff
    using Parameters

    include("pg_struct.jl")
    export PGGeneration
    
    include("line_parameters.jl")
    export get_effective_distances
    
    include("PhaseAmplitudeOscillator.jl")
    export parameter_DroopControlledInverterApprox

    include("pg_generation.jl")
    export random_PD_grid

    include("validators.jl")

    include("operationpoint.jl")

    include("active_power.jl")

    include("line_dynamics.jl")
    
    include("nodal_dynamics.jl")

    include("probabilistic_capacity_expansion.jl")
end
