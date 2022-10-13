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
    using StatsBase

    include("pg_struct.jl")
    export PGGeneration

    include("line_parameters.jl")
    export get_effective_distances
    
    include("PhaseAmplitudeOscillator.jl")
    export parameter_DroopControlledInverterApprox

    include("pg_generation.jl")
    export generate_powergrid_dynamics

    include("validators.jl")

    include("operationpoint.jl")

    include("fluc_node.jl")

    include("active_power.jl")

    include("line_dynamics.jl")
    
    include("nodal_dynamics.jl")
    export get_DroopControlledInverterApprox, get_ThirdOrderMachineApprox, get_PQ, get_dVOCapprox, get_normalform 

    include("probabilistic_capacity_expansion.jl")

end
