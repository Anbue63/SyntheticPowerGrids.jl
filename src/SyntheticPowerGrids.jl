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
    using KernelDensity

    include("pg_struct.jl")
    export PGGeneration3
    
    include("line_parameters.jl")
    export get_effective_distances
    
    include("PhaseAmplitudeOscillator.jl")
    export parameter_schiffer

    include("pg_generation.jl")
    export random_PD_grid

    include("validators.jl")

    include("utils.jl")

    include("power_distribution.jl")

    include("line_dynamics.jl")
    
    include("nodal_dynamics.jl")
end
