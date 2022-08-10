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
    export PGGeneration2

    include("line_parameters.jl")

    include("PhaseAmplitudeOscillator.jl")

    include("pg_generation.jl")
    export random_PD_grid

    include("tests.jl")

    include("utils.jl")

    include("power_distribution.jl")

    include("line_dynamics.jl")
    
    include("nodal_dynamics.jl")
end
