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

    include("line_parameters.jl")

    include("PhaseAmplitudeOscillator.jl")

    include("pg_generation.jl")
    export random_PD_grid

    include("tests.jl")

    include("utils.jl")
end
