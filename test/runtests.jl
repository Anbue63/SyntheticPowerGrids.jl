using SyntheticPowerGrids
using SyntheticNetworks
using Graphs
using Test
using Random 
using PowerDynamics
import SyntheticPowerGrids.get_line_admittance_matrix

@testset "Mixed_Network" begin
    mixed_network_path = joinpath(@__DIR__, "..", "examples", "mixed_network.jl")
    include(mixed_network_path)

    @test length(pg.nodes) == num_nodes                                           # Correct number of nodes has been created
    @test unique(typeof.(pg.lines)) == [PiModelLine]                              # Correct line type
    @test sort(Symbol.(unique(typeof.(pg.nodes)))) == [:NormalForm, :PQAlgebraic] # Correct Node types
end

@testset "Own_Topology" begin
    own_topology_path = joinpath(@__DIR__, "..", "examples", "own_topology.jl")
    include(own_topology_path)

    @test length(pg.nodes) == num_nodes                      # Correct number of nodes has been created
    @test length(edges(own_graph.graph)) == length(pg.lines) # Correct number of edges has been created

    @test isapprox(op[:, :p], P_vec, rtol = 0.1) # Power distribution used correctly, rtol = 0.1 is the the bottleneck from power dynamics, power model uses 0.9 * P and 1.1 * P as the bounds for the power set points

    @test sort(Symbol.(unique(typeof.(pg.nodes)))) == [:NormalForm, :SlackAlgebraic] # Correct Node types 
    @test unique(typeof.(pg.lines)) == [StaticLine] # Correct Line type
end 

@testset "Own_Nodal_Dynamics" begin
    own_dynamics_path = joinpath(@__DIR__, "..", "examples", "own_nodal_dynamics.jl")
    include(own_dynamics_path)

    @test unique(typeof.(pg.nodes)) == [SwingEqLVS] # Only SwingEqLVS should be created
    @test unique(typeof.(pg.lines)) == [StaticLine] # Correct Line type
end