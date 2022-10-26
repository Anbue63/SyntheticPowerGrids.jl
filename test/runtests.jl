using SyntheticPowerGrids
using SyntheticNetworks
using Graphs
using Test
using Random 
using PowerDynamics
import SyntheticPowerGrids.get_line_admittance_matrix

@testset "Mixed_Network" begin
    path = joinpath(@__DIR__, "..", "examples", "mixed_network.jl")
    include(path)

    @test length(pg.nodes) == num_nodes                                           # Correct number of nodes has been created
    @test unique(typeof.(pg.lines)) == [PiModelLine]                              # Correct line type
    @test Set(unique(typeof.(pg.nodes))) == Set([NormalForm{1}, PQAlgebraic, SlackAlgebraic]) # Correct Node types

    file = joinpath(@__DIR__,"grid.json")
    write_powergrid(pg, file, Json)
    pg_read = read_powergrid(file, Json)
    @test pg_read.nodes == pg.nodes # Parsing is working correctly
    rm(file)
end

@testset "Own_Topology" begin
    path = joinpath(@__DIR__, "..", "examples", "own_topology.jl")
    include(path)

    @test length(pg.nodes) == num_nodes                      # Correct number of nodes has been created
    @test length(edges(own_graph.graph)) == length(pg.lines) # Correct number of edges has been created

    @test isapprox(op[:, :p], P_vec, rtol = 0.1) # Power distribution used correctly, rtol = 0.1 is the the bottleneck from power dynamics, power model uses 0.9 * P and 1.1 * P as the bounds for the power set points

    @test Set(unique(typeof.(pg.nodes))) == Set([NormalForm{1}, SlackAlgebraic]) # Correct Node types 
    @test unique(typeof.(pg.lines)) == [StaticLine] # Correct Line type
end 

@testset "Own_Nodal_Dynamics" begin
    path = joinpath(@__DIR__, "..", "examples", "own_nodal_dynamics.jl")
    include(path)

    @test unique(typeof.(pg.nodes)) == [SwingEqLVS] # Only SwingEqLVS should be created
    @test unique(typeof.(pg.lines)) == [StaticLine] # Correct Line type
end

@testset "Capacity_Expansion" begin
    path = joinpath(@__DIR__, "..", "examples", "probabilistic_capacity_expansion.jl")
    include(path)
    
    @test all(unique(pg_struct_updated.cables_vec) .> 0)
    @test length(unique(pg_struct_updated.cables_vec)) > 1 
end