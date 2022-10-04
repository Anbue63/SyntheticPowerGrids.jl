using SyntheticPowerGrids
using SyntheticNetworks
using Graphs
using Test
using Random 
using PowerDynamics
import SyntheticPowerGrids.get_line_admittance_matrix

@testset "MixedNetwork" begin
    # Generate a grid with a predefined topology!
    num_nodes = 100

    nodal_parameters = Dict(:X => 1.0, :γ => 0.2, :α => 2.0, :τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) 
    nodal_shares = Dict(:load_share => 0.5, :ThirdOrderMachineApprox_share => 0.25, :DroopControlledInverterApprox_share => 0.25)

    x = PGGeneration(num_nodes = num_nodes, nodal_parameters = nodal_parameters, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :Mixed, nodal_shares = nodal_shares)
    pg, op, embedded_graph, rejections = random_PD_grid(x)

    @test length(pg.nodes) == num_nodes                                           # Correct number of nodes has been created
    @test unique(typeof.(pg.lines)) == [PiModelLine]                              # Correct line type
    @test sort(Symbol.(unique(typeof.(pg.nodes)))) == [:NormalForm, :PQAlgebraic] # Correct Node types
end

@testset "OwnTopology" begin
    # Generate a grid with a predefined topology!
    num_nodes = 10
    own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph
    P_vec = rand(MersenneTwister(42), num_nodes) # Random Power Distribution
    P_vec .-= sum(P_vec) / (num_nodes)  # Assure power balance
    
    e = edges(own_graph.graph)
    cables_vec = 3 * ones(length(e))
    
    L = get_effective_distances(own_graph; mean_len_km = 42, shortest_line_km = 0.06) # Effective spacial distances

    nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0])
    nodal_shares = Dict(:DroopControlledInverterApprox_share => 1.0, :load_share => 0.0)

    x = PGGeneration(num_nodes = num_nodes, cables_vec = cables_vec, nodal_parameters = nodal_parameters, nodal_shares = nodal_shares, P_vec = P_vec, embedded_graph = own_graph, coupling = :predefined, lines = :StaticLine, slack = false)
    Y, Y_shunt = get_line_admittance_matrix(x, L)                             # Admittances and shunts
    x.edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 
    
    
    pg, op, embedded_graph, rejections = random_PD_grid(x)

    @test length(pg.nodes) == num_nodes                      # Correct number of nodes has been created
    @test length(edges(own_graph.graph)) == length(pg.lines) # Correct number of edges has been created

    @test isapprox(op[1:num_nodes - 1, :p], P_vec[1:end - 1], atol = 10^-6) # Power distribution used correctly
    @test isapprox(op[num_nodes, :p], P_vec[end], atol = 0.01)              # Higher tolerance because of the compensation of losses

    @test unique(typeof.(pg.nodes)) == [NormalForm] # Only NormalForms should be created
    @test unique(typeof.(pg.lines)) == [StaticLine] # Correct Line type
end 
