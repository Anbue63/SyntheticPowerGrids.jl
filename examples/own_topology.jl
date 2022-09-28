#using Pkg
#Pkg.activate(@__DIR__)

using Revise
using SyntheticPowerGrids
using SyntheticNetworks
import SyntheticPowerGrids.get_line_admittance_matrix
using Random

##
# Generate grid with own topology!
seed = MersenneTwister(42)
num_nodes = 10
own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...))
P_vec = rand(seed, num_nodes)
P_vec .-= sum(P_vec) / (num_nodes)  # Assure power balance

L = get_effective_distances(own_graph; mean_len_km = 42, shortest_line_km = 0.06)
Y, Y_shunt = get_line_admittance_matrix(own_graph, L)

edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt)
nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0])
nodal_shares = Dict(:DroopControlledInverterApprox_share => 1.0, :load_share => 0.0)

x = PGGeneration(num_nodes = num_nodes, nodal_parameters = nodal_parameters, nodal_shares = nodal_shares, P_vec = P_vec, embedded_graph = own_graph, coupling = :predefined, edge_parameters = edge_parameters, lines = :StaticLine)
pg, op, embedded_graph, rejections = random_PD_grid(x)
