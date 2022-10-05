using Pkg
Pkg.activate(@__DIR__)

using SyntheticPowerGrids
using SyntheticNetworks
import SyntheticPowerGrids.get_line_admittance_matrix
using Random
using Graphs

##
# Generate grid with own topology!
seed = MersenneTwister(42)

num_nodes = 10
own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph
P_vec = rand(seed, num_nodes) # Random Power Distribution
P_vec .-= sum(P_vec) / (num_nodes)  # Assure power balance

e = edges(own_graph.graph)
cables_vec = 3 * ones(Int, length(e))

L = get_effective_distances(own_graph; mean_len_km = 42, shortest_line_km = 0.06) # Effective spacial distances

nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 1.0)
nodal_dynamics = [(1.0, get_DroopControlledInverterApprox, nodal_parameters)]

x = PGGeneration(num_nodes = num_nodes, cables_vec = cables_vec, nodal_dynamics = nodal_dynamics, P_vec = P_vec, embedded_graph = own_graph, coupling = :predefined, lines = :StaticLine, slack = false)
Y, Y_shunt = get_line_admittance_matrix(x, L)                             # Admittances and shunts
x.edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 

pg, op, pg_struct_new, rejections = random_PD_grid(x)