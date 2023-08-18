using SyntheticPowerGrids

# In this example we will show how to use your own predefined topology, in the form of an [EmbeddedGraph](https://github.com/PIK-ICoNe/EmbeddedGraphs.jl), and set-points for the grid generation.
# We begin by loading all relevant packages:

using SyntheticNetworks
import SyntheticPowerGrids.get_line_admittance_matrix
using Random
using Graphs

##
# Generate grid with own topology!
seed = MersenneTwister(42)

# Then we generate our topology using [SyntheticNetworks.jl](https://github.com/PIK-ICoNe/SyntheticNetworks.jl):
num_nodes = 100
own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph

# And sample random power set-points and assure power balance:
P_vec = rand(seed, num_nodes) 
P_vec .-= sum(P_vec) / (num_nodes)  
Q_vec = zeros(num_nodes)
nodes_load_flow = fill(:PVAlgebraic, num_nodes)
nodes_load_flow[10] = :PQAlgebraic

# Then we assign the typical number of cables to each transmission line and calculate the admittances and shunts.
# Here we use the default algorithms given by the package but generating it another way is possible as well.
e = edges(own_graph.graph)
cables_vec = 3 * ones(Int, length(e))
L = get_effective_distances(own_graph; mean_len_km = 37.12856121212121, shortest_line_km = 0.06) # Effective spacial distances
Y, Y_shunt = get_line_admittance_matrix(own_graph; L_matrix = L, cables_vec = cables_vec, num_nodes = num_nodes) # Admittances and shunts

# Then we define the nodal parameters again, but this time we also use the `edge_parameters` option to give the system the admittances of each transmission line:
edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 
nodal_parameters = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5, :τ_P => 5.0)
nodal_dynamics = [(1.0, get_DroopControlledInverterApprox, nodal_parameters)]

# Then we generate the struct again but use multiple new optional arguments such as `P_vec` for the active power set-points and `embedded_graph` to hand over the topology:
pg_struct = PGGeneration(num_nodes = num_nodes, cables_vec = cables_vec, nodal_dynamics = nodal_dynamics, Q_vec = Q_vec, P_vec = P_vec, embedded_graph = own_graph, coupling = :predefined, lines = :StaticLine, slack = true, edge_parameters = edge_parameters, node_types_ancillary = nodes_load_flow)

# Finally we generate the synthetic network:
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)