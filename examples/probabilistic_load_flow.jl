#using Pkg
#Pkg.activate(@__DIR__)

using Revise
using SyntheticPowerGrids

##
#
num_nodes = 10
own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph

e = edges(own_graph.graph)
cables_vec = ones(Int, length(e))

L = get_effective_distances(own_graph; mean_len_km = 37.12856121212121, shortest_line_km = 0.06) # Effective spacial distances

Y, Y_shunt = get_line_admittance_matrix(own_graph; L_matrix = L, cables_vec = cables_vec, num_nodes = num_nodes) # Admittances and shunts
edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 

nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)
nodal_dynamics =  [(1/2, get_DroopControlledInverterApprox, nodal_parameters), (1/2, get_PQ, nothing)]

##
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, probabilistic_capacity_expansion = true, validators = false, num_tries = 20, edge_parameters = edge_parameters, cables_vec = cables_vec, P0 = 10)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)

unique(pg_struct_updated.cables_vec)