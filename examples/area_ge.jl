using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))

using Revise
using SyntheticPowerGrids
using PowerDynamics

using OrdinaryDiffEq
using Statistics
import SyntheticPowerGrids.get_effective_distances
import SyntheticPowerGrids.get_line_admittance_matrix
using SyntheticNetworks
using Graphs
using Plots

##
area_ge = 357111 # [km²]
a_ge = sqrt(area_ge) 
dist_to_km = a_ge / 2
num_nodes = 400

##
eg = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph
dist_nodes = weights(eg) # Euclidean distance of the edges in EmbeddedGraphs

# Remove all "unconnected" distances!
dist_nodes_connected = vcat(dist_nodes...)
unconnected_idx = findall(iszero, dist_nodes_connected) # Unconnected nodes have a length of d = 0.0
deleteat!(dist_nodes_connected, unconnected_idx) 

dist_nodes_mean = mean(dist_nodes_connected) # Mean length with respect to connected nodes (otherwise we skew the mean!)

mean_len_km = dist_nodes_mean * dist_to_km

e = edges(eg.graph)
cables_vec = 3 * ones(Int, length(e))
L = get_effective_distances(eg; mean_len_km = mean_len_km, shortest_line_km = 0.06) # Effective spacial distances
Y, Y_shunt = get_line_admittance_matrix(eg; L_matrix = L, cables_vec = cables_vec, num_nodes = num_nodes) # Admittances and shunts

edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 

## Area matches area of Germany
parameters_third_order = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
parameters_droop_controlled = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5.0, :τ_P => 1.0) 

nodal_dynamics = [(1/3, get_ThirdOrderMachineApprox, parameters_third_order), (1/3, get_DroopControlledInverterApprox, parameters_droop_controlled), (1/3, get_PQ, nothing)]

pg_struct = PGGeneration(num_nodes = num_nodes, cables_vec = cables_vec, nodal_dynamics = nodal_dynamics, P0 = 1.3, mean_len_km = mean_len_km, maxiters = 10, embedded_graph = eg, coupling = :predefined, lines = :StaticLine, slack = true, edge_parameters = edge_parameters)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

## Generate grid the usual way!
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, maxiters = 10, lines = :StaticLine)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

mean_len_km2 = pg_struct_new.mean_len_km
eg2 = pg_struct_new.embedded_graph
weights(eg2)

dist_nodes2 = weights(eg2) # Euclidean distance of the edges in EmbeddedGraphs

# Remove all "unconnected" distances!
dist_nodes_connected2 = vcat(dist_nodes2...)
unconnected_idx2 = findall(iszero, dist_nodes_connected2) # Unconnected nodes have a length of d = 0.0
deleteat!(dist_nodes_connected2, unconnected_idx2) 

dist_nodes_mean2 = mean(dist_nodes_connected2)

a = (mean_len_km2 / dist_nodes_mean2) * 2

a^2

a^2/area_ge