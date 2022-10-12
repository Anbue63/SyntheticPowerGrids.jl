# Predefined Topology and Set-Points

In this example we will show how to use your own predefined topology and set-points for the grid generation.

We begin by loading all relevant packages:
```@julia
using SyntheticNetworks
import SyntheticPowerGrids.get_line_admittance_matrix
using Random
using Graphs

seed = MersenneTwister(42)
```

Then we generate our topology using synthetic networks:
num_nodes = 100
own_graph = generate_graph(RandomPowerGrid(num_nodes, [1, 1/5, 3/10, 1/3, 1/10, 0.0]...)) # Generate embedded graph

And sample random power set-points and assure power balance

```@julia
P_vec = rand(seed, num_nodes) 
P_vec .-= sum(P_vec) / (num_nodes) 
```


e = edges(own_graph.graph)
cables_vec = 3 * ones(Int, length(e))

L = get_effective_distances(own_graph; mean_len_km = 37.12856121212121, shortest_line_km = 0.06) # Effective spacial distances

nodal_parameters = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5, :τ_P => 5.0)
nodal_dynamics = [(1.0, get_DroopControlledInverterApprox, nodal_parameters)]

pg_struct = PGGeneration(num_nodes = num_nodes, cables_vec = cables_vec, nodal_dynamics = nodal_dynamics, P_vec = P_vec, embedded_graph = own_graph, coupling = :predefined, lines = :StaticLine, slack = true)
Y, Y_shunt = get_line_admittance_matrix(pg_struct, L) # Admittances and shunts
x.edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 

pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)