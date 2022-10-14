using SyntheticPowerGrids
using Revise

x_dims = 2 # x_dims > 1 
crand() = complex(-rand(),rand())

parameters_nf = Dict(:Bᵤ => 1im*rand(1,x_dims), :Cᵤ => crand(), :Gᵤ => crand(), :Hᵤ => crand(), 
                     :Bₓ => -rand(x_dims,x_dims), :Cₓ => -rand(x_dims), :Gₓ => -rand(x_dims), :Hₓ => -rand(x_dims),
                     :x_dims => x_dims)  
nodal_dynamics = [(1.0, get_normalform, parameters_nf)]

pg_struct = PGGeneration(num_nodes = 100, nodal_dynamics = nodal_dynamics)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)