#using Pkg
#Pkg.activate(@__DIR__)

using Revise
using SyntheticPowerGrids

##
#
nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)
nodal_dynamics =  [(1/2, get_DroopControlledInverterApprox, nodal_parameters), (1/2, get_PQ, nothing)]

pg_struct = PGGeneration(num_nodes = 10, nodal_dynamics = nodal_dynamics, probabilistic_capacity_expansion = true, P0 = 1000, validators = false, num_tries = 20)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)

unique(pg_struct_updated.cables_vec)