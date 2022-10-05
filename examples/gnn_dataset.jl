using Pkg
Pkg.activate(@__DIR__)
using SyntheticPowerGrids

##
nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 5.0) 
nodal_parameters_b = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 1.0) 
nodal_parameters_c = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)

nodal_dynamics = [(1/6, get_DroopControlledInverterApprox, nodal_parameters_a), (1/6, get_DroopControlledInverterApprox, nodal_parameters_b), (1/6, get_DroopControlledInverterApprox, nodal_parameters_c), (0.5, get_PQ, nothing)]

##
a = PGGeneration(num_nodes = 101, nodal_dynamics = nodal_dynamics, lines = :StaticLine)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(a)