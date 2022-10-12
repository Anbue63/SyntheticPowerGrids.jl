#using Pkg
#Pkg.activate(@__DIR__)

using SyntheticPowerGrids

##
parameters_third_order = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
parameters_droop_controlled = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 1.0) 
#parameters_normalform = Dict(:Bᵤ = Bᵤ, :Cᵤ => Cᵤ, :Gᵤ => Gᵤ, :Hᵤ => Hᵤ, :Bₓ = Bₓ, :Cₓ = [], :Gₓ = [], :Hₓ = [], x_dims = 0)

nodal_dynamics = [(1/3, get_ThirdOrderMachineApprox, parameters_third_order), (1/3, get_DroopControlledInverterApprox, parameters_droop_controlled), (1/3, get_PQ, nothing)]
num_nodes = 100

pg_mixed = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics)

pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_mixed)