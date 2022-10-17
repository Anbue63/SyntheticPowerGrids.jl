using SyntheticPowerGrids

##
# Generating networks with different types of dynamics on the nodes is also possible. In this example we will generate a network consisting of Droop controlled inverters, Third order machines and PQ Busses. 
# We begin again by defining the parameters for the nodes, this time each type of node gets its own Dict: 
parameters_third_order = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
parameters_droop_controlled = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5.0, :τ_P => 5.0) 

# Then we write the vector containing the dynamics, each node should show up with a share of 1/3:
nodal_dynamics = [(1/3, get_ThirdOrderMachineApprox, parameters_third_order), (1/3, get_DroopControlledInverterApprox, parameters_droop_controlled), (1/3, get_PQ, nothing)]

# Again we generate the struct and simply call `generate_powergrid_dynamics`
num_nodes = 100
pg_mixed = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, slack = true)

##
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_mixed)