using Revise
using SyntheticPowerGrids

nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0])

a = PGGeneration2(num_nodes = 100, nodal_parameters = nodal_parameters)

pg, op, rejections = random_PD_grid(a)



