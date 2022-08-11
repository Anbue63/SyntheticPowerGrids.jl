using Revise
using SyntheticPowerGrids

#nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0], :Nps => 0.6603, :Npt => 3.2308, :Nqs => -2.2439, :Nqt => 18.3881, :Tp => 0.0135, :Tq => 0.1017, :V0 => 1.0)

nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0])

a = PGGeneration(num_nodes = 20, nodal_parameters = nodal_parameters, loads = :PQAlgebraic)

pg, op, rejections = random_PD_grid(a)