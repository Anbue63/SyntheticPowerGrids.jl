using Revise
using SyntheticPowerGrids

#nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0], :Nps => 0.6603, :Npt => 3.2308, :Nqs => -2.2439, :Nqt => 18.3881, :Tp => 0.0135, :Tq => 0.1017, :V0 => 1.0)

nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0], :Nps => 0.6789, :Npt => 2.6579, :Nqs => 1.1585, :Nqt => 23.5035, :Tp => 0.0140, :Tq => 0.0348, :V0 => 1.0)

a = PGGeneration3(num_nodes = 100, nodal_parameters = nodal_parameters, loads = :ExponentialRecovery)

pg, op, rejections = random_PD_grid(a)
