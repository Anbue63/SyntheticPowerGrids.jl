using Revise
using SyntheticPowerGrids

#nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0], :Nps => 0.6603, :Npt => 3.2308, :Nqs => -2.2439, :Nqt => 18.3881, :Tp => 0.0135, :Tq => 0.1017, :V0 => 1.0)

nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) # This option recover the GNN dataset grids, please do not delete
nodal_parameters_b = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 

a = PGGeneration(num_nodes = 100, nodal_parameters = nodal_parameters_a, loads = :PQAlgebraic, lines = :PiModelLine)
b = PGGeneration(num_nodes = 100, nodal_parameters = nodal_parameters_b, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :Schmietendorf)

pg, op, rejections = random_PD_grid(b)

##
using OrdinaryDiffEq
using PowerDynamics
using Plots

rpg = rhs(pg)

ω_idx = findall(map(x -> 'ω' == string(x)[1], rpg.syms))

x0 = copy(op)

x0.vec[ω_idx[1]] = 1

prob = ODEProblem(rpg, x0.vec, (0.0, 100.0), nothing)
sol = solve(prob, Rodas4())

plot(sol, idxs = ω_idx)