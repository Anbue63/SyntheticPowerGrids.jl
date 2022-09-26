#using Pkg
#Pkg.activate(@__DIR__)

using Revise
using SyntheticPowerGrids

nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) # This option recover the GNN dataset grids, please do not delete it
nodal_parameters_b = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
nodal_parameters_c = Dict(:X => 1.0, :γ => 0.2, :α => 2.0, :τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) 
nodal_parameters_d = Dict(:D => 10.0, :H => 5, :Γ => 10, :V => 1.0, :Ω => 2π * 50)  # for Jakob and Mehrnaz networks

nodal_shares_a = Dict(:load_share => 0.5, :DroopControlledInverterApprox_share => 0.5)
nodal_shares_b = Dict(:load_share => 0.5, :ThirdOrderMachineApprox_share => 0.5)
nodal_shares_c = Dict(:load_share => 0.5, :ThirdOrderMachineApprox_share => 0.25, :DroopControlledInverterApprox_share => 0.25)
nodal_shares_d = Dict(:load_share => 0.0, :swingLVS_share => 1.0)
edge_parameters = Dict(:K => -10im)

a = PGGeneration1(num_nodes = 100, nodal_parameters = nodal_parameters_a, loads = :PQAlgebraic, lines = :StaticLine, nodal_shares = nodal_shares_a)
b = PGGeneration1(num_nodes = 100, nodal_parameters = nodal_parameters_b, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :ThirdOrderMachineApprox, nodal_shares = nodal_shares_b)
c = PGGeneration1(num_nodes = 100, nodal_parameters = nodal_parameters_c, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :Mixed, nodal_shares = nodal_shares_c)
d = PGGeneration1(num_nodes = 100, power_distribution = :Plus_Minus_1, nodal_parameters = nodal_parameters_d, loads = :PQAlgebraic, lines = :StaticLine, generation_dynamics = :SwingEqLVS, coupling = :homogenous, edge_parameters = edge_parameters, nodal_shares = nodal_shares_d, slack = false)

for_mehrnaz = PGGeneration1(num_nodes = 100, nodal_parameters = nodal_parameters_d, loads = :PQAlgebraic, lines = :StaticLine, generation_dynamics = :SwingEqLVS, nodal_shares = nodal_shares_d)

##
pg, op, embedded_graph, rejections = random_PD_grid(d)
op[:, :p] # vector containing all nodal active powers
# get_effective_distances(embedded_graph, mean_len_km = c.mean_len_km, shortest_line_km = c.shortest_line_km)

##
using OrdinaryDiffEq
using PowerDynamics
using Plots

rpg = rhs(pg)
x0 = copy(op)

x0.vec[1] = 2

prob = ODEProblem(rpg, x0.vec, (0.0, 100.0), nothing)
sol = solve(prob, Rodas4())

sol = PowerGridSolution(sol, pg)

plot(sol, :, :v, legend=false)