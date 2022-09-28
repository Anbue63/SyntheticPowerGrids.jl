using Pkg
#Pkg.activate(@__DIR__)
using Revise
using SyntheticPowerGrids

nodal_parameters = Dict(:η => 3 * 10^(-3), :α => 5.0,  :κ => π/2, :Ω => 0)
nodal_shares = Dict(:load_share => 0.0, :dVOC_share => 1.0)
pg_struct = PGGeneration(num_nodes = 100, nodal_parameters = nodal_parameters, nodal_shares = nodal_shares,  lines = :StaticLine, generation_dynamics = :dVOCapprox)

##
pg, op, embedded_graph, rejections = random_PD_grid(pg_struct)

##
using OrdinaryDiffEq
using PowerDynamics
using Plots

rpg = rhs(pg)
x0 = copy(op)

x0[1, :v] = 1.5

prob = ODEProblem(rpg, x0.vec, (0.0, 100.0), nothing)
sol = solve(prob, Rodas4())
solution1 = PowerGridSolution(sol, pg)

##
plot(solution1, :, :v, legend=false)