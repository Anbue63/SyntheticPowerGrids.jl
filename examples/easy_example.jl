#using Pkg
#Pkg.activate(@__DIR__)
using SyntheticPowerGrids

##
nodal_parameters = Dict(:η => 3 * 10^(-3), :α => 5.0,  :κ => π/2)
nodal_dynamics = [(1.0, get_dVOCapprox, nodal_parameters)]

pg_struct = PGGeneration(num_nodes = 10, nodal_dynamics = nodal_dynamics,  lines = :StaticLine)

##
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

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