using Pkg
Pkg.activate(@__DIR__)
using SyntheticPowerGrids
using PowerDynamics
using OrdinaryDiffEq
using Plots
using PowerGridNoise
using Interpolations

default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=1.5, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)
##
nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 5.0) 
nodal_parameters_b = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 1.0) 
nodal_parameters_c = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)

nodal_dynamics = [(1/6, get_DroopControlledInverterApprox, nodal_parameters_a), (1/6, get_DroopControlledInverterApprox, nodal_parameters_b), (1/6, get_DroopControlledInverterApprox, nodal_parameters_c), (0.5, get_PQ, nothing)]

num_nodes = 100

##
a = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(a)

##
#
nodes = pg.nodes
fluc_node_idx = findfirst(typeof.(pg.nodes) .== PQAlgebraic)
P_set = nodes[fluc_node_idx].P
Q_set = nodes[fluc_node_idx].Q

##
# Easy Sine Example
fluc_func(t) = 0.1 * sin(t)

# add fluctuation to the load (exchange PQAlgebraic with FluctuationNode)
nodes[fluc_node_idx] = FluctuationNode(t -> P_set + fluc_func(t), t -> Q_set)
pg_fluc = PowerGrid(nodes, pg.lines)

##
#
tspan = (0.0, 50.0)
ode = ODEProblem(rhs(pg_fluc), op.vec, tspan)
sol = solve(ode, Rodas4())

solution1 = PowerGridSolution(sol, pg_fluc)
plot(solution1, [fluc_node_idx], :p, legend=false)

##
# Wind Model

D = 0.1 # Intermittence strength
p = 0.2 # Penetration parameter
x, t = wind_power_model(tspan, D = D)

x_intpol = linear_interpolation(t, x)

plot(t, x, idxs = 1, xlabel = "t[s]", ylabel = "x(t)", label = "Time series", lw = 3)
plot!(t, x_intpol(t), idxs = 1,label = "Interpolated time series", line_style = :dash, xlabel = "t[s]", ylabel = "x(t)")

##
# Add fluctuation to the load (exchange PQAlgebraic with FluctuationNode)
nodes[fluc_node_idx] = FluctuationNode(t -> P_set + p * x_intpol(t), t -> Q_set)
pg_fluc_wind = PowerGrid(nodes, pg.lines)

##
#
ode = ODEProblem(rhs(pg_fluc_wind), op.vec, tspan)
sol = solve(ode, Rodas4())

solution2 = PowerGridSolution(sol, pg_fluc_wind)
hline([P_set], label = "Set Point", alpha = 0.3, c = :black)
plot!(solution2, [fluc_node_idx], label = "Active Power",:p, lw = 3, ylabel = "P[p.u.]", xlabel = "t[s]")

##