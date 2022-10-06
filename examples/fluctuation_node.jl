using Pkg
Pkg.activate(@__DIR__)
using SyntheticPowerGrids

##
nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 5.0) 
nodal_parameters_b = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 1.0) 
nodal_parameters_c = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)

nodal_dynamics = [(1/6, get_DroopControlledInverterApprox, nodal_parameters_a), (1/6, get_DroopControlledInverterApprox, nodal_parameters_b), (1/6, get_DroopControlledInverterApprox, nodal_parameters_c), (0.5, get_PQ, nothing)]

num_nodes = 100

##
a = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine, validators = false)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(a)

##
fluc_func(t) = 0.1 * sin(t)

##
using PowerDynamics

nodes = pg.nodes
fluc_node_idx = findfirst(typeof.(pg.nodes) .== PQAlgebraic)
P_set = nodes[fluc_node_idx].P
Q_set = nodes[fluc_node_idx].Q

# add fluctuation to the load (exchange PQAlgebraic with FluctuationNode)
nodes[fluc_node_idx] = FluctuationNode(t -> P_set + fluc_func(t), t -> Q_set)
pg_fluc = PowerGrid(nodes, pg.lines)

##
using OrdinaryDiffEq
timespan = (-5.,50.)
ode = ODEProblem(rhs(pg_fluc), op.vec, timespan)
dqsol = solve(ode, Rodas4())

solution1 = PowerGridSolution(dqsol, pg_fluc)

##
using Plots
plot(solution1, :, :p, legend=false)
