# Nodal Fluctuations
In this example we will show how to add fluctuations to the synthetic grids.

We begin by loading all relevant packages:
```@julia
using SyntheticPowerGrids
using PowerDynamics
using OrdinaryDiffEq
using Plots
using Interpolations

default(grid = false, foreground_color_legend = nothing, bar_edges = false,  lw=1.5, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 15, tickfontsize = 10)
```

Then we generate a synthetic grid as usual:

```@julia
nodal_parameters = Dict(:τ_Q => 5.0, :K_P => 5, :K_Q => 5.0, :τ_P => 0.5)
nodal_dynamics = [(1/2, get_DroopControlledInverterApprox, nodal_parameters), (0.5, get_PQ, nothing)]
num_nodes = 100
pg_struct = PGGeneration(num_nodes = num_nodes, nodal_dynamics = nodal_dynamics, lines = :StaticLine)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)
```

We load the power grid nodes and their set points and choose a node we were we want to simulate a fluctuation.
```@julia
nodes = copy(pg.nodes)
fluc_node_idx = findfirst(typeof.(pg.nodes) .== PQAlgebraic)
P_set = nodes[fluc_node_idx].P
Q_set = nodes[fluc_node_idx].Q
```

We simply use a sine as the fluctuation in this example:
```@julia
fluc_func(t) = 0.1 * sin(t)
```

We add the fluctuation to the load by exchanging a load bus (PQAlgebraic) with a FluctuationNode. 
Then we generate a new power grid with the updated nodes:
```@julia
nodes[fluc_node_idx] = FluctuationNode(t -> P_set + fluc_func(t), t -> Q_set)
pg_fluc = PowerGrid(nodes, pg.lines)
```

We can simulate this system and see that the power is fluctuating around its set-point.
```@julia
tspan = (0.0, 50.0)
ode = ODEProblem(rhs(pg_fluc), op.vec, tspan)
sol = solve(ode, Rodas4())

solution1 = PowerGridSolution(sol, pg_fluc)
plot(solution1, [fluc_node_idx], :p, legend=false)
```