# Getting Started
Lets start by generating a simple network with only one type of dynamics actor on the nodes. In this case the dispatchable virtual oscillator `dVOC`. We begin by saving its control parameters η, α, κ and its frequency set point Ω
in a Dict. Note that `PowerDynamics.jl` uses a co-rotating reference frame and hence Ω = 0.
```@julia
parameters = Dict(:η => 3 * 10^(-3), :α => 5.0,  :κ => π/2, :Ω => 0)
```

The dynamics of the grid are defined in a vector. Each entry contains a tuple that holds first the nodal share then the function that generates the dynamics and finally the parameters of the nodes in the form of a Dict. As there is only one dynamic in this example our vector has only one entry:
```@julia
dVOC_share = 1.0
dynamics = [(dVOC_share, get_dVOCapprox, nodal_parameters)]
```

Using the nodal dynamics we can generate the struct that holds all relevant data for the synthetic grids. The struct uses predefined values for e.g. the Per Unit system unless the user explicitly defines them.  
```@julia
pg_struct = PGGeneration(num_nodes = 10, nodal_dynamics = dynamics,  lines = :StaticLine)
```

Then we can simply call the `generate_powergrid_dynamics` function to generate a synthetic grid!
```@julia
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)
```
# TODO: Sources, Normalform