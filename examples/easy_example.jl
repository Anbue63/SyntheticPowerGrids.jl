using SyntheticPowerGrids

##
# Lets start by generating a simple network with only one type of dynamics actor on the nodes. 
# In this case the dispatchable virtual oscillator `dVOC`. 
# We begin by saving its control parameters η, α, κ and its frequency set point Ω
nodal_parameters = Dict(:η => 3 * 10^(-3), :α => 5.0,  :κ => π/2, :Ω => 0)

# The dynamics of the grid are defined in a vector. 
# Each entry contains a tuple that holds first the nodal share then the function that generates the dynamics and finally the parameters of the nodes in the form of a Dict.
# As there is only one dynamic in this example our vector has only one entry:
nodal_dynamics = [(1.0, get_dVOCapprox, nodal_parameters)]


# Using the nodal dynamics we can generate the struct that holds all relevant data for the synthetic grids.
# The struct uses predefined values for e.g. the Per Unit system unless the user explicitly defines them.  
pg_struct = PGGeneration(num_nodes = 10, nodal_dynamics = nodal_dynamics,  lines = :StaticLine)

# Then we can simply call the `generate_powergrid_dynamics` function to generate a synthetic grid!
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

##
# We can use PowerDynamics.jl and the OrdinaryDiffEq.jl package to perform simulations of these grids.
using OrdinaryDiffEq
using PowerDynamics
using Plots

rpg = rhs(pg)
x0 = copy(op)

x0[1, :v] = 1.5

prob = ODEProblem(rpg, x0.vec, (0.0, 100.0), nothing)
sol = solve(prob, Rodas4())
solution1 = PowerGridSolution(sol, pg)

plot(solution1, :, :v, legend=false)