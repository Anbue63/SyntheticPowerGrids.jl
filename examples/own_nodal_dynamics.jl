# Adding new nodal dynamics
using Pkg
Pkg.activate(@__DIR__)
using SyntheticPowerGrids
##
#It is easy to add new node types for the synthetic grids, with out having to extend the package.
# We begin by loading the `PowerDynamics.jl` package to make use of its library of node dynamics.
using PowerDynamics

# In general new nodal dynamics have to be defined in the form of:
# `get_node_type(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters)`
# where P_set, Q_set and V_set are the set-points.

# As example we define a function which generates a `SwingEqLVS` node from the set-points and the nodal parameters:
function get_swingLVS(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    D = nodal_parameters[:D] # Damping Coefficient
    Γ = nodal_parameters[:Γ] # Voltage stability Coefficient

    SwingEqLVS(H = H, P = P_set, D = D, Ω = Ω, Γ = Γ, V = V_set)
end

# Then we can use the function `get_swingLVS` like any other function for the nodal dynamics.
# Again we define the parameters and the `PGGeneration` struct and generate the synthetic power grid.
parameters = Dict(:D => 10.0, :H => 5, :Γ => 10, :V => 1.0, :Ω => 2π * 50)  
dynamics = [(1.0, get_swingLVS, parameters_swing)]

pg_struct = PGGeneration(num_nodes = 100, nodal_dynamics = dynamics, lines = :StaticLine)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)