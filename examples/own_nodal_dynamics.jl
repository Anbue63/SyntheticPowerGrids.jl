#using Pkg
#Pkg.activate(@__DIR__)

using SyntheticPowerGrids
using PowerDynamics

##
# Define your own function for a node type
# Here we use simply the swing equation!

function get_swingLVS(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    D = nodal_parameters[:D] # Damping Coefficient
    Γ = nodal_parameters[:Γ] # Voltage stability Coefficient

    SwingEqLVS(H = H, P = P_set, D = D, Ω = Ω, Γ = Γ, V = V_set)
end

parameters_swing = Dict(:D => 10.0, :H => 5, :Γ => 10, :V => 1.0, :Ω => 2π * 50)  # for Jakob and Mehrnaz networks
nodal_dynamics = [(1.0, get_swingLVS, parameters_swing)]

##
# mehrnaz
pg_struct = PGGeneration(num_nodes = 100, nodal_dynamics = nodal_dynamics, lines = :StaticLine)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)

##
# Jakob
edge_parameters = Dict(:K => -10im)
d = PGGeneration(num_nodes = 100, power_distribution = :Plus_Minus_1, nodal_dynamics = nodal_dynamics, lines = :StaticLine, coupling = :homogenous, edge_parameters = edge_parameters)

pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(d)