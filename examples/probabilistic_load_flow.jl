#using Pkg
#Pkg.activate(@__DIR__)

using Revise
using SyntheticPowerGrids

##
#
nodal_parameters = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)
nodal_shares = Dict(:load_share => 0.5, :DroopControlledInverterApprox_share => 0.5)
pg_struct = PGGeneration(num_nodes = 10, nodal_parameters = nodal_parameters, loads = :PQAlgebraic, lines = :StaticLine, nodal_shares = nodal_shares, probabilistic_capacity_expansion = true, P0 = 50, validators = false, num_tries = 20)

pg, op, pg_struct_updated, rejections = random_PD_grid(pg_struct)

unique(pg_struct_updated.cables_vec)