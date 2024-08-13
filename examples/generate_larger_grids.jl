using Pkg
Pkg.activate(".")

using SyntheticPowerGrids
using PowerDynamics


const N_lower_bound = 10
const N_upper_bound = 10
const num_grids = 10000
const grid_index_start = 1
const grid_index_end = 1
const P0_offpeak = 1.3
const P0_onpeak = 3.18

nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 5.0)
nodal_parameters_b = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 1.0)
nodal_parameters_c = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :τ_P => 0.5)

nodal_dynamics = [
    (1 / 6, get_DroopControlledInverterApprox, nodal_parameters_a),
    (1 / 6, get_DroopControlledInverterApprox, nodal_parameters_b),
    (1 / 6, get_DroopControlledInverterApprox, nodal_parameters_c),
    (0.5, get_PQ, nothing),
]

# get a list containing the grid size using uniform sampling
random_N = rand(N_lower_bound:N_upper_bound, num_grids)
random_P0 = rand(P0_offpeak:P0_onpeak, num_grids)

path_storage = joinpath(@__DIR__, "./data")

c = SyntheticPowerGrids.PGGeneration(
        num_nodes=random_N[1],
        nodal_dynamics=nodal_dynamics,
        lines=:PiModelLine,
        slack=true,
        P0=random_P0[1],
        use_static_lines_for_helper_grid = true,
    )
pg, op, grid_pg_struct, rejections = generate_powergrid_dynamics(c)