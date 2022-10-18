using SyntheticPowerGrids

# Example: droop-controlled inverter with additional (exponentially declining) variable

V_r = 1.0
τ_P = 1.0
τ_Q = 8.0
k_P = 5.0
k_Q = 0.1
τ_X = 10.0

parameters_nf = Dict(:Bᵤ => [1im 0],
                     :Cᵤ => - 1/(2*τ_Q*V_r^2),
                     :Gᵤ => 0,
                     :Hᵤ => -k_Q/(τ_Q*V_r), 
                     :Bₓ => [-1/τ_P 0; 0 -1/τ_X],
                     :Cₓ => [0, 0],
                     :Gₓ => [-k_P/τ_P, 0],
                     :Hₓ => [0, 0],
                     :x_dims => 2)

nodal_dynamics = [(1/2, get_normalform, parameters_nf), (1/2, get_PQ, nothing)]
pg_struct = PGGeneration(num_nodes = 100, nodal_dynamics = nodal_dynamics)
pg, op, pg_struct_updated, rejections = generate_powergrid_dynamics(pg_struct)