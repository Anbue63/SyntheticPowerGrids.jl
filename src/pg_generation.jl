"""
    random_PD_grid(N::Int)

Generates a random power grid using SyntheticNetworks and then turns it into a PowerDynamics.PowerGrid type.
"""
function random_PD_grid(pg_struct::PGGeneration2; τ_Q = 8.0, K_P = 5, K_Q = 0.1)
    test_struct(pg_struct) # Test if all options given by the user are valid

    # Accessing the data from the struct
    N = pg_struct.num_nodes                              # Number of nodes
    n0, p, q, r, s, u = pg_struct.SyntheticNetworksParas # Parameters for SyntheticNetworks
    Y_base = pg_struct.P_base / (pg_struct.V_base)^2     # Base admittance
    power_dist = get_power_distribution(pg_struct)       # Distribution for the active power

    rejections = 0
    for i in 1:pg_struct.maxiters # maxiters until a stable grid is found 
        pg = generate_graph(RandomPowerGrid(N, n0, p, q, r, s, u)) # Random power grid topology
        
        P_vec = rand(power_dist, N - 1) # Power Generation / Consumption of the nodes
        P_vec .-= sum(P_vec) / (N - 1)  # Assure power balance

        Y = get_line_admittance_matrix(pg) # Line admittance matrix, Entry's are in Ohm
                
        lines = get_lines(pg, Y, Y_base)                     # Line dynamics
        op_ancillary = get_ancillary_grid(pg, P_vec, lines)  # Operation point of Ancillary power grid
        nodes = get_nodes(pg, op_ancillary, pg_struct)       # Nodal dynamics

        pg = PowerGrid(nodes, lines)
        rpg = rhs(pg)

        ic_guess = get_initial_guess(rpg, op_ancillary)                # Initial guess for rootfind
        op = find_operationpoint(pg, ic_guess, sol_method = :rootfind) #, solve_powerflow = true) # find operation point of the full power grid

        if test_power_grid(pg, op) == true
            return pg, op, rejections
        end
        rejections += 1
    end
end 




