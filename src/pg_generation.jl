"""
    random_PD_grid(N::Int)

Generates a random power grid using SyntheticNetworks and then turns it into a PowerDynamics.PowerGrid type.
"""
function random_PD_grid(pg_struct::PGGeneration)
    validate_struct(pg_struct) # Test if options given by the user are valid

    # Accessing the data from the struct
    N = pg_struct.num_nodes                              # Number of nodes
    n0, p, q, r, s, u = pg_struct.SyntheticNetworksParas # Parameters for SyntheticNetworks
    Y_base = pg_struct.P_base / (pg_struct.V_base)^2     # Base admittance
    rejections = 0
    
    for i in 1:pg_struct.maxiters # maxiters until a stable grid is found 
        if pg_struct.embedded_graph === nothing
            pg_struct.embedded_graph = generate_graph(RandomPowerGrid(N, n0, p, q, r, s, u)) # Random power grid topology
        end
        
        if typeof(pg_struct.P_vec) == Vector{Nothing} 
            pg_struct.P_vec = get_power_distribution(pg_struct)
        end

        lines = get_lines(pg_struct, Y_base) # Line dynamics

        S_vec = get_ancillary_apparent_power(pg_struct, lines)
        pg_struct.P_vec = real.(S_vec)

        if typeof(pg_struct.Q_vec) == Vector{Nothing} 
            pg_struct.Q_vec = imag.(S_vec) # Reactive Power of the ancillary power grid
        end

        nodes = get_nodes(pg_struct) # Nodal dynamics
        pg = PowerGrid(nodes, lines)

        op = find_operationpoint(pg, sol_method = :rootfind) #, solve_powerflow = true) # find operation point of the full power grid

        if pg_struct.validators == true              # Sanity checks before returning
            if validate_power_grid(pg, op, pg_struct) == true
                return pg, op, pg_struct.embedded_graph, rejections
            end
        else
            return pg, op, embedded_graph, rejections
        end
        rejections += 1
    end
    return nothing, nothing, nothing, rejections
end 