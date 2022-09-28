"""
    random_PD_grid(N::Int)

Generates a random power grid using SyntheticNetworks and then turns it into a PowerDynamics.PowerGrid type.
"""
function random_PD_grid(pg_struct::PGGeneration)
    validate_struct(pg_struct) # Test if options given by the user are valid
    
    # Accessing the data from the struct
    N = pg_struct.num_nodes                              # Number of nodes
    rejections = 0
    
    for i in 1:pg_struct.maxiters # maxiters until a stable grid is found 
        if pg_struct.embedded_graph === nothing
            n0, p, q, r, s, u = pg_struct.SyntheticNetworksParas # Parameters for SyntheticNetworks
            pg_struct.embedded_graph = generate_graph(RandomPowerGrid(N, n0, p, q, r, s, u)) # Random power grid topology
        end
        
        if typeof(pg_struct.P_vec) == Vector{Nothing} 
            pg_struct.P_vec = get_power_distribution(pg_struct)
        end

        if pg_struct.cables_vec === nothing
            e = edges(pg_struct.embedded_graph.graph)
            pg_struct.cables_vec = 3 * ones(length(e))
        end

        lines = get_lines(pg_struct) # Line dynamics

        op_ancillary = get_ancillary_operationpoint(pg_struct, lines)
        pg_struct.P_vec = op_ancillary[:, :p]
        pg_struct.V_vec = op_ancillary[:, :v]

        if typeof(pg_struct.Q_vec) == Vector{Nothing} 
            pg_struct.Q_vec = op_ancillary[:, :q] # Reactive Power of the ancillary power grid
        end

        nodes = get_nodes(pg_struct) # Nodal dynamics
        pg = PowerGrid(nodes, lines)
        ic_guess = get_initial_guess(pg, op_ancillary)

        op = find_operationpoint(pg, ic_guess, sol_method = :rootfind) #, solve_powerflow = true) # find operation point of the full power grid

        if pg_struct.validators == true # Sanity checks before returning
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