"""
    generate_powergrid_dynamics(pg_struct::PGGeneration)

Generates a random power grid.
"""
function generate_powergrid_dynamics(pg_struct::PGGeneration)
    validate_struct(pg_struct) # Test if options given by the user are valid
    
    # Accessing the data from the struct
    N = pg_struct.num_nodes                              # Number of nodes
    rejections = Rejections()
    
    for i in 1:pg_struct.maxiters # maxiters until a stable grid is found 
        pg_struct_updated = deepcopy(pg_struct) # Generate container for new info

        if pg_struct.embedded_graph === nothing
            n0, p, q, r, s, u = pg_struct.SyntheticNetworksParas # Parameters for SyntheticNetworks
            pg_struct_updated.embedded_graph = generate_graph(RandomPowerGrid(N, n0, p, q, r, s, u)) # Random power grid topology
        end
        
        if typeof(pg_struct.P_vec) == Vector{Nothing} 
            pg_struct_updated.P_vec = get_power_distribution(pg_struct)
        end

        if pg_struct.slack_idx === nothing
            pg_struct_updated.slack_idx = findmax(pg_struct_updated.P_vec)[2] 
        end

        if pg_struct.cables_vec === nothing           # If there is no predefined number of cables use typical number for all lines
            e = edges(pg_struct_updated.embedded_graph.graph)
            pg_struct_updated.cables_vec = 3 * ones(length(e)) 
        end

        if pg_struct.probabilistic_capacity_expansion == true # Use a probabilistic load flow to expand the capacity such that it full fills all scenarios
            if pg_struct.dist_load === nothing
                pg_struct_updated.dist_load = get_bimodal_distribution # Default leads to 
                k = Vector{Any}(undef, 2)
                k[1] = pg_struct_updated.num_nodes
                k[2] = pg_struct_updated.P0
                pg_struct_updated.dist_args = k
            end
            pg_struct_updated, lines = probabilistic_capacity_expansion(pg_struct_updated, pg_struct_updated.dist_load, pg_struct_updated.dist_args)
        else
            lines = get_lines(pg_struct_updated) # Line dynamics
        end

        op_ancillary = get_ancillary_operationpoint(pg_struct_updated.P_vec, pg_struct_updated.Q_vec, pg_struct_updated.V_vec, pg_struct_updated.node_types_ancillary, pg_struct_updated.slack_idx, lines)

        pg_struct_updated.P_vec = op_ancillary[:, :p]
        pg_struct_updated.V_vec = op_ancillary[:, :v]

        if typeof(pg_struct.Q_vec) == Vector{Nothing} 
            pg_struct_updated.Q_vec = op_ancillary[:, :q] # Reactive Power of the ancillary power grid
        end

        nodes = get_nodes(pg_struct_updated.P_vec, pg_struct_updated.Q_vec, pg_struct_updated.V_vec, pg_struct_updated.slack, pg_struct_updated.slack_idx, pg_struct_updated.nodal_dynamics)
        pg = PowerGrid(nodes, lines)

        ic_guess = get_initial_guess(pg, op_ancillary)

        println("Searching for the operation point...")
        op = find_operationpoint(pg, ic_guess, sol_method=:rootfind)
            
        if pg_struct.validators == true # Sanity checks before returning
            if validate_power_grid(pg, op, pg_struct_updated, rejections) == true
                return pg, op, pg_struct_updated, rejections
            end
        else
            return pg, op, pg_struct_updated, rejections
        end
    end
    return nothing, nothing, nothing, rejections
end 