function probabilistic_capacity_expansion(pg_struct::PGGeneration, dist_load)
    num_nodes = pg_struct.num_nodes
    lines = get_lines(pg_struct)

    # Run over all scenarios
    for k in 1:pg_struct.num_tries
        P_set_prob = map(node -> rand(dist_load[node]), 1:num_nodes) # Sample new set points from the Distributions
        P_set_prob .-= sum(P_set_prob) / (num_nodes)                       # Assure power balance

        op_prob = get_ancillary_operationpoint(P_set_prob, pg_struct.V_vec, num_nodes, pg_struct.slack_idx, lines)
        save_network, save_flow = validate_power_flow_on_lines(op_prob, pg_struct)

        if save_network == false 
            overload_idx = findall(save_flow .== false) # Find overloaded line(s)
            cables_vec = pg_struct.cables_vec
            cables_vec[overload_idx] += 1               # Add a new cables to the overloaded line(s)

            pg_struct.cables_vec = cables_vec

            Y, Y_shunt = get_line_admittance_matrix(pg_struct, L) # Update Admittances and shunts
            pg_struct.edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 

            lines = get_lines(pg_struct) # Update Lines
        end
    end
    return pg_struct, lines
end

function nodal_power_normal_distribution(P_set, num_nodes)
    μ = P_set
    σ = abs.(P_set)

    dist_load = Vector{Any}(undef, num_nodes)

    for n in 1:num_nodes
        dist_load[n] = Normal(μ[n], σ[n]) 
    end

    return dist_load
end