"""
    probabilistic_capacity_expansion(pg_struct::PGGeneration, dist_load, dist_args)

Performs a probabilistic capacity expansion using the following steps:

1. Sample new Power Set Points from the distribution `dist_load`
2. Calculate the load flow using PowerDynamics.jl
3. Check the line limit using `validate_power_flow_on_lines` -> increase the capacity by adding a new cable to any overloaded line
4. Repeat steps 1.-3. N times
"""
function probabilistic_capacity_expansion(pg_struct::PGGeneration, dist_load, dist_args)
    num_nodes = pg_struct.num_nodes
    lines = get_lines(pg_struct)

    # Run over all scenarios
    for k in 1:pg_struct.num_tries
        P_set_prob = dist_load(dist_args...) # Sample new set points from the Distributions

        op_prob = get_ancillary_operationpoint(P_set_prob, pg_struct.V_vec, num_nodes, pg_struct.slack_idx, lines)
        save_network, save_flow = validate_power_flow_on_lines(op_prob, pg_struct)

        if save_network == false 
            overload_idx = findall(save_flow .== false) # Find overloaded line(s)
            cables_vec = pg_struct.cables_vec
            cables_vec[overload_idx] .= 3 .+ cables_vec[overload_idx] # Add a new cables to the overloaded line(s)

            pg_struct.cables_vec = cables_vec

            L = get_effective_distances(pg_struct.embedded_graph; mean_len_km = pg_struct.mean_len_km, shortest_line_km = pg_struct.shortest_line_km) # Effective spacial distances
            Y, Y_shunt = get_line_admittance_matrix(pg_struct.embedded_graph; L_matrix = L, cables_vec = pg_struct.cables_vec, num_nodes = pg_struct.num_nodes) # Update Admittances and shunts
            pg_struct.edge_parameters = Dict(:Y => Y, :Y_shunt => Y_shunt) 

            lines = get_lines(pg_struct) # Update Lines
        end
    end
    return pg_struct, lines
end

"""
    consumer_producer_nodal_power(P_set, num_nodes)

Generates a new load flow from the power stet points.
"""
function consumer_producer_nodal_power(P_set, num_nodes)
    μ = 1
    σ = 1/6
    
    net_producer = findall(P_set .> 0.0) 
    net_consumer = findall(P_set .< 0.0)

    on_peak_off_peak_factor = rand(Normal(μ, σ))
    P_set_new = Vector{Float64}(undef, num_nodes)

    for n in 1:num_nodes
        if n ∈ net_consumer
            P_set_new[n] = on_peak_off_peak_factor * P_set[n]
        elseif n ∈ net_producer
            P_set_new[n] = rand(Normal(μ, σ)) * P_set[n]
        end
    end

    P_set_new .-= sum(P_set_new) / (num_nodes)  # Assure power balance  
    
    return P_set_new
end