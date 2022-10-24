"""
    get_power_distribution(pg_struct)

Generates the active power set points for the nodes. The default option is the bimodal distribution.
"""
function get_power_distribution(pg_struct)
    if pg_struct.power_distribution == :Bimodal # Bimodal Distribution for the active power
        P_vec = get_bimodal_distribution(pg_struct)
    elseif pg_struct.power_distribution ==  :Plus_Minus_1
        P_vec = get_plus_minus_1_distribution(pg_struct)
    else 
        error("This option for the power distribution is not supported.")
    end
    return P_vec
end

"""
    get_bimodal_distribution(pg_struct)

Generates the active power set points for the nodes. Samples the set-points from a bimodal distribution as suggested in [1]

[1] H. Taher, S. Olmi and E. Schöll, 2019, Physical Review E, 100 062306
"""
function get_bimodal_distribution(pg_struct)
    num_nodes = pg_struct.num_nodes                                        # Number of nodes
    P0 = pg_struct.P0
    power_dist = MixtureModel(Normal[Normal(P0, P0/2), Normal(-P0, P0/2)]) # Distribution for the active power
    P_vec = rand(power_dist, num_nodes)                                    # Power Generation / Consumption of the nodes
    P_vec .-= sum(P_vec) / (num_nodes)                                     # Assure power balance

    return P_vec
end

"""
get_bimodal_distribution(pg_struct)

Generates the active power set points for the nodes.
Turns half of them into producers with P = +1 and the other half into consumers with P = -1. Often used in the theoretical physics community.
"""
function get_plus_minus_1_distribution(pg_struct)
    num_nodes = pg_struct.num_nodes # Number of nodes
    P_vec = ones(num_nodes)
    consumers = sample(1:num_nodes, Int(num_nodes/2); replace=false)
    P_vec[consumers] .= -1

    return P_vec
end