function get_power_distribution(pg_struct)
    N = pg_struct.num_nodes                              # Number of nodes

    if pg_struct.power_distribution == :Bimodal # Bimodal Distribution for the active power
        P0 = pg_struct.P0
        power_dist = MixtureModel(Normal[Normal(P0, P0/2), Normal(-P0, P0/2)]) # Distribution for the active power
        P_vec = rand(power_dist, N - 1) # Power Generation / Consumption of the nodes
        P_vec .-= sum(P_vec) / (N - 1)  # Assure power balance
    elseif pg_struct.power_distribution ==  :Plus_Minus_1
        P_vec = ones(N)
        consumers = sample(1:N, Int(N/2); replace=false)
        P_vec[consumers] .= -1
    else 
        error("This option for the power distribution is not supported.")
    end

    return P_vec
end