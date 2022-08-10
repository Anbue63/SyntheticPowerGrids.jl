function get_power_distribution(pg_struct)
    if pg_struct.power_distribution == :Bimodal # Bimodal Distribution for the active power
        P0 = pg_struct.P0
        return MixtureModel(Normal[Normal(P0, P0/2), Normal(-P0, P0/2)])
    else 
        error("This option for the power distribution is not supported.")
    end
end