@with_kw struct PGGeneration
    P_base::Float64 = 400 * 10^6
    V_base::Float64 = 380 * 10^3
    loads::Symbol = :PQAlgebraic
    lines::Symbol = :PiModelLine
    fluctuations::Symbol = :off
    generation_dynamics::Symbol = :SchifferApprox
    power_distribution::Symbol = :Bimodal
    maxiters::Int64 = 1000
    P0::Float64 = 1.0
    V_ref::Float64 = 1.0
    tests::Bool = true
    SyntheticNetworksParas::Vector{Float64} = [1, 1/5, 3/10, 1/3, 1/10, 0.0]
    shortest_line_km::Float64 = 0.06
    mean_len_km::Float64 = 42.872746445497626
    num_nodes::Int64
    nodal_parameters::Dict
    nodal_shares::Dict
end

function test_struct(pg_struct::PGGeneration)
    if sum(values(pg_struct.nodal_shares)) != 1.0    
        error("The sum of all nodal share has to equal 1.0!")
    end
    if pg_struct.V_base != 380 * 10^3
        error("This voltage level is not supported. Please use V_base = 380 * 10^3 instead.")
    end
    
    if pg_struct.fluctuations != :off
        error("This fluctuation process is not supported. Please use fluctuations = :off instead.")
    end

    if pg_struct.loads != :PQAlgebraic
        error("This option for the loads is not supported. Please use loads = :PQAlgebraic instead.")
    end

    if pg_struct.power_distribution != :Bimodal
        error("This option for the power distribution is not supported. Please use power_distribution = :Bimodal instead.")
    end

    if pg_struct.generation_dynamics != :SchifferApprox && pg_struct.generation_dynamics != :Schmietendorf && pg_struct.generation_dynamics != :Mixed
        error("This option for the nodal dynamics is not supported. Please use generation_dynamics = :SchifferApprox or :Schmietendorf or :Mixed instead.")
    end

    try 
        pg_struct.nodal_shares[:load_share]
    catch err
        error("Please define the share of loads in the network.")
    end

    if pg_struct.generation_dynamics == :Schmietendorf
        try 
            pg_struct.nodal_shares[:schmietendorf_share]
        catch err
            error("Please define the share of Schmietendorf nodes in the network.")
        end
    end

    if pg_struct.generation_dynamics == :SchifferApprox
        try 
            pg_struct.nodal_shares[:schiffer_share]
        catch err
            error("Please define the share of Schiffer nodes in the network.")
        end
    end

    if pg_struct.generation_dynamics == :Mixed
        try 
            pg_struct.nodal_shares[:schmietendorf_share]
        catch err
            error("Please define the share of Schmietendorf nodes in the network.")
        end

        try 
            pg_struct.nodal_shares[:schiffer_share]
        catch err
            error("Please define the share of Schiffer nodes in the network.")
        end
    end

    if pg_struct.lines != :StaticLine && pg_struct.lines != :PiModelLine
        error("This option for the line dynamics is not supported. Please use lines = :StaticLine or :PiModelLine instead.")
    end

    if pg_struct.tests == false
        warning("The tests have been turned off. This option is not advised unless it is for debugging purposes.")
    end
end