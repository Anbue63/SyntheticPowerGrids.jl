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
    num_nodes::Int64
    nodal_parameters::Dict
end

function test_struct(pg_struct::PGGeneration)
    if pg_struct.V_base != 380 * 10^3
        error("This voltage level is not supported. Please use V_base = 380 * 10^3 instead.")
    end
    
    if pg_struct.fluctuations != :off
        error("This fluctuation process is not supported. Please use fluctuations = :off instead.")
    end

    if pg_struct.loads != :PQAlgebraic && pg_struct.loads != :ExponentialRecovery && pg_struct.loads != :PQDynamic
        error("This option for the loads is not supported. Please use loads = :PQAlgebraic or :PQDynamic or :ExponentialRecovery instead.")
    end

    if pg_struct.power_distribution != :Bimodal
        error("This option for the power distribution is not supported. Please use power_distribution = :Bimodal instead.")
    end

    if pg_struct.generation_dynamics != :SchifferApprox && pg_struct.generation_dynamics != :Schmietendorf
        error("This option for the nodal dynamics is not supported. Please use generation_dynamics = :SchifferApprox or :Schmietendorf instead.")
    end

    if pg_struct.lines != :StaticLine && pg_struct.lines != :PiModelLine
        error("This option for the line dynamics is not supported. Please use lines = :StaticLine or :PiModelLine instead.")
    end

    if pg_struct.tests == false
        warning("The tests have been turned off. This option is not advised unless it is for debugging purposes.")
    end
end