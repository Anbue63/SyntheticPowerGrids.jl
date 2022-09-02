@with_kw struct PGGeneration
    P_base::Float64 = 100 * 10^6; @assert P_base > 0.0 "Base Power has to be positive."
    V_base::Float64 = 380 * 10^3; @assert V_base > 0.0 "Base Voltage has to be positive."
    loads::Symbol = :PQAlgebraic;
    coupling::Symbol = :line_lengths
    lines::Symbol = :PiModelLine;
    generation_dynamics::Symbol = :SchifferApprox;
    power_distribution::Symbol = :Bimodal;
    maxiters::Int64 = 1000; @assert maxiters > 0.0 "Maxiters has to be positive."
    P0::Float64 = 1.31; @assert maxiters > 0.0 "Reference power for power distribution has to be positive."
    V_ref::Float64 = 1.0; @assert V_ref > 0.0 "Reference voltage magnitude has to be positive."
    validators::Bool = true;
    SyntheticNetworksParas::Vector{Float64} = [1, 1/5, 3/10, 1/3, 1/10, 0.0];
    shortest_line_km::Float64 = 0.06; @assert mean_len_km >= 0.0 "The shortest line length has to be not be negative."
    mean_len_km::Float64 = 37.12856121212121; @assert mean_len_km > 0.0 "The mean line length has to be bigger than 0.0."
    num_nodes::Int64; @assert num_nodes > 0.0 "Number of nodes can not be negative."
    nodal_parameters::Dict;
    edge_parameters::Dict = Dict();
    nodal_shares::Dict; @assert sum(values(nodal_shares)) == 1.0 "The sum of all nodal share has to equal 1.0!"
    slack::Bool = true
end

function validate_struct(pg_struct::PGGeneration)
    if pg_struct.V_base != 380 * 10^3
        error("This voltage level is not supported. Please use V_base = 380 * 10^3 instead.")
    end

    if pg_struct.loads != :PQAlgebraic
        error("This option for the loads is not supported. Please use loads = :PQAlgebraic instead.")
    end

    if pg_struct.power_distribution != :Bimodal
        error("This option for the power distribution is not supported. Please use power_distribution = :Bimodal instead.")
    end

    if pg_struct.coupling != :line_lengths && pg_struct.coupling != :homogenous
        error("This option for the coupling is not supported. Please use coupling = :line_lengths or :homogenous instead.")
    end

    if pg_struct.generation_dynamics != :SchifferApprox && pg_struct.generation_dynamics != :Schmietendorf && pg_struct.generation_dynamics != :Mixed  && pg_struct.generation_dynamics != :SwingEqLVS #&& pg_struct.generation_dynamics != :SwingEq
        error("This option for the nodal dynamics is not supported. Please use generation_dynamics = :SchifferApprox, :Schmietendorf, :SwingEqLVS, :SwingEq or :Mixed instead.")
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

    if pg_struct.validators == false
       @warn "The validators have been turned off. This option is not advised unless it is for debugging purposes."
    end
end