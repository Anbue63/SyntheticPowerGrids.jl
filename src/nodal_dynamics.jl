get_load(P,Q,V, pars) = PQAlgebraic(P = P, Q = Q)
get_normal_form(P,Q,V, pars) = SwingEqLVS(H = 1., P = P, D = 0.1, Ω = 50., Γ = 3., V = V)

node_types_default = [(0.5, get_load, nothing), (0.5, get_normal_form, nothing)]

function get_nodes2(power_flow_solution; node_types = node_types_default)
    num_nodes = length(power_flow_solution)
    node_idxs = collect(1:num_nodes)
    nt_idxs = []
    for nt in node_types
        n_nodes = round(Int, nt[1] * num_nodes)
        s = sample(node_idxs, n_nodes, replace=false)
        symdiff!(node_idxs, s)
        push!(nt_idxs, s)
    end

    if ! (node_idxs == []) # Just to make sure all nodes are actually distributed to the different types
        append!(nt_idxs[end], node_idxs)
    end

    nodes = Array{Any}(undef, num_nodes)
    @assert isapprox(sum(x -> x[1], node_types), 1)
    for (nt, idxs) in zip(node_types, nt_idxs)
        for i in idxs
            nodes[i] = nt[2](power_flow_solution[i]..., nt[3]) # parameters
        end
    end
    nodes
end

"""
    get_nodes(pg_struct)

Get the nodal dynamics of the full power grid.
- `pg_struct`: Struct containing all data of the power grid
"""
function get_nodes(pg_struct::PGGeneration)
    nodal_dyn_prob = rand(pg_struct.num_nodes) # Turns nodes into loads or PhaseAmplitudeOscillators
    nodes = Array{Any}(undef, pg_struct.num_nodes)
    load_share = pg_struct.nodal_shares[:load_share]

    if pg_struct.generation_dynamics == :DroopControlledInverterApprox 
        DroopControlledInverterApprox_share = pg_struct.nodal_shares[:DroopControlledInverterApprox_share]

        nodes = get_nodes_DroopControlledInverterApprox(pg_struct, nodes, nodal_dyn_prob, [load_share, load_share + DroopControlledInverterApprox_share])
    elseif pg_struct.generation_dynamics == :ThirdOrderMachineApprox
        ThirdOrderMachineApprox_share = pg_struct.nodal_shares[:ThirdOrderMachineApprox_share]

        nodes = get_nodes_ThirdOrderMachineApprox(pg_struct, nodes, nodal_dyn_prob, [load_share, load_share + ThirdOrderMachineApprox_share])
    elseif pg_struct.generation_dynamics == :Mixed
        ThirdOrderMachineApprox_share = pg_struct.nodal_shares[:ThirdOrderMachineApprox_share]
        DroopControlledInverterApprox_share = pg_struct.nodal_shares[:DroopControlledInverterApprox_share]

        nodes = get_nodes_DroopControlledInverterApprox(pg_struct, nodes, nodal_dyn_prob, [load_share, load_share + DroopControlledInverterApprox_share])
        nodes = get_nodes_ThirdOrderMachineApprox(pg_struct, nodes, nodal_dyn_prob, [load_share + DroopControlledInverterApprox_share, load_share + DroopControlledInverterApprox_share + ThirdOrderMachineApprox_share])
    elseif pg_struct.generation_dynamics == :SwingEqLVS
        swingLVS_share =  pg_struct.nodal_shares[:swingLVS_share]

        nodes = get_nodes_swingLVS(pg_struct, nodes, nodal_dyn_prob, [load_share, load_share + swingLVS_share])
    elseif pg_struct.generation_dynamics == :dVOCapprox 
        dVOCapprox_share = pg_struct.nodal_shares[:dVOC_share]

        nodes = get_nodes_dVOCapprox(pg_struct, nodes, nodal_dyn_prob, [load_share, load_share + dVOCapprox_share])
    end
    if pg_struct.loads == :PQAlgebraic
        nodes = get_nodes_PQ(pg_struct, nodes, nodal_dyn_prob, load_share)
    end

    if pg_struct.slack == true
        nodes[pg_struct.slack_idx] = SlackAlgebraic(U = complex(pg_struct.V_vec[pg_struct.slack_idx]))
    end
    return nodes
end

function get_nodes_DroopControlledInverterApprox(pg_struct::PGGeneration, nodes, nodal_dyn_prob, threshold)
    embedded_graph = pg_struct.embedded_graph
    nodal_parameters = pg_struct.nodal_parameters

    P_vec = pg_struct.P_vec
    Q_vec = pg_struct.Q_vec 
    V_vec = pg_struct.V_vec 

    τ_P = nodal_parameters[:τ_P] # Time constant low pass filter measuring the active power
    τ_Q = nodal_parameters[:τ_Q] # Time constant low pass filter measuring the reactive power
    K_P = nodal_parameters[:K_P] # Gain constant low pass filter measuring the active power
    K_Q = nodal_parameters[:K_Q] # Gain constant low pass filter measuring the reactive power
    
    for n in 1:nv(embedded_graph.graph)
        if nodal_dyn_prob[n] > threshold[1]  # Grid Forming
            if threshold[2] >= nodal_dyn_prob[n] 
                β = rand(1:length(nodal_parameters[:τ_P])) # Randomly chooses one of three possible time constant for the low pass filter measuring the active power
                τ_P_node = τ_P[β]
                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, x_dims = parameter_DroopControlledInverterApprox(P_set = P_vec[n], Q_set = Q_vec[n], τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_vec[n], τ_P = τ_P_node, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
            end
        end
    end
    return nodes
end

function get_nodes_ThirdOrderMachineApprox(pg_struct::PGGeneration, nodes, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    P_vec = pg_struct.P_vec
    Q_vec = pg_struct.Q_vec 
    V_vec = pg_struct.V_vec

    X = nodal_parameters[:X] # Reactance
    α = nodal_parameters[:α] # Voltage dynamics time constant
    γ = nodal_parameters[:γ] # Damping Coefficient

    for n in 1:pg_struct.num_nodes
        if nodal_dyn_prob[n] > threshold[1] # Third Order Machine
            if threshold[2] >= nodal_dyn_prob[n] 
                E_f = V_vec[n] + (X * Q_vec[n] / (V_vec[n])) # Electric field voltage that results in the correct nodal voltage magnitude

                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, x_dims = parameter_ThirdOrderMachineApprox(P_m = P_vec[n], E_f = E_f, E_set = V_vec[n], X = X, α = α, γ = γ, Y_n = 0.0)

                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
            end
        end
    end
    return nodes
end

function get_nodes_PQ(pg_struct::PGGeneration, nodes, nodal_dyn_prob, threshold)
    P_vec = pg_struct.P_vec
    Q_vec = pg_struct.Q_vec 

    for n in 1:pg_struct.num_nodes
        if nodal_dyn_prob[n] < threshold # Grid following / loads
            nodes[n] = PQAlgebraic(P = P_vec[n], Q = Q_vec[n]) 
        end
    end
    return nodes
end

function get_nodes_swingLVS(pg_struct::PGGeneration, nodes, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    P_vec = pg_struct.P_vec
    V_vec = pg_struct.V_vec

    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    V = nodal_parameters[:V] # Reference Voltage, typically 1 [p.u.]
    D = nodal_parameters[:D] # Damping Coefficient
    Γ = nodal_parameters[:Γ] # Voltage stability Coefficient

    for n in 1:pg_struct.num_nodes
        if nodal_dyn_prob[n] > threshold[1]
            if threshold[2] >= nodal_dyn_prob[n] 
                nodes[n] = SwingEqLVS(H = H, P = P_vec[n], D = D, Ω = Ω, Γ = Γ, V = V_vec[n])
            end
        end
    end
    return nodes
end

function get_nodes_dVOCapprox(pg_struct::PGGeneration, nodes, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    P_vec = pg_struct.P_vec
    Q_vec = pg_struct.Q_vec 
    V_vec = pg_struct.V_vec

    Ω = nodal_parameters[:Ω] # Rated Frequency
    η = nodal_parameters[:η] # positive control parameter
    α = nodal_parameters[:α] # positive control parameter
    κ = nodal_parameters[:κ] # uniform complex phase

    for n in 1:pg_struct.num_nodes
        if nodal_dyn_prob[n] > threshold[1]
            if threshold[2] >= nodal_dyn_prob[n]
                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, x_dims = parameter_dVOC(P_set = P_vec[n], Q_set = Q_vec[n], V_set = V_vec[n], Ω_set = Ω, η = η, α = α, κ = κ, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
            end
        end
    end
    return nodes
end
