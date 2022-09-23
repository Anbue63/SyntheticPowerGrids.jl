"""
    get_nodes(pg, op_ancillary, pg_struct)

get the nodal dynamics of the full power grid.
- `pg`: Graph structure of the power grid
- `op_ancillary`: Operation point of the helper grid
"""
function get_nodes(pg, op_ancillary, pg_struct)
    nodal_dyn_prob = rand(nv(pg.graph)) # Turns nodes into loads or PhaseAmplitudeOscillators
    nodes = Array{Any}(undef, nv(pg.graph))
    
    load_share = pg_struct.nodal_shares[:load_share]

    if pg_struct.generation_dynamics == :SchifferApprox 
        schiffer_share = pg_struct.nodal_shares[:schiffer_share]

        nodes = get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share, load_share + schiffer_share])
    elseif pg_struct.generation_dynamics == :Schmietendorf
        schmietendorf_share = pg_struct.nodal_shares[:schmietendorf_share]

        nodes = get_nodes_schmietendorf(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share, load_share + schmietendorf_share])
    elseif pg_struct.generation_dynamics == :Mixed
        schmietendorf_share = pg_struct.nodal_shares[:schmietendorf_share]
        schiffer_share = pg_struct.nodal_shares[:schiffer_share]

        nodes = get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share, load_share + schiffer_share])
        nodes = get_nodes_schmietendorf(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share + schiffer_share, load_share + schiffer_share + schmietendorf_share])
    elseif pg_struct.generation_dynamics == :SwingEqLVS
        swingLVS_share =  pg_struct.nodal_shares[:swingLVS_share]

        nodes = get_nodes_swingLVS(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share, load_share + swingLVS_share])
    elseif pg_struct.generation_dynamics == :dVOCapprox 
        dVOCapprox_share = pg_struct.nodal_shares[:dVOC_share]

        nodes = get_nodes_dVOCapprox(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share, load_share + dVOCapprox_share])
    #elseif pg_struct.generation_dynamics == :SwingEq
    #    swing_share =  pg_struct.nodal_shares[:swing_share]
    
    #    nodes = get_nodes_swing(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [load_share, load_share + swing_share])
    end
    if pg_struct.loads == :PQAlgebraic
        nodes = get_nodes_PQ(pg, op_ancillary, nodes, nodal_dyn_prob, load_share)
    end

    if pg_struct.slack == true
        nodes[end] = SlackAlgebraic(U = complex(pg_struct.V_ref))
    end
    return nodes
end

function get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    τ_P = nodal_parameters[:τ_P] # Time constant low pass filter measuring the active power
    τ_Q = nodal_parameters[:τ_Q] # Time constant low pass filter measuring the reactive power
    V_r = nodal_parameters[:V_r] # Reference voltage magnitude. In [p.u.] systems the default is V_r = 1.0
    K_P = nodal_parameters[:K_P] # Gain constant low pass filter measuring the active power
    K_Q = nodal_parameters[:K_Q] # Gain constant low pass filter measuring the reactive power
    
    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] > threshold[1]  # Grid Forming
            if threshold[2] >= nodal_dyn_prob[n] 
                β = rand(1:length(nodal_parameters[:τ_P])) # Randomly chooses one of three possible time constant for the low pass filter measuring the active power
                τ_P_node = τ_P[β]
                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims = parameter_schiffer(P_set = op_ancillary[n, :p], Q_set = op_ancillary[n, :q], τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_r, τ_P = τ_P_node, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, xdims = xdims)
            end
        end
    end
    return nodes
end

function get_nodes_schmietendorf(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    X = nodal_parameters[:X] # Reactance
    α = nodal_parameters[:α] # Voltage dynamics time constant
    γ = nodal_parameters[:γ] # Damping Coefficient
    
    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] > threshold[1] # Third Order Machine
            if threshold[2] >= nodal_dyn_prob[n] 
                E_set = op_ancillary[n, :v]
                Q_n = op_ancillary[n, :q]
                E_f = E_set + (X * Q_n / (E_set)) # Electric field voltage that results in the correct nodal voltage magnitude

                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims = parameter_schmietendorf(P_m = op_ancillary[n, :p], E_f = E_f, E_set = E_set, X = X, α = α, γ = γ, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, xdims = xdims)
            end
        end
    end
    return nodes
end

function get_nodes_PQ(pg, op_ancillary, nodes, nodal_dyn_prob, threshold)
    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] < threshold # Grid following / loads
            nodes[n] = PQAlgebraic(P = op_ancillary[n, :p], Q = op_ancillary[n, :q]) 
        end
    end
    return nodes
end

function get_nodes_swingLVS(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    V = nodal_parameters[:V] # Reference Voltage, typically 1 [p.u.]
    D = nodal_parameters[:D] # Damping Coefficient
    Γ = nodal_parameters[:Γ] # Voltage stability Coefficient

    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] > threshold[1]
            if threshold[2] >= nodal_dyn_prob[n] 
                nodes[n] = SwingEqLVS(H = H, P = op_ancillary[n, :p], D = D, Ω = Ω, Γ = Γ, V = V)
            end
        end
    end
    return nodes
end

function get_nodes_dVOCapprox(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters

    Ω = nodal_parameters[:Ω] # Rated Frequency
    η = nodal_parameters[:η] # positive control parameter
    α = nodal_parameters[:α] # positive control parameter
    κ = nodal_parameters[:κ] # uniform complex phase

    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] > threshold[1]
            if threshold[2] >= nodal_dyn_prob[n]
                P_set = op_ancillary[n, :p] # active power set point
                Q_set = op_ancillary[n, :q] # reactive power set point
                V_set = op_ancillary[n, :v] # voltage magnitude set point
                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims = parameter_dVOC(P_set = P_set, Q_set = Q_set, V_set = V_set, Ω_set = Ω, η = η, α = α, κ = κ, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, xdims = xdims)
            end
        end
    end
    return nodes
end


#=
function get_nodes_swing(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, threshold)
    nodal_parameters = pg_struct.nodal_parameters
    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    D = nodal_parameters[:D] # Damping Coefficient
    
    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] > threshold[1]
            if threshold[2] >= nodal_dyn_prob[n] 
                nodes[n] = SwingEq(H = H, P = op_ancillary[n, :p], D = D, Ω = Ω)
            end
        end
    end
    return nodes
end
=#