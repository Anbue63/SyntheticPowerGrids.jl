"""
    get_nodes(pg, op_ancillary, pg_struct)

get the nodal dynamics of the full power grid.
- `pg`: Graph structure of the power grid
- `op_ancillary`: Operation point of the helper grid
"""
function get_nodes(pg, op_ancillary, pg_struct)
    nodal_dyn_prob = rand(nv(pg.graph)) # Turns nodes into loads or PhaseAmplitudeOscillators
    nodes = Array{Any}(undef, nv(pg.graph))
    if pg_struct.generation_dynamics == :SchifferApprox 
        nodes = get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [0.5, 1.0])
    elseif pg_struct.generation_dynamics == :Schmietendorf
        nodes = get_nodes_schmietendorf(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [0.5, 1.0])
    elseif pg_struct.generation_dynamics == :Mixed
        nodes = get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [0.5, 0.75])
        nodes = get_nodes_schmietendorf(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob, [0.75, 1.0])
    end
    if pg_struct.loads == :PQAlgebraic
        nodes = get_nodes_PQ(pg, op_ancillary, nodes, nodal_dyn_prob)
    elseif pg_struct.loads == :ExponentialRecovery
        nodes = get_nodes_ExponentialRecovery(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob)

    elseif pg_struct.loads == :PQDynamic
        nodes = get_nodes_PQDynamic(pg, op_ancillary, nodes, pg_struct, nodal_dyn_prob)
    end
    nodes[end] = SlackAlgebraic(U = complex(pg_struct.V_ref))
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
                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n = parameter_schiffer(P_set = op_ancillary[n, :p], Q_set = op_ancillary[n, :q], τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_r, τ_P = τ_P_node, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ)
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

                Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n = parameter_schmietendorf(P_m = op_ancillary[n, :p], E_f = E_f, E_set = E_set, X = X, α = α, γ = γ, Y_n = 0.0)
                nodes[n] = NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ)
            end
        end
    end
    return nodes
end

function get_nodes_PQ(pg, op_ancillary, nodes, nodal_dyn_prob)
    for n in 1:nv(pg.graph)
        if nodal_dyn_prob[n] < 0.5 # Grid following / loads
            nodes[n] = PQAlgebraic(P = op_ancillary[n, :p], Q = op_ancillary[n, :q]) 
        end
    end
    return nodes
end