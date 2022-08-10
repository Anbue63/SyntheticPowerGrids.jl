"""
    get_nodes(pg, op_ancillary, pg_struct)

get the nodal dynamics of the full power grid.
- `pg`: Graph structure of the power grid
- `op_ancillary`: Operation point of the helper grid
"""
function get_nodes(pg, op_ancillary, pg_struct)
    α = rand(nv(pg.graph)) # Turns nodes into loads or PhaseAmplitudeOscillators
    nodes = Array{Any}(undef, nv(pg.graph))
    if pg_struct.generation_dynamics == :SchifferApprox 
        nodes = get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, α)
    end
    if pg_struct.loads == :PQ
        nodes = get_nodes_PQ(pg, op_ancillary, nodes, α)
    end
    nodes[end] = SlackAlgebraic(U = complex(pg_struct.V_ref))
    return nodes
end

function get_nodes_schiffer(pg, op_ancillary, nodes, pg_struct, α)
    nodal_parameters = pg_struct.nodal_parameters

    τ_P = nodal_parameters[:τ_P] # Time constant low pass filter measuring the active power
    τ_Q = nodal_parameters[:τ_Q] # Time constant low pass filter measuring the reactive power
    V_r = nodal_parameters[:V_r] # Reference voltage magnitude. In [p.u.] systems the default is V_r = 1.0
    K_P = nodal_parameters[:K_P] # Gain constant low pass filter measuring the active power
    K_Q = nodal_parameters[:K_Q] # Gain constant low pass filter measuring the reactive power

    for n in 1:nv(pg.graph)
        if α[n] > 0.5 # Grid Forming
            β = rand(1:length(nodal_parameters[:τ_P])) # Randomly chooses one of three possible time constant for the low pass filter measuring the active power
            τ_P_node = τ_P[β]
            nodes[n] = SchifferApprox(τ_P = τ_P_node, τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_r, P = op_ancillary[n, :p], Q = op_ancillary[n, :q], Y_n = 0)
        end
    end
    return nodes
end

function get_nodes_PQ(pg, op_ancillary, nodes, α)
    for n in 1:nv(pg.graph)
        if α[n] < 0.5 # Grid following / loads
            nodes[n] = PQAlgebraic(P = op_ancillary[n, :p], Q = op_ancillary[n, :q]) 
        end
    end
    return nodes
end