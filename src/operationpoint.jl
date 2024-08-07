"""
    get_ancillary_operationpoint(P_vec, Q_vec, V_vec, node_types, slack_idx, lines)

Get an Ancillary Power Grid with the same power flow on the lines as the full power grid.
We use PV-Nodes to find the reactive power at the nodes with results in power grids where the voltage magnitude V is equal to the Reference voltage magnitude V_r.

- `P_vec`: Vector containing the nodal active powers
- `Q_vec`: Vector containing the nodal reactive powers
- `V_vec`: Vector containing the nodal voltage magnitudes
- `node_types`: Vector congaing the different node type for the ancillary grid
- `slack_idx`: Index of the slack bus in the ancillary grid
- `lines`: Line dynamics of the power grid
"""
function get_ancillary_operationpoint(P_vec, Q_vec, V_vec, node_types, slack_idx, lines, rejections)
    num_nodes = length(node_types)
    nodes = Array{Any}(undef, num_nodes)

    # Get Ancillary Power Grid
    for n in 1:num_nodes
        if node_types[n] == :PVAlgebraic
            nodes[n] = PVAlgebraic(P = P_vec[n], V = V_vec[n])
        elseif node_types[n] == :PQAlgebraic
            nodes[n] = PQAlgebraic(P = P_vec[n], Q = Q_vec[n])
        end
    end
    nodes[slack_idx] = SlackAlgebraic(U = complex(V_vec[slack_idx]))
    
    pg_cons = PowerGrid(nodes, lines)

    try 
        op = find_operationpoint(pg_cons, sol_method=:rootfind, solve_powerflow = true)
        return op, rejections
    catch
        println("Could not find the operation point. Problem unfeasible.")
        rejections.unfeasible += 1
        rejections.total += 1
    
        return nothing, rejections
    end
end

"""
    get_initial_guess(rpg, op_ancillary)

Initial guess for the operation point search of PowerDynamics.jl for the Operation point of the helper grid with the correct power flows and voltage magnitudes.
- `rpg`: Right-Hand side function of the power grid
- `op_ancillary`: Operation point of the helper grid
"""
function get_initial_guess(pg, op_ancillary)
    rpg = rhs(pg)

    ic_guess = zeros(length(rpg.syms))
    u_r_idx = findall(map(x -> occursin("u_r", string(x)), rpg.syms))
    u_i_idx = findall(map(x -> occursin("u_i", string(x)), rpg.syms))
    ic_guess[u_r_idx] = op_ancillary[:, :u_r]
    ic_guess[u_i_idx] = op_ancillary[:, :u_i]

    return ic_guess
end

