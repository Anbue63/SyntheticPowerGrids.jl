"""
    get_ancillary_operationpoint(P_vec, V_vec, num_nodes, slack_idx, lines)

Get an Ancillary Power Grid with the same power flow on the lines as the full power grid.
We use PV-Nodes to find the reactive power at the nodes with results in power grids where the voltage magnitude V is equal to the Reference voltage magnitude V_r.

- `pg_struct`: Struct containing all data of the power grid
- `lines`: Line dynamics of the power grid
"""
function get_ancillary_operationpoint(P_vec, V_vec, num_nodes, slack_idx, lines)
    nodes = Array{Any}(undef, num_nodes)

    # Get Ancillary Power Grid
    for n in 1:num_nodes
        nodes[n] = PVAlgebraic(P = P_vec[n], V = V_vec[n])
    end
    nodes[slack_idx] = SlackAlgebraic(U = complex(V_vec[slack_idx]))
    
    pg_cons = PowerGrid(nodes, lines)
    op = find_operationpoint(pg_cons, sol_method=:rootfind, solve_powerflow = true)
    return op
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

