"""
    get_ancillary_operationpoint(pg, P_vec, lines)

get an Ancillary Power Grid with the same power flow on the lines as the full power grid.
We use PV-Nodes to find the reactive power at the nodes with results in power grids where the voltage magnitude V is equal to the Reference voltage magnitude V_r.

- `pg`: Graph structure of the power grid
- `P_vec`: Vector containing the power consumption / generation at each node
- `lines`: Line dynamics of the power grid
- `V_r`: Reference voltage magnitude. In [p.u.] system the default is V_r = 1.0
"""
function get_ancillary_operationpoint(pg, P_vec, lines, V_r = 1.0)
    nodes = Array{Any}(undef, nv(pg.graph))

    # Get Ancillary Power Grid
    for n in 1:nv(pg.graph)-1
        nodes[n] = PVAlgebraic(P = P_vec[n], V = V_r)
    end
    nodes[end] = SlackAlgebraic(U = complex(1.0))
    
    pg_cons = PowerGrid(nodes, lines)
    # Find the operation point of Ancillary grid -> we need the reactive power at the nodes
    op = find_operationpoint(pg_cons, sol_method=:rootfind)#, solve_powerflow = true)
    return op
end

"""
    get_initial_guess(rpg, op_ancillary)

Initial guess for the operation point search of PowerDynamics.jl for the Operation point of the helper grid with the correct power flows and voltage magnitudes.
- `rpg`: Right-Hand side function of the power grid
- `op_ancillary`: Operation point of the helper grid
"""
function get_initial_guess(rpg, op_ancillary)
    
    ic_guess = zeros(length(rpg.syms))
    u_r_idx = findall(map(x -> occursin("u_r", string(x)), rpg.syms))
    u_i_idx = findall(map(x -> occursin("u_i", string(x)), rpg.syms))
    ic_guess[u_r_idx] = op_ancillary[:, :u_r]
    ic_guess[u_i_idx] = op_ancillary[:, :u_i]

    return ic_guess
end

