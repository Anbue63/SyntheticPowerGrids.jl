"""
    random_PD_grid(N::Int)

Generates a random power grid using SyntheticNetworks and then turns it into a PowerDynamics.PowerGrid type.
"""
function random_PD_grid(N::Int; τ_Q = 8.0, K_P = 5, K_Q = 0.1, P0 = 1.0, P_base = 400 * 10^6, V_base = 380 * 10^3)
    rejections = 0
    Y_base = P_base / (V_base)^2 # Base Admittance

    power_dist = MixtureModel(Normal[Normal(P0, P0/2), Normal(-P0, P0/2)])  # Bimodal Distribution for the active power

    for i in  1:1000 # maxiters until a stable grid is found 
        pg = generate_graph(RandomPowerGrid(N, 1, 1/5, 3/10, 1/3, 1/10, 0.0)) # Generates the random power grid topology
        
        P_vec = rand(power_dist, N - 1) # Generate Power Generation / Consumption at the nodes
        P_vec .-= sum(P_vec) / (N - 1)  # Assure power balance

        Y = generate_line_admittance_matrix(pg) # Line admittance matrix, Entry's are in [Ohm]
                
        lines = generate_lines(pg, Y, Y_base)                     # Generate line dynamics
        op_ancillary = generate_ancillary_grid(pg, P_vec, lines)  # Operation point of Ancillary power grid
        
        if test_power_flow_on_lines(op_ancillary) == true             # Check if power flow on the lines is on a save limit
            nodes = generate_nodes(pg, op_ancillary, τ_Q, K_P, K_Q)   # Generate nodal dynamics

            pg = PowerGrid(nodes, lines)
            rpg = rhs(pg)

            ic_guess = generate_initial_guess(rpg, op_ancillary) # Generate an initial guess for rootfind
            op = find_operationpoint(pg, ic_guess, sol_method = :rootfind)#, solve_powerflow = true) # find operation point of the full power grid

            if test_power_grid(pg, op) == true
                return pg, op, rejections
            end
        end
        rejections += 1
    end
end 


"""
   generate_lines(pg, Y, Y_base)

Generate the lines dynamics of the power grid.
- `pg`: Graph structure of the power grid
- `Y`:  Line admittance matrix, Entry's are in [Ohm] and depend on the line length
- `Y_base`: Base admittance of the power system for [p.u.] conversion
"""
function generate_lines(pg, Y, Y_base)
    e = edges(pg.graph)
    from_vec = src.(e)
    to_vec = dst.(e)
    
    lines = Array{Any}(undef, length(e))

    for l in 1:length(e)
        Y_pu = Y[from_vec[l], to_vec[l]] / Y_base # Convert admittance [ohm] to p.u. system
        lines[l] = StaticLine(from = from_vec[l], to = to_vec[l], Y = Y_pu) 
    end

    return lines
end


"""
    generate_nodes(pg, op_ancillary, τ_Q, K_P, K_Q, V_r = 1.0)

Generate the nodal dynamics of the full power grid.
- `pg`: Graph structure of the power grid
- `op_ancillary`: Operation point of the helper grid
- `τ_Q`: Time constant low pass filter measuring the reactive power
- `K_P`: Gain constant low pass filter measuring the active power
- `K_Q`: Gain constant low pass filter measuring the reactive power
- `V_r`: Reference voltage magnitude. In [p.u.] system the default is V_r = 1.0
"""
function generate_nodes(pg, op_ancillary, τ_Q, K_P, K_Q, V_r = 1.0)
    nodes = Array{Any}(undef, nv(pg.graph))
    τ_P_options = [0.5, 1.0 , 5.0]
        for n in 1:nv(pg.graph)
            α = rand() # Randomly turns nodes into PQAlgebraic or PhaseAmplitudeOscillators
            if α > 0.5 # Grid Forming
                β = rand(1:3) # Randomly chooses one of three possible time constant for the low pass filter measuring the active power
                τ_P = τ_P_options[β]

                nodes[n] = SchifferApprox(τ_P = τ_P, τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_r, P = op_ancillary[n, :p], Q = op_ancillary[n, :q], Y_n = 0)
            else       # Grid Following
                nodes[n] = PQAlgebraic(P = op_ancillary[n, :p], Q = op_ancillary[n, :q]) 
            end
        end
        nodes[end] = SlackAlgebraic(U = complex(V_r))
    
    return nodes
end