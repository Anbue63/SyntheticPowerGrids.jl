"""
    test_power_grid(pg, op)

Performs a variety of test to assure that the dynamical system can represent a meaningful power grid.
"""
function test_power_grid(pg, op, pg_struct)
    if test_voltage(op) == true                             # Voltage magnitudes in the operation point
        if test_power_flow_on_lines(op, pg_struct) == true  # Power flow on the lines
            _, stable = test_eigenvalues(pg, op)            # Operation point is linearly stable
            if stable == true
                if test_slack_bus(pg, op) == true           # Power consumption of the slack bus
                    return true                             # All test have been passed. The power grid can be returned
                end
            end
        end
    end
    println("Restarting power grid generation algorithm.")
    println("")
    return false
end

"""
    test_power_flow_on_lines(state::State)

Calculates the power flow on the transmission lines of a grid. Checks if it is below the threshold of 70% of the physical limit.
"""
function test_power_flow_on_lines(state::State, pg_struct)
    lines = state.grid.lines
    save_flow = Vector{Bool}(undef, length(lines))

    for j in 1:length(lines)
        l = lines[j]
        m = l.from
        k = l.to
        
        if pg_struct.lines == :StaticLine
            g_mk = real(l.Y) # Line Conductance
            b_mk = imag(l.Y) # Line Susceptance
            
        elseif pg_struct.lines == :PiModelLine
            g_mk = real(l.y) # Line Conductance
            b_mk = imag(l.y) # Line Susceptance
        end

        v_m = state[m, :v]  # Voltage magnitude node m 
        v_k = state[k, :v]  # Voltage magnitude node k

        φ_mk = state[m, :φ] - state[k, :φ] # Phase difference

        P_mk = abs(v_k * v_m * b_mk * sin(φ_mk)) # Flow on the line connecting m, k
        P_save = abs(b_mk * 0.7)                 # Save flow on the line, 70% of the physically possible level

        save_flow[j] = P_mk < P_save
    end 
    save_network = all(save_flow)

    if save_network == false
        println("The power lines are overloaded.")
    end
    return save_network
end

"""
    test_voltage(op)    

Checks if all voltage magnitude are close to the correct operation voltage magnitude.
"""
function test_voltage(op)
    V = op[:, :v] # Voltage magnitudes in the operation point
    if all(isapprox.(V, 1.0, atol = 0.1))
        return true
    else 
        println("The voltage conditions for the power grid could not be met.")
        return false
    end
end

"""
    test_eigenvalues(pg::PowerGrid, s::State)

Calculates the eigenvalues λ of the jacobian of the right hand side of
the power grid at the state s.
The sign of the real part of the eigenvalues will decide whether a fixed point
will be attracting (λ_max < 0) or repelling (λ_max > 0).
"""
function test_eigenvalues(pg::PowerGrid, s::State)
    rpg = rhs(pg)
    M = Array(rpg.mass_matrix)
    f!(dx, x) = rpg(dx, x, nothing, 0.0)
    j(x) = (dx = similar(x); ForwardDiff.jacobian(f!, dx, x))
    λ = eigvals(j(s.vec) * pinv(M) * M) .|> real |> extrema
    stable = isapprox(last(λ), 0, atol=1e-8)

    if stable == false
        println("The eigenvalues of the jacobian have positive real parts.")
    end
    return λ, stable
end

"""
test_slack_bus(pg, op, threshold_power = 5.0)


Assures that the power consumed by slack bus is below a threshold.
"""
function test_slack_bus(pg, op, threshold_power = 5.0)
    slack_idx = findfirst(typeof.(pg.nodes) .== SlackAlgebraic)
    slack_test = abs(op[slack_idx, :p]) < threshold_power

    if slack_test == false
        println("The slack bus consumes to much power.")
    end
    return slack_test
end