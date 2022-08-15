"""
    test_power_grid(pg, op)

Performs a variety of test to assure that the dynamical system can represent a meaningful power grid.
"""
function test_power_grid(pg, op, pg_struct)
    if test_voltage(op) == true                             # Voltage magnitudes in the operation point
        if test_power_flow_on_lines(op, pg_struct) == true  # Power flow on the lines
            #stable = test_eigenvalues(pg, op)             
            stable = small_signal_stability_analysis(rhs(pg), op.vec) # Operation point is linearly stable
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

    for j in eachindex(lines)
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
    return stable
end


"""
    small_signal_stability_analysis(h::ODEFunction, eq_point, p = nothing)

Performs a small signal stability analysis according to: 
[1] "Power System Modelling and Scripting", F. Milano, Chapter 7.2.
We find the reduced Jacobian (or State Matrix A_s) and calculate its eigenvalues.
If the eigenvalues have positive real parts we classify the grid as unstable. 

- `h`: Full DAE system
- `eq_point`: Equilibrium point / fixed point of h. h(eq_point) = 0.0
"""
function small_signal_stability_analysis(h::ODEFunction, eq_point, p = nothing)
    h!(dx, x) = h(dx, x, p, 0.0)
    
    j(x) = (dx = similar(x); ForwardDiff.jacobian(h!, dx, x)) # Jacobian
    J = j(eq_point) # Full Jacobian at the equilibrium point

    c_idx, d_idx = separate_differential_constraint_eqs(h, p) # Constraint and differential indices 

    f_x = J[d_idx, d_idx] # Differential equations evaluated at the differential variables
    f_y = J[d_idx, c_idx] # Differential equations evaluated at the constrained variables

    g_x = J[c_idx, d_idx] # Constrained equations evaluated at the differential variables
    g_y = J[c_idx, c_idx] # Constrained equations evaluated at the constrained variables

    D = f_y * pinv(g_y) * g_x # Degradation matrix
    A_s = f_x - D             # State matrix / Reduced Jacobian (eq. 7.16 in [1])

    λ = eigvals(A_s) .|> real |> extrema # Eigenvalues of the state matrix

    if all(λ .< 0.0)
        stable = true
    else
        stable = isapprox(last(λ), 0, atol=1e-8)
    end

    if stable == false
        println("The eigenvalues of the reduced jacobian have positive real parts.")
    end

    return stable
end

"""
    separate_differential_constraint_eqs(DAE::ODEFunction, p=nothing)
Returns the constraint equations and differential equations indices from an ODEFunction h(x) used in DifferentialEquations.jl.
The ODE h must be in Mass Matrix form meaning: M ẋ = h(x), with M diagonal. h should be inplace.
"""
function separate_differential_constraint_eqs(h::ODEFunction, p=nothing)
    M = h.mass_matrix
    M == I && error("There are no constraints in the system!")
    M != Diagonal(M) && error("The constraints are not diagonal.")
    
    c_idx = findall(diag(M) .== 0)
    d_idx = findall(diag(M) .== 1)

    return c_idx, d_idx  
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