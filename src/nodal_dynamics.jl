"""
    get_nodes(P_set::Vector{Float64}, Q_set::Vector{Float64}, V_set::Vector{Float64}, slack::Bool, slack_idx::Int64, node_dynamics)

Generates the nodal dynamics of the network.

- `P_set::Vector{Float64}`: Active Power Set-Points
- `Q_set::Vector{Float64}`: Reactive Power Set-Points
- `V_set::Vector{Float64}`: Voltage Set-Points
- `slack::Bool`: Activates or deactivates the slack bus
- `slack_idx::Int64`: Gives the node index of the slack bus if one is used
- `node_dynamics`: Stores all the different node dynamics. The first arg gives the share the second gives the node type and the third stores the nodal parameters
"""
function get_nodes(P_set::Vector{Float64}, Q_set::Vector{Float64}, V_set::Vector{Float64}, slack::Bool, slack_idx::Int64, node_dynamics)
    num_nodes = length(P_set)        # Number of all nodes
    node_idxs = collect(1:num_nodes) 
    node_dynamics_idxs = [] 

    n_nodes = [round(Int, num_nodes * nd[1]) for nd in node_dynamics] # Fraction of nodes of type `nd`

    # Catch rounding artifacts!
    min_max_normalize_int_array!(num_nodes, n_nodes) # Note that this can introduce nodes with weight 0. into the system.

    for n_node in n_nodes # Run over node types to distribute them
        s = sample(node_idxs, n_node, replace = false) # Sample nodes (without replacement) with should have type nd
        symdiff!(node_idxs, s) # Remove nodes `s` which were already sampled
        push!(node_dynamics_idxs, s)
    end

    nodes = Array{Any}(undef, num_nodes)
    for (nd, idxs) in zip(node_dynamics, node_dynamics_idxs)
        for i in idxs
            nodal_parameters = nd[3]
            nodes[i] = nd[2](P_set[i], Q_set[i], V_set[i], nodal_parameters)
        end
    end

    if slack == true # Adds a slack bus to the grid, the default option is not to use a slack
        nodes[slack_idx] = SlackAlgebraic(U = complex(V_set[slack_idx]))
    end
    return nodes
end

"""
    min_max_normalize_int_array!(sum_arr, arr)

During the rounding the sum_arr might be unequal to sum(arr). 
This function finds the largest / smallest entry in arr and removes / adds +1 until sum_arr = sum(arr).
If sum(arr) == sum_arr this function does nothing
"""
function min_max_normalize_int_array!(sum_arr, arr)
    while sum(arr) > sum_arr
        _, i = findmax(arr)
        arr[i] -= 1
    end

    while sum(arr) < sum_arr
        _, i = findmin(arr)
        arr[i] += 1
    end
    nothing
end

function parameter_DroopControlledInverterApprox(;τ_P, τ_Q, K_P, K_Q, V_r, Y_n)
    Bᵤ = 1im 
    Cᵤ = - 1 / (2 * τ_Q * V_r^2)
    Gᵤ = 0
    Hᵤ = - K_Q / (τ_Q * V_r)
    Bₓ = - 1 / τ_P
    Cₓ = 0
    Gₓ = - K_P / τ_P
    Hₓ = 0 

    return [Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n]
end

"""
   get_DroopControlledInverterApprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)

Normal Form Approximation of a Droop Controlled Inverter.   
"""
function get_DroopControlledInverterApprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    τ_P = nodal_parameters[:τ_P] # Time constant low pass filter measuring the active power
    τ_Q = nodal_parameters[:τ_Q] # Time constant low pass filter measuring the reactive power
    K_P = nodal_parameters[:K_P] # Gain constant low pass filter measuring the active power
    K_Q = nodal_parameters[:K_Q] # Gain constant low pass filter measuring the reactive power
    
    Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n = parameter_DroopControlledInverterApprox(τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_set, τ_P = τ_P, Y_n = 0.0)
    NormalForm(P = P_set, Q = Q_set, V = V_set, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ)
end

function parameter_ThirdOrderMachineApprox(;E_f, E_set, X, α, γ, Y_n)
    Bᵤ = 1im
    Cᵤ = (- E_f / (2E_set^3)) / α
    Gᵤ = 0.0
    Hᵤ = - X / (α * (E_set)^2)
    Bₓ = - γ
    Cₓ = 0.0 
    Gₓ = -1
    Hₓ = 0

    return [Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n]
end

"""
    get_ThirdOrderMachineApprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)

Normal Form Approximation of a Third Order Machine.   
"""
function get_ThirdOrderMachineApprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    X = nodal_parameters[:X] # Reactance
    α = nodal_parameters[:α] # Voltage dynamics time constant
    γ = nodal_parameters[:γ] # Damping Coefficient

    E_f = V_set + (X * Q_set / (V_set)) # Electric field voltage that results in the correct nodal voltage magnitude
    Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n = parameter_ThirdOrderMachineApprox(E_f = E_f, E_set = V_set, X = X, α = α, γ = γ, Y_n = 0.0)
    NormalForm(P = P_set, Q = Q_set, V = V_set, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ)
end

function get_PQ(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters)
    PQAlgebraic(P = P_set, Q = Q_set)
end

function parameter_dVOC(;P_set, Q_set, V_set, η, α, κ, Y_n)
    Bᵤ = nothing
    Cᵤ = - α * η / (V_set^2) + η * exp(κ * 1im) * complex(P_set, - Q_set) / (V_set^4)
    Gᵤ = - η * exp(κ * 1im) / (V_set^2)
    Hᵤ = 1im * η * exp(κ * 1im) / (V_set^2)
    Bₓ = nothing
    Cₓ = nothing
    Gₓ = nothing
    Hₓ = nothing

    return [Cᵤ, Gᵤ, Hᵤ, Y_n]
end

"""
    get_dVOCapprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)

Normal Form Approximation of the dispatchable virtual oscillator control.
"""
function get_dVOCapprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    η = nodal_parameters[:η] # positive control parameter
    α = nodal_parameters[:α] # positive control parameter
    κ = nodal_parameters[:κ] # uniform complex phase

    Cᵤ, Gᵤ, Hᵤ, Y_n = parameter_dVOC(P_set = P_set, Q_set = Q_set, V_set = V_set, η = η, α = α, κ = κ, Y_n = 0.0)
    NormalForm(P = P_set, Q = Q_set, V = V_set, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ)
end

"""
    get_normalform(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)

General Normal Form.
"""
function get_normalform(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)

    Bᵤ = nodal_parameters[:Bᵤ] 
    Cᵤ = nodal_parameters[:Cᵤ] 
    Gᵤ = nodal_parameters[:Gᵤ] 
    Hᵤ = nodal_parameters[:Hᵤ] 

    Bₓ = nodal_parameters[:Bₓ]  
    Cₓ = nodal_parameters[:Cₓ]  
    Gₓ = nodal_parameters[:Gₓ]  
    Hₓ = nodal_parameters[:Hₓ] 

    NormalForm(P = P_set, Q = Q_set, V = V_set, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ)
end