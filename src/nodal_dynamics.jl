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

function get_nodes(P_set::Vector{Float64}, Q_set::Vector{Float64}, V_set::Vector{Float64}, slack::Bool, slack_idx::Int64, node_dynamics)
    num_nodes = length(P_set)        # Number of all nodes
    node_idxs = collect(1:num_nodes) 
    node_dynamics_idxs = [] 

    n_nodes = [round(Int, num_nodes * nd[1]) for nd in node_dynamics] # Fraction of nodes of type `nd`

    min_max_normalize_int_array!(num_nodes, n_nodes)  # Note that this can introduce nodes with weight 0. into the system.

    for n_node in n_nodes # Run over node types to distribute them
        s = sample(node_idxs, n_node, replace = false) # Sample nodes (without replacement) with should have type nd
        symdiff!(node_idxs, s) # Remove nodes `s` which were already sampled
        push!(node_dynamics_idxs, s)
    end

    nodes = Array{Any}(undef, num_nodes)
    # @assert isapprox(sum(x -> x[1], node_dynamics), 1) should go in validate pg not here!
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

function get_DroopControlledInverterApprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    τ_P = nodal_parameters[:τ_P] # Time constant low pass filter measuring the active power
    τ_Q = nodal_parameters[:τ_Q] # Time constant low pass filter measuring the reactive power
    K_P = nodal_parameters[:K_P] # Gain constant low pass filter measuring the active power
    K_Q = nodal_parameters[:K_Q] # Gain constant low pass filter measuring the reactive power
    
    P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n, x_dims = parameter_DroopControlledInverterApprox(P_set = P_set, Q_set = Q_set, τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_set, τ_P = τ_P, Y_n = 0.0)
    NormalForm(P = P, Q = Q, V = V, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, x_dims = x_dims)
end

function get_ThirdOrderMachineApprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    X = nodal_parameters[:X] # Reactance
    α = nodal_parameters[:α] # Voltage dynamics time constant
    γ = nodal_parameters[:γ] # Damping Coefficient

    E_f = V_set + (X * Q_set / (V_set)) # Electric field voltage that results in the correct nodal voltage magnitude
    P, Q, V, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n, x_dims = parameter_ThirdOrderMachineApprox(P_set = P_set, Q_set = Q_set, E_f = E_f, E_set = V_set, X = X, α = α, γ = γ, Y_n = 0.0)
    NormalForm(P = P, Q = Q, V = V, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, x_dims = x_dims)
end

function get_PQ(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters)
    PQAlgebraic(P = P_set, Q = Q_set)
end

function get_dVOCapprox(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    η = nodal_parameters[:η] # positive control parameter
    α = nodal_parameters[:α] # positive control parameter
    κ = nodal_parameters[:κ] # uniform complex phase

    Bᵤ, Cᵤ, Gᵤ, Hᵤ, Bₓ, Cₓ, Gₓ, Hₓ, Y_n, x_dims = parameter_dVOC(P_set = P_set, Q_set = Q_set, V_set = V_set, η = η, α = α, κ = κ, Y_n = 0.0)
    NormalForm(P = P, Q = Q, V = V, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, x_dims = x_dims)
end

function get_normalform(P_set::Float64, Q_set::Float64, V_set::Float64, nodal_parameters::Dict)
    x_dims = nodal_parameters[:x_dims] # Number of internal variables

    Bᵤ = nodal_parameters[:Bᵤ] 
    Cᵤ = nodal_parameters[:Cᵤ] 
    Gᵤ = nodal_parameters[:Gᵤ] 
    Hᵤ = nodal_parameters[:Hᵤ] 

    Bₓ = nodal_parameters[:Bₓ]  
    Cₓ = nodal_parameters[:Cₓ]  
    Gₓ = nodal_parameters[:Gₓ]  
    Hₓ = nodal_parameters[:Hₓ] 

    NormalForm(P = P_set, Q = Q_set, V = V_set, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, x_dims = x_dims)
end