
node_types_default = [(1/3, get_PQ, nothing), (1/3, get_DroopControlledInverterApprox, Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5)), (1/3, get_DroopControlledInverterApprox, Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => 0.5))]

function get_nodes2(P_set::Vector{Float64}, Q_set::Vector{Float64}, V_set::Vector{Float64}, slack::Bool, slack_idx::Int64, node_types = node_types_default)
    num_nodes = length(P_set)        # Number of all nodes
    node_idxs = collect(1:num_nodes) 
    node_types_idxs = [] 

    for nt in node_types 
        n_nodes = round(Int, nt[1] * num_nodes)
        s = sample(node_idxs, n_nodes, replace = false)
        symdiff!(node_idxs, s)
        push!(node_types_idxs, s)
    end

    if ! (node_idxs == []) # Check if all nodes were distributed to the different types
        append!(node_types_idxs[end], node_idxs)
    end

    nodes = Array{Any}(undef, num_nodes)
    # @assert isapprox(sum(x -> x[1], node_types), 1) should go in validate pg not here!
    for (nt, idxs) in zip(node_types, node_types_idxs)
        for i in idxs
            nodal_parameters = nt[3]
            nodes[i] = nt[2](P_set[i], Q_set[i], V_set[i], nodal_parameters)
        end
    end

    if slack == true
        nodes[slack_idx] = SlackAlgebraic(U = complex(V_set[slack_idx]))
    end
    return nodes
end

function get_DroopControlledInverterApprox(P_set, Q_set, V_set, nodal_parameters)
    τ_P = nodal_parameters[:τ_P] # Time constant low pass filter measuring the active power
    τ_Q = nodal_parameters[:τ_Q] # Time constant low pass filter measuring the reactive power
    K_P = nodal_parameters[:K_P] # Gain constant low pass filter measuring the active power
    K_Q = nodal_parameters[:K_Q] # Gain constant low pass filter measuring the reactive power
    
    Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, x_dims = parameter_DroopControlledInverterApprox(P_set = P_set, Q_set = Q_set, τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = V_set, τ_P = τ_P, Y_n = 0.0)
    NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
end

function get_ThirdOrderMachineApprox(P_set, Q_set, V_set, nodal_parameters)
    X = nodal_parameters[:X] # Reactance
    α = nodal_parameters[:α] # Voltage dynamics time constant
    γ = nodal_parameters[:γ] # Damping Coefficient

    E_f = V_set + (X * Q_set / (V_set)) # Electric field voltage that results in the correct nodal voltage magnitude
    Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, x_dims = parameter_ThirdOrderMachineApprox(P_m = P_set, E_f = E_f, E_set = V_set, X = X, α = α, γ = γ, Y_n = 0.0)
    NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
end

function get_PQ(P_set, Q_set, V_set, nodal_parameters)
    PQAlgebraic(P = P_set, Q = Q_set)
end

function get_swingLVS(P_set, Q_set, V_set, nodal_parameters)
    H = nodal_parameters[:H] # Inertia Constant
    Ω = nodal_parameters[:Ω] # Rated Frequency
    D = nodal_parameters[:D] # Damping Coefficient
    Γ = nodal_parameters[:Γ] # Voltage stability Coefficient

    SwingEqLVS(H = H, P = P_set, D = D, Ω = Ω, Γ = Γ, V = V_set)
end

function get_dVOCapprox(P_set, Q_set, V_set, nodal_parameters)
    Ω = nodal_parameters[:Ω] # Rated Frequency
    η = nodal_parameters[:η] # positive control parameter
    α = nodal_parameters[:α] # positive control parameter
    κ = nodal_parameters[:κ] # uniform complex phase

    Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, x_dims = parameter_dVOC(P_set = P_set, Q_set = Q_set, V_set = V_set, Ω_set = Ω, η = η, α = α, κ = κ, Y_n = 0.0)
    NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
end

function get_normalform(P_set, Q_set, V_set, nodal_parameters)
    x_dims = nodal_parameters[:x_dims] # Number of internal variables

    Aᵤ = nodal_parameters[:Aᵤ] 
    Bᵤ = nodal_parameters[:Bᵤ] 
    Cᵤ = nodal_parameters[:Cᵤ] 
    Gᵤ = nodal_parameters[:Gᵤ] 
    Hᵤ = nodal_parameters[:Hᵤ] 

    Aₓ = nodal_parameters[:Aₓ] 
    Bₓ = nodal_parameters[:Bₓ]  
    Cₓ = nodal_parameters[:Cₓ]  
    Gₓ = nodal_parameters[:Gₓ]  
    Hₓ = nodal_parameters[:Hₓ] 
    Mₓ = nodal_parameters[:Mₓ] 

    NormalForm(Aᵤ = Aᵤ, Bᵤ = Bᵤ, Cᵤ = Cᵤ, Gᵤ = Gᵤ, Hᵤ = Hᵤ, Aₓ = Aₓ, Bₓ = Bₓ, Cₓ = Cₓ, Gₓ = Gₓ, Hₓ = Hₓ, Mₓ = Mₓ, x_dims = x_dims)
end