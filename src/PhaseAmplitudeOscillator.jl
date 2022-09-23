using NetworkDynamics: ODEVertex
import PowerDynamics: dimension, symbolsof, construct_vertex 
import PowerDynamics: showdefinition

function parameter_schiffer(;P_set, Q_set, τ_P, τ_Q, K_P, K_Q, V_r, Y_n)
    Aᵤ = (V_r + 2 * K_Q * Q_set) / (2 * τ_Q * V_r)
    Bᵤ = 1im 
    Cᵤ = - 1 / (2 * τ_Q * V_r^2)
    Gᵤ = 0
    Hᵤ = - K_Q / (τ_Q * V_r)
    Aₓ = K_P * P_set 
    Bₓ = - 1 
    Cₓ = 0
    Gₓ = - K_P 
    Hₓ = 0 
    Mₓ = τ_P
    xdims = 1

    return [Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims]
end

function parameter_schmietendorf(;P_m, E_f, E_set, X, α, γ, Y_n) 
    Aᵤ = ((3E_f / 2E_set) - 1) / α
    Bᵤ = 1im 
    Cᵤ = (- E_f / (2E_set^3)) / α
    Gᵤ = 0.0
    Hᵤ = - X / (α * (E_set)^2)
    Aₓ = P_m 
    Bₓ = - γ
    Cₓ = 0.0 
    Gₓ = -1
    Hₓ = 0
    Mₓ = 1.0
    xdims = 1

    return [Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims]
end

function parameter_dVOC(;P_set, Q_set, V_set, Ω_set, η, α, κ, Y_n)
    #ToDo: Is Ω_set always 0 ?
    Cᵤ = - α * η / (V_set^2) + η * exp(κ * 1im) * complex(P_set, - Q_set) / (V_set^4)
    Gᵤ = - η * exp(κ * 1im) / (V_set^2)
    Hᵤ = 1im * η * exp(κ * 1im) / (V_set^2)
    Aᵤ = 1im * Ω_set  - Cᵤ * (V_set)^2 - Gᵤ * P_set - Hᵤ * Q_set

    Bᵤ = 0
    
    Aₓ = []
    Bₓ = []
    Cₓ = []
    Gₓ = [] 
    Hₓ = [] 
    Mₓ = []
    xdims = 0

    return [Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims]
end

struct NormalForm <: AbstractNode
    Aᵤ
    Bᵤ
    Cᵤ
    Gᵤ
    Hᵤ
    Aₓ
    Bₓ
    Cₓ
    Gₓ
    Hₓ
    Mₓ
    Y_n # Shunt admittance for power dynamics
    xdims
end

NormalForm(; Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n = 0, xdims) = NormalForm(Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n, xdims)

function construct_vertex(nf::NormalForm)
    sym = symbolsof(nf)
    dim = dimension(nf)
    mass_matrix = ones(Int64,dim,dim) |> Diagonal
    
    Aᵤ = nf.Aᵤ
    Bᵤ = nf.Bᵤ
    Cᵤ = nf.Cᵤ
    Gᵤ = nf.Gᵤ
    Hᵤ = nf.Hᵤ
    
    Aₓ = nf.Aₓ
    Bₓ = nf.Bₓ
    Cₓ = nf.Cₓ
    Gₓ = nf.Gₓ
    Hₓ = nf.Hₓ
    Mₓ = nf.Mₓ

    Y_n = nf.Y_n

    if xdims > 0

        @assert length(Aₓ) == nf.xdims
        @assert length(Bₓ) == nf.xdims
        @assert length(Cₓ) == nf.xdims
        @assert length(Gₓ) == nf.xdims
        @assert length(Hₓ) == nf.xdims
        @assert length(Mₓ) == nf.xdims

        function rhs!(dz, z, edges, p, t)
            i = total_current(edges) + Y_n * (z[1] + z[2] * 1im) # Current, couples the NF to the rest of the network, Y_n Shunt admittance
            u = z[1] + z[2] * im  # Complex Voltage
            x = z[3:dim]          # Internal Variables
            s = u * conj(i)       # Apparent Power S = P + iQ
            v2 = abs2(u)          # Absolute squared voltage
        
            dx = (Aₓ + Bₓ .* x + Cₓ .* v2 + Gₓ .* real(s) + Hₓ .* imag(s)) ./ Mₓ
            du = (Aᵤ + Bᵤ  * x + Cᵤ  * v2 + Gᵤ  * real(s) + Hᵤ  * imag(s)) * u
            
            # Splitting the with the complex parameters
            dz[1] = real(du)  
            dz[2] = imag(du)
            dz[3:dim] = real(dx)
            return nothing
        end

    else

        function rhs!(dz, z, edges, p, t)
            i = total_current(edges) + Y_n * (z[1] + z[2] * 1im) # Current, couples the NF to the rest of the network, Y_n Shunt admittance
            u = z[1] + z[2] * im  # Complex Voltage
            s = u * conj(i)       # Apparent Power S = P + iQ
            v2 = abs2(u)          # Absolute squared voltage
        
            du = (Aᵤ + Cᵤ  * v2 + Gᵤ  * real(s) + Hᵤ  * imag(s)) * u
            
            # Splitting the with the complex parameters
            dz[1] = real(du)  
            dz[2] = imag(du)
            return nothing
        end

    end

    ODEVertex(rhs!, dim, mass_matrix, sym)

end

function symbolsof(nf::NormalForm)
    symbols = [:u_r, :u_i]
    append!(symbols, [Symbol("x_$i") for i in 1:nf.xdims])
end

dimension(nf::NormalForm) = 2 + nf.xdims

export NormalForm