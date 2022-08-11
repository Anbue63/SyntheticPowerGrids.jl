using NetworkDynamics: ODEVertex
import PowerDynamics: dimension, symbolsof, construct_vertex 
import PowerDynamics: showdefinition

function parameter_schiffer(;P_set, Q_set, τ_P = 5.0, τ_Q = 8.0, K_P = 5.0, K_Q = 0.1, V_r = 1.0, Y_n = 0.0)
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

    return [Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n]
end

function parameter_schmietendorf(;P_m, Y_n, E_f = 1.0, E_c = 1.0, X = 1.0, Q_c = 0, α = 2, γ = 0.2) 
    Aᵤ = ((3E_f / 2E_c) - 1 + X * Q_c / (E_c)^2 ) / α
    Bᵤ = 1im 
    Cᵤ = (- E_f / (2E_c^3) + X * Q_c / (E_c)^4 ) / α
    Gᵤ = 0.0
    Hᵤ = - X / (α * (E_c)^2)
    Mₓ = 1.0
    Aₓ = P_m 
    Bₓ = - γ
    Cₓ = 0.0 
    Gₓ = -1
    Hₓ = 0

    return [Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n]
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
end

NormalForm(; Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n = 0) = NormalForm(Aᵤ, Bᵤ, Cᵤ, Gᵤ, Hᵤ, Aₓ, Bₓ, Cₓ, Gₓ, Hₓ, Mₓ, Y_n)

function construct_vertex(nf::NormalForm)
    sym = symbolsof(nf)
    dim = dimension(nf)
    mass_matrix = ones(Int64,3,3) |> Diagonal
    
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

    function rhs!(dx, x, edges, p, t)
        i = total_current(edges) + Y_n * (x[1] + x[2] * 1im) # Current, couples the NF to the rest of the network, Y_n Shunt admittance
        u = x[1] + x[2] * im  # Complex Voltage
        ω = x[3]              # Frequency (or any other internal variable x)
        s = u * conj(i)       # Apparent Power S = P + iQ
        v2 = abs2(u)          # Absolute squared voltage
    
        dω = (Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s)) / Mₓ
        du = (Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s)) * u
        
        # Splitting the with the complex parameters
        dx[1] = real(du)  
        dx[2] = imag(du)
        dx[3] = real(dω)
        return nothing
    end

    ODEVertex(f! = rhs!, dim = dim, mass_matrix = mass_matrix, sym = sym)

end

symbolsof(nf::NormalForm) = [:u_r, :u_i, :ω]

dimension(nf::NormalForm) = 3

export NormalForm