using NetworkDynamics: ODEVertex
import PowerDynamics: dimension, symbolsof, construct_vertex 
import PowerDynamics: showdefinition

@DynamicNode PQDynamic(P_0, Q_0, τ) begin
    @assert τ > 0 "Time constant (τ) should be > 0"
end [] begin
    s = u * conj(i)
    du = (1/τ) * (complex(P_0, Q_0) - s)
end

export PQDynamic