using Revise
using SyntheticPowerGrids

nodal_parameters_a = Dict(:τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) # This option recover the GNN dataset grids, please do not delete it
nodal_parameters_b = Dict(:X => 1.0, :γ => 0.2, :α => 2.0) 
nodal_parameters_c = Dict(:X => 1.0, :γ => 0.2, :α => 2.0, :τ_Q => 8.0, :K_P => 5, :K_Q => 0.1, :V_r => 1.0, :τ_P => [0.5, 1.0 , 5.0]) 

nodal_shares_a = Dict(:load_share => 0.5, :schiffer_share => 0.5)
nodal_shares_b = Dict(:load_share => 0.5, :schmietendorf_share => 0.5)
nodal_shares_c = Dict(:load_share => 0.5, :schmietendorf_share => 0.25, :schiffer_share => 0.25)

a = PGGeneration(num_nodes = 100, nodal_parameters = nodal_parameters_a, loads = :PQAlgebraic, lines = :StaticLine, nodal_shares = nodal_shares_a)
b = PGGeneration(num_nodes = 100, nodal_parameters = nodal_parameters_b, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :Schmietendorf, nodal_shares = nodal_shares_b)
c = PGGeneration(num_nodes = 100, nodal_parameters = nodal_parameters_c, loads = :PQAlgebraic, lines = :PiModelLine, generation_dynamics = :Mixed, nodal_shares = nodal_shares_c)

##
pg, op, embedded_graph, rejections = random_PD_grid(c)
get_effective_distances(embedded_graph, mean_len_km = c.mean_len_km, shortest_line_km = c.shortest_line_km)

##
using OrdinaryDiffEq
using PowerDynamics
using Plots

rpg = rhs(pg)

ω_idx = findall(map(x -> 'ω' == string(x)[1], rpg.syms))

x0 = copy(op)

x0.vec[ω_idx[1]] = 1

prob = ODEProblem(rpg, x0.vec, (0.0, 100.0), nothing)
sol = solve(prob, Rodas4())

plot(sol, idxs = ω_idx)