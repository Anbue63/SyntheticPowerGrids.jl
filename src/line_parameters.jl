"""
    impedance_380kV(length, ω = 2π * 50)

Tabelle 5.6 Standardfreileitungen in der 220-kV- und 380-kV-Ebene 
https://www.dena.de/newsroom/publikationsdetailansicht/pub/dena-verteilnetzstudie-ausbau-und-innovationsbedarf-der-stromverteilnetze-in-deutschland-bis-2030/
"""
function impedance_380kV(length, ω = 2π * 50)
    R = 0.025 * length
    X = 0.25  * length
    C_shunt = 13.7 * 10^-9 * length
    Y_shunt = 1im * C_shunt * ω
    Z = R + 1im * X
    return Z, Y_shunt
end

function get_line_admittance_matrix(g::EmbeddedGraph{Int64}, mean_len_km = 42.872746445497626, shortest_line_km = 0.06)
    dist_nodes = EmbeddedGraphs.weights(g) # Euclidean distance of the edges in EmbeddedGraphs

    # Remove all "unconnected" distances!
    dist_nodes_connected = vcat(dist_nodes...)
    unconnected_idx = findall(iszero, dist_nodes_connected) # Unconnected nodes have a length of d = 0.0
    deleteat!(dist_nodes_connected, unconnected_idx) 

    dist_nodes_mean = mean(dist_nodes_connected) # Mean length with respect to connected nodes (otherwise we skew the mean!)
    dist_to_km = mean_len_km / dist_nodes_mean   # Conversion factor from euclidean distance to [km]
    
    N = nv(g)                                 # Number of nodes
    Y = zeros(Complex{Float64}, N, N)         # Admittance Matrix
    Y_shunt = zeros(Complex{Float64}, N, N)   # Admittance Matrix
    for i in 1:N
        for j in 1:i-1
            if dist_nodes[i, j] > 0 && dist_nodes[j, i] > 0 # If nodes are connected -> Calculate Admittance
                len_line = dist_nodes[i,j] * dist_to_km     # Euclidean distance conversion to km
                if len_line < shortest_line_km              # SyntheticNetworks can generate very short lines. We fix this by adding a threshold to the admittance. The shortest line length is taken from SciGrids
                    impedance_line, shunt_admittance_line = impedance_380kV(shortest_line_km)  
                else
                    impedance_line, shunt_admittance_line = impedance_380kV(len_line)  # Total impedance of the line in Ohm
                end
                # Conversion to Admittance
                Y[i, j] = 1 / impedance_line 
                Y[j, i] = 1 / impedance_line

                Y_shunt[i, j] = shunt_admittance_line / 2 # Shunt is distributed: 1/2 at the beginning, 1/2 at end of the line
                Y_shunt[j, i] = shunt_admittance_line / 2
            else # Unconnected Edges don't need an admittance
                Y[i, j] = 0 + 0im
                Y[j, i] = 0 + 0im

                Y_shunt[i, j] = 0 + 0im
                Y_shunt[j, i] = 0 + 0im
            end
        end
    end
    return Y, Y_shunt
end