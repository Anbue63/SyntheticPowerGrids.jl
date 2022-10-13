"""
    line_properties_380kV(length, ω = 2π * 50)

Tabelle 5.6 Standardfreileitungen in der 220-kV- und 380-kV-Ebene 
https://www.dena.de/newsroom/publikationsdetailansicht/pub/dena-verteilnetzstudie-ausbau-und-innovationsbedarf-der-stromverteilnetze-in-deutschland-bis-2030/
"""
function line_properties_380kV(length, cables, ω = 2π * 50)
    c_wire = 1 #wires / wires_typical
    c_cable = cables / 3

    R = 0.025 * length / (c_cable * c_wire)
    X = 0.25  * length / (c_cable * c_wire)
    C_shunt = 13.7 * 10^-9 * length * (c_cable * c_wire)

    Y_shunt = 1im * C_shunt * ω
    Z = R + 1im * X
    return Z, Y_shunt
end

"""
    get_line_admittance_matrix(embedded_graph::EmbeddedGraph{Int64}; L_matrix::Matrix{Float64}, cables_vec, num_nodes::Int64)

Calculates the line admittances and the shunts with respect to the line length.
"""
function get_line_admittance_matrix(embedded_graph::EmbeddedGraph{Int64}; L_matrix::Matrix{Float64}, cables_vec, num_nodes::Int64)
    Y = zeros(Complex{Float64}, num_nodes, num_nodes)       # Admittance Matrix
    Y_shunt = zeros(Complex{Float64}, num_nodes, num_nodes) # Shunt Admittance Matrix

    e = edges(embedded_graph.graph)
    dest_vec = dst.(e)
    source_vec = src.(e)

    for i in eachindex(source_vec) # If nodes are connected -> Calculate Admittance
        source = source_vec[i]
        dest = dest_vec[i]

        impedance_line, shunt_admittance_line = line_properties_380kV(L_matrix[source, dest], cables_vec[i])  # Total impedance of the line in Ohm
                
        # Conversion to Admittance
        Y[source, dest] = 1 / impedance_line 
        Y[dest, source] = 1 / impedance_line

        Y_shunt[source, dest] = shunt_admittance_line / 2 # Shunt is distributed: 1/2 at the beginning, 1/2 at end of the line
        Y_shunt[dest, source] = shunt_admittance_line / 2
    end
    return Y, Y_shunt
end

"""
    get_effective_distances(g::EmbeddedGraph{Int64}; mean_len_km, shortest_line_km)

Calculates the geographic distances in [km] from an embedded graph.
"""
function get_effective_distances(g::EmbeddedGraph{Int64}; mean_len_km, shortest_line_km)
    dist_nodes = EmbeddedGraphs.weights(g) # Euclidean distance of the edges in EmbeddedGraphs

    # Remove all "unconnected" distances!
    dist_nodes_connected = vcat(dist_nodes...)
    unconnected_idx = findall(iszero, dist_nodes_connected) # Unconnected nodes have a length of d = 0.0
    deleteat!(dist_nodes_connected, unconnected_idx) 

    dist_nodes_mean = mean(dist_nodes_connected) # Mean length with respect to connected nodes (otherwise we skew the mean!)
    dist_to_km = mean_len_km / dist_nodes_mean   # Conversion factor from euclidean distance to [km]
    
    N = nv(g) 
    L_matrix = zeros(Float64, N, N)  # Line Length Matrix

    for i in 1:N
        for j in 1:i-1
            len_line = dist_nodes[i,j] * dist_to_km # Euclidean distance conversion to km
            if len_line < shortest_line_km          # SyntheticNetworks can generate very short lines. We fix this by adding a threshold to the admittance. The shortest line length is taken from SciGrids
                L_matrix[i, j] = shortest_line_km
                L_matrix[j, i] = shortest_line_km
            else
                L_matrix[i, j] = len_line
                L_matrix[j, i] = len_line
            end
        end
    end
    return L_matrix
end