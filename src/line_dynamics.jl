"""
    get_lines(pg_struct)

get the lines dynamics of the power grid.
- `pg_struct`: Struct containing all data of the power grid
"""
function get_lines(pg_struct::PGGeneration)
    embedded_graph = pg_struct.embedded_graph

    # Different possible coupling types
    if pg_struct.coupling == :line_lengths
        L_matrix = get_effective_distances(embedded_graph, mean_len_km = pg_struct.mean_len_km, shortest_line_km = pg_struct.shortest_line_km) # Matrix containing the line lengths in km
        Y, Y_shunt = get_line_admittance_matrix(pg_struct, L_matrix) # Line admittance matrix and Shunts, Entry's are in Ohm

    elseif pg_struct.coupling == :predefined
        Y = pg_struct.edge_parameters[:Y]
        Y_shunt = pg_struct.edge_parameters[:Y_shunt]
        
    elseif pg_struct.coupling == :homogenous
        K = pg_struct.edge_parameters[:K]
        lines = get_lines_homogenous(pg_struct, K)

        return lines
    end

    if pg_struct.lines == :StaticLine
        lines = get_lines_static(pg_struct, Y)
        
    elseif pg_struct.lines == :PiModelLine
        lines = get_lines_Pi(pg_struct, Y, Y_shunt)
    end

    return lines
end

"""
   get_lines_static(pg_struct::PGGeneration, Y, Y_base)

- `pg_struct`: Struct containing all data of the power grid
"""
function get_lines_static(pg_struct::PGGeneration, Y)
    Y_base = pg_struct.Y_base
    e, from_vec, to_vec, lines = line_data(pg_struct)


    for l in 1:length(e)
        Y_pu = Y[from_vec[l], to_vec[l]] / Y_base # Convert admittance [ohm] to p.u. system
        lines[l] = StaticLine(from = from_vec[l], to = to_vec[l], Y = Y_pu) 
    end
    return lines
end

"""
   get_lines_Pi(pg_struct::PGGeneration, Y, Y_shunt)

- `pg_struct`: Struct containing all data of the power grid
- `Y`:  Line admittance matrix, Entry's are in [Ohm] and depend on the line length
- `Y_shunt`: Shunt Admittance of the lines. Used in the PI Model
"""
function get_lines_Pi(pg_struct::PGGeneration, Y, Y_shunt)
    Y_base = pg_struct.Y_base
    e, from_vec, to_vec, lines = line_data(pg_struct)

    for l in 1:length(e)
        Y_pu = Y[from_vec[l], to_vec[l]] / Y_base # Convert admittance [ohm] to p.u. system
        Y_shunt_pu = Y_shunt[from_vec[l], to_vec[l]] / Y_base # Convert shunt admittance [ohm] to p.u. system
        lines[l] = PiModelLine(from = from_vec[l], to = to_vec[l], y = Y_pu, y_shunt_km = Y_shunt_pu, y_shunt_mk = Y_shunt_pu) 
    end
    return lines
end

"""
   get_lines_homogenous(pg_struct::PGGeneration, K)

- `pg_struct`: Struct containing all data of the power grid
- `K`:  Coupling strength of all lines in [p.u.]. 
"""
function get_lines_homogenous(pg_struct::PGGeneration, K)
    embedded_graph = pg_struct.embedded_graph
    e, from_vec, to_vec, lines = line_data(pg_struct)

    for l in 1:length(e)
        lines[l] = StaticLine(from = from_vec[l], to = to_vec[l], Y = K) 
    end
    return lines
end

function line_data(pg_struct::PGGeneration)
    embedded_graph = pg_struct.embedded_graph

    e = edges(embedded_graph.graph)
    from_vec = src.(e)
    to_vec = dst.(e)
    lines = Array{Any}(undef, length(e))

    return e, from_vec, to_vec, lines 
end