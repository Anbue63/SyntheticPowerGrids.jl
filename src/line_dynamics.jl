"""
   get_lines(pg, Y, Y_base)

get the lines dynamics of the power grid.
- `pg`: Graph structure of the power grid
- `Y`:  Line admittance matrix, Entry's are in [Ohm] and depend on the line length
- `Y_shunt`: Shunt Admittance of the lines. Used in the PI Model
- `Y_base`: Base admittance of the power system for [p.u.] conversion
"""
function get_lines(pg, pg_struct, Y, Y_shunt, Y_base)
    if pg_struct.lines == :StaticLine
        lines = get_lines_static(pg, Y, Y_base)
    elseif pg_struct.lines == :PiModelLine
        lines = get_lines_Pi(pg, Y, Y_shunt, Y_base)
    end
    return lines
end

function get_lines_static(pg, Y, Y_base)
    e = edges(pg.graph)
    from_vec = src.(e)
    to_vec = dst.(e)
    
    lines = Array{Any}(undef, length(e))

    for l in 1:length(e)
        Y_pu = Y[from_vec[l], to_vec[l]] / Y_base # Convert admittance [ohm] to p.u. system
        lines[l] = StaticLine(from = from_vec[l], to = to_vec[l], Y = Y_pu) 
    end
    return lines
end

function get_lines_Pi(pg, Y, Y_shunt, Y_base)
    e = edges(pg.graph)
    from_vec = src.(e)
    to_vec = dst.(e)
    
    lines = Array{Any}(undef, length(e))

    for l in 1:length(e)
        Y_pu = Y[from_vec[l], to_vec[l]] / Y_base # Convert admittance [ohm] to p.u. system
        Y_shunt_pu = Y_shunt[from_vec[l], to_vec[l]] / Y_base # Convert shunt admittance [ohm] to p.u. system
        lines[l] = PiModelLine(from = from_vec[l], to = to_vec[l], y = Y_pu, y_shunt_km = Y_shunt_pu, y_shunt_mk = Y_shunt_pu) 
    end
    return lines
end