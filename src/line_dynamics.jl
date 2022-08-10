"""
   get_lines(pg, Y, Y_base)

get the lines dynamics of the power grid.
- `pg`: Graph structure of the power grid
- `Y`:  Line admittance matrix, Entry's are in [Ohm] and depend on the line length
- `Y_base`: Base admittance of the power system for [p.u.] conversion
"""
function get_lines(pg, Y, Y_base)
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