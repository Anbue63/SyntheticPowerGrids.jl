import PowerDynamics: make_generator!, _make_generator_header
import PowerDynamics: make_bus_ac!, _make_bus_ac_header

function make_generator!(dict::Dict{String,Any}, key_b::Int, node::DroopControlledInverterApprox)
    _make_generator_header(dict, key_b)
    key = length(dict["gen"])
    ((dict["gen"])[string(key)])["pg"] = node.P
    ((dict["gen"])[string(key)])["pmin"] = 0.9 * node.P
    ((dict["gen"])[string(key)])["pmax"] = 1.1 * node.P
    ((dict["gen"])[string(key)])["vg"] = node.V_r
end

function make_bus_ac!(data::Dict{String,Any}, node::DroopControlledInverterApprox)
    bus_dict = _make_bus_ac_header(data)
    bus_dict["bus_type"] = 2
    bus_dict["vm"] = abs(node.V_r)  # assumed p.u.
    bus_dict["vmin"] = 0.9 * bus_dict["vm"]
    bus_dict["vmax"] = 1.1 * bus_dict["vm"]
    bus_dict["va"] = angle(node.V_r)
    make_generator!(data, bus_dict["index"], node)
end