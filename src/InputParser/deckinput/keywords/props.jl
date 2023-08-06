
function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:RPTPROPS})
    read_record(f)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:SWOF}, Val{:SGOF}})
    k = unpack_val(v)
    sat_tab = parse_saturation_table(f, outer_data)
    for tab in sat_tab
        swap_unit_system_axes!(tab, units, (:identity, :identity, :identity, :pressure))
    end
    data["$k"] = sat_tab
end


function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVDG})
    pvdg = parse_dead_pvt_table(f, outer_data)
    for tab in pvdg
        swap_unit_system_axes!(tab, units, (:pressure, :gas_formation_volume_factor, :viscosity))
    end
    data["PVDG"] = pvdg
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVTO})
    pvto = parse_live_pvt_table(f, outer_data)
    for tab in pvto
        swap_unit_system_axes!(tab["data"], units, (:pressure, :liquid_formation_volume_factor, :viscosity))
        swap_unit_system!(tab["key"], units, :u_rs)
    end
    data["PVTO"] = pvto
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVTW})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN]
    utypes = (:pressure, :liquid_formation_volume_factor, :compressibility, :viscosity, :compressibility)
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        t = parse_defaulted_line(rec, tdims)
        swap_unit_system_axes!(t, units, utypes)
        @assert all(isfinite, t) "PVTW cannot be defaulted, found defaulted record in region $i"
        push!(out, t)
    end
    data["PVTW"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PVCDO})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN]
    utypes = (:pressure, :liquid_formation_volume_factor, :compressibility, :viscosity, :compressibility)
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        t = parse_defaulted_line(rec, tdims)
        swap_unit_system_axes!(t, units, utypes)
        @assert all(isfinite, t) "PVCDO cannot be defaulted, found defaulted record in region $i"
        push!(out, t)
    end
    data["PVCDO"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ROCK})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN, NaN, NaN, NaN]
    utypes = [:pressure, :compressibility, :compressibility, :compressibility, :id, :id]
    out = []
    nreg = number_of_tables(outer_data, :pvt)
    for i = 1:nreg
        l = parse_defaulted_line(rec, tdims)
        swap_unit_system_axes!(l, units, utypes)
        push!(out, l)
    end
    data["ROCK"] = out
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:DENSITY})
    rec = read_record(f)
    tdims = [NaN, NaN, NaN]
    nreg = number_of_tables(outer_data, :pvt)
    out = []
    for i = 1:nreg
        t = parse_defaulted_line(rec, tdims)
        swap_unit_system!(t, units, :density)
        push!(out, t)
    end
    data["DENSITY"] = out
end
