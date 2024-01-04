function parse_and_set_grid_data!(data, outer_data, units, cfg, f, k; unit = :id, T = Float64)
    bdims = get_boxdims(outer_data)
    cdims = get_cartdims(outer_data)
    vals = parse_grid_vector(f, bdims, T)
    if unit != :id
        vals = swap_unit_system!(vals, units, Val(unit))
    end
    if bdims == cdims
        data["$k"] = vals
    else
        error("Not implemented")
    end
end

function finish_current_section!(data, units, cfg, outer_data, ::Val{:GRID})
    if !haskey(data, "MINPV")
        io = IOBuffer("1e-6\n/\n")
        parse_keyword!(data, outer_data, units, cfg, io, Val(:MINPV))
    end
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:GRIDFILE})
    rec = read_record(f)
    tdims = [0, 1];
    data["GRIDFILE"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Union{Val{:MINPORV}, Val{:MINPV}})
    rec = read_record(f)
    tdims = [1e-6];
    rec = parse_defaulted_line(rec, tdims)
    zcorn = swap_unit_system!(rec, units, :volume)
    data["MINPV"] = only(rec)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:INIT})
    data["INIT"] = true
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COORDSYS})
    read_record(f)
    parser_message(cfg, outer_data, "COORDSYS", PARSER_MISSING_SUPPORT)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MAPUNITS})
    # TODO: This needs to be handled
    partial_parse!(data, outer_data, units, cfg, f, :MAPUNITS)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:GRIDUNIT})
    # TODO: This needs to be handled
    partial_parse!(data, outer_data, units, cfg, f, :GRIDUNIT)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FILEUNIT})
    # TODO: This needs to be handled
    partial_parse!(data, outer_data, units, cfg, f, :FILEUNIT)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:GDORIENT})
    # TODO: This needs to be handled
    partial_parse!(data, outer_data, units, cfg, f, :GDORIENT)
end

function partial_parse!(data, outer_data, units, cfg, f, k::Symbol)
    rec = read_record(f)
    parser_message(cfg, outer_data, "$k", PARSER_MISSING_SUPPORT)
    data["$k"] = rec
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MAPAXES})
    rec = parse_deck_vector(f, Float64)
    data["MAPAXES"] = rec
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:COORD})
    coord = parse_deck_vector(f, Float64)
    coord = swap_unit_system!(coord, units, Val(:length))
    data["COORD"] = coord
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:ZCORN})
    zcorn = parse_deck_vector(f, Float64)
    zcorn = swap_unit_system!(zcorn, units, Val(:length))
    data["ZCORN"] = zcorn
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:PERMX}, Val{:PERMY}, Val{:PERMZ}})
    k = unpack_val(v)
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, k, unit = :permeability)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Val{:MULTPV})
    k = unpack_val(v)
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, k)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:PRATIO}, Val{:BIOTCOEF}})
    k = unpack_val(v)
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, k)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:YMODULE}})
    k = unpack_val(v)
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, k, unit = :gigapascal)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:POELCOEF}, Val{:THELCOEF}, Val{:THERMEXR}, Val{:THCONR}})
    k = unpack_val(v)
    vals = parse_grid_vector(f, get_cartdims(outer_data), Float64)
    parser_message(cfg, outer_data, "$k", PARSER_PARTIAL_SUPPORT)
    data["$k"] = vals
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:FIPNUM}, Val{:PVTNUM}, Val{:SATNUM}, Val{:EQLNUM}, Val{:ROCKNUM}})
    k = unpack_val(v)
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, k, T = Int)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PORO})
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, :PORO)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:NTG})
    parse_and_set_grid_data!(data, outer_data, units, cfg, f, :NTG)
end

function parse_keyword!(data, outer_data, units, cfg, f, v::Union{Val{:DX}, Val{:DY}, Val{:DZ}})
    k = unpack_val(v)
    data["$k"] = parse_grid_vector(f, get_cartdims(outer_data), Float64)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:TOPS})
    tops = parse_deck_vector(f, Float64)
    data["TOPS"] = swap_unit_system!(tops, units, Val(:length))
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:DIMENS})
    rec = read_record(f)
    to_int = x -> Parsers.parse(Int, x)
    d = to_int.(filter!(x -> length(x) > 0, split(only(rec), DECK_SPLIT_REGEX)))
    data["DIMENS"] = d
    set_cartdims!(outer_data, d)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:SPECGRID})
    rec = read_record(f)
    tdims = [1, 1, 1, 1, "F"]
    data["SPECGRID"] = parse_defaulted_line(rec, tdims)
    set_cartdims!(outer_data, data["SPECGRID"][1:3])
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:PINCH})
    rec = read_record(f)
    tdims = [0.001, "GAP", Inf, "TOPBOT", "TOP"]
    parser_message(cfg, outer_data, "PINCH", PARSER_JUTULDARCY_PARTIAL_SUPPORT)
    data["PINCH"] = parse_defaulted_line(rec, tdims)
end

function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:FAULTS})
    read_record
    tdims = ["NAME", -1, -1, -1, -1, -1, -1, "XYZ_IJK"]
    if !haskey(data, "FAULTS")
        data["FAULTS"] = Dict{String, Any}()
    end
    faults = data["FAULTS"]
    while true
        rec = read_record(f)
        if length(rec) == 0
            break
        end
        parsed = parse_defaulted_line(rec, tdims, required_num = length(tdims), keyword = "FAULTS")
        name = parsed[1]
        flt = (
            i = parsed[2]:parsed[3],
            j = parsed[4]:parsed[5],
            k = parsed[6]:parsed[7],
            direction = parsed[8]
        )
        if haskey(faults, name)
            push!(faults[name], flt)
        else
            faults[name] = [flt]
        end
    end
end


function parse_keyword!(data, outer_data, units, cfg, f, ::Val{:MULTFLT})
    read_record
    tdims = ["NAME", 1.0, 1.0]
    faults = outer_data["GRID"]["FAULTS"]
    out = []
    while true
        rec = read_record(f)
        if length(rec) == 0
            break
        end
        parsed = parse_defaulted_line(rec, tdims, required_num = length(tdims), keyword = "FAULTS")
        name = parsed[1]
        @assert haskey(faults, name) "Fault $name used in MULTFLT, but it is not in the list of declared faults: $(keys(faults))"
        push!(out, parsed)
    end
    data["MULTFLT"] = out
end

