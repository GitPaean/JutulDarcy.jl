
function equilibriate_state(model, contacts,
        datum_depth = missing,
        datum_pressure = JutulDarcy.DEFAULT_MINIMUM_PRESSURE;
        cells = missing,
        rs = missing,
        rv = missing,
        kwarg...
    )
    model = reservoir_model(model)
    D = model.data_domain
    G = physical_representation(D)
    cc = D[:cell_centroids][3, :]
    if ismissing(cells)
        cells = 1:number_of_cells(G)
        pts = cc
    else
        pts = view(cc, cells)
    end

    if ismissing(datum_depth)
        datum_depth = pmin
    end
    sys = flow_system(model.system)

    init = Dict{Symbol, Any}()
    init = equilibriate_state!(init, pts, model, sys, contacts, datum_depth, datum_pressure;
        cells = cells, rv = rv, rs = rs, kwarg...
    )

    is_blackoil = sys isa StandardBlackOilSystem
    if is_blackoil
        sat = init[:Saturations]
        pressure = init[:Pressure]
        nph, nc = size(sat)

        if has_disgas(sys)
            if ismissing(rs)
                Rs = zeros(nc)
            else
                Rs = rs.(pts)
            end
            init[:Rs] = Rs
        end
        if has_vapoil(sys)
            if ismissing(rv)
                Rv = zeros(nc)
            else
                Rv = rv.(pts)
            end
            init[:Rv] = Rv
        end
    end
    return init
end

function equilibriate_state!(init, depths, model, sys, contacts, depth, datum_pressure;
        cells = 1:length(depths),
        rs = missing,
        rv = missing,
        T_z = missing,
        s_min = missing,
        contacts_pc = missing,
        pvtnum = 1,
        kwarg...
    )
    if ismissing(contacts_pc)
        contacts_pc = zeros(number_of_phases(sys)-1)
    end
    zmin = minimum(depths)
    zmax = maximum(depths)

    nph = number_of_phases(sys)

    @assert length(contacts) == nph-1

    rho_s = JutulDarcy.reference_densities(sys)
    phases = JutulDarcy.get_phases(sys)
    disgas = JutulDarcy.has_disgas(sys)
    vapoil = JutulDarcy.has_vapoil(sys)

    if disgas || vapoil
        if JutulDarcy.has_other_phase(sys)
            _, rhoOS, rhoGS = rho_s
        else
            rhoOS, rhoGS = rho_s
        end
    end
    rho = model.secondary_variables[:PhaseMassDensities]
    if rho isa Pair
        rho = last(rho)
    end
    relperm = model.secondary_variables[:RelativePermeabilities]
    if relperm isa Pair
        relperm = last(relperm)
    end

    pvt = rho.pvt
    reg = Int[pvtnum]
    function density_f(p, z, ph)
        pvt_i = pvt[ph]
        if phases[ph] == LiquidPhase() && disgas
            rs_max = table_by_region(sys.rs_max, pvtnum)
            Rs = min(rs(z), rs_max(p))
            b = JutulDarcy.shrinkage(pvt_i, reg, p, Rs, 1)
            rho = b*(rhoOS + Rs*rhoGS)
        elseif phases[ph] == VaporPhase() && vapoil
            rv_max = table_by_region(sys.rv_max, pvtnum)
            Rv = min(rv(z), rv_max(p))
            b = JutulDarcy.shrinkage(pvt_i, reg, p, Rv, 1)
            rho = b*(rhoGS + Rv*rhoOS)
        else
            rho = rho_s[ph]*JutulDarcy.shrinkage(pvt_i, reg, p, 1)
        end
        return rho
    end
    pressures = determine_hydrostatic_pressures(depths, depth, zmin, zmax, contacts, datum_pressure, density_f, contacts_pc)
    if nph > 1
        s, pc = determine_saturations(depths, contacts, pressures; s_min = s_min, kwarg...)

        nc_total = number_of_cells(model.domain)
        kr = zeros(nph, nc_total)
        s_eval = zeros(nph, nc_total)
        s_eval[:, cells] .= s
        phases = get_phases(sys)
        if length(phases) == 3 && AqueousPhase() in phases
            swcon = zeros(nc_total)
            if !ismissing(s_min)
                swcon[cells] .= s_min[1]
            end
            JutulDarcy.update_kr!(kr, relperm, model, s_eval, swcon, cells)
        else
            JutulDarcy.update_kr!(kr, relperm, model, s_eval, cells)
        end
        kr = kr[:, cells]
        init[:Saturations] = s
        init[:Pressure] = init_reference_pressure(pressures, contacts, kr, pc, 2)
    else
        p = copy(vec(pressures))
        init[:Pressure] = p
        init[:Saturations] = ones(1, length(p))
    end
    if !ismissing(T_z)
        init[:Temperature] = T_z.(depths)
    end

    return init
end


function parse_state0_equil(model, datafile)
    has_water = haskey(datafile["RUNSPEC"], "WATER")
    has_oil = haskey(datafile["RUNSPEC"], "OIL")
    has_gas = haskey(datafile["RUNSPEC"], "GAS")

    is_single_phase = (has_water + has_oil + has_gas) == 1

    sys = model.system
    d = model.data_domain
    has_sat_reg = haskey(d, :satnum)
    ncells = number_of_cells(d)
    if has_sat_reg
        satnum = d[:satnum]
    else
        satnum = ones(Int, ncells)
    end
    if haskey(d, :pvtnum)
        pvtnum = d[:pvtnum]
    else
        pvtnum = ones(Int, nc)
    end
    eqlnum = model.data_domain[:eqlnum]

    if is_single_phase
        has_pc = false
        pc_functions = missing
        krw_fn = missing
    else
        has_pc = haskey(model.secondary_variables, :CapillaryPressure)
        if has_pc
            pcvar = model.secondary_variables[:CapillaryPressure]
            pc_functions = pcvar.pc
            if has_sat_reg
                @assert !isnothing(pcvar.regions)
                @assert pcvar.regions == satnum
            end
        else
            pc_functions = missing
        end

        if has_water
            kr_var = unwrap_reservoir_variable(model.secondary_variables[:RelativePermeabilities])
            krw_fn = kr_var.krw
            if has_sat_reg
                @assert !isnothing(kr_var.regions)
                @assert kr_var.regions == satnum
            end
        else
            krw_fn = missing
        end
    end

    init = Dict{Symbol, Any}()
    sol = datafile["SOLUTION"]
    G = physical_representation(model.data_domain)
    nc = number_of_cells(G)
    nph = number_of_phases(model.system)
    actnum_ix = G.cell_map
    is_blackoil = sys isa StandardBlackOilSystem
    disgas = JutulDarcy.has_disgas(model.system)
    vapoil = JutulDarcy.has_vapoil(model.system)

    equil = sol["EQUIL"]
    nequil = GeoEnergyIO.InputParser.number_of_tables(datafile, :eqlnum)
    npvt = GeoEnergyIO.InputParser.number_of_tables(datafile, :pvtnum)
    nsat = GeoEnergyIO.InputParser.number_of_tables(datafile, :satnum)

    if haskey(sol, "RTEMP")
        Ti = convert_to_si(only(sol["RTEMP"]), :Celsius)
        T_z = z -> Ti
    elseif haskey(sol, "TEMPVD")
        z = vec(sol["TEMPVD"][:, 1])
        Tvd = vec(sol["TEMPVD"][:, 2] + 273.15)
        T_z = get_1d_interpolator(z, Tvd)
    else
        T_z = missing
    end

    @assert length(equil) == nequil
    inits = []
    inits_cells = []
    for ereg in 1:nequil
        eq = equil[ereg]
        cells_eqlnum = findall(isequal(ereg), eqlnum)
        for sreg in 1:nsat
            cells_satnum = findall(isequal(sreg), satnum)
            cells_sat_and_pvt = intersect_sorted(cells_satnum, cells_eqlnum)
            for preg in 1:npvt
                cells_pvtnum = findall(isequal(preg), pvtnum)
                cells = intersect_sorted(cells_pvtnum, cells_sat_and_pvt)
                ncells_reg = length(cells)
                if ncells_reg == 0
                    continue
                end
                actnum_ix_for_reg = actnum_ix[cells]
                datum_depth = eq[1]
                datum_pressure = eq[2]

                woc = eq[3]
                woc_pc = eq[4]
                goc = eq[5]
                goc_pc = eq[6]
                # Contact depths
                s_max = 1.0
                s_min = 0.0

                non_connate = ones(ncells_reg)
                s_max = Vector{Float64}[]
                s_min = Vector{Float64}[]
                if has_pc
                    pc = []
                    for (i, f) in enumerate(pc_functions)
                        f = table_by_region(f, sreg)
                        s = f.X
                        cap = f.F
                        if s[1] < 0
                            s = s[2:end]
                            cap = cap[2:end]
                        end
                        ix = unique(i -> cap[i], 1:length(cap))
                        if i == 1 && get_phases(model.system)[1] isa AqueousPhase
                            @. cap *= -1
                        end
                        s = s[ix]
                        cap = cap[ix]
                        if length(s) == 1
                            push!(s, s[end])
                            push!(cap, cap[end]+1.0)
                        end
                        push!(pc, (s = s, pc = cap))
                    end
                else
                    pc = nothing
                end
                if is_single_phase
                    push!(s_min, zeros(ncells_reg))
                    push!(s_max, ones(ncells_reg))
                else
                    if has_water
                        krw = table_by_region(krw_fn, sreg)
                        if haskey(datafile, "PROPS") && haskey(datafile["PROPS"], "SWL")
                            swl = vec(datafile["PROPS"]["SWL"])
                            swcon = swl[actnum_ix_for_reg]
                        else
                            swcon = fill(krw.connate, ncells_reg)
                        end
                        push!(s_min, swcon)
                        push!(s_max, ones(ncells_reg))
                        @. non_connate -= swcon
                    end
                    if has_oil
                        push!(s_min, zeros(ncells_reg))
                        push!(s_max, non_connate)
                    end
                    if has_gas
                        push!(s_min, zeros(ncells_reg))
                        push!(s_max, non_connate)
                    end
                end

                if nph == 1
                    contacts = []
                    contacts_pc = []
                elseif nph == 2
                    if has_oil && has_gas
                        contacts = (goc, )
                        contacts_pc = (goc_pc, )
                    else
                        contacts = (woc, )
                        contacts_pc = (-woc_pc, )
                    end
                else
                    contacts = (woc, goc)
                    contacts_pc = (-woc_pc, goc_pc)
                end

                rhoGS = map(x -> x[3], datafile["PROPS"]["DENSITY"])
                rhoOS = map(x -> x[1], datafile["PROPS"]["DENSITY"])

                Rs_scale = (rhoGS[preg]/rhoGS[1])*(rhoOS[preg]/rhoOS[1])
                Rv_scale = 1.0/Rs_scale
                if disgas
                    if haskey(sol, "RSVD")
                        rsvd = sol["RSVD"][ereg]
                        z = rsvd[:, 1]
                        Rs = rsvd[:, 2]
                    else
                        @assert haskey(sol, "PBVD")
                        pbvd = sol["PBVD"][ereg]
                        z = pbvd[:, 1]
                        pb = vec(pbvd[:, 2])
                        Rs = sys.rs_max[preg].(pb)
                    end
                    rs = Jutul.LinearInterpolant(z, Rs_scale.*Rs)
                else
                    rs = missing
                end
                if vapoil
                    @assert haskey(sol, "RVVD")
                    rvvd = sol["RVVD"][ereg]
                    z = rvvd[:, 1]
                    Rv = rvvd[:, 2]
                    rv = Jutul.LinearInterpolant(z, Rv_scale.*Rv)
                else
                    rv = missing
                end

                subinit = equilibriate_state(
                        model, contacts, datum_depth, datum_pressure,
                        cells = cells,
                        pvtnum = preg,
                        contacts_pc = contacts_pc,
                        s_min = s_min,
                        s_max = s_max,
                        T_z = T_z,
                        rs = rs,
                        rv = rv,
                        pc = pc
                    )
                push!(inits, subinit)
                push!(inits_cells, cells)
            end
        end
    end
    if length(inits) == 1
        init = only(inits)
    else
        # Handle multiple regions by merging each init
        init = Dict{Symbol, Any}()
        nc = number_of_cells(model.domain)
        touched = [false for i in 1:nc]
        for (k, v) in first(inits)
            if v isa AbstractVector
                init[k] = zeros(nc)
            else
                @assert v isa AbstractMatrix
                init[k] = zeros(size(v, 1), nc)
            end
        end
        for (subinit, cells) in zip(inits, inits_cells)
            for c in cells
                if touched[c]
                    @warn "Equils overlap for cell $c?"
                end
                touched[c] = true
            end

            for (k, v) in subinit
                fill_subinit!(init[k], cells, v)
            end
        end
        @assert all(touched) "Some cells are not initialized by equil: $(findall(!, touched))"
    end
    return init
end

function intersect_sorted(a::Vector{T}, b::Vector{T}) where T
    na = length(a)
    nb = length(b)
    c = Vector{T}()
    ptr_a = ptr_b = 1
    while ptr_a <= na && ptr_b <= nb
        a_val = a[ptr_a]
        b_val = b[ptr_b]
        if a_val == b_val
            push!(c, a_val)
            ptr_a += 1
            ptr_b += 1
        elseif a_val < b_val
            ptr_a += 1
        else
            ptr_b += 1
        end
    end
    return c
end

function fill_subinit!(x::Vector, cells, v::Vector)
    @assert length(v) == length(cells)
    for (i, c) in enumerate(cells)
        x[c] = v[i]
    end
end

function fill_subinit!(x::Matrix, cells, v::Matrix)
    @assert size(x, 1) == size(v, 1)
    @assert size(v, 2) == length(cells)
    for (i, c) in enumerate(cells)
        for j in axes(x, 1)
            x[j, c] = v[j, i]
        end
    end
end

function init_reference_pressure(pressures, contacts, kr, pc, ref_ix = 2)
    nph, nc = size(kr)
    p = zeros(nc)
    ϵ = 1e-12
    for i in eachindex(p)
        p[i] = pressures[ref_ix, i]
        kr_ref = kr[ref_ix, i]
        @assert kr_ref >= -ϵ "Evaluated rel. perm. was $kr_ref for phase reference phase (index $ref_ix)."
        if kr[ref_ix, i] <= ϵ
            for ph in 1:nph
                if kr[ph, i] > ϵ
                    p[i] = pressures[ph, i]
                end
            end
        end
    end
    return p
end

function determine_hydrostatic_pressures(depths, depth, zmin, zmax, contacts, datum_pressure, density_f, contacts_pc)
    nc = length(depths)
    nph = length(contacts) + 1
    ref_ix = min(2, nph)
    I_ref = phase_pressure_depth_table(depth, zmin, zmax, datum_pressure, density_f, ref_ix)
    pressures = zeros(nph, nc)
    pos = 1
    for ph in 1:nph
        if ph == ref_ix
            I = I_ref
        else
            contact = contacts[pos]
            datum_pressure_ph = I_ref(contact) + contacts_pc[pos]
            I = phase_pressure_depth_table(contact, zmin, zmax, datum_pressure_ph, density_f, ph)
            pos += 1
        end
        for (i, z) in enumerate(depths)
            pressures[ph, i] = I(z)
        end
    end
    return pressures
end



function phase_pressure_depth_table(depth, zmin, zmax, datum_pressure, density_f, phase)
    if zmin > depth
        z_up = Float64[]
        p_up = Float64[]
    else
        z_up, p_up = integrate_phase_density(depth, zmin, datum_pressure, density_f, phase)
        # Remove top point
        z_up = reverse!(z_up[2:end])
        p_up = reverse!(p_up[2:end])
    end
    if zmax < depth
        z_down = Float64[]
        p_down = Float64[]
    else
        z_down, p_down = integrate_phase_density(depth, zmax, datum_pressure, density_f, phase)
    end
    z = vcat(z_up, z_down)
    p = vcat(p_up, p_down)
    return Jutul.LinearInterpolant(z, p)
end

function integrate_phase_density(z_datum, z_end, p0, density_f, phase; n = 1000, g = Jutul.gravity_constant)
    dz = (z_end - z_datum)/n
    pressure = zeros(n+1)
    z = zeros(n+1)
    pressure[1] = p0
    z[1] = z_datum

    for i in 2:(n+1)
        p = pressure[i-1]
        depth = z[i-1] + dz
        pressure[i] = p + dz*density_f(p, depth, phase)*g
        z[i] = depth
    end
    return (z, pressure)
end

function determine_saturations(depths, contacts, pressures; s_min = missing, s_max = missing, pc = nothing)
    nc = length(depths)
    nph = length(contacts) + 1
    if ismissing(s_min)
        s_min = [zeros(nc) for i in 1:nph]
    end
    if ismissing(s_max)
        s_max = [ones(nc) for i in 1:nph]
    end
    sat = zeros(nph, nc)
    sat_pc = similar(sat)
    if isnothing(pc)
        for i in eachindex(depths)
            z = depths[i]
            ph = current_phase_index(z, contacts)
            for j in axes(sat, 1)
                is_main = ph == j
                s = is_main*s_max[j][i] + !is_main*s_min[j][i]
                sat[j, i] = s
            end
        end
    else
        ref_ix = 2
        offset = 1
        for ph in 1:nph
            if ph != ref_ix
                s, pc_pair = pc[offset]
                pc_max = maximum(pc_pair)
                pc_min = minimum(pc_pair)
                I = get_1d_interpolator(pc_pair, s, constant_dx = false)
                I_pc = get_1d_interpolator(s, pc_pair, constant_dx = false)
                for i in eachindex(depths)
                    z = depths[i]

                    dp = pressures[ph, i] - pressures[ref_ix, i]
                    if dp > pc_max
                        s_eff = s_max[ph][i]
                    elseif dp < pc_min
                        s_eff = s_min[ph][i]
                    else
                        s_eff = I(dp)
                    end
                    s_eff = clamp(s_eff, s_min[ph][i], s_max[ph][i])
                    sat[ph, i] = s_eff
                    sat_pc[ph, i] = I_pc(s_eff)
                end
                offset += 1
            end
        end
        bad_cells = Int[]
        for i in eachindex(depths)
            sat_i = view(sat, :, i)
            sat_tot = sum(sat_i)
            s_fill = 1 - sat_tot
            if s_fill < 0
                push!(bad_cells, i)
                sat[ref_ix, i] = 0.0
                for j in axes(sat, 1)
                    sat[j, i] /= sat_tot
                end
            else
                sat[ref_ix, i] = s_fill
            end
        end
        if length(bad_cells) > 0
            @warn "Negative saturation in $(length(bad_cells)) cells for phase $ref_ix. Normalizing."
        end
    end
    return (sat, sat_pc)
end

function current_phase_index(z, depths; reverse = true)
    n = length(depths)+1
    out = -1
    if reverse
        i = current_phase_index(z, Base.reverse(depths), reverse = false)
        out = n - i + 1
    else
        if z < depths[1]
            out = 1
        elseif z > depths[end]
            out = n
        else
            for (i, d) in enumerate(depths)
                if d >= z
                    out = i
                    break
                end
            end
        end
    end
    @assert out > 0
    return out
end
