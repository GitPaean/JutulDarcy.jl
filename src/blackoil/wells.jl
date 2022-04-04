function apply_well_reservoir_sources!(sys::BlackOilSystem, res_q, well_q, state_res, state_well, param_res, param_well, perforations, sgn)
    p_res = state_res.Pressure
    p_well = state_well.Pressure

    val = x -> local_ad(x, nothing)

    μ = state_res.PhaseViscosities
    kr = state_res.RelativePermeabilities
    ρ = state_res.PhaseMassDensities
    rs = state_res.Rs
    b = state_res.ShrinkageFactors

    ρ_w = state_well.PhaseMassDensities
    s_w = state_well.Saturations
    rs_w = state_well.Rs
    b_w = state_well.ShrinkageFactors

    rhoS = tuple(param_res[:reference_densities]...)

    perforation_sources_blackoil!(well_q, perforations, val(p_res),         p_well,  val(kr), val(μ), val(ρ), val(b), val(rs),    ρ_w,      b_w,      s_w,      rs_w,  sgn, rhoS)
    perforation_sources_blackoil!(res_q,  perforations,     p_res,      val(p_well),     kr,      μ,      ρ,      b,      rs, val(ρ_w), val(b_w), val(s_w), val(rs_w), sgn, rhoS)
end

function perforation_sources_blackoil!(target, perf, p_res, p_well, kr, μ, ρ, b, rs, ρ_w, b_w, s_w, rs_w, sgn, rhoS)
    # (self -> local cells, reservoir -> reservoir cells, WI -> connection factor)
    nc = size(ρ, 1)
    nph = size(μ, 1)
    a, l, v = 1, 2, 3

    rhoOS = rhoS[l]
    rhoGS = rhoS[v]
    @inbounds for i in eachindex(perf.self)
        si, ri, wi, gdz = unpack_perf(perf, i)
        if gdz != 0
            ρ_mix = @views mix_by_saturations(s_w[:, si], ρ_w[:, si])
            ρgdz = gdz*ρ_mix
        else
            ρgdz = 0
        end
        @inbounds dp = wi*(p_well[si] - p_res[ri] + ρgdz)
        λ_a = kr[a, ri]/μ[a, ri]
        λ_l = kr[l, ri]/μ[l, ri]
        λ_v = kr[v, ri]/μ[v, ri]
        if dp > 0
            # Injection
            λ_t = λ_a + λ_l + λ_v
            Q = sgn*λ_t*dp

            bO = b_w[l, si]
            bG = b_w[v, si]
            rs_i = rs_w[si]

            sO = s_w[l, si]
            sG = s_w[v, si]

            target[a, i] = s_w[a, si]*ρ_w[a, si]*Q
            target[l, i] = sO*rhoOS*bO*Q
            target[v, i] = rhoGS*(sG*bG + rs_i*sO*bO)*Q
        else
            # Production
            Q = sgn*dp
            target[a, i] = Q*ρ[a, ri]*λ_a

            bO = b[l, ri]
            bG = b[v, ri]

            α_l = bO*λ_l
            α_v = bG*λ_v

            q_o = Q*bO*λ_l*rhoOS
            q_g = Q*(bO*λ_l*rs[ri] + bG*λ_v)*rhoGS

            target[l, i] = q_o
            target[v, i] = q_g
        end
    end
end


function flash_wellstream_at_surface(well_model::SimulationModel{D, S}, well_state, rhoS) where {D, S<:BlackOilSystem}
    vol = well_state.TotalMasses[:, 1]./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end
