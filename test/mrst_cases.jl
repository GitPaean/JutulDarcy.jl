using JutulDarcy, Test

function test_mrst_case(casename; tol = 0.01, wtol = tol, otol = tol, gtol = tol, bhptol = tol)
    mrst_data_folder = joinpath(pathof(JutulDarcy), "..", "..", "test", "mrst")
    states, reports, output_path, setup = simulate_mrst_case(joinpath(mrst_data_folder, "$(casename).mat"), info_level = -1, verbose = false)
    model = setup.sim.model
    forces = setup.case.forces
    dt = setup.case.dt
    ws = full_well_outputs(model, states, forces, shortname = true)
    ref = setup.mrst["extra"][1]["mrst_solution"]
    @assert isdir(mrst_data_folder)
    for k in keys(ws)
        @testset "$k" begin
            ix = findfirst(isequal("$k"), vec(ref["names"]))
            w = (abs.(ws[k][:mass]) .> 1e-10).*dt
            for wfield in [:bhp, :wrat, :orat, :grat]
                if haskey(ws[k], wfield)
                    is_bhp = false
                    if wfield == :bhp
                        mrst = ref["bhp"]
                        is_bhp = true
                        ϵ = bhptol
                    elseif wfield == :wrat
                        mrst = ref["qWs"]
                        ϵ = wtol
                    elseif wfield == :orat
                        mrst = ref["qOs"]
                        ϵ = otol
                    else
                        @assert wfield == :grat
                        mrst = ref["qGs"]
                        ϵ = gtol
                    end
                    jutul = ws[k][wfield]
                    m = mrst[:, ix]
                    ref_norm = norm(m.*w, 1)
                    # Don't test pure noise. MRST and Jutul do not use exactly the same wells.
                    if ref_norm > 1e-8*sum(dt)
                        @testset "$wfield" begin
                            err = norm(jutul.*w - m.*w, 1)/ref_norm
                            if err > ϵ
                                @info "$k: $wfield" err ref_norm #jutul m ws[k][:mass] w
                                println("Compare")
                                display(hcat(jutul, m, w))
                            end
                            @test err < ϵ
                        end
                    end
                end
            end
        end
    end
    return (ws, ref)
end

@testset "SPE1" begin
    ws, ref = test_mrst_case("spe1")
end
@testset "SPE3" begin
    ws, ref = test_mrst_case("spe3", wtol = 0.2, gtol = 0.02)
end
@testset "SPE9" begin
    ws, ref = test_mrst_case("spe9", wtol = Inf, gtol = 0.1, otol = 0.05)
end
@testset "Egg" begin
    ws, ref = test_mrst_case("egg")
end