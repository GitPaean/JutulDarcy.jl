# # Validation of equation-of-state compositiona simulator
# This example solves a 1D two-phase, three component miscible displacement
# problem and compares against existing simulators (E300, AD-GPRS) to verify
# correctness.
#
# The case is loaded from an input file that can be run in other simulators. For
# convenience, we provide solutions from the other simulators as a binary file
# to perform a comparison without having to run and convert results from other
# the simulators.
#
# This case is a small compositional problem inspired by the examples in Voskov
# et al (JPSE, 2012). A 1D reservoir of 1,000 meters length is discretized into
# 1,000 cells. The model initially contains a mixture made up of 0.6 parts C10,
# 0.1 parts CO2, and 0.3 parts C1 by moles at 150 degrees C and 75 bar pressure.
# Wells are placed in the leftmost and rightmost cells of the domain, with the
# leftmost well injecting pure CO$_2$ at a fixed bottom-hole pressure of 100 bar
# and the other well producing at 50 bar. The model is isothermal and contains a
# phase transition from the initial two-phase mixture to single-phase gas as
# injected CO$_2$ eventually displaces the resident fluids. For further details
# on this setup, see Møyner and Tchelepi (SPE J. 2018)
# [moyner_tchelepi_2018](@cite).
using JutulDarcy
using Jutul
using GLMakie
# using Infiltrator
dpth = @__DIR__
data_path = joinpath(dpth, "temp_data/SIMPLE_COMP_SMALLZ.DATA");
case = setup_case_from_data_file(data_path)
result = simulate_reservoir(case, info_level = 1)
ws, states = result;
# ## Plot solutions and compare
# The 1D displacement can be plotted as a line plot. We pick a step midway
# through the simulation and plot compositions, saturations and pressure.
cmap = :tableau_hue_circle
# TODO: to add the reference data here
# ref_path = joinpath(dpth, "reference.jld2")
# ref = Jutul.JLD2.load(ref_path)

ref = Dict(
    :ZMF => [
        [0.99000001E+00, 0.98999697E+00, 0.98996443E+00, 0.98977447E+00,
         0.98896033E+00, 0.98596948E+00, 0.97544080E+00, 0.93774980E+00,
         0.87621981E+00, 0.83818489E+00, 0.81011033E+00, 0.78687298E+00,
         0.76531929E+00, 0.74298692E+00, 0.71850479E+00, 0.69188869E+00,
         0.66415775E+00, 0.63668633E+00, 0.61071956E+00, 0.58713597E+00,
         0.56638771E+00, 0.54854435E+00, 0.53338718E+00, 0.52051598E+00,
         0.50944388E+00, 0.49966964E+00, 0.49072978E+00, 0.48224518E+00,
         0.47400683E+00, 0.50000000E+00],

        [0.89999996E-02, 0.90016052E-02, 0.90176007E-02, 0.91000479E-02,
         0.93920352E-02, 0.10183870E-01, 0.11883356E-01, 0.14489659E-01,
         0.18011266E-01, 0.24292562E-01, 0.33787258E-01, 0.47048371E-01,
         0.64401381E-01, 0.85644074E-01, 0.10994513E+00, 0.13601145E+00,
         0.16237332E+00, 0.18763749E+00, 0.21066682E+00, 0.23068388E+00,
         0.24729614E+00, 0.26045233E+00, 0.27035528E+00, 0.27736199E+00,
         0.28189558E+00, 0.28438190E+00, 0.28521711E+00, 0.28477025E+00,
         0.28344235E+00, 0.30000001E+00],

        [0.10000000E-02, 0.10014048E-02, 0.10179788E-02, 0.11254736E-02,
         0.16476334E-02, 0.38466749E-02, 0.12675846E-01, 0.47760539E-01,
         0.10576894E+00, 0.13752255E+00, 0.15610243E+00, 0.16607863E+00,
         0.17027932E+00, 0.17136902E+00, 0.17155010E+00, 0.17209983E+00,
         0.17346892E+00, 0.17567620E+00, 0.17861365E+00, 0.18218017E+00,
         0.18631616E+00, 0.19100332E+00, 0.19625755E+00, 0.20212200E+00,
         0.20866056E+00, 0.21594846E+00, 0.22405311E+00, 0.23298459E+00,
         0.24255082E+00, 0.19999997E+00]
    ],
    :SGAS => [
        0.10000000E+01, 0.10000000E+01, 0.10000000E+01, 0.10000000E+01,
        0.10000000E+01, 0.10000000E+01, 0.10000000E+01, 0.98710263E+00,
        0.83166999E+00, 0.75620270E+00, 0.72579676E+00, 0.72351509E+00,
        0.73700809E+00, 0.75643218E+00, 0.77562833E+00, 0.79210407E+00,
        0.80557418E+00, 0.81655985E+00, 0.82568872E+00, 0.83348954E+00,
        0.84037882E+00, 0.84668940E+00, 0.85269928E+00, 0.85865647E+00,
        0.86480319E+00, 0.87140530E+00, 0.87879455E+00, 0.88743806E+00,
        0.89806920E+00, 0.93658108E+00
    ]
)

step_to_plot = 100
fig = with_theme(theme_latexfonts()) do
    x = reservoir_domain(case)[:cell_centroids][1, :]
    mz = 20
    ix = step_to_plot
    mt = :circle
    fig = Figure(size = (800, 400))
    ax = Axis(fig[2, 1], xlabel = "Cell center / m")
    cnames = ["CO2", "METHANE", "DECANE"]
    cnames = ["CO₂", "C₁", "C₁₀"]
    cnames = [L"\text{CO}_2", L"\text{C}_1", L"\text{C}_{10}"]
    lineh = []
    lnames = []
    crange = (1, 4)
    for i in range(crange...)
        if i == 4
            cname = L"\text{S}_g"
            ecl = ref[:SGAS][:]
            #@infiltrate
            ju = states[ix][:Saturations][2, :]
        else
            @assert i <= 4
            #@infiltrate
            ecl = ref[:ZMF][i]
            ju = states[ix][:OverallMoleFractions][i, :]
            cname = cnames[i]
        end
        h = lines!(ax, x, ju, colormap = cmap, color = i, colorrange = crange, label = cname)
        push!(lnames, cname)
        push!(lineh, h)
        scatter!(ax, x, ecl, markersize = mz, colormap = cmap, color = i, colorrange = crange)
    end

    l_ju = LineElement(color = :black, linestyle = nothing)
    l_ecl = MarkerElement(color = :black, markersize = mz, marker = mt)

    Legend(
        fig[1, 1],
        [[l_ju, l_ecl], lineh],
        [[L"\text{JutulDarcy}", L"\text{E300}"], lnames],
        ["Simulator", "Result"],
        tellwidth = false,
        orientation = :horizontal,
    )
    fig
end