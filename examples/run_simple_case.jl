using JutulDarcy
using Jutul
using GLMakie
## run the simulation
# data_path = joinpath(@__DIR__, "data/SIMPLE_COMP.DATA")
dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("SIMPLE_COMP")
data_path = joinpath(dpth, "SIMPLE_COMP.DATA")
case = setup_case_from_data_file(data_path)
result = simulate_reservoir(case, info_level = 1)
ws, states = result;

## plotting the results
plot_reservoir(case, states, step = length(case.dt), key = :Saturations)
