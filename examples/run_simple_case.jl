using JutulDarcy
using Jutul
using GLMakie
using Logging


# using PrettyPrint:pretty_print
# import PrettyPrint

# # global_logger(Logging.ConsoleLogger(stderr, Logging.Debug));

# run the simulation
data_path = joinpath(@__DIR__, "data/SIMPLE_COMP_SMALL.DATA")
# dpth = JutulDarcy.GeoEnergyIO.test_input_file_path("SIMPLE_COMP")
# data_path = joinpath(dpth, "SIMPLE_COMP.DATA")
case = setup_case_from_data_file(data_path)
# println("dump usage")
# dump(case)
# println("dump usage end")
# println("show")
# @show case
# println("show usage end")
# println("pretty print usage end")
# PrettyPrint.pprint(case)
# println("pretty print usage end")

result = simulate_reservoir(case, info_level = 1)
ws, states = result;

## plotting the results
#plot_reservoir(case, states, step = length(case.dt), key = :Saturations);

## plot the well results
plot_well_results(ws);