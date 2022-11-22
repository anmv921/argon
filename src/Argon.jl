module Argon

# todo combine the autocorrelation functions with the  force calculation for efficiency
# todo cell lists
# todo use array operations (ex .*) and dot products
# https://www.ovito.org/docs/current/reference/file_formats/input/xyz.html#file-formats-input-xyz-extended-format
# todo 3d
# todo predictor-corrector
# todo write velocities to a file and stresses for the correlation
# for each particle - more samples...

# switch to cbd & srk algorithm,
# add division
# voronoi tesselation
# self-intermediate scattering

# test execution time vs t
# do the simulations for each rho
# switch 

# todo initial velocities gaussian
# initial positions read from a file - organize files by N and rho

# organize this into files and folders, 
# write each particle velocity for correlation
# write each particle shear stress
# average over every timestep and realizations
# diffusion...
# other integration schemes
# active particles, brownian/langevin dynamics
# precompile
# lookup table
# plots



import JSON
import Plots

using DataFrames
using Dates
using DelimitedFiles
using GLM
using Query

include("./Workflow.jl")
using .Workflow

open("params.json", "r") do f
    global params
    dicttxt = read(f, String)
    params = JSON.parse(dicttxt)
end

# Initial spacing between particles. determines the density
spacing0 = parse(Float64, params["inicial_spacing"])
tmax = parse(Float64, params["tmax"])
npart = parse(Int64, params["npart"]) # Must have integer sqrt
sampleStep = parse(Int64, params["sampleStep"])

datetoday = Dates.today()
rightnow = Dates.Time(Dates.now())
global folderName = "data\\" * string(datetoday) *
 "_" * replace(string(rightnow), ":" => ".")

testing = false
global tempo = @elapsed Workflow.md(testing,
 spacing0, tmax, npart, sampleStep, folderName, params)

println("The execution took " * string(tempo) * " seconds.")
timeFile = open(joinpath(folderName, "exec_time_sec.txt"), "w");
write(timeFile, string(tempo)); # in seconds
close(timeFile);

end

# folderName = "2022-05-08-22.19.34.042"
# data = readdlm(joinpath(folderName, "diff.dat"), '\t', Float64, '\n', header=false,
#     comments=true, comment_char='#');

# tempos = data[:,1];
# corrsv = data[:, 2];
# msds = data[:,3];

# Plots.plot(tempos, corrsv)

# df = DataFrame(X=tempos, Y=msds);

# q1 = @from i in df begin
#     @where i.X>2
#     @select {X=i.X, i.Y}
#     @collect DataFrame
# end;

# ols = lm(@formula(Y ~ X), q1)

# intercep, declive = coef(ols)
# err_int, err_dec = stderror(ols)

# Plots.plot(tempos, msds)
# Plots.plot!(q1.X, q1.X .* declive .+ intercep)