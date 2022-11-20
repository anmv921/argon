module Argon

# todo combine the autocorrelation functions with the  force calculation for efficiency
# todo cell lists
# todo use array operations (ex .*) and dot products
# todo xyz file 
# https://www.ovito.org/docs/current/reference/file_formats/input/xyz.html#file-formats-input-xyz-extended-format
# todo 3d
# todo leapfrog algorithm...
# todo write velocities to a file and stresses for the correlation with fft for all t-t0 values
# for each particle - more samples...

# switch to cbd & srk algorithm,
# add division
# voronoi tesselation
# self-intermediate scattering

# test execution time vs t
# do the simulations for each rho
# switch 

# todo initial velocities gaussian

# organize this into files and folders, 
# write each particle velocity for correlation
# write each particle shear stress
# average over every timestep and realizations
# diffusion...
# other integration schemes
# c 3d following thijjsen & rapaport
# active particles, brownian/langevin dynamics
# ants
# gui
# precompile
# see jupyter nb for calculations



    using DelimitedFiles
    import Plots
    using GLM, DataFrames
    using Query
    using Dates

    include("./Workflow.jl")
    using .Workflow

    spacing0 = 1 # initial spacing between particles. determines the density
    tmax = 1.0
    npart = 10^2
    sampleStep = 500
    testing = false
    datetoday = Dates.today()
    rightnow = Dates.Time(Dates.now())
    global folderName = string(datetoday) * "-" * replace(string(rightnow), ":" => ".")

    global tempo = @elapsed Workflow.md(testing, spacing0, tmax, npart, sampleStep, folderName)

    println("The execution took " * string(tempo) * " seconds.")
    timeFile = open(joinpath(folderName, "exec_time_sec.txt"), "w");
    write(timeFile, string(tempo)); # in seconds
    close(timeFile);

end

