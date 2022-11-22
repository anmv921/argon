module Workflow

    function dif(switch#=, nsamp=#) # switch=0 init, 1 sample, 2 results; 
        # Diffusion
        
        if switch==0 # initialization
            global vacf = 0.0
            global r2t = 0.0
            global vt0 = zeros(Float64, npart, 2) # velocidade @t0
            global xt0 = zeros(Float64, npart, 2)
            global sigmaxy = zeros(Float64, npart, 1)
            global sigmaxy0 = 0.0
            global sacf = 0.0
            
            
            for i=1:npart
                global vt0[i, 1], vt0[i, 2] = v[i, 1],  v[i, 2]
                global xt0[i, 1], xt0[i, 2] = x[i, 1], x[i, 2]
                for j=i+1:npart # loop over all pairs
                    global sigmaxy0 += v[i,1]*v[i,2] + (x[i,1]-x[j,1])*f[i,2]
                end
            end
            
        
        elseif (switch==1) # sample
            
            global vacf = 0.0
            global r2t = 0.0
            global sigmaxy = 0.0
            global sacf = 0.0
            
            for i=1:npart
                global vacf += v[i, 1]*vt0[i, 1] + v[i, 2]*vt0[i, 2] # update velocity autocorr.
                global r2t += (x[i, 1] - xt0[i, 1] )^2 + (x[i, 2] - xt0[i, 2] )^2
                # update mean-squared displ.
                for j=i+1:npart # loop over all pairs
                    global sigmaxy += v[i,1]*v[i,2] + (x[i,1]-x[j,1])*f[i,2]
                end
            end
            
            global vacf /= npart
            global r2t /= npart
            
            global sacf = sigmaxy0 * sigmaxy / npart
            
            write(diffFile,
                string(t - step0diff*delt) * "\t" * string(vacf) *
                "\t" * string(r2t) * "\t" * string(sacf) *  "\n");
        end
        
    end

    function gr(switch) # TODO combine this with force calculation for efficiency
        # radial distribution function
        # switch 0 initialization
        # 1 sample
        # 2 results
        
        if (switch == 0) # initialization
            global nhis = 200 # total number of bins
            global ngr = 0
            global delg = box / (2 * nhis) # bin size
            global g = zeros(Float64, nhis, 1)     
        
        elseif (switch == 1) # sample
            global ngr += 1
            for i=1:npart-1
                for j=i+1:npart # loop over all pairs
                    xr = x[i, :] - x[j, :]
                    xr = xr - box * [round(XR/box) for XR in xr] # periodic bounds
                    r = sqrt(sum(xr.*xr))
                    
                    #raioCorte = box / 2
                    raioCorte = 3.2
                    if (r < box / 2) # half of box length
                        ig = floor(Int, r/delg)
                        global g[ig] = g[ig] + 2 # contribution for particle i and j
                    end
                end
            end

        elseif (switch == 2) # determine g(r)
            
            write(grfile, "# r" * "\t" * "g" * "\n")
            
            for i=1:nhis
                r = delg * (i + 0.5) # distance r
                #vb = ( (i + 1)^3 - i^3) * delg^3 # volume between bin i+1 and i
                vb =  delg * r # r dr
                #nid = (4/3) * pi * vb * rho # ideal gas 
                # http://www.physics.emory.edu/faculty/weeks//idl/gofr2.html
                nid = 2 * pi * vb * rho # ideal gas part in vb
                global g[i] = g[i] / (ngr * npart * nid) # normalize gr
                
                write(grfile, string(r) * "\t" * string(g[i]) * "\n")
            end
            
            
        end

    end

    function term()
        gr(2)
        
        close(io);
        close(ovito);
        close(grfile);
        close(diffFile);
        println("Program execution done.")
    end

    function lattice_pos(index)
        # Index begins at 1
        ix = mod(index - 1, width)
        iy = (index - 1) ÷ width # Integer division
        xPos = poffset + ix*dx
        yPos = poffset + iy*dy
        return xPos, yPos
    end

    function init(in_Testing, spacing0, in_Npart,
         in_SampleStep, in_FolderName, tmax, params)
        # Initialization of md program
        global t = 0.0;

        global testing = in_Testing

        #rightnow = Dates.Time(Dates.now())
        #datetoday = Dates.today()
        #folderName = string(datetoday) * "-" * replace(string(rightnow), ":" => ".")

        folderName = in_FolderName
        mkdir(folderName)

        enFilePath = joinpath(folderName, "e.dat")
        global io = open(enFilePath, "w");
        ovFilePath = joinpath(folderName, "ovito.xyz")
        global ovito = open(ovFilePath, "w")
        global diffFile = open(joinpath(folderName, "diff.dat"), "w")

        global dt = parse(Float64, params["dt"])

        global grfile = open(joinpath(folderName, "gr.dat"), "w")

        global npart = in_Npart

        global rad = parse(Float64, params["mol_radius"])

        # Initial conditions parameters

        # Distance between the particle centers; determines the density
        global dx = 2*rad + spacing0
        global dy = dx
        global width = sqrt(npart) # Number of particles (laterally)
        global height = width
        global poffset = rad + spacing0/2
        
        ########################################################################
        global box = sqrt(npart) * 2 * rad + spacing0 * (sqrt(npart) - 1) +
         spacing0
        # Simulation cell size!
        ########################################################################
        
        global rho = npart / (box * box)
        global temp = parse(Float64, params["temperatura"])
        
        global rc = parse(Float64, params["cutoff_radius"])
        global rc2 = rc*rc
        global ecut = 4.0 * ( 1.0 / (rc^12) - 1.0 / (rc^6) )

        global Ndim = parse(Int64, params["Ndim"])

        global x = zeros(Float64, npart, Ndim)
        global v = zeros(Float64, npart, Ndim)
        global xm = zeros(Float64, npart, Ndim)
        global f = zeros(Float64, npart, Ndim)
        global press = zeros(Float64, npart, 1) # Pressure

        global sumv = zeros(Float64, Ndim)
        global sumv2 = 0

        
        
        global en = 0.0

        global delt = dt
        global sampleStep = in_SampleStep
        
        global passo = 0
        global passoMax = tmax / delt
        global particulaTipo = "8"
        
        global etot = 0.0
        global ecin = 0.0
        
        global virial = 0.0
        
        
        #global step0diff = floor(Int, 0.1 * tmax / delt) # initial step for vacf & diff @ tmax*...
        global step0diff = 0.0
        
        for i=1:npart

            global x[i, 1], x[i, 2] = lattice_pos(i) # place the particles on a lattice
            # section 3.2.2
            # lies... there is nothing in section 3.2.2. NOTHING!
            global v[i, 1], v[i, 2] = ( rand() - 0.5 ), ( rand() - 0.5 )
            global sumv[1], sumv[2] = sumv[1] + v[i, 1], sumv[2] + v[i, 2]
            global sumv2 = sumv2 + ( v[i, 1]*v[i, 1] + v[i, 2]*v[i, 2] )
        end

        global sumv[1], sumv[2] = sumv[1] / npart, sumv[2] / npart # velocity center of mass

        global sumv2 = sumv2 / npart # mean squared velocity
        global fs = sqrt(Ndim * temp / sumv2) # scale factor of the velocities

        for i=1:npart # set the kinetic energy and
            global v[i, 1],  v[i, 2] = (v[i, 1] - sumv[1]) * fs, (v[i, 2] - sumv[1]) * fs 
            # set veclocity center of mass to zero
            global xm[i, 1], xm[i, 2]  = x[i, 1] - v[i, 1]*dt, x[i, 2] - v[i, 2]*dt # position previous time step
        end

        write(ovito, string(npart) * "\n");
        write(ovito, "Lattice=\"" * string(box) * " 0.0 0.0 0.0 "* string(box) * " 0.0 0.0 0.0 0.0\" ")
        write(ovito, "Properties=pos:R:2:Radius:R:1:Type:S:1:Velocity:R:1:Pressure:R:1")
        write(ovito, " Origin=\"0.0 0.0 0.0\"")
        write(ovito, " pbc=\"T T F\"")
        write(ovito, " t=" * string(t) * "\n")
        
        for i=1:npart
            vel = sqrt(v[i, 1] * v[i, 1] + v[i, 2] * v[i, 2])
            ostr = string(x[i, 1]) * "\t" * string(x[i, 2]) * "\t" * string(rad) * "\t" * string(i) * "\t" * string(vel) 
            ostr = ostr * "\t" * string(press[i]) * "\n"
            write(ovito, ostr)
        end
        
        write(io, "# t" * "\t" * "u" * "\t" * "ec" * "\t" * "e" * "\t" * "T" * "\t" * "P" * "\n");
        
        write(diffFile, "# t-t0" * "\t" * "Cv" * "\t" * "<r^2>" * "\t" * "Csigmaxy" * "\n");
        #dif(0)
        #dif(1)
        
        gr(0)
    end

    function sample()
        outData = string(t) * "\t" * string(en/npart) * "\t" * string(ecin) * "\t" * string(etot) * "\t" * string(temp)
        outData = outData * "\t" * string(virial) * "\n"
        write(io, outData);
        
        write(ovito, string(npart) * "\n");
        write(ovito, "Lattice=\"" * string(box) * " 0.0 0.0 0.0 "* string(box) * " 0.0 0.0 0.0 0.0\" ")
        write(ovito, "Properties=pos:R:2:Radius:R:1:Type:S:1:Velocity:R:1:Pressure:R:1")
        write(ovito, " Origin=\"0.0 0.0 0.0\"")
        write(ovito, " pbc=\"T T F\"")
        write(ovito, " t=" * string(t) * "\n")
        for i=1:npart
            
            vel = sqrt(v[i, 1] * v[i, 1] + v[i, 2] * v[i, 2])
            
            outData = string(mod(x[i, 1], box)) * "\t" * string(mod(x[i, 2], box))
            outData = outData * "\t" * "0.5" * "\t" * string(i) * "\t" * string(vel) * "\t" * string(press[i]) * "\n"
            write(ovito, outData)
        end
        
        
        if (passo > step0diff) #
            dif(1)
        end
    end

    function force(f, en) # determine the force
        global en = 0.0 # and energy
        global etot = 0.0
        global ecin = 0.0
        global virial = 0.0
        for i=1:npart
            global f[i, 1], f[i, 2] = 0.0, 0.0 # set forces to 0
            global press[i] = 0
        end

        for i=1:npart-1
            for j=i+1:npart # loop over all pairs
                xr = Array{Float64}(undef, 2)
                xr[1], xr[2] = x[i, 1] - x[j, 1], x[i, 2] - x[j, 2]
                xr[1], xr[2] = xr[1] - box*round(xr[1]/box), xr[2] - box*round(xr[2]/box)
                # periodic boundary conditions # section 3.2.2. nothing... lies again...
                # box - diameter of the periodic box

                r2 = xr[1]*xr[1] + xr[2]*xr[2]
                if (r2 < rc2)
                    r2i = 1 / r2
                    r6i = r2i * r2i * r2i
                    ff = 48 * r2i * r6i * (r6i - 0.5) # Lennard-Jones potential
                    global f[i, 1], f[i, 2] = f[i, 1] + ff*xr[1], f[i, 2] + ff*xr[2]
                    global f[j, 1], f[j, 2] = f[j, 1] - ff*xr[1], f[j, 2] - ff*xr[2]
                    global en = en + 4 * r6i * (r6i - 1) - ecut
                    global virial += (f[i, 1]*xr[1] + f[i, 2]*xr[2])
                    global press[i] += (f[i, 1]*xr[1] + f[i, 2]*xr[2])
                end
            end
        end
    end

    function integrate(f, en)
        global sumv[1], sumv[2] = 0.0, 0.0
        global sumv2 = 0.0

        xx = Array{Float64}(undef, 2)
        
        for i=1:npart
            xx[1], xx[2] = 2*x[i, 1] - xm[i, 1] + delt*delt * f[i, 1], 2*x[i, 2] - xm[i, 2] + delt*delt * f[i, 2]
            global v[i, 1], v[i, 2] = (xx[1] - xm[i, 1]) / (2*delt), (xx[2] - xm[i, 2]) / (2*delt)
            global sumv[1], sumv[2]  = sumv[1] + v[i, 1], sumv[2] + v[i, 2]
            global sumv2 = sumv2 + (v[i, 1]*v[i, 1] + v[i, 2]*v[i, 2])
            global xm[i, 1], xm[i, 2] = x[i, 1], x[i, 2]
            global x[i, 1], x[i, 2] = xx[1], xx[2]
            
            global press[i] = press[i] / (npart * Ndim * box * box)
            
        end
        
        global temp = sumv2 / (Ndim * npart)
        
        global ecin = 0.5 * sumv2 / npart
        global etot = (en + 0.5*sumv2) / npart
        
        global virial = (virial / (npart * Ndim * box * box)) + rho * temp
    end

    function md(in_Testing, spacing0, tmax, in_Npart,
         in_SampleStep, in_FolderName, params)
        # simple md program
        println("Starting the initialization...")
        init(in_Testing, spacing0, in_Npart,
         in_SampleStep, in_FolderName, tmax, params)
        println("Initialization finished. Starting main loop...")
        
        println("densidade="*string(rho))
        # todo throw an error if these values are different
        #println(1 / (spacing0 + 2*rad)^2)
        
        while t < tmax # md loop
            force(f, en) # determine the forces # appendix F for tricks 
            
            #if (passo == step0diff) #
                #println("VACF initialized.")
                #dif(0)
            #end
            dif(0)
            
            integrate(f, en) # integrate the equations of motion
            
            global t = t + delt
            global passo += 1
            
            if mod(passo, sampleStep ) == 0
                sample() # sample averages
            end
        end
        println("Main loop ended.")
        
        gr(1)
        
        
        term()
    end

end