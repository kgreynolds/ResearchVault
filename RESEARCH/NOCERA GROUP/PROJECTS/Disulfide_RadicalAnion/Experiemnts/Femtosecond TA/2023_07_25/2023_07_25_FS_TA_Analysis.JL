using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

Phssph_20mM_DMSO_360nm = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_25_FS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_DMSO/Ph_SS_Ph_DMSO_360pump_3ns.csv", header = true, DataFrame)
Phssph_20mM_KHCO2_20mM_DMSO_360nmpump = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_25_FS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_20mM_DMSO/Ph_SS_Ph_formate_DMSO_360pump_3ns.csv", header = true, DataFrame)

# file3 = CSV.read("path", header = true, DataFrame)


### Make Figures 

function Fig1(f = Figure())
    Data = Phssph_20mM_DMSO_360nm
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Spectra of Ph-SS-Ph (20mM) in DMSO", titlesize = 20,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA (mOD)", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(360,690)
    CairoMakie.ylims!(-0.05,nothing)
    width = 4
    x = Data[:,1]

    lines!(x,Data[:,3]*1000,linewidth = width,label = "1 ps",color = cm[n])
    lines!(x,Data[:,4]*1000,linewidth = width,label = "2 ps",color = cm[2*n])
    lines!(x,Data[:,5]*1000,linewidth = width,label = "5 ps",color = cm[3*n])
    lines!(x,Data[:,6]*1000,linewidth = width,label = "10 ps",color = cm[4*n])
    lines!(x,Data[:,7]*1000,linewidth = width,label = "25 ps",color = cm[5*n])
    lines!(x,Data[:,8]*1000,linewidth = width,label = "50 ps",color = cm[6*n])
    lines!(x,Data[:,9]*1000,linewidth = width,label = "100 ps",color = cm[7*n])
    lines!(x,Data[:,10]*1000,linewidth = width,label = "500 ps",color = cm[8*n])
    lines!(x,Data[:,11]*1000,linewidth = width,label = "1000 ps",color = cm[9*n])
    lines!(x,Data[:,12]*1000,linewidth = width,label = "3000 ps",color = cm[9*n])

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

Fig1()

save("JULIA/Outputs/Fig1.png",Fig1())

function Fig2(f = Figure())
    Data = Phssph_20mM_KHCO2_20mM_DMSO_360nmpump
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Spectra of Ph-SS-Ph (20mM) + KHCO2 (20mM) in DMSO", titlesize = 20,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA (mOD)", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(360,700)
    CairoMakie.ylims!(-0.05,nothing)
    width = 4
    x = Data[:,1]

    lines!(x,Data[:,2]*1000,linewidth = width,label = "1 ps",color = cm[n])
    lines!(x,Data[:,3]*1000,linewidth = width,label = "2 ps",color = cm[2*n])
    lines!(x,Data[:,4]*1000,linewidth = width,label = "5 ps",color = cm[3*n])
    lines!(x,Data[:,5]*1000,linewidth = width,label = "10 ps",color = cm[4*n])
    lines!(x,Data[:,6]*1000,linewidth = width,label = "25 ps",color = cm[5*n])
    lines!(x,Data[:,7]*1000,linewidth = width,label = "50 ps",color = cm[6*n])
    lines!(x,Data[:,8]*1000,linewidth = width,label = "100 ps",color = cm[7*n])
    lines!(x,Data[:,9]*1000,linewidth = width,label = "500 ps",color = cm[8*n])
    lines!(x,Data[:,10]*1000,linewidth = width,label = "1000 ps",color = cm[9*n])
    lines!(x,Data[:,11]*1000,linewidth = width,label = "3000 ps",color = cm[9*n])

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :lt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

Fig2()

save("JULIA/Outputs/Fig2.png",Fig2())

### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = "/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_450nm_Kinetics_180us_COnverted.csv"

    Data = (CSV.read(file,header = true,DataFrame))
    t = Data[:,1]
    I = Data[:,2]

    # Set start and end time for dataset
    uppertimebound = 180
    t0 = 0.1

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]


    ## FITTING FUNCTIONS

    # fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
    # p0 = [0.01, 11.0, 0.0]

    fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
    p0 = [0.01, 9.0, 0.02, 15.0, 0.0]

    Plots.plot(time, KineticTrace)
    Plots.plot!(time, fitfunction(p0))

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,BFGS(),Optim.Options(time_limit = 2000))
    values = Optim.minimizer(res) # Kinetics results
    mexp(values)
    Plots.plot!(time, fitfunction(values), linewidth = 3, color = :black)

    # Fit Plots  
    Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
    width = 2

    ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
    xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # xlims!(ax1,350,650)
    # ylims!(nothing,)

    lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
    lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
    lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
    axislegend(ax1, position = :rt, framevisible = false)

    ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
    xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

    lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)

    colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
    rowsize!(Fit_Fig.layout, 1, 300)
    rowsize!(Fit_Fig.layout,2,100)
    resize_to_layout!(Fit_Fig)

    Fit_Fig

    save("JULIA/Outputs/Fit_Fig.png",Fit_Fig, pix_per_unit = 8)

# SOME FUNCTIONS

    #Normalization to x=z
    function Norm(x,y,z)
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/x_interp(z) #normalize to desired value
        return y_norm
    end
    # Normalize to Maximum
    function MaxNorm(x,y)
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/maximum(y_interp) #normalize to desired value
        return y_norm
    end