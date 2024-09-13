using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

PhSSPh_20mM_KHCO2_5mM_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_DMSO_TA_Spectra.csv", header = true, DataFrame)


## Flowing (3 mL/min)
PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_455nm_45us_TA_Kinetics_100avg_Converted.csv", header = true, DataFrame)
Norm_PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics = MaxNorm(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1],PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,2])
Plots.plot(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1],PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,2] , xlimits = (-1,40))

PhSSPh_20mM_KHCO2_5mM_365nm_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_365nm_2us_TA_Kinetics_100avg_Converted.csv", header = true, DataFrame)


## Not Flowing
PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_455nm_45us_TA_Kinetics_100avg_NoFlow_Conveted.csv", header = true, DataFrame)
Norm_PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow = MaxNorm(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,1],PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,2])
Plots.plot!(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,1],Norm_PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow, xlimits = (-1,40), ylimits = (-0.1, 1.1))

PhSSPh_20mM_KHCO2_5mM_365nm_Kinetics_NoFLow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_365nm_2us_TA_Kinetics_100avg_Noflow_Converted.csv", header = true, DataFrame)


# Data from 2023_07_20 Ph-SS-Ph 10mM in DMSO by itslef flowing (but slower)
Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_DMSO_450nmKinetics_100avg_900us_flow_Converted.csv", header = true, DataFrame)
Norm_Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,2])
Plots.plot(Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow , xlimits = (-1,40), ylimits = (-0.1, 1.1))

### Make Figures 

function Fig1(f = Figure())
    Data1 = MaxNorm(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1],PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,2])
    Data2 = MaxNorm(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,1],PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,2])
    #Set color Scheme (n increments over the colorscheme depends on dataframe size)

    #Make Cairo Mackie figure
    
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Kinetics of Ph-SS-Ph (20mM) + KHCO2 (5 mM) in DMSO", titlesize = 20,
    xlabel = "time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Norm ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(0,5)
    CairoMakie.ylims!(-0.05,nothing)
    width = 3
    x1 = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1]
    x2 = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,1]
    lines!(x1,Data1,linewidth = width,label = "Flow",color = :blue)
    lines!(x2,Data2,linewidth = width,label = "No Flow",color = :red)

    axislegend("455 nm", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    ax2 = CairoMakie.Axis(f[2,1], titlefont = :bold,
    xlabel = "time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Norm ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(-0,40)
    CairoMakie.ylims!(-0.05,nothing)
    width = 3
    x1 = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1]
    x2 = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics_NoFLow[:,1]
    lines!(ax2,x1,Data1,linewidth = width,label = "Flow",color = :blue)
    lines!(ax2,x2,Data2,linewidth = width,label = "No Flow",color = :red)

    axislegend("455 nm", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    
    # Label(f[1,1,TopLeft()],"(a)")
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end


Fig1()

save("JULIA/Outputs/Fig1.png",Fig1())

function Fig2(f = Figure())
    Data1 = Norm(Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,2], 0.5)
    Data2 = Norm(PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1],PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,2],0.5)
    #Set color Scheme (n increments over the colorscheme depends on dataframe size)

    #Make Cairo Mackie figure
    
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "Quenching Kinetics of Ph-SS-Ph (20mM) + KHCO2 (5 mM) in DMSO", titlesize = 20,
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Norm ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(-1,40)
    CairoMakie.ylims!(-0.05,nothing)
    width = 3
    x1 = Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1]
    x2 = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics[:,1]


    lines!(x1,Data1,linewidth = width,label = "Ph-SS-Ph",color = :blue)
    lines!(x2,Data2,linewidth = width,label = "Ph-SS-Ph + KHCO₂",color = :red)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("450 nm", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

Fig2()

save("JULIA/Outputs/Fig2.png",Fig2())

function Fig3(f = Figure())
    Data = PhSSPh_20mM_KHCO2_5mM_TA_Spectra
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Spectra of Ph-SS-Ph (20 mM) + KHCO2 (5 mM) in DMSO Flow (3 ml/min)", titlesize = 20,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA (mOD)", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(375,600)
    # CairoMakie.ylims!(-0.05,nothing)
    width = 4
    x = Data[:,1]

    lines!(x,(Data[:,2] * 4),linewidth = width,label = "100 ns",color = cm[n])
    # lines!(x,Data[:,3],linewidth = width,label = "200 ns",color = cm[2*n])
    # lines!(x,Data[:,4],linewidth = width,label = "300 ns",color = cm[3*n])
    lines!(x,Data[:,5] * 1.65,linewidth = width,label = "500 ns",color = cm[4*n])
    lines!(x,Data[:,6] * 2.4,linewidth = width,label = "1500 ns",color = cm[5*n])
    lines!(x,Data[:,7] * 4.8,linewidth = width,label = "5000 ns",color = cm[6*n])


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :lb,nbanks = 3,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

Fig3()

save("JULIA/Outputs/Fig3.png",Fig3())

### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = "/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_455nm_45us_TA_Kinetics_100avg_Converted.csv"
    Data = (CSV.read(file,header = true,DataFrame))

    Data = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics

    t = Data[:,1]
    I = Data[:,2]

    # Set start and end time for dataset
    uppertimebound = 41
    t0 = 0.5

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]


    ## FITTING FUNCTIONS

    # fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
    # p0 = [0.01, 11.0, 0.0]

    fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
    p0 = [0.01, 9.0, 0.02, 15.0, 0.0]

    # fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
    # p0 = [0.05, 1.19, 0.04, 8.7, 0.03, 39.0, 0.011]

    Plots.plot(time, KineticTrace, xaxis = :log)
    Plots.plot!(time, fitfunction(p0))

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,BFGS(),Optim.Options(time_limit = 2000))
    values = Optim.minimizer(res) # Kinetics results
    mexp(values)
    Plots.plot!(time, fitfunction(values), linewidth = 3, color = :black)

    # Fit Plots  
    Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
    width = 2

    ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Ph-SS-Ph (20 mM) + KHCO₂ (5 mM) TA Kinetics", #xscale = log10,
    xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # xlims!(ax1,350,650)
    # ylims!(nothing,)

    lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "455 nm", color = :blue)
    lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "τ₁ = 2.9 μs, τ₂ = 9.7 μs", color = :red)
    # lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
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

    save("JULIA/Outputs/Fit_Fig.png",Fit_Fig)

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

    function GatherAllignedPoints(coordinate_sets::Vector{Vector{Vector{Float64}}},x_interval)
        interps::Vector{Any} = []
        out::Vector{Vector{Float64}} = []
    
        for v in coordinate_sets
            x,y = v
            push!(interps,LinearInterpolation(x,y))
        end
    
        for i in interps
            push!(out,i(collect(x_interval)))
        end
    
        return out
    end
    