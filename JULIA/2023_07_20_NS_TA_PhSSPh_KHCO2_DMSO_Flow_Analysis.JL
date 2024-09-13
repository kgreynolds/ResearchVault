using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

## TA Spectra
Ph_SS_Ph_10mM_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_DMSO_1cm_TA_Spectra_100avg.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_KHCO2_20mM_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_KHCO2_20mM_DMSO_1cm_TA_Spectra_100avg.csv", header = true, DataFrame)


## TA Kinetics

Ph_SS_Ph_10mM_370nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_DMSO_370nmKinetics.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_DMSO_450nmKinetics_100avg_900us_noflow_Converted.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_DMSO_450nmKinetics_100avg_900us_flow_Converted.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_DMSO_600nmKinetics_100avg_900us_flow_Converted.csv", header = true, DataFrame)

Ph_SS_Ph_10mM_KHCO2_20mM_370nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_KHCO2_20mM_DMSO_370nmKinetics.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_KHCO2_20mM_430nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_KHCO2_20mM_DMSO_430nmKinetics_100avg_900us_flow.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_KHCO2_20mM_DMSO_450nmKinetics_100avg_900us_flow_Converted.csv", header = true, DataFrame)
Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_20_NS_TA_PhSSPh_KHCO2_DMSO_Flow/Ph_SS_Ph_10mM_KHCO2_20mM_DMSO_550nmKinetics_100avg_900us_flow.csv", header = true, DataFrame)

Norm_Ph_SS_Ph_10mM_370nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_370nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_370nm_Kinetics_900us_flow[:,2])
Norm_Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,2])
Norm_Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow[:,2])

Norm_Ph_SS_Ph_10mM_KHCO2_20mM_370nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_KHCO2_20mM_370nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_KHCO2_20mM_370nm_Kinetics_900us_flow[:,2])
Norm_Ph_SS_Ph_10mM_KHCO2_20mM_430nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_KHCO2_20mM_430nm_Kinetics_900us_flow[:,1], Ph_SS_Ph_10mM_KHCO2_20mM_430nm_Kinetics_900us_flow[:,2])
Norm_Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow[:,2])
Norm_Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow = MaxNorm(Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow[:,1],Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow[:,2])

### Make Figures 

function PhSSPh_TA_Spectra(f = Figure())
    Data = Ph_SS_Ph_10mM_TA_Spectra
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = font, title = "PhSSPh TA Spectra", titlesize = 20,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # CairoMakie.xlims!(375,600)
    # ylims!(nothing,)
    width = 4
    x = Data[:,1]

    lines!(x,Data[:,2],linewidth = width,label = "50 ns",color = cm[n])
    # lines!(x,Data[:,3],linewidth = width,label = "100 ns",color = cm[2*n])
    # lines!(x,Data[:,4],linewidth = width,label = "150 ns",color = cm[3*n])
    # lines!(x,Data[:,5],linewidth = width,label = "1000 ns",color = cm[4*n])
    lines!(x,Data[:,6],linewidth = width,label = "2000 ns",color = cm[3*n])
    lines!(x,Data[:,7],linewidth = width,label = "5000 ns",color = cm[6*n])
    lines!(x,Data[:,8],linewidth = width,label = "10000 ns",color = cm[9*n])

    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

PhSSPh_TA_Spectra()
save("JULIA/Outputs/PhSSPh_TA_Spectra.png",PhSSPh_TA_Spectra(), pix_per_unit = 8)

function PhSSPh_Formate_TA_Spectra(f = Figure())
    Data = Ph_SS_Ph_10mM_KHCO2_20mM_TA_Spectra
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = font, title = "PhSSPh + Formate TA Spectra", titlesize = 20,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(375,600)
    # ylims!(nothing,)
    width = 4
    x = Data[:,1]

    lines!(x,Data[:,2],linewidth = width,label = "50 ns",color = cm[n])
    lines!(x,Data[:,3],linewidth = width,label = "100 ns",color = cm[2*n])
    lines!(x,Data[:,4],linewidth = width,label = "150 ns",color = cm[3*n])
    lines!(x,Data[:,5],linewidth = width,label = "200 ns",color = cm[4*n])
    lines!(x,Data[:,6],linewidth = width,label = "500 ns",color = cm[5*n])
    lines!(x,Data[:,7],linewidth = width,label = "1000 ns",color = cm[6*n])
    lines!(x,Data[:,8],linewidth = width,label = "2000 ns",color = cm[7*n])
    lines!(x,Data[:,9],linewidth = width,label = "5000 ns",color = cm[8*n])
    lines!(x,Data[:,10],linewidth = width,label = "10000 ns",color = cm[9*n])

    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

PhSSPh_Formate_TA_Spectra()
save("JULIA/Outputs/PhSSPh_Formate_TA_Spectra.png",PhSSPh_Formate_TA_Spectra(), pix_per_unit = 8)

function PhSSPh_TA_Kinetics(f = Figure())

    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = font, title = "PhSSPh TA Kinetics", titlesize = 20,
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(-5,900)
    CairoMakie.ylims!(-0.01, nothing)
    width = 4

    # lines!(Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow[:,1],Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow[:,2],linewidth = width,label = "450 nm NoFlow",color = :blue)
    lines!(Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow,linewidth = width,label = "600 nm",color = :red)
    lines!(Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow,linewidth = width,label = "450 nm",color = :blue)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

PhSSPh_TA_Kinetics()
save("JULIA/Outputs/PhSSPh_TA_Kinetics.png",PhSSPh_TA_Kinetics(), pix_per_unit = 32)


function PhSSPh_Formate_TA_Kinetics(f = Figure())

    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = font, title = "PhSSPh Formate TA Kinetics", titlesize = 20,
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(-0,100)
    CairoMakie.ylims!(-0.01,nothing)
    width = 4

    # lines!(Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow[:,1],Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow[:,2],linewidth = width,label = "450 nm NoFlow",color = :blue)
    lines!(Ph_SS_Ph_10mM_KHCO2_20mM_430nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_KHCO2_20mM_430nm_Kinetics_900us_flow,linewidth = width,label = "430 nm",color = :blue)
    lines!(Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow,linewidth = width,label = "450 nm",color = :red)
    lines!(Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow,linewidth = width,label = "550 nm",color = :black)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

PhSSPh_Formate_TA_Kinetics()
save("JULIA/Outputs/PhSSPh_Formate_TA_Kinetics.png",PhSSPh_Formate_TA_Kinetics(), pix_per_unit = 8)

function PhSSPh_Formate_Kinetics_Comparison(f = Figure())

    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = font, title = "PhSSPh and Formate Kinetics Comparison", titlesize = 20,
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(-0,100)
    CairoMakie.ylims!(-0.01,nothing)
    width = 3

    lines!(Ph_SS_Ph_10mM_370nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_370nm_Kinetics_900us_flow,linewidth = width,label = "370 nm",color = :black)
    # lines!(Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow[:,1],Ph_SS_Ph_10mM_450nm_Kinetics_900us_noflow[:,2],linewidth = width,label = "450 nm NoFlow",color = :blue)
    # lines!(Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow,linewidth = width,label = "450 nm",color = :blue)
    # lines!(Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_600nm_Kinetics_900us_flow,linewidth = width,label = "600 nm",color = :black)

    lines!(Ph_SS_Ph_10mM_KHCO2_20mM_370nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_KHCO2_20mM_370nm_Kinetics_900us_flow,linewidth = width,label = "370 nm + Form",color = :green)
    # lines!(Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_KHCO2_20mM_450nm_Kinetics_900us_flow,linewidth = width,label = "450 nm +Formate",color = :red)
    # lines!(Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow[:,1],Norm_Ph_SS_Ph_10mM_KHCO2_20mM_550nm_Kinetics_900us_flow,linewidth = width,label = "550 nm + Formate",color = :purple)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

PhSSPh_Formate_Kinetics_Comparison()
save("JULIA/Outputs/PhSSPh_Formate_Kinetics_Comparison.png",PhSSPh_Formate_Kinetics_Comparison(), pix_per_unit = 32)

### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim


   Data = Ph_SS_Ph_10mM_450nm_Kinetics_900us_flow
    t = Data[:,1]
    I = Data[:,2]

    # Set start and end time for dataset
    uppertimebound = 10
    t0 = 0.2


    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]


    ## FITTING FUNCTIONS

    # fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
    # p0 = [0.01, 11.0, 0.0]

    fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
    p0 = [0.05, 1.19, 0.04, 8.7, 0.02]

    # fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
    # p0 = [0.05, 1.19, 0.04, 8.7, 0.03, 39.0, 0.011]

    Plots.plot(time, KineticTrace)
    Plots.plot!(time, fitfunction(p0))

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,BFGS(),Optim.Options(time_limit = 2000))
    values = Optim.minimizer(res) # Kinetics results
    fitfunction(values)
    Plots.plot!(time, fitfunction(values), linewidth = 3, color = :black)

    # Fit Plots  
    Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
    width = 2

    ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Ph-SS-Ph 370nm Triexp Fit",
    xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # xlims!(ax1,350,650)
    # ylims!(nothing,)

    lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
    lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "τ₁ = 1.2 μs, τ₂ = 6.9 μs, τ₃ = 39.2 μs", color = :red)
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

    save("JULIA/Outputs/Ph-SS-Ph 370nm Triexp Fit.png",Fit_Fig, pix_per_unit = 32)

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