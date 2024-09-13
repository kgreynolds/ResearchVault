using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes, Statistics ,CairoMakie; CairoMakie.activate!()
plotly()

IrK_PhssPh_450TAKin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/RAW_Data/NS/2023_08_24/IrK_PhSSPh_450nm_TAKinetics_500avg_18us.trc", skipto = 6,drop = [3], header = ["Time", "Ampl"], DataFrame)
Plots.plot(IrK_PhssPh_450TAKin[:,1]*1000000,IrK_PhssPh_450TAKin[:,2])

IrK_PhSSPh_480nm_EmissionKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_PhSSPh_480nm_EmissionKinetics_25avg_900ns.csv", header = true, DataFrame)
IrK_PhSSPh_480nm_EmissionKinetics = Resize(IrK_PhSSPh_480nm_EmissionKinetics,1000)
Plots.plot!(IrK_PhSSPh_480nm_EmissionKinetics[:,1],(IrK_PhSSPh_480nm_EmissionKinetics[:,2] .* 0.265).- 0.0768)

(IrK_PhssPh_450TAKin[:,2] .- ((IrK_PhSSPh_480nm_EmissionKinetics[:,2] .* 0.265).- 0.0768))

# UvVis

IrK_MeCN_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/RAW_Data/UV-Vis/2023_09_19/Ir(III)K_MeCN_UvVis.csv",header = true,DataFrame)
IrK_MeCN_UvVis = InterpOnRange(IrK_MeCN_UvVis,300,650,1)
IrK_exCoeff = ((IrK_MeCN_UvVis[:,2] ./ IrK_MeCN_UvVis[81,2]) .* 6165.95) ./ 1000 # based on Exctinction coeff from Chem. Mater., Vol. 17, No. 23, 2005 5715


IrK_MeCN_UvVis[:,1]

# TA Spectra

IrK = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/TA_Spectra/IrK_TA_Spectra_1000avg.csv",header = true,DataFrame)
IrK_PhSSPh = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/TA_Spectra/IrK_PhSSPH_TA_Spectra_500avg.csv",header = true,DataFrame)


# Emission Kinetics
IrK_480nm_Emission = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_480nm_EmissionKinetics_500avg_900ns.csv", header = true, DataFrame)
IrK_480nm_Emission = Resize(IrK_480nm_Emission,1000)

Plots.plot(IrK_480nm_Emission[:,1],IrK_480nm_Emission[:,2])
IrK_480nm_Emission_Fit = MonoFitIRF(IrK_480nm_Emission,2,0.03,0.8,[-0.052, 0.036, 0.011, 0.489, 0.186, -0.029],true,false,1.0) ## τ = 186 ns

IrK_450nm_EmissionKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_450nm_EmissionKinetics_500avg_900ns.csv",header = true,DataFrame)
IrK_450nm_EmissionKinetics = Resize(IrK_450nm_EmissionKinetics,1000)
Plots.plot(IrK_450nm_EmissionKinetics[:,1],IrK_450nm_EmissionKinetics[:,2])

IrK_PhSSPh_480nm_EmissionKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_PhSSPh_480nm_EmissionKinetics_25avg_900ns.csv", header = true, DataFrame)
IrK_PhSSPh_480nm_EmissionKinetics = Resize(IrK_PhSSPh_480nm_EmissionKinetics,1000)
Plots.plot!(IrK_PhSSPh_480nm_EmissionKinetics[:,1],IrK_PhSSPh_480nm_EmissionKinetics[:,2])
IrK_PhSSPh_480nm_EmissionKinetics_Fit = MonoFitIRF(IrK_PhSSPh_480nm_EmissionKinetics,2,0.03,0.8,[-0.062, 0.031, 0.019, 0.491, 0.082, 0.0001],true,true,1) # τ = 82 ns



# TA Kinetics
IrK_450nm_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_450nm_TAKinetics_500avg_900ns_Converted.csv",header = true,DataFrame)
IrK_450nm_TAKinetics = Resize(IrK_450nm_TAKinetics,1000)
Plots.plot!(IrK_450nm_TAKinetics[:,1],IrK_450nm_TAKinetics[:,2])

IrK_400nm_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_400nm_TAKinetics_500avg_450ns_Converted.csv",header = true,DataFrame)
IrK_400nm_TAKinetics = Resize(IrK_400nm_TAKinetics,1000)
Plots.plot(IrK_400nm_TAKinetics[:,1],IrK_400nm_TAKinetics[:,2])
IrK_400nm_TAKinetics_Fit = MonoFitIRF(IrK_400nm_TAKinetics,2,0.0445,0.4,[0.034, 0.042, 0.014, 0.488, 0.174, -0.029],true,true,1.0) # τ = 174 ns

IrK_PhSSPh_450nm_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_PhSSPh_450nm_TAKinetics_500avg_18us_Converted.csv",header = true,DataFrame)
IrK_PhSSPh_450nm_TAKinetics = Resize(IrK_PhSSPh_450nm_TAKinetics,1000)
Plots.plot!(IrK_PhSSPh_450nm_TAKinetics[:,1],IrK_PhSSPh_450nm_TAKinetics[:,2])
IrK_PhSSPh_450nm_TAKinetics_FitIRF = BiFitIRF(IrK_PhSSPh_450nm_TAKinetics,2,0.5,17,[0.013, 0.699, 0.0012, 0.999, 6.666, 0.860, 1.119, 0.229],true,true,1) # τ₁ = 1.12 μs, τ₂ = 6.67 μs


function IrK_UvVis(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrK TA Uv-Vis",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 1, xminortickalign =1, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ϵ (10³ M⁻¹ cm⁻¹)", ylabelsize = 20, ytickalign = 1,yminortickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    CairoMakie.xlims!(300, 625)
    ax.yticks = 0:2:27
    # CairoMakie.ylims!(-0.025,nothing)
    width = 3
    
    # lines!(IrK[:,1],IrK[:,2],linewidth = width,label = "100 ns")
    # lines!(IrK_MeCN_UvVis[:,1],IrK_MeCN_UvVis[:,2],linewidth = width,label = "Ir(III)K MeCN")
    lines!(IrK_MeCN_UvVis[:,1],IrK_exCoeff,linewidth = width,label = "Ir(III)K MeCN ex Coeff")

    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend(position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrK_UvVis()
save("JULIA/Outputs/IrK_UvVis.eps",IrK_UvVis())


function IrK_Spectra(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrK TA Spectra",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 1, xminortickalign =1, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 1,yminortickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    CairoMakie.xlims!(375, 625)
    # CairoMakie.ylims!(-0.025,nothing)
    width = 3
    

    # lines!(IrK[:,1],IrK[:,2],linewidth = width,label = "100 ns")
    lines!(IrK[:,1],IrK[:,3],linewidth = width,label = "300 ns")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrK_Spectra()
save("JULIA/Outputs/IrK_Spectra.eps",IrK_Spectra())

function IrK_kinetics(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrK Kinetics",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 1, xminortickalign =1, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "I", ylabelsize = 20, ytickalign = 1,yminortickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    # CairoMakie.xlims!(375, 625)
    # CairoMakie.ylims!(-0.025,nothing)
    width = 3
    

    lines!(IrK_480nm_Emission_Fit[:,1],IrK_480nm_Emission_Fit[:,2],linewidth = width,label = "480nm emission")
    # lines!(IrK_450nm_TAKinetics[:,1],IrK_450nm_TAKinetics[:,2],linewidth = width,label = "450nm TA Kinetics")


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rb,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrK_kinetics()
save("JULIA/Outputs/IrK_kinetics.eps",IrK_kinetics())

function IrK_PhSSPh_Spectra(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrK + PhSSPh TA Spectra",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 1, xminortickalign =1, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 1,yminortickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.winter, 4),))
    CairoMakie.xlims!(375, 625)
    CairoMakie.ylims!(-0.008,0.031)
    
    width = 3
    
    CairoMakie.scatter!(IrK[:,1],IrK[:,3]*2.5,markersize = 5, color = :red) #IrK TA Spectrum for reference


    lines!(IrK_PhSSPh[:,1],IrK_PhSSPh[:,2],linewidth = width,label = "100 ns")
    lines!(IrK_PhSSPh[:,1],IrK_PhSSPh[:,3],linewidth = width,label = "200 ns")
    lines!(IrK_PhSSPh[:,1],IrK_PhSSPh[:,7],linewidth = width,label = "2000 ns")
    lines!(IrK_PhSSPh[:,1],IrK_PhSSPh[:,8],linewidth = width,label = "4000 ns")


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)", textsize = 25)
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrK_PhSSPh_Spectra()
save("JULIA/Outputs/IrK_PhSSPh_TASpectra.png",IrK_PhSSPh_Spectra())

function IrK_PhSSPh_TA_Kinetics(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrK + PhSSPh TA Spectra", 
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 1, xminortickalign =1, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 1,yminortickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    # ax.xticks= 0:1:17
    CairoMakie.xlims!(-0.1,10)
    # CairoMakie.ylims!(-0.025,nothing)
    width = 4
    

    lines!(IrK_450nm_TAKinetics[:,1],MaxNorm(IrK_450nm_TAKinetics[:,2]),linewidth = width,label = "IrK")
    lines!(IrK_PhSSPh_450nm_TAKinetics[:,1],MaxNorm(IrK_PhSSPh_450nm_TAKinetics[:,2]),linewidth = width,label = "IrK + PhSSPh")
    # lines!(IrK_PhSSPh_450nm_TAKinetics_FitIRF[:,1],IrK_PhSSPh_450nm_TAKinetics_FitIRF[:,3],linewidth = 3,label = "IrK + PhSSPh Fit")


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend(position = :rb,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrK_PhSSPh_TA_Kinetics()
save("JULIA/Outputs/IrK_PhSSPh_TA_Kinetics.png",IrK_PhSSPh_TA_Kinetics())


function Emission_Quenching(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrK + PhSSPh 480nm Emission Kinetics",palette = (color = palette(ColorSchemes.hsv, 4),),
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 1, xminortickalign =1, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized Intensity", ylabelsize = 20, ytickalign = 1,yminortickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0))
    # CairoMakie.xlims!(1250,1790)
    # CairoMakie.ylims!(-0.025,nothing)
    width = 3
    

    # lines!(IrK[:,1],IrK[:,2],linewidth = width,label = "100 ns")
    lines!(IrK_480nm_Emission[:,1],MaxNorm(-IrK_480nm_Emission[:,2]),linewidth = width,label = "IrK Emission", color = :red)
    lines!(IrK_PhSSPh_480nm_EmissionKinetics[:,1],MaxNorm(-IrK_PhSSPh_480nm_EmissionKinetics[:,2]),linewidth = width,label = "IrK + PhSSPh Emission", color = :green)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend(position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    text!(0.3,0.25,text = "τ = 186 ns", color = :red)
    text!(0.1,0.05,text = "τ = 82 ns", color = :green)

    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
Emission_Quenching()
save("JULIA/Outputs/Emission_Quenching.eps",Emission_Quenching())









# SOME FUNCTIONS
function MonoFit(file,column,t0,uppertimebound, residuals, ReturnFit, coeffs::Vector{Float64})

    Data = file
    t = Data[:,1]
    I = Data[:,column]

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
    p0 = coeffs

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual =  KineticTrace .- fitfunction(values)
    print("Fit Coeff:", values)

    if ReturnFit == true
        return hcat(time,KineticTrace,fitfunction(values), residual)
    

    elseif residuals == true
        Fit_Fig = Figure(figure_padding = 25,fontsize =20, resolution = (600,600))
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

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        Fit_Fig
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
    else

        Plots.plot(time, KineticTrace, color = :blue, label = "Data")
        Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
        Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")

    end

    # return fitfunction(values)
end

function BiFit(file,column,t0,uppertimebound, residuals, ReturnFit, coeffs::Vector{Float64})

    Data = file
    t = Data[:,1]
    I = Data[:,column]

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
    p0 = coeffs

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual =  KineticTrace .- fitfunction(values)
    print("Fit Coeff:", values)

    if ReturnFit == true
        return hcat(time,KineticTrace,fitfunction(values), residual)
    

    elseif residuals == true
        Fit_Fig = Figure( figure_padding = 25,fontsize =20, resolution = (600,600))
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

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        Fit_Fig
        # save("JULIA/Outputs/BiFit_Fig.png",Fit_Fig)
    else

        Plots.plot(time, KineticTrace, color = :blue, label = "Data")
        Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
        Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")

    end

    # return fitfunction(values)
end

function TriFit(file,column,t0,uppertimebound, residuals, ReturnFit, coeffs::Vector{Float64})

    Data = file
    t = Data[:,1]
    I = Data[:,column]

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
    p0 = coeffs

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual =  KineticTrace .- fitfunction(values)
    print("Fit Coeff:", values)

    if ReturnFit == true
        return hcat(time,fitfunction(values), residual)
    

    elseif residuals == true
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

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        Fit_Fig
        # save("JULIA/Outputs/TriFit_Fig.png",Fit_Fig)
    else

        Plots.plot(time, KineticTrace, color = :blue, label = "Data")
        Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
        Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")

    end

    # return fitfunction(values)
end

function MonoFitIRF(file,column,t0,uppertimebound,coeffs, residuals, ReturnFit, reduce)

    function Fitfunc(x,coeffs)
    
        Gaussian(p) =  exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
        mexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ p[6])
    
        IRF = Gaussian(coeffs)
        Exp = mexp(coeffs)
    
        function DirectLinearConvolution(f,g)
            N = length(f) # for This to work both f and g need to be the same size
            M = length(g)
        
            Conv = zeros(N)
            for n=1:N
                for m=1:M
                    if n-m+1 > 0
                        Conv[n] = Conv[n] + f[m] * g[n-m+1]
                    end
                    # n+1 <= m
                end
            end
            return Conv
        end
    
        Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))
    
        return Convolution
    end

    Data = file

    rebin = round(length(Data[:,1])/reduce)
    stepsize = (Data[end,1]-Data[1,1])/rebin
    InterpData = LinearInterpolation(Data[:,1],Data[:,2])
    
    t = collect(Data[1,1]:stepsize:Data[end,1])
    I = InterpData(t)

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)


    res = optimize(peakfit0, coeffs,LBFGS(),Optim.Options(time_limit = 200))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual = KineticTrace .- Fitfunc(time,values)
    print("Fit Coeff:", values)

    df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)

    if ReturnFit == true
        return df #hcat(time,KineticTrace,Fitfunc(time,values),residual) DataFrame
    
    elseif residuals == true
        Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
        width = 2

        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
        xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
        # xlims!(ax1,350,650)
        # ylims!(nothing,)

        lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
        lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = 2, label = "Fit", color = :red)
        lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = 2, label = "Guess", color = :green)
        axislegend(ax1, position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
        xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        Fit_Fig
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
    else

        Plots.plot(time, KineticTrace, color = :blue, label = "Data")
        Plots.plot!(time, Fitfunc(time,coeffs), color = :green, label = "Guess")
        Plots.plot!(time, Fitfunc(time,values), linewidth = 3, color = :red, label = "Fit")

    end

    # return fitfunction(values)
end

function BiFitIRF(file,column,t0,uppertimebound,coeffs, residuals, ReturnFit, reduce)

    function Fitfunc(x,coeffs)
    
        Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
        Biexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ (p[6] .* exp.((-(x) ./ p[7]))) .+ p[8])
    
        IRF = Gaussian(coeffs)
        Exp = Biexp(coeffs)
    
        function DirectLinearConvolution(f,g)
            N = length(f) # for This to work both f and g need to be the same size
            M = length(g)
        
            Conv = zeros(N)
            for n=1:N
                for m=1:M
                    if n-m+1 > 0
                        Conv[n] = Conv[n] + f[m] * g[n-m+1]
                    end
                    # n+1 <= m
                end
            end
            return Conv
        end
    
        Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))
    
        return Convolution
    end

    Data = file

    rebin = round(length(Data[:,1])/reduce)
    stepsize = (Data[end,1]-Data[1,1])/rebin
    InterpData = LinearInterpolation(Data[:,1],Data[:,2])
    
    t = collect(Data[1,1]:stepsize:Data[end,1])
    I = InterpData(t)

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)


    res = optimize(peakfit0, coeffs,LBFGS(),Optim.Options(time_limit = 200))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual = KineticTrace .- Fitfunc(time,values)
    print("Fit Coeff:", values)

    if ReturnFit == true
        return hcat(time,KineticTrace,Fitfunc(time,values),residual)
    
    elseif residuals == true
        Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
        width = 2

        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
        xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
        # xlims!(ax1,350,650)
        # ylims!(nothing,)

        lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
        lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = 2, label = "Fit", color = :red)
        lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = 2, label = "Guess", color = :green)
        axislegend(ax1, position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
        xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        Fit_Fig
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
    else

        Plots.plot(time, KineticTrace, color = :blue, label = "Data")
        Plots.plot!(time, Fitfunc(time,coeffs), color = :green, label = "Guess")
        Plots.plot!(time, Fitfunc(time,values), linewidth = 3, color = :red, label = "Fit")

    end

    # return fitfunction(values)
end


#Normalization to x=z
function Norm(x,y,z)
    x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
    y_interp = x_interp(x) # compute y values for the x vector
    y_norm = y_interp/x_interp(z) #normalize to desired value
    return y_norm
end
# Normalize to Maximum
function MaxNorm(y)
    y_norm = y/maximum(abs.(y)) #normalize to desired value
    return y_norm
end

function SpectralCorrection(Wavelength,Spectrum,Time,Kinetic,CorrWL,CorrTime)
    SpectrumInterp = LinearInterpolation(Wavelength,Spectrum)
    Spectrum_Y = SpectrumInterp(Wavelength)

    KineticInterp = LinearInterpolation(Time,Kinetic)
    
    CorrFactor = KineticInterp[CorrTime]/SpectrumInterp[CorrWL]

    CorrSpectrum = Spectrum_Y .* CorrFactor

    print("factor is:", CorrFactor)
    return CorrSpectrum
end


function Resize(Data, length)

    x_out = collect(range(Data[1,1],Data[end,1],length))
    Interp  = LinearInterpolation(Data[:,1],Data[:,2])
    y_out = Interp(x_out)

    df = DataFrame(x = x_out,y = y_out)
    return df
end

function InterpOnRange(Data,start,stop, interval)

    x_out = collect(start:interval:stop)
    Interp  = LinearInterpolation(Data[:,1],Data[:,2])
    y_out = Interp(x_out)

    df = DataFrame(x = x_out,y = y_out)
    return df
end
