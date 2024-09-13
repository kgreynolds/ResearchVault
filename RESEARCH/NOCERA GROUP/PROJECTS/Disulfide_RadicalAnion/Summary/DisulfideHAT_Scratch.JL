using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes, Statistics ,CairoMakie; CairoMakie.activate!()
plotly()

IrC = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrC_MeCN/IrC_TA_Spectra_1000avg.csv",header = true,DataFrame)
IrC_PhSSPh = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrC_MeCN/IrC_PhSSPH_TA_Spectra_500avg.csv",header = true,DataFrame)

Ph_SS_Ph_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/Spectra/20mM_PhSSPh_MeCN_TA_Spectrum.csv",header = true,DataFrame)
Ph_SS_Ph_Spectra_150ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,2],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,0.150)
Ph_SS_Ph_Spectra_500ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,3],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,0.500)
Ph_SS_Ph_Spectra_1000ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,4],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,1.0)
Ph_SS_Ph_Spectra_2000ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,5],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,2.0)

Ph_SS_Ph_TBAHCO2_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/Spectra/20mM_PhSSPh_10mM_TBAHCO2_MeCN_TA_Spectrum.csv",header = true,DataFrame)
Ph_SS_Ph_TBAHCO2_Spectra_150ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,2],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,0.150)
Ph_SS_Ph_TBAHCO2_Spectra_500ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,3],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,0.500)
Ph_SS_Ph_TBAHCO2_Spectra_1000ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,4],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,1.0)
Ph_SS_Ph_TBAHCO2_Spectra_2000ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,5],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,2.0)



function IrC_Spectra(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrC TA Spectra",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    # CairoMakie.xlims!(1250,1790)
    CairoMakie.ylims!(-0.025,nothing)
    width = 2
    

    lines!(IrC[:,1],IrC[:,2],linewidth = width,label = "100 ns")
    lines!(IrC[:,1],IrC[:,3],linewidth = width,label = "300 ns")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrC_Spectra()
save("JULIA/Outputs/IrC_Spectra.png",IrC_Spectra())


function IrC_PhSSPh_Spectra(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrC + PhSSPh TA Spectra",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    # CairoMakie.xlims!(1250,1790)
    CairoMakie.ylims!(-0.025,nothing)
    width = 2
    

    lines!(IrC_PhSSPh[:,1],IrC_PhSSPh[:,2],linewidth = width,label = "100 ns")
    lines!(IrC_PhSSPh[:,1],IrC_PhSSPh[:,3],linewidth = width,label = "200 ns")
    lines!(IrC_PhSSPh[:,1],IrC_PhSSPh[:,7],linewidth = width,label = "2000 ns")
    lines!(IrC_PhSSPh[:,1],IrC_PhSSPh[:,8],linewidth = width,label = "4000 ns")


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
IrC_PhSSPh_Spectra()
save("JULIA/Outputs/IrC_PhSSPh_Spectra.png",IrC_PhSSPh_Spectra())

function PhSSPh_Spectra_Comparison(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "IrC + PhSSPh TA Spectra",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    CairoMakie.xlims!(400,600)
    # CairoMakie.ylims!(-0.025,nothing)
    width = 3
    


    lines!(IrC_PhSSPh[:,1],MaxNorm(IrC_PhSSPh[:,7]),linewidth = width,label = "IrC + PhSSPh")



    lines!(Ph_SS_Ph_Spectra[:,1],MaxNorm(Ph_SS_Ph_Spectra_2000ns),linewidth = width,label = "PhSSPh")
    lines!(Ph_SS_Ph_TBAHCO2_Spectra[:,1],MaxNorm(Ph_SS_Ph_TBAHCO2_Spectra_1000ns),linewidth = width,label = "PhSSPh + TBAHCO₂")


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :lb,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
PhSSPh_Spectra_Comparison()
save("JULIA/Outputs/PhSSPh_Spectra_Comparison.png",PhSSPh_Spectra_Comparison())



function SpectralCorrection(Wavelength,Spectrum,Time,Kinetic,CorrWL,CorrTime)
    SpectrumInterp = LinearInterpolation(Wavelength,Spectrum)
    Spectrum_Y = SpectrumInterp(Wavelength)

    KineticInterp = LinearInterpolation(Time,Kinetic)
    
    CorrFactor = KineticInterp[CorrTime]/SpectrumInterp[CorrWL]

    CorrSpectrum = Spectrum_Y .* CorrFactor

    print("factor is:", CorrFactor)
    return CorrSpectrum
 end

 function MaxNorm(y)

    y_norm = y/maximum(y) #normalize to desired value
    return y_norm
end
