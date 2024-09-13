using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes, Statistics ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

## TA Spectra
Ph_SS_Ph_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/Spectra/20mM_PhSSPh_MeCN_TA_Spectrum.csv",header = true,DataFrame)
Ph_SS_Ph_Spectra_150ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,2],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,0.150)
Ph_SS_Ph_Spectra_500ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,3],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,0.500)
Ph_SS_Ph_Spectra_1000ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,4],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,1.0)
Ph_SS_Ph_Spectra_2000ns = SpectralCorrection(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra[:,5],Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],450,2.0)

Ph_SS_Ph_TBAHCO2_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/Spectra/20mM_PhSSPh_10mM_TBAHCO2_MeCN_TA_Spectrum.csv",header = true,DataFrame)
Ph_SS_Ph_TBAHCO2_Spectra_150ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,2],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,0.150)
Ph_SS_Ph_TBAHCO2_Spectra_500ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,3],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,0.500)
Ph_SS_Ph_TBAHCO2_Spectra_1000ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,4],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,1.0)
Ph_SS_Ph_TBAHCO2_Spectra_2000ns = SpectralCorrection(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra[:,5],Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],450,2.0)

## TA Kinetics
Ph_SS_Ph_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/20mM_PhSSPh_MeCN_ArPurge_TAKinetics_500avg.csv", header = true, DataFrame)


Ph_SS_Ph_TBAHCO2_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/20mM_PhSSPh_10mM_TBAHCO2_MeCN_ArPurge_420nm_TAKinetics_500avg.csv", header = true, DataFrame)


Ph_SS_Ph_Kinetics_420nm = BiFit(Resize(Ph_SS_Ph_Kinetics,2,1000),2,0.13,15,true,true,[0.02, 0.4, 0.02, 3.0, 0.0]) # t1 = 494ns t2 = 3648 ns
Ph_SS_Ph_Kinetics_450nm = BiFit(Resize(Ph_SS_Ph_Kinetics,3,1000),2,0.13,15,true,true,[0.01, 0.4, 0.02, 3.0, 0.0]) #t1 = 485ns t2 = 3490 ns
Ph_SS_Ph_Kinetics_470nm = BiFit(Resize(Ph_SS_Ph_Kinetics,4,1000),2,0.13,15,true,true,[0.01, 0.4, 0.02, 3.0, 0.0])# t1 = 504 ns t2 = 3627 ns

Ph_SS_Ph_TBAHCO2_Kinetics_420 = BiFit(Resize(Ph_SS_Ph_TBAHCO2_Kinetics,2,1000),2,0.13,15,true,true,[0.05, 0.9, 0.02, 3.5, 0.0]) # t1 = 1155 ns t2 = 4323 ns
Ph_SS_Ph_TBAHCO2_Kinetics_450 = BiFit(Resize(Ph_SS_Ph_TBAHCO2_Kinetics,3,1000),2,0.13,15,true,true,[0.05, 0.9, 0.02, 3.5, 0.0]) # t1 = 802 ns t2 = 3138 ns
Ph_SS_Ph_TBAHCO2_Kinetics_550 = BiFit(Resize(Ph_SS_Ph_TBAHCO2_Kinetics,4,1000),2,0.13,15,true,true,[0.05, 0.9, 0.02, 3.5, 0.0]) # t1  = 1212 ns t2 = 4518 ns

## FTIR Photolysis
Ph_SS_Ph_Photolysis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/FTIR/FTIR_Photolysis_PhSSPh_20mM_MeCN_390um_Ar.csv", header = true, DataFrame)
Ph_SS_Ph_TBAHCO2_Photolysis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/FTIR/FTIR_Photolysis_PhSSPh_20mM_TBAHCO2_10mM_MeCN_390um_Ar.csv", header = true, DataFrame)



### Make Figures 

function PhSSPh_Spectra(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "Ph-SS-Ph (20 mM) TA Spectra in MeCN",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.winter, 3),))
    # CairoMakie.xlims!(1250,1790)
    CairoMakie.ylims!(-0.025,nothing)
    width = 3
    

    # lines!(x,Data[:,2],linewidth = width,label = "420 nm",color = cm[n])
    # lines!(x,Data[:,3],linewidth = width,label = "450 nm",color = cm[2*n])
    # lines!(x,Data[:,4],linewidth = width,label = "550 nm",color = cm[3*n])

    lines!(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra_150ns,linewidth = width,label = "150 ns")
    lines!(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra_1000ns,linewidth = width,label = "1000 ns")
    lines!(Ph_SS_Ph_Spectra[:,1],Ph_SS_Ph_Spectra_2000ns,linewidth = width,label = "2000 ns")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 1,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
PhSSPh_Spectra()
save("JULIA/Outputs/PhSSPh_Spectra_MeCN.png",PhSSPh_Spectra())

function PhSSPh_TBAHCO2_Spectra(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "Ph-SS-Ph (20 mM) + TBAHCO₂ (10 mM) TA Spectra in MeCN",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,xautolimitmargin = (0.0, 0.0),
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, yautolimitmargin = (0.0, 0.2),palette = (color = palette(ColorSchemes.winter, 3),))
    # CairoMakie.xlims!(1250,1790)
    CairoMakie.ylims!(-0.025,nothing)
    width = 3
    

    # lines!(x,Data[:,2],linewidth = width,label = "420 nm",color = cm[n])
    # lines!(x,Data[:,3],linewidth = width,label = "450 nm",color = cm[2*n])
    # lines!(x,Data[:,4],linewidth = width,label = "550 nm",color = cm[3*n])

    lines!(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra_150ns,linewidth = width,label = "150 ns")
    lines!(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra_500ns,linewidth = width,label = "500 ns")
    lines!(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra_1000ns,linewidth = width,label = "1000 ns")
    # lines!(Ph_SS_Ph_TBAHCO2_Spectra[:,1],Ph_SS_Ph_TBAHCO2_Spectra_2000ns,linewidth = width,label = "2000 ns")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 1,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
PhSSPh_TBAHCO2_Spectra()
save("JULIA/Outputs/PhSSPh_TBAHCO2_Spectra_MeCN.png",PhSSPh_TBAHCO2_Spectra())


function PhSSPh_Kinetics(f = Figure())
    Data = Ph_SS_Ph_Kinetics

    ax = CairoMakie.Axis(f[1,1], title = "Ph-SS-Ph (20 mM) TA Kinetics",
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 5),))
    CairoMakie.xlims!(-1,10)
    CairoMakie.ylims!(nothing,0.02)
    width = 2
    x = Data[:,1]

    # lines!(x,Data[:,2],linewidth = width,label = "420 nm",color = cm[n])
    # lines!(x,Data[:,3],linewidth = width,label = "450 nm",color = cm[2*n])
    # lines!(x,Data[:,4],linewidth = width,label = "550 nm",color = cm[3*n])

    # lines!(x,MaxNorm(Data[:,2]),linewidth = width,label = "420 nm")
    lines!(x,Data[:,3],linewidth = width,label = "450 nm")
    # lines!(x,MaxNorm(Data[:,4]),linewidth = width,label = "470 nm")

    lines!(PhSSPh_TA_Spectra[2:end,1]/1000,PhSSPh_TA_Spectra[2:end,12],linewidth = width,label = "450 nm Kinetics")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Wavelength", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
PhSSPh_Kinetics()
save("JULIA/Outputs/PhSSPh_Kinetics_old2.png",PhSSPh_Kinetics(),pt_per_unit = 20)

function PhSSPh_TBAHCO2_Kinetics(f = Figure())
    Data = Ph_SS_Ph_TBAHCO2_Kinetics

    ax = CairoMakie.Axis(f[1,1], title = "Ph-SS-Ph (20 mM) + TBAHCO₂ (10 mM) TA Kinetics",
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv,5; rev = false),))
    CairoMakie.xlims!(nothing,15)
    CairoMakie.ylims!(-0.05,nothing)
    width = 2
    x = Data[:,1]

    # lines!(x,Data[:,2],linewidth = width,label = "420 nm")
    # lines!(x,Data[:,3],linewidth = width,label = "450 nm")
    # lines!(x,Data[:,4],linewidth = width,label = "470 nm")

    lines!(x,MaxNorm(Data[:,2]),linewidth = width,label = "420 nm")
    lines!(x,MaxNorm(Data[:,3]),linewidth = width,label = "450 nm")
    lines!(x,MaxNorm(Data[:,4]),linewidth = width,label = "550 nm")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Wavelength", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
PhSSPh_TBAHCO2_Kinetics()
save("JULIA/Outputs/PhSSPh_TBAHCO2_Kinetics.png",PhSSPh_TBAHCO2_Kinetics())

function KineticsComparison_420nm(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "420nm Kinetics Comparison",
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.viridis, 5),))
    
    CairoMakie.ylims!(-0.1,nothing)

    width = 3
    
    lines!(Ph_SS_Ph_Kinetics[:,1],MaxNorm(Ph_SS_Ph_Kinetics[:,2]),linewidth = width,label = "PhSSPh", color = :blue)
    # lines!(Ph_SS_Ph_Kinetics_420nm[:,1],Ph_SS_Ph_Kinetics_420nm[:,3],linewidth = width,label = "PhSSPh Fit")

    lines!(Ph_SS_Ph_TBAHCO2_Kinetics[:,1],MaxNorm(Ph_SS_Ph_TBAHCO2_Kinetics[:,2]),linewidth = width,label = "PhSSPh + TBAHCO2", color = :darkorange)
    # lines!(Ph_SS_Ph_TBAHCO2_Kinetics_420[:,1],Ph_SS_Ph_TBAHCO2_Kinetics_420[:,3],linewidth = width,label = "PhSSPh + TBAHCO2 Fit")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    text!(5,0.7,text = "τ₁ = 494 ns, τ₂ = 3648 ns", color = :blue)
    text!(5,0.5,text = "τ₁ = 1155 ns, τ₂ = 4323 ns", color = :darkorange)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
KineticsComparison_420nm()
save("JULIA/Outputs/KineticsComparison_420nm.png",KineticsComparison_420nm())

function KineticsComparison_450nm(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "450nm Kinetics Comparison",
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.rainbow, 5),))
    
    # CairoMakie.ylims!(-0.1,nothing)

    width = 3
    
    ## Normalized Kinetics
    lines!(Ph_SS_Ph_Kinetics[:,1],MaxNorm(Ph_SS_Ph_Kinetics[:,3]),linewidth = width,label = "PhSSPh", color = :blue)
    lines!(Ph_SS_Ph_TBAHCO2_Kinetics[:,1],MaxNorm(Ph_SS_Ph_TBAHCO2_Kinetics[:,3]),linewidth = width,label = "PhSSPh + TBAHCO2", color = :orange)

    ### Fits 
    # Makie.scatter!(Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,3],linewidth = width,label = "PhSSPh")
    # lines!(Ph_SS_Ph_Kinetics_450nm[:,1],Ph_SS_Ph_Kinetics_450nm[:,3],linewidth = width,label = "PhSSPh Fit")

    # lines!(Ph_SS_Ph_TBAHCO2_Kinetics[:,1],Ph_SS_Ph_TBAHCO2_Kinetics[:,3],linewidth = width,label = "PhSSPh + TBAHCO2")
    # lines!(Ph_SS_Ph_TBAHCO2_Kinetics_450[:,1],Ph_SS_Ph_TBAHCO2_Kinetics_450[:,3],linewidth = width,label = "PhSSPh + TBAHCO2 Fit")

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    text!(5,0.7,text = "τ₁ = 485 ns, τ₂ = 3490 ns", color = :blue)
    text!(5,0.5,text = "τ₁ = 802 ns, τ₂ = 3138 ns", color = :orange)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
KineticsComparison_450nm()
save("JULIA/Outputs/KineticsComparison_450nm.png",KineticsComparison_450nm())

function KineticsComparison(f = Figure())

    ax = CairoMakie.Axis(f[1,1], title = "Kinetics Comparison in MeCN",
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0))
    
    CairoMakie.xlims!(-0.5,15)
    CairoMakie.ylims!(-0.05,nothing)

    width = 3
    
    ## Normalized Kinetics
    lines!(Ph_SS_Ph_Kinetics[:,1],MaxNorm(Ph_SS_Ph_Kinetics[:,3]),linewidth = width,label = "PhSSPh 450nm", color = :blue)
    # lines!(Ph_SS_Ph_Kinetics[:,1],MaxNorm(Ph_SS_Ph_Kinetics[:,2]),linewidth = width,label = "PhSSPh 420nm", color = :darkblue)
    # lines!(Ph_SS_Ph_Kinetics[:,1],MaxNorm(Ph_SS_Ph_Kinetics[:,4]),linewidth = width,label = "PhSSPh 470nm", color = :darkblue)


    lines!(Ph_SS_Ph_TBAHCO2_Kinetics[:,1],MaxNorm(Ph_SS_Ph_TBAHCO2_Kinetics[:,3]),linewidth = width,label = "PhSSPh + TBAHCO2 450nm", color = :red)
    # lines!(Ph_SS_Ph_TBAHCO2_Kinetics[:,1],MaxNorm(Ph_SS_Ph_TBAHCO2_Kinetics[:,2]),linewidth = width,label = "PhSSPh + TBAHCO2 420nm", color = :darkorange)
    # lines!(Ph_SS_Ph_TBAHCO2_Kinetics[:,1],MaxNorm(Ph_SS_Ph_TBAHCO2_Kinetics[:,4]),linewidth = width,label = "PhSSPh + TBAHCO2 550nm", color = :green)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    text!(5,0.7,text = "τ1 = 485 ns τ2 = 3490 ns", color = :blue)
    text!(5,0.5,text = "τ1 = 802 ns τ2 = 3138 ns", color = :red)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
end
KineticsComparison()
save("JULIA/Outputs/KineticsComparison.png",KineticsComparison())

function FTIR_Photolysis(f = Figure(fontsize = 15))

    ax1 = CairoMakie.Axis(f[1,1],
     xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,xticks = 1600:25:2000,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 3),))
        CairoMakie.xlims!(ax1,1590,1699)
        # CairoMakie.ylims!(ax1,-0.025,0.6)

        width = 4

        # lines!(Ph_SS_Ph_Photolysis[:,1],Ph_SS_Ph_Photolysis[:,2],linewidth = width,label = "pre")
        # lines!(Ph_SS_Ph_Photolysis[:,1],Ph_SS_Ph_Photolysis[:,3],linewidth = width,label = "post")
        lines!(Ph_SS_Ph_Photolysis[:,1],Ph_SS_Ph_Photolysis[:,4],linewidth = width,label = "Ph-SS-Ph", color = :blue)

        # lines!(Ph_SS_Ph_TBAHCO2_Photolysis[:,1],Ph_SS_Ph_TBAHCO2_Photolysis[:,2],linewidth = width,label = "pre")
        # lines!(Ph_SS_Ph_TBAHCO2_Photolysis[:,1],Ph_SS_Ph_TBAHCO2_Photolysis[:,3],linewidth = width,label = "post")
        lines!(Ph_SS_Ph_TBAHCO2_Photolysis[:,1],Ph_SS_Ph_TBAHCO2_Photolysis[:,4],linewidth = width,label = "Ph-SS-Ph + TBAHCO₂", color = :red)
        axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    ax2 = CairoMakie.Axis(f[1,2],
     xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,xticks = 2300:25:2400,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0),palette = (color = palette(ColorSchemes.hsv, 3),))
        hideydecorations!(ax2)
        CairoMakie.xlims!(ax2,2311,2362)
        # CairoMakie.ylims!(ax2,-0.025,0.6)

        # lines!(ax2,Ph_SS_Ph_Photolysis[:,1],Ph_SS_Ph_Photolysis[:,2],linewidth = width,label = "pre")
        # lines!(ax2,Ph_SS_Ph_Photolysis[:,1],Ph_SS_Ph_Photolysis[:,3],linewidth = width,label = "post")
        lines!(Ph_SS_Ph_Photolysis[:,1],Ph_SS_Ph_Photolysis[:,4],linewidth = width,label = "diff", color = :blue)

        # lines!(ax2,Ph_SS_Ph_TBAHCO2_Photolysis[:,1],Ph_SS_Ph_TBAHCO2_Photolysis[:,2],linewidth = width,label = "pre")
        # lines!(ax2,Ph_SS_Ph_TBAHCO2_Photolysis[:,1],Ph_SS_Ph_TBAHCO2_Photolysis[:,3],linewidth = width,label = "post")
        lines!(Ph_SS_Ph_TBAHCO2_Photolysis[:,1],Ph_SS_Ph_TBAHCO2_Photolysis[:,4],linewidth = width,label = "diff", color = :red)


    Label(f[2,1:2],"Wavenumber (cm⁻¹)", fontsize = 20)
    rowsize!(f.layout,2,20) 
    
    colgap!(f.layout,1,10)
    colsize!(f.layout,1,Relative(3/4)) #Set aspect ration of the y vs. x axis)
    # colsize!(f.layout,2,200) #Set aspect ration of the y vs. x axis)
   

    resize_to_layout!(f)
    f
end

FTIR_Photolysis()
save("JULIA/Outputs/FTIR_Photolysis_MeCN.png",FTIR_Photolysis())


### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = "/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/20mM_PhSSPh_MeCN_ArPurge_420nm_TAKinetics_500avg.csv"

    ## FITTING FUNCTIONS

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
    


# SOME FUNCTIONS

    #Normalization to x=z
    function Norm(x,y,z)
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/x_interp(z) #normalize to desired value
        return y_norm
    end
    # Normalize to Maximum

    function MaxNorm(y)

        y_norm = y/maximum(y) #normalize to desired value
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

     function Resize(Data,column, length)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(Data[1,1],Data[end,1],length))
        Interp  = LinearInterpolation(Data[:,1],Data[:,column])
        y_out = Interp(x_out)
    
        df = DataFrame(x = x_out,y = y_out)
        return df
    end

    MonoFit(Ph_SS_Ph_Kinetics,2,0.1,15.0,true,false)
    BiFit(Ph_SS_Ph_Kinetics,2,0.12,15.0,true,false)
    BiFit(Ph_SS_Ph_TBAHCO2_Kinetics,2,0.12,15.0,true,false)