using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

Scheme = ColorSchemes.hsv
## IMPORT DATA


## FTIR Photolysis
PhSSPh_20mM_FTIR_Photolysis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/FTIR_Photolysis/PhSSPh_20mM_DMSO_Ar_Purged_390um_FTIR_Photolysis.csv", header = true, DataFrame)
PhSSPh_20mM_KHCO2_10mM_FTIR_Photolysis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/FTIR_Photolysis/PhSSPh_20mM_KHCO2_10mM_DMSO_Ar_Purged_390um_FTIR_Photolysis.csv", header = true, DataFrame)

## NS TA Spectra
PhSSPh_20mM_TA_Spectra= CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/PhSSPh_20mM_DMSO_Ar_Purged_TA_Spectra.csv", header = true, DataFrame)
# Plots.plot(PhSSPh_20mM_TA_Spectra[:,1],PhSSPh_20mM_TA_Spectra[:,7])

PhSSPh_20mM_KHCO2_10mM_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/PhSSPh_20mM_KHCO2_10mM_DMSO_Ar_Purged_TA_Spectra_2mlpmin.csv", header = true, DataFrame)
# Plots.plot(PhSSPh_20mM_KHCO2_10mM_TA_Spectra[:,1],PhSSPh_20mM_KHCO2_10mM_TA_Spectra[:,6])

## NS TA Kinetics
PhSSPh_20mM_TA_Kinetics_2mlpmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/Kinetics/PhSSPh_20mM_DMSO_Ar_Purged_TAKinetics_flow_2mlpmin.csv", header = true, DataFrame)
# Plots.plot(PhSSPh_20mM_TA_Kinetics_2mlpmin[:,1], PhSSPh_20mM_TA_Kinetics_2mlpmin[:,4])
PhSSPh_20mM_TA_Kinetics_2mlpmin_470nm = Resize(PhSSPh_20mM_TA_Kinetics_2mlpmin,4,1000)

PhSSPh_20mM_TA_Kinetics_2mlpmin_470nm = BiFit_new(Resize(PhSSPh_20mM_TA_Kinetics_2mlpmin,4,1000),2,0.13,15,true,false,[0.01, 0.4, 0.02, 3.0, 0.0])# t1 = 504 ns t2 = 3627 ns
PhSSPh_MonoFitIRF = MonoFitIRF(PhSSPh_20mM_TA_Kinetics_2mlpmin_470nm, 2, 0.06, 15, [0.0301, 0.071, 0.045, 1.279, 4.138, 0.550], true, true, 1)
# PhSSPh_BiFitIRF = BiFitIRF(PhSSPh_20mM_TA_Kinetics_2mlpmin_470nm,2,0.08,15,[0.0326, 0.11, 0.2414, 1.0193, 6.6790, 0.8262, 1.1017, 0.4393],true, false,2)


PhSSPh_20mM_TA_Kinetics_3mlpmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/Kinetics/PhSSPh_20mM_DMSO_Ar_Purged_TAKinetics_flow_3mlpmin.csv", header = true, DataFrame)

PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/Kinetics/PhSSPh_20mM_KHCO2_10mM_DMSO_Ar_Purged_TAKinetics_flow_2mlpmin.csv", header = true, DataFrame)
# Plots.plot(PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin[:,1],PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin[:,4])
PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_470 = Resize(PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin,4,1000)
PhSSPh_HCO2_BiFitIRF = BiFitIRF(PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_470,2,0.08,15,[0.098, -15.666, 0.990, 3.482, 11.251, 3.581, 2.295, 0.594],true, true,1)


PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_3mlpmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/Kinetics/PhSSPh_20mM_KHCO2_10mM_DMSO_Ar_Purged_TAKinetics_flow_3mlpmin.csv", header = true, DataFrame)

PhSSPh_20mM_KHCO2_10mM_420nm_TA_Kinetics_Flowrates = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_01_FTIR_and_NS_TA_PhSSPh_KHCO2_Flow_Ar_purged/NS_TA/Kinetics/PhSSPh_20mM_KHCO2_10mM_DMSO_Ar_Purged_420nm_TAKinetics_flowrates.csv", header = true, DataFrame)


### Make Figures 

function FTIR_Photolysis(f = Figure(figure_padding = 25, fontsize = 20))
    Data1 = PhSSPh_20mM_FTIR_Photolysis
    Data2 = PhSSPh_20mM_KHCO2_10mM_FTIR_Photolysis

    Scheme = ColorSchemes.winter

    #Make Cairo Mackie figure
    
    ax1 = CairoMakie.Axis(f[1,1],  titlesize = 20, palette = (color = palette(Scheme, 4),),
    xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(1550,1651)
    # CairoMakie.ylims!(-0.1,0.8)
    width = 4
    x1 = Data1[:,1]
    x2 = Data2[:,1]

    # CairoMakie.lines!(ax1,x1,Data1[:,2],linewidth = width,label = "0 sec")
    # CairoMakie.lines!(ax1,x1,Data1[:,3],linewidth = width,label = "60 sec")

    CairoMakie.lines!(ax1,x1,Data1[:,3] .- Data1[:,2],linewidth = width,label = "delta 60 sec")


    ## Second Axis
    ax2 = CairoMakie.Axis(f[1,2],  titlesize = 20, palette = (color = palette(Scheme, 4),),
     xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    hideydecorations!(ax2) # hides labels and tick marks

    CairoMakie.xlims!(2311,2362)
    # CairoMakie.ylims!(-0.1,0.8)

    # CairoMakie.lines!(ax2,x1,Data1[:,2],linewidth = width,label = "0 sec")
    # CairoMakie.lines!(ax2,x1,Data1[:,3],linewidth = width,label = "60 sec")

    CairoMakie.lines!(ax2,x1,Data1[:,3].- Data1[:,2],linewidth = width,label = "60 sec")


    axislegend("Ph-SS-Ph",position = :ct,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)


    
    ax3 = CairoMakie.Axis(f[3,1],  titlesize = 20, palette = (color = palette(Scheme, 4),),
    xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(1550,1651)
    # CairoMakie.ylims!(-0.1,0.8)
    # CairoMakie.lines!(ax3,x2,Data2[:,2],linewidth = width,label = "0 sec")
    # CairoMakie.lines!(ax3,x2,Data2[:,3],linewidth = width,label = "30 sec")
    # CairoMakie.lines!(ax3,x2,Data2[:,4],linewidth = width,label = "60 sec")
    # CairoMakie.lines!(ax3,x2,Data2[:,5],linewidth = width,label = "120 sec")

    CairoMakie.lines!(ax3,x2,Data2[:,2].- Data2[:,2],linewidth = width,label = "0 sec")
    CairoMakie.lines!(ax3,x2,Data2[:,3].- Data2[:,2],linewidth = width,label = "30 sec")
    CairoMakie.lines!(ax3,x2,Data2[:,4].- Data2[:,2],linewidth = width,label = "60 sec")
    CairoMakie.lines!(ax3,x2,Data2[:,5].- Data2[:,2],linewidth = width,label = "120 sec")

    ## Second Axis
    ax4 = CairoMakie.Axis(f[3,2],  titlesize = 20, palette = (color = palette(Scheme, 4),),
     xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    hideydecorations!(ax4) # hides labels and tick marks

    CairoMakie.xlims!(2311,2362)
    # CairoMakie.ylims!(-0.1,0.8)

    # CairoMakie.lines!(ax4,x2,Data2[:,2],linewidth = width,label = "0 sec")
    # CairoMakie.lines!(ax4,x2,Data2[:,3],linewidth = width,label = "30 sec")
    # CairoMakie.lines!(ax4,x2,Data2[:,4],linewidth = width,label = "60 sec")
    # CairoMakie.lines!(ax4,x2,Data2[:,5],linewidth = width,label = "120 sec")

    CairoMakie.lines!(ax4,x2,Data2[:,2].- Data2[:,2],linewidth = width,label = "0 sec")
    CairoMakie.lines!(ax4,x2,Data2[:,3].- Data2[:,2],linewidth = width,label = "30 sec")
    CairoMakie.lines!(ax4,x2,Data2[:,4].- Data2[:,2],linewidth = width,label = "60 sec")
    CairoMakie.lines!(ax4,x2,Data2[:,5].- Data2[:,2],linewidth = width,label = "120 sec")


    axislegend("Ph-SS-Ph+KHCO₂",position = :cb,nbanks = 2,framevisible = false)

    Label(f[2,(1:2)],"Wavenumber (cm⁻¹)",halign = :center,valign = :top)
    Label(f[4,(1:2)],"Wavenumber (cm⁻¹)",halign = :center,valign = :top)
    # colsize!(f.layout,1,Aspect(1, 1))
    # colsize!(f.layout,2,Aspect(1, 1)) #Set aspect ration of the y vs. x axis)
    colgap!(f.layout,1,5)
    rowgap!(f.layout,1,5)
    rowgap!(f.layout,2,20)
    rowgap!(f.layout,3,5)
    rowsize!(f.layout,1,300)
    rowsize!(f.layout,3,300)
    resize_to_layout!(f)

    f
end
FTIR_Photolysis()
save("JULIA/Outputs/FTIR_Photolysis.eps",FTIR_Photolysis())

function PhSSPh_TASpectra(f = Figure(figure_padding = 25, fontsize = 20, font = :bold))
    Data1 = PhSSPh_20mM_TA_Spectra

    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Spectra of Ph-SS-Ph (20mM) in DMSO", titlesize = 20, palette = (color = palette(Scheme, 5),),
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(375,600)
    # CairoMakie.ylims!(-0.05,nothing)
    width = 4
    x1 = Data1[:,1]

    # lines!(x1,Data1[:,2],linewidth = width,label = "20 ns")
    # lines!(x1,Data1[:,3],linewidth = width,label = "100 ns")
    # lines!(x1,Data1[:,4],linewidth = width,label = "500 ns")
    # lines!(x1,Data1[:,5],linewidth = width,label = "1000 ns")
    # lines!(x1,Data1[:,6],linewidth = width,label = "2000 ns")
    # lines!(x1,Data1[:,7],linewidth = width,label = "10000 ns")

    # lines!(x1,Data1[:,3] * 1.686,linewidth = width,label = "100 ns")
    lines!(x1,Data1[:,4] * 4.714,linewidth = width,label = "500 ns")
    # lines!(x1,Data1[:,5]* 1.308,linewidth = width,label = "1000 ns")
    lines!(x1,Data1[:,6] * 2.909,linewidth = width,label = "2000 ns")
    lines!(x1,Data1[:,7]* 4.034,linewidth = width,label = "10000 ns")

    axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    f
end
PhSSPh_TASpectra()
save("JULIA/Outputs/PhSSPh_TASpectra.png",PhSSPh_TASpectra())


function PhSSPh_KHCO2_TASpectra(f = Figure(figure_padding = 25, fontsize = 20, font = :bold))
    Data1 = PhSSPh_20mM_KHCO2_10mM_TA_Spectra

    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Spectra of Ph-SS-Ph (20 mM) + KHCO₂ (10 mM) in DMSO", titlesize = 20,  palette = (color = palette(Scheme, 5),),
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(375,600)
    CairoMakie.ylims!(-0.001,nothing)
    width = 4
    x1 = Data1[:,1]

    # lines!(x1,Data1[:,2],linewidth = width,label = "100 ns")
    # lines!(x1,Data1[:,3],linewidth = width,label = "500 ns")
    # lines!(x1,Data1[:,4],linewidth = width,label = "1000 ns")
    # lines!(x1,Data1[:,5],linewidth = width,label = "2000 ns")
    # lines!(x1,Data1[:,6],linewidth = width,label = "10000 ns")

    # lines!(x1,Data1[:,2]*5.535,linewidth = width,label = "100 ns")
    lines!(x1,Data1[:,3]*4.106,linewidth = width,label = "500 ns")
    # lines!(x1,Data1[:,4]*2.199,linewidth = width,label = "1000 ns")
    lines!(x1,Data1[:,5]*3.440,linewidth = width,label = "2000 ns")
    lines!(x1,Data1[:,6]*1.736,linewidth = width,label = "10000 ns")

    axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    f
end
PhSSPh_KHCO2_TASpectra()
save("JULIA/Outputs/PhSSPh_KHCO2_TASpectra.png",PhSSPh_KHCO2_TASpectra())

function PhSSPh_Kinetics(f = Figure(figure_padding = 25, fontsize = 20, font = :bold))
    Data1 = PhSSPh_20mM_TA_Kinetics_2mlpmin
    Data2 = PhSSPh_20mM_TA_Kinetics_3mlpmin

    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Kinetics of Ph-SS-Ph (20 mM) in DMSO", titlesize = 20, xautolimitmargin = (0.01, 0.0),
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # CairoMakie.xlims!(0,15)
    # CairoMakie.ylims!(-0.005,nothing)
    width = 3
    x1 = Data1[:,1]
    x2 = Data2[:,1]

    # lines!(x1,MaxNorm(x1,Data1[:,2]),linewidth = width,label = "380 nm")
    # lines!(x1,MaxNorm(x1,Data1[:,3]),linewidth = width,label = "420 nm")
    # lines!(x1,MaxNorm(x1,Data1[:,4]),linewidth = width,label = "470 nm")

    lines!(PhSSPh_MonoFitIRF[:,1],PhSSPh_MonoFitIRF[:,2],linewidth = width,label = "470 nm", color=:blue)
    lines!(PhSSPh_MonoFitIRF[:,1],PhSSPh_MonoFitIRF[:,3],linewidth = 2,label = "fit",color=:red)



    # lines!(x2,MaxNorm(x2,Data2[:,2]),linewidth = width,label = "380 nm 3")
    # lines!(x2,MaxNorm(x2,Data2[:,3]),linewidth = width,label = "420 nm 3")
    # lines!(x2,MaxNorm(x2,Data2[:,4]),linewidth = width,label = "470 nm 3")

    axislegend(position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    f
end
PhSSPh_Kinetics()
save("JULIA/Outputs/PhSSPh_Kinetics.png",PhSSPh_Kinetics())


function PhSSPh_KHCO2_Kinetics(f = Figure(figure_padding = 25, fontsize = 20, font = :bold))
    Data1 = PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin

    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Kinetics of Ph-SS-Ph (20 mM) + HCO2 (10 mM) in DMSO", titlesize = 20, xautolimitmargin = (0.01,0),
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # CairoMakie.xlims!(0,5)
    # CairoMakie.ylims!(-0.005,nothing)
    width = 3
    x1 = Data1[:,1]

    # lines!(x1,MaxNorm(x1,Data1[:,2]),linewidth = width,label = "380 nm")
    # lines!(x1,Norm(x1,Data1[:,3], 0.1),linewidth = width,label = "420 nm")
    # lines!(x1,Norm(x1,Data1[:,4],0.1),linewidth = width,label = "470 nm")

    lines!(PhSSPh_HCO2_BiFitIRF[:,1],PhSSPh_HCO2_BiFitIRF[:,2],linewidth = width, label = "470 nm", color = :blue)
    lines!(PhSSPh_HCO2_BiFitIRF[:,1],PhSSPh_HCO2_BiFitIRF[:,3],linewidth = 2, label = "fit",color = :red)

    axislegend(position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    f
end
PhSSPh_KHCO2_Kinetics()
save("JULIA/Outputs/PhSSPh_KHCO2_Kinetics.png",PhSSPh_KHCO2_Kinetics())


function PhSSPh_KHCO2_420nm_FlowrateKinetics(f = Figure(figure_padding = 25, fontsize = 20, font = :bold))
    Data1 = PhSSPh_20mM_KHCO2_10mM_420nm_TA_Kinetics_Flowrates

    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Kinetics of Ph-SS-Ph (20 mM) in DMSO", titlesize = 20,
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(0,15)
    CairoMakie.ylims!(-0.005,nothing)
    width = 3
    x1 = Data1[:,1]

    lines!(x1,Data1[:,2],linewidth = width,label = "2")
    lines!(x1,Data1[:,3],linewidth = width,label = "3")
    lines!(x1,Data1[:,4],linewidth = width,label = "4")
    lines!(x1,Data1[:,5],linewidth = width,label = "6")
    lines!(x1,Data1[:,6],linewidth = width,label = "10")



    # lines!(x1,MaxNorm(x1,Data1[:,2]),linewidth = width,label = "2")
    # lines!(x1,MaxNorm(x1,Data1[:,3]),linewidth = width,label = "3")
    # lines!(x1,MaxNorm(x1,Data1[:,4]),linewidth = width,label = "4")
    # lines!(x1,MaxNorm(x1,Data1[:,5]),linewidth = width,label = "6")
    # lines!(x1,MaxNorm(x1,Data1[:,6]),linewidth = width,label = "10")


    axislegend("Wavelength", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    f
end
PhSSPh_KHCO2_420nm_FlowrateKinetics()
 save("JULIA/Outputs/PhSSPh_KHCO2_420nm_FlowrateKinetics.png",PhSSPh_KHCO2_420nm_FlowrateKinetics())

function Kinetics_Comparison(f = Figure(figure_padding = 25, fontsize = 20, font = :bold))
    Data1 = PhSSPh_20mM_TA_Kinetics_2mlpmin_470nm_MonoFit = MonoFit(PhSSPh_20mM_TA_Kinetics_2mlpmin,4,0.5,15,false,true)
    Data2 = PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin_470nm_BiexpFit = BiFit(PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin,4,0.5,15,false,true)
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Kinetics Comparison in DMSO", titlesize = 20,
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true, xscale = identity,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(nothing,15)
    CairoMakie.ylims!(-0.005,nothing)
    width = 4
    x1 = Data1[:,1]
    x2 = Data2[:,1]
    # lines!(x1,MaxNorm(x1,Data1[:,2]),linewidth = width,label = "380 nm")
    lines!(x1,Norm(Data1[:,1],Data1[:,2],0.51),linewidth = width,label = "470 nm", color = :blue)
    # lines!(x1,Data1[:,3],linewidth = width, color =:black)

    lines!(x2,Norm(Data2[:,1],Data2[:,2],0.51),linewidth = width,label = "470 nm + KHCO2", color = :red)
    # lines!(x2,Norm(Data2[:,1],Data2[:,3],0.61),linewidth = width, color = :black)


    text!(7.5, 0.65, text = "τ₁ = 4.3 μs", color = :blue)
    text!(2.5, 0.2, text = "τ₁ = 1.4 μs, τ₂ = 6.6 μs", color = :red)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)

    f
end
Kinetics_Comparison()
save("JULIA/Outputs/Kinetics_Comparison.png",Kinetics_Comparison())


PhSSPh_20mM_TA_Kinetics_2mlpmin_470nm_MonoFit = MonoFit(PhSSPh_20mM_TA_Kinetics_2mlpmin,4,0.5,15,false,true)
PhSSPh_20mM_TA_Kinetics_2mlpmin_420nm_MonoFit = MonoFit(PhSSPh_20mM_TA_Kinetics_2mlpmin,3,0.5,15,false,false)

PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin_470nm_BiexpFit = BiFit(PhSSPh_20mM_KHCO2_10mM_TA_Kinetics_2mlpmin,4,0.5,15,false,false)

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

    function MonoFit(file,column,t0,uppertimebound, residuals, ReturnFit)

        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
        p0 = [0.01, 11.0, 0.0]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values))
        
    
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
    
            lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)
    
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
    
    function BiFit(file,column,t0,uppertimebound, residuals, ReturnFit)

        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
        p0 = [0.01, 9.0, 0.02, 15.0, 0.0]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values))
        
    
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
    
            lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)
    
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
    
    function TriFit(file,column,t0,uppertimebound, residuals, ReturnFit)
    
        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
        p0 = [0.01, 2.19, 0.01, 8.7, 0.01, 19.0, 0.01]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,fitfunction(values))
        
    
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
    
            lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)
    
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
    
    function BiFit_new(file,column,t0,uppertimebound, residuals, ReturnFit, coeffs::Vector{Float64})

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
        InterpData = LinearInterpolation(Data[:,1],Data[:,column])
        
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
            return df
        
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

    function MonoFitIRF(file,column,t0,uppertimebound,coeffs, residuals, ReturnFit, reduce)

        function Fitfunc(x,coeffs)
        
            Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
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
        InterpData = LinearInterpolation(Data[:,1],Data[:,column])
        
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
            return df
        
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