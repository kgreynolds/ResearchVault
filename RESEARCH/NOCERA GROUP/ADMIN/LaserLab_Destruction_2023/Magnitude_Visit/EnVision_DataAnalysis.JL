using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,Statistics,MAT ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

matfile = matopen("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/Sample Data Structure copy/ADMIN/LaserLab_Destruction_2023/Magnitude_Visit/mAGNITUDE/PhotoB2 only emission 575 nm.mat")
read(matfile,"data")

Sample_UvVis_Cary = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/LASERS/LaserLab_Destruction_2023/Magnitude_Visit/Magnitude_Visit_UvVis.csv", header = true, DataFrame)


IrK_Emission = CSV.read("/Users/kris/Desktop/mAGNITUDE/Ir(III)K emission.csv", header = false, DataFrame)

Photobeta_alpha_575nm_emission = CSV.read("/Users/kris/Desktop/mAGNITUDE/PhotoB2+wt alpha emission 575 nm.csv", header = true,DataFrame)
Photobeta_alpha_fit = MonoFitIRF(Photobeta_alpha_575nm_emission,2,20,1500,[1.857, -14.033, 2.819, 52.25, 307.286, 5.3],true,false,1) ## τ = 186 ns
Photobeta_alpha_fit = BiFitIRF(Photobeta_alpha_575nm_emission,2,20,1500,[1.857, -14.033, 2.819, 52.25, 307.286, 50,100, 5.3],true,false,1) ## τ = 186 ns


Photobeta_575nm_emission = CSV.read("/Users/kris/Desktop/mAGNITUDE/PhotoB2 only emission 575 nm.csv", header = true, DataFrame)
Photobeta_fit = MonoFitIRF(Photobeta_575nm_emission,2,-100,1500,[1.857, -14.033, 2.819, 52.25, 307.286, 5.3],true,false,1) ## τ = 186 ns
Photobeta_monofit = MonoFit(Photobeta_575nm_emission,2,70,1500,true, false,[1.25, 337.286, 0.06])
Photobeta_Bifit = BiFitIRF(Photobeta_575nm_emission,2,20,1500,[1.857, -14.033, 2.819, 52.25, 307.286, 50,100, 5.3],true,false,1) ## τ = 186 ns


PhSSPh_455nm_kinetics = CSV.read("/Users/kris/Desktop/mAGNITUDE/PhSSPh 455nm kinetics.csv", header = true, DataFrame)
PhSSPh_TA_Spectra= CSV.read("/Users/kris/Desktop/mAGNITUDE/Ph-SS-Ph NS TA Spectra.csv", header = false, DataFrame)

PtTFP_665nm_PL_Airfree = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 665nm PL Kinetics 262us airfree FPT.csv", header = true, DataFrame)
PtTFP_740nm_PL_Airfree = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 740nm PL Kinetics 262us airfree FPT.csv", header = true, DataFrame)
PtTFP_665nm_PLKinetics_Air = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 665nm PL Kinetics 65us opened to air.csv", header = true, DataFrame)
PtTFP_665nm_PL_Airfree_fit = MonoFitIRF(PtTFP_665nm_PL_Airfree,2,10,40000,[0.5436, -12.357, 4.165, -1.257, 456.4331, 25.8714],true,false,10) #τ₁ = 1.195 μs


PtTFP_740nm_PLKinetics_Air = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 740nm PL Kinetics 65us opened to air.csv", header =true, DataFrame)
PtTFP_N2_TA_spectra = Matrix(CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP N2 bubbled TA Spectrum 420-600nm, 5nm, 2nm step- photodegraded.csv", header = false, DataFrame))
PtTFP_N2_TA_spectra_T = smoothMat(transpose(PtTFP_N2_TA_spectra),3)
PtTFP_air_TA_Spectra = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP air TA Spectrum 420-600 nm, 5 nm resolution, 35 mins.csv", header = false,DataFrame)
PtTFP_air_TA_Spectra_Kinetics = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP air TA Spectrum 420-600 nm, 5 nm resolution, 35 mins.csv", header = true,DataFrame)

PtTFP_air_455nmTA_Kinetics = MonoFitIRF(PtTFP_air_TA_Spectra_Kinetics,9,-100,3400,[0.029, -12.246, 4.4328, 67.434, 1195.628, -0.7811],true,true,1) #τ₁ = 1.184 μs
PtTFP_air_515nmTA_Kinetics = MonoFitIRF(PtTFP_air_TA_Spectra_Kinetics,21,-100,3400,[-0.006, -12.469, 3.740, 67.442, 1164.530, -0.241],true,true,1) #τ₁ = 1.157 μs
PtTFP_air_560nmTA_Kinetics = MonoFitIRF(PtTFP_air_TA_Spectra_Kinetics,32,-100,3400,[0.004, -12.398, 4.178, 67.538, 1135.799, 0.1697],true,true,1) #τ₁ = 1.132 μs


### Make Figures 

function UVvis(f = Figure())
    Data = CSV.read("RESEARCH/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/Paper/Data/IrB_C12Bicarb_TRIR.csv",DataFrame)
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], 
            palette = (color = palette(ColorSchemes.hsv, 6),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "A"
        )
        ax.yticks= 0:0.5:10
        ax.xticks= 200:50:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(Sample_UvVis_Cary[:,1],Sample_UvVis_Cary[:,2],linewidth = width,label = "Ir(III)K 75 μM MeCN")
        # lines!(Sample_UvVis_Cary[:,1],Sample_UvVis_Cary[:,3],linewidth = width,label = "Ph-SS-Ph 20 mM MeCN")
        # lines!(Sample_UvVis_Cary[:,1],Sample_UvVis_Cary[:,4],linewidth = width,label = "Ru(bpy)₃ MeCN")
        # lines!(Sample_UvVis_Cary[:,1],Sample_UvVis_Cary[:,5],linewidth = width,label = "Ru(bpy)₃ 2mm MeCN")
        # lines!(Sample_UvVis_Cary[:,1],Sample_UvVis_Cary[:,6],linewidth = width,label = "Zn-TPP Concentrated MeCN")
        lines!(Sample_UvVis_Cary[:,1],Sample_UvVis_Cary[:,7],linewidth = width,label = "Zn-TPP MeCN")

        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Sample", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
UVvis()
save("JULIA/Outputs/Zn-TPP MeCN.png",UVvis())

function PhotoRNR(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Photo RNR 575nm emission kinetics",
            palette = (color = palette(ColorSchemes.hsv, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "A.U."
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 200:50:800

        # CairoMakie.xlims!(nothing,10)
        # CairoMakie.ylims!(nothing,0.2)
        
        lines!(Photobeta_575nm_emission[:,1],MaxNorm(Photobeta_575nm_emission[:,2]),linewidth = width,label = "Photo-β")
        lines!(Photobeta_alpha_575nm_emission[:,1],MaxNorm(Photobeta_alpha_575nm_emission[:,2]),linewidth = width,label = "Photo-β + α")


        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Sample", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PhotoRNR()
save("JULIA/Outputs/PhotoRNR.png",PhotoRNR())

function PtTFP_PLKinetics(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pt-TFP emission kinetics",
            palette = (color = palette(ColorSchemes.hsv, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (μs)", 
            ylabel = "Normalized Intensity"
        )
        # ax.yticks= 0:0.5:10
        ax.xticks= 0:50:800

        # CairoMakie.xlims!(nothing,10)
        # CairoMakie.ylims!(nothing,0.2)
        
        lines!(PtTFP_665nm_PL_Airfree[230:end,1]/1000,MaxNorm(PtTFP_665nm_PL_Airfree[230:end,2]),linewidth = width,label = "PtTFP 665nm Airfree")
        # lines!(PtTFP_740nm_PL_Airfree[:,1]/1000,MaxNorm(PtTFP_740nm_PL_Airfree[:,2]),linewidth = width,label = "PtTFP 740nm Airfree")

        # lines!(PtTFP_665nm_PLKinetics_Air[700:end,1]/1000,MaxNorm(PtTFP_665nm_PLKinetics_Air[700:end,2]),linewidth = width,label = "PtTFP 665nm Air")
        # lines!(PtTFP_740nm_PLKinetics_Air[:,1]/1000,MaxNorm(PtTFP_740nm_PLKinetics_Air[:,2]),linewidth = width,label = "PtTFP 740nm Air")

        # lines!(PtTFP_air_TA_Spectra[2:end,1]/1000,MaxNorm(PtTFP_air_TA_Spectra[2:end,8]),linewidth = width,label = "PtTFP 665nm Air")


        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend(position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_PLKinetics()
save("JULIA/Outputs/PtTFP_PLKinetics.png",PtTFP_PLKinetics())

function PtTFP_air_TA_spectra_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pt-TFP TA Spectra opened to air 3pt smooth",
            palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 225:25:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(PtTFP_N2_TA_spectra[:,1],PtTFP_N2_TA_spectra[:,10],linewidth = width,label = "Zn-TPP MeCN")

        # lines!(collect(PtTFP_air_TA_Spectra[1,2:end]),collect(PtTFP_air_TA_Spectra[352,2:end]),linewidth = width,label = "100 ns")
        # lines!(collect(PtTFP_air_TA_Spectra[1,2:end]),collect(PtTFP_air_TA_Spectra[402,2:end]),linewidth = width,label = "200 ns")
        # lines!(collect(PtTFP_air_TA_Spectra[1,2:end]),collect(PtTFP_air_TA_Spectra[552,2:end]),linewidth = width,label = "500 ns")
        # lines!(collect(PtTFP_air_TA_Spectra[1,2:end]),collect(PtTFP_air_TA_Spectra[802,2:end]),linewidth = width,label = "1000 ns")
        # lines!(collect(PtTFP_air_TA_Spectra[1,2:end]),collect(PtTFP_air_TA_Spectra[1302,2:end]),linewidth = width,label = "2000 ns")

        lines!(smooth(collect(PtTFP_air_TA_Spectra[1,2:end]),3),smooth(collect(PtTFP_air_TA_Spectra[352,2:end]),3),linewidth = width,label = "100 ns")
        lines!(smooth(collect(PtTFP_air_TA_Spectra[1,2:end]),3),smooth(collect(PtTFP_air_TA_Spectra[402,2:end]),3),linewidth = width,label = "200 ns")
        lines!(smooth(collect(PtTFP_air_TA_Spectra[1,2:end]),3),smooth(collect(PtTFP_air_TA_Spectra[552,2:end]),3),linewidth = width,label = "500 ns")
        lines!(smooth(collect(PtTFP_air_TA_Spectra[1,2:end]),3),smooth(collect(PtTFP_air_TA_Spectra[802,2:end]),3),linewidth = width,label = "1000 ns")
        lines!(smooth(collect(PtTFP_air_TA_Spectra[1,2:end]),3),smooth(collect(PtTFP_air_TA_Spectra[1302,2:end]),3),linewidth = width,label = "2000 ns")

        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_air_TA_spectra_plot()
save("JULIA/Outputs/PtTFP_air_TA_spectra_plot_smoothed.png",PtTFP_air_TA_spectra_plot())

function PtTFP_air_TA_Kinetics_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pt-TFP TA Kinetics opened to air",
            palette = (color = palette(ColorSchemes.hsv, 6),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 225:25:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)


        # lines!(PtTFP_air_TA_Spectra[2:end,1],MaxNorm(PtTFP_air_TA_Spectra[2:end,9]),linewidth = width,label = "455 nm ")
        # lines!(PtTFP_air_TA_Spectra[2:end,1],MaxNorm(PtTFP_air_TA_Spectra[2:end,21]),linewidth = width,label = "515 nm ")
        # lines!(PtTFP_air_TA_Spectra[2:end,1],MaxNorm(PtTFP_air_TA_Spectra[2:end,32]),linewidth = width,label = "560 nm ")

        lines!(PtTFP_air_TA_Spectra[2:end,1],PtTFP_air_TA_Spectra[2:end,9],linewidth = width,label = "455 nm ")
        lines!(PtTFP_air_455nmTA_Kinetics[:,1],PtTFP_air_455nmTA_Kinetics[:,3],linewidth = 2,label = "455 nm fit, 1.19 μs ")


        lines!(PtTFP_air_TA_Spectra[2:end,1],PtTFP_air_TA_Spectra[2:end,21],linewidth = width,label = "515 nm ")
        lines!(PtTFP_air_515nmTA_Kinetics[:,1],PtTFP_air_515nmTA_Kinetics[:,3],linewidth = 2,label = "515 nm fit, 1.16 μs ")


        lines!(PtTFP_air_TA_Spectra[2:end,1],PtTFP_air_TA_Spectra[2:end,32],linewidth = width,label = "560 nm ")
        lines!(PtTFP_air_560nmTA_Kinetics[:,1],PtTFP_air_560nmTA_Kinetics[:,3],linewidth = 2,label = "560 nm fit, 1.13 μs ")


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_air_TA_Kinetics_plot()
save("JULIA/Outputs/PtTFP_air_TA_Kinetics_plot.png",PtTFP_air_TA_Kinetics_plot())


function PtTFP_N2_TA_spectra_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pt-TFP N2 purged TA Spectra",
            palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 200:50:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(PtTFP_N2_TA_spectra[:,1],PtTFP_N2_TA_spectra[:,10],linewidth = width,label = "Zn-TPP MeCN")

        lines!(PtTFP_N2_TA_spectra_T[2:end,1],PtTFP_N2_TA_spectra_T[2:end,215],linewidth = width,label = "100 ns")
        lines!(PtTFP_N2_TA_spectra_T[2:end,1],PtTFP_N2_TA_spectra_T[2:end,827],linewidth = width,label = "5000 ns")
        lines!(PtTFP_N2_TA_spectra_T[2:end,1],PtTFP_N2_TA_spectra_T[2:end,1452],linewidth = width,label = "10000 ns")
        lines!(PtTFP_N2_TA_spectra_T[2:end,1],PtTFP_N2_TA_spectra_T[2:end,5452],linewidth = width,label = "50000 ns")
        lines!(PtTFP_N2_TA_spectra_T[2:end,1],PtTFP_N2_TA_spectra_T[2:end,10202],linewidth = width,label = "80000 ns")

        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_N2_TA_spectra_plot()
save("JULIA/Outputs/PtTFP_N2_TA_spectra_plot.png",PtTFP_N2_TA_spectra_plot())

function IrK_emission_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ir(III)K Emission Spectra",
            palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "A.U."
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 200:50:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(PtTFP_N2_TA_spectra[:,1],PtTFP_N2_TA_spectra[:,10],linewidth = width,label = "Zn-TPP MeCN")

        lines!(collect(IrK_Emission[1,2:end]),collect(IrK_Emission[207,2:end]),linewidth = width,label = "10 ns")
        lines!(collect(IrK_Emission[1,2:end]),collect(IrK_Emission[227,2:end]),linewidth = width,label = "50 ns")
        lines!(collect(IrK_Emission[1,2:end]),collect(IrK_Emission[252,2:end]),linewidth = width,label = "100 ns")
        lines!(collect(IrK_Emission[1,2:end]),collect(IrK_Emission[302,2:end]),linewidth = width,label = "200 ns")
        lines!(collect(IrK_Emission[1,2:end]),collect(IrK_Emission[402,2:end]),linewidth = width,label = "400 ns")

        # lines!(smooth(collect(IrK_Emission[1,2:end]),2),smooth(collect(IrK_Emission[207,2:end]),2),linewidth = width,label = "10 ns")
        # lines!(smooth(collect(IrK_Emission[1,2:end]),2),smooth(collect(IrK_Emission[227,2:end]),2),linewidth = width,label = "50 ns")
        # lines!(smooth(collect(IrK_Emission[1,2:end]),2),smooth(collect(IrK_Emission[252,2:end]),2),linewidth = width,label = "100 ns")
        # lines!(smooth(collect(IrK_Emission[1,2:end]),2),smooth(collect(IrK_Emission[302,2:end]),2),linewidth = width,label = "200 ns")
        # lines!(smooth(collect(IrK_Emission[1,2:end]),2),smooth(collect(IrK_Emission[402,2:end]),2),linewidth = width,label = "400 ns")

        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
IrK_emission_plot()
save("JULIA/Outputs/IrK_emission_plot.png",IrK_emission_plot())

function IrK_emissionKinetics_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ir(III)K Emission Kinetics",
            palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "A.U."
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 200:50:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(PtTFP_N2_TA_spectra[:,1],PtTFP_N2_TA_spectra[:,10],linewidth = width,label = "Zn-TPP MeCN")

        lines!(collect(IrK_Emission[2:end,1]),MaxNorm(collect(IrK_Emission[2:end,5])),linewidth = width,label = "475 nm")

        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
IrK_emissionKinetics_plot()
save("JULIA/Outputs/IrK_emissionKinetics_plot.png",IrK_emissionKinetics_plot())

function PhSSPh_Spectra(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ph-SS-Ph TA Spectra 3pt smooth",
            palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 200:50:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(PtTFP_N2_TA_spectra[:,1],PtTFP_N2_TA_spectra[:,10],linewidth = width,label = "Zn-TPP MeCN")

        # lines!(collect(PhSSPh_TA_Spectra[1,2:end]),collect(PhSSPh_TA_Spectra[240,2:end]),linewidth = width,label = "100 ns")
        # lines!(collect(PhSSPh_TA_Spectra[1,2:end]),collect(PhSSPh_TA_Spectra[1477,2:end]),linewidth = width,label = "10000 ns")
        # lines!(collect(PhSSPh_TA_Spectra[1,2:end]),collect(PhSSPh_TA_Spectra[2727,2:end]),linewidth = width,label = "20000 ns")
        # lines!(collect(PhSSPh_TA_Spectra[1,2:end]),collect(PhSSPh_TA_Spectra[6477,2:end]),linewidth = width,label = "50000 ns")
        # lines!(collect(PhSSPh_TA_Spectra[1,2:end]),collect(PhSSPh_TA_Spectra[12727,2:end]),linewidth = width,label = "100000 ns")

        lines!(smooth(collect(PhSSPh_TA_Spectra[1,2:end]),3),smooth(collect(PhSSPh_TA_Spectra[240,2:end]),3),linewidth = width,label = "10 ns")
        lines!(smooth(collect(PhSSPh_TA_Spectra[1,2:end]),3),smooth(collect(PhSSPh_TA_Spectra[1477,2:end]),3),linewidth = width,label = "10000 ns")
        lines!(smooth(collect(PhSSPh_TA_Spectra[1,2:end]),3),smooth(collect(PhSSPh_TA_Spectra[2727,2:end]),3),linewidth = width,label = "20000 ns")
        lines!(smooth(collect(PhSSPh_TA_Spectra[1,2:end]),3),smooth(collect(PhSSPh_TA_Spectra[6477,2:end]),3),linewidth = width,label = "50000 ns")
        lines!(smooth(collect(PhSSPh_TA_Spectra[1,2:end]),3),smooth(collect(PhSSPh_TA_Spectra[12727,2:end]),3),linewidth = width,label = "100000 ns")

        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PhSSPh_Spectra()
save("JULIA/Outputs/PhSSPh_Spectra_smoothed.png",PhSSPh_Spectra())

function PhSSPh_TA_Kinetics(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ph-SS-Ph TA Kinetics",
            palette = (color = palette(ColorSchemes.hsv, 2),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (μs)", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 200:50:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        

        lines!(PhSSPh_TA_Spectra[2:end,1]/1000,PhSSPh_TA_Spectra[2:end,12],linewidth = width,label = "450 nm Kinetics")

        # lines!(PhSSPh_455nm_kinetics[:,1]/1000,PhSSPh_455nm_kinetics[:,2],linewidth = width,label = "455 nm Kinetics")


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PhSSPh_TA_Kinetics()
save("JULIA/Outputs/PhSSPh_TA_Kinetics_2.png",PhSSPh_TA_Kinetics())

### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_450nm_Kinetics_180us_COnverted.csv",header = true,DataFrame)

    ## FITTING FUNCTIONS

    function LinearFit(x,y,residuals, ReturnFit)
    
        time = x
        KineticTrace = y
    
        ## FITTING FUNCTION
        fitfunction = mexp(p) = p[1] .* x .+ p[2]
        p0 = [5188.825, 0.00768]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        residual = KineticTrace .- fitfunction(values)
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values),residual)
        
    
        elseif residuals == true
            Fit_Fig = Figure(figure_padding = 25,fontsize =20, resolution = (600,600))
            width = 2
    
            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
            xlabel = "x",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)
    
            lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
            lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
            lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
            axislegend(ax1, position = :rt, framevisible = false)
    
            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "x",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
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
        residual = KineticTrace .- fitfunction(values)
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
        """Normalize Data to a value Z """
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/x_interp(z) #normalize to desired value
        return y_norm
    end

    # Normalize to Maximum
    function MaxNorm(y)
        """Normalize vector to its absolute largest value"""
        y_norm = y/(maximum(abs.(y))) #normalize to desired value
        return y_norm
    end

    function Resize(Data,column, length)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(Data[1,1],Data[end,1],length))
        Interp  = LinearInterpolation(Data[:,1],Data[:,column])
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
    
    function SpectralCorrection(Wavelength,Spectrum,Time,Kinetic,CorrWL,CorrTime)
        """Correct a TA spectrum to match with a given kinetic trace by scaling the spectrum by a factor"""
        SpectrumInterp = LinearInterpolation(Wavelength,Spectrum)
        Spectrum_Y = SpectrumInterp(Wavelength)
    
        KineticInterp = LinearInterpolation(Time,Kinetic)
        
        CorrFactor = KineticInterp[CorrTime]/SpectrumInterp[CorrWL]
    
        CorrSpectrum = Spectrum_Y .* CorrFactor
    
        print("factor is:", CorrFactor)
        return CorrSpectrum
     end

    function ReadFilesinDir(directory::String)
        """ Read files in a directory, recurisley, of same file structure and read the data into one matrix"""
        files = []
        for (root, dirs, filenames) in walkdir(directory)
        for filename in filenames
            if isfile(joinpath(root, filename))
            push!(files, readdlm(joinpath(root, filename),skipstart = 60))
            end
        end
        end
    
        return files
    end

    function smoothMat(Mat,SmoothPts)
        """This function takes and input matrix Mar and generates a new matrix sMat with the moving average over
         some number of points defined by the SmoothPts argument"""
        sMat = zeros(size(Mat,1)-SmoothPts,size(Mat,2))
    
        sMat[1,:] = Mat[1,:]

        for i in 2:size(Mat,1)-SmoothPts
            for j in 1:size(Mat,2)
                sMat[i,j] = mean(skipmissing(Mat[i:i+SmoothPts,j])) # This works! calculates the mean of the i to i+2th row in the jth columns and replaces the E[i,j]th value with that mean
            end
        end
        return sMat
    
    end 
