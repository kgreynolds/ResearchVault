using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,MAT,LineSearches ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

PhSSPh_450nm_TA_400Hz = EnVisionKinetic(2,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Disulfide_RadicalAnion/Experiemnts/Nanosecond TA/2025_05_30/PhSSPh_TAKinetics_450nm_400Hz_3mLFlow_2ms_Delay.mat")

PhSSPh_450nm_TA_120Hz = EnVisionKinetic(2,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Disulfide_RadicalAnion/Experiemnts/Nanosecond TA/2025_05_30/PhSSPh_TAKinetics_450nm_120Hz_3mLFlow.mat")
PhSSPh_450nm_TA_500Hz = EnVisionKinetic(2,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Disulfide_RadicalAnion/Experiemnts/Nanosecond TA/2025_05_30/PhSSPh_TAKinetics_450nm_500Hz_3mLFlow.mat")
PhSSPh_450nm_TA_500Hz_6mL = EnVisionKinetic(2,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Disulfide_RadicalAnion/Experiemnts/Nanosecond TA/2025_05_30/PhSSPh_TAKinetics_450nm_500Hz_6mLFlow.mat")

PhSSPh_450nm_TA_1000Hz = EnVisionKinetic(2,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Disulfide_RadicalAnion/Experiemnts/Nanosecond TA/2025_05_30/PhSSPh_TAKinetics_450nm_1000Hz_3mLFlow.mat")
Rubpy_Flow = EnVisionKinetic(2,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Disulfide_RadicalAnion/Experiemnts/Nanosecond TA/2025_05_30/Ru(bpy)_450nm_TaKinetics_355nm_3mLFlow.mat")

# Data Workup
####################################################################################################

Plots.plot(PhSSPh_450nm_TA_120Hz[:,1],PhSSPh_450nm_TA_120Hz[:,2])
Plots.plot!(PhSSPh_450nm_TA_500Hz[:,1],PhSSPh_450nm_TA_500Hz[:,2])
Plots.plot!(PhSSPh_450nm_TA_500Hz_6mL[:,1],PhSSPh_450nm_TA_500Hz_6mL[:,2])
Plots.plot!(PhSSPh_450nm_TA_1000Hz[:,1],PhSSPh_450nm_TA_1000Hz[:,2])

Plots.plot(Rubpy_Flow[:,1],Rubpy_Flow[:,2])



PhSSPh_450nm_TA_400Hz_BifitIRF = BiFitIRF(PhSSPh_450nm_TA_400Hz,2,-10.0,1500000.0,[4.539, 1.571, 9.047, 1716.938, 210721.211, 3621.499, 24264.975, 38.294, 43.471],1,Inf,false,1000)
PhSSPh_450nm_TA_500Hz_BifitIRF = BiFitIRF(PhSSPh_450nm_TA_500Hz,2,-10.0,522272.0,[4.777, 1.571, 9.04, 2007.634, 126152.315, 2877.244, 17258.239, 188.943, 44.091],0,0,false,1000)
PhSSPh_450nm_TA_1000Hz_BifitIRF = BiFitIRF(PhSSPh_450nm_TA_1000Hz,2,-10.0,522272.0,[4.699, 1.571, 9.047, 177.997, 123589.058, 257.93, 17261.395, 20.452, 43.749],0,0,false,1000)

Rubpy_BifitIRF = BiFitIRF(Rubpy_Flow,2,-10.0,3488.0,[-9.476, 5.272, 9.542, 200, 165038.892, 834064.318, 394.081, 0.0, 6.509],0,0,false,1000)


### Make Figures 
####################################################################################################

function Fig1(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], 
            palette = (color = palette(ColorSchemes.hsv, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavenumber", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        lines!(x,Data[:,2],linewidth = width,label = "10 ns")


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Fig1()
save("JULIA/Outputs/Fig1.pdf",Fig1())






