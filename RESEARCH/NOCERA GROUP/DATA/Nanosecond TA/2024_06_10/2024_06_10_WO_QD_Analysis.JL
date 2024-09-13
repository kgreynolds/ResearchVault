## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

# Open and read MatLab Files (For Magnitude)

WO_PL = matopen("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/DATA/Nanosecond TA/2024_06_10/WO_PL_355nm_Pump_370-700nm_2ns.mat")
    newmat = read(WO_PL,"data")

PL_spectra = transpose(newmat)
smoothed_PL_Spectra = smoothMat(transpose(newmat),3)

WO_400nm_TA_Kinetics = EnVisionKinetic(0,"/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/DATA/Nanosecond TA/2024_06_10/WO_TA_400nm_355nm_Pump.mat")



Plots.plot(newmat[2:end,1],newmat[2:end, 10])

Plots.plot(PL_spectra[2:end,1],PL_spectra[2:end,104])
Plots.plot(smoothed_PL_Spectra[2:end,1],smoothed_PL_Spectra[2:end,104])


Plots.plot(WO_400nm_TA_Kinetics[:,1],WO_400nm_TA_Kinetics[:,2])

## IMPORT DATA
file1 = CSV.read("path", header = true, DataFrame)




# Data Workup
####################################################################################################

WO_400nm_TA_Kinetics_MonoFit = MonoFitIRF_New(WO_400nm_TA_Kinetics,2,-10.0,100,[1.0,3.0,-100.0,2.0,0.1],Inf,Inf,false,100)




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



