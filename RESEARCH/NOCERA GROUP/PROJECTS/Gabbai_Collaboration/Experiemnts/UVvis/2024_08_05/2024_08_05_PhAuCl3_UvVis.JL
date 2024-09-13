## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

# Open and read MatLab Files (For Magnitude)
File = EnVisionKinetic(2,"path")

## IMPORT DATA
Ph_Au_Cl3_MeCN_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/UVvis/Ph_Au_Cl3_MeCN_UvVis.csv", header = true, DataFrame)
Ph_Au_Cl3_MeCN_UvVis_post25nmTA = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/UVvis/Ph_Au_Cl3_MeCN_UvVis_postTA_425nm.csv", header = true, DataFrame)




# Data Workup
####################################################################################################





### Make Figures 
####################################################################################################

function UvVis(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "PhAuCl3 TA in MeCN",
            palette = (color = palette(ColorSchemes.hsv, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "A"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(270,550)
        # CairoMakie.ylims!(nothing,)
        
        lines!(Ph_Au_Cl3_MeCN_UvVis[:,1],Ph_Au_Cl3_MeCN_UvVis[:,2],linewidth = width,label = "Pre-TA")
        lines!(Ph_Au_Cl3_MeCN_UvVis_post25nmTA[:,1],Ph_Au_Cl3_MeCN_UvVis_post25nmTA[:,2],linewidth = width,label = "Pre-TA")


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
UvVis()
save("JULIA/Outputs/UvVis.pdf",UvVis())



