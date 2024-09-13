## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Outside(),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)

# Open and read MatLab Files (For Magnitude)
File = EnVisionKinetic(2,"path")

## IMPORT DATA
file1 = CSV.read("path", header = true, DataFrame)

# Data Workup
####################################################################################################
Data = MonoFitIRF_New(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017, 0.0], 10,Inf,false,10000)
Data = BiFitIRF_new(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017,10.0, 100, 0.0], 10,Inf,false,10000)


### Make Figures 
####################################################################################################


function Fig(f = Figure())

    Data = PhAuCl3_Photolyzed_TA_Spectra

    ### First Panel ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                title = "Title",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

            for i in 2:ncol(Data) ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, Data[:,1],Data[:,i],linewidth = width,label = " $(names(Data)[i])")
            
            end

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        

    # ### Second Panel ##############################################################################

        ax2 =    CairoMakie.Axis(
                f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "Title",
                xlabel = "Time (ns)", 
                ylabel = "ΔA (OD)", 
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

            lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,2],linewidth = width,label = "nm")
            lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

            lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,2],linewidth = width,label = "nm")
            lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

            lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,2],linewidth = width,label = "nm")
            lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

            axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
Fig()
save("JULIA/Outputs/Fig.pdf",Fig())

function UvVis_Fig(f = Figure())

    data = MesAcr_PhosphineOxide_UvVis
    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.hsv, 3),), ## or replace ncol with Integer 

                title = "Title",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(350,600)
            CairoMakie.ylims!(-0.02,0.7)

        ### PlotUvVis  Data ##############################################################################


            Makie.lines!(ax, data[:, 1], data[:, 2], # This is the raw kinetic Trace
                linewidth = 3, 
                label = "pre TA"
            )
            Makie.lines!(ax, data[:, 1], data[:, 3], # This is the raw kinetic Trace
                linewidth = 3, 
                label = "Post TA"
            )
  

            axislegend(
                position = :rt, 
                nbanks = 3, 
                framevisible = false, 
                labelsize = 15
            )

            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        
    resize_to_layout!(f)
    f
end
UvVis_Fig()
save("JULIA/Outputs/UvVis_Fig.pdf",UvVis_Fig())


function TA_Fig(f = Figure())

    TA_Spectra = MesAcr_PhosphineOxide_TASpectra
    TA_Kinetics = [
        (MesAcr_PhosphineOxide_TAKinetics_391nm, " 391 nm"),
        (MesAcr_PhosphineOxide_TAKinetics_560nm, "560 nm"),
        (MesAcr_PhosphineOxide_TAKinetics_595nm, "595 nm")
    ]
    
    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 

                title = "Title",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

            for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
            
            end

            axislegend("Delay Time", 
                position = :rt, 
                nbanks = 3, 
                framevisible = false, 
                labelsize = 15
            )

            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        

    # ### Second Panel TA KINETICS ##############################################################################

        ax2 =    CairoMakie.Axis(
                f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "Title",
                xlabel = "Time (ns)", 
                ylabel = "ΔA (OD)", 
                # xscale = Makie.Symlog10(1.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                xminorticks = IntervalsBetween(2)
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(-10, 2000)
            CairoMakie.ylims!(-0.001,0.02)

        ### Plot Kinetic Data ##############################################################################

        for (data, label) in TA_Kinetics
            Makie.lines!(ax2, data[:, 1], data[:, 2], # This is the raw kinetic Trace
                linewidth = 3, 
                label = label
            )
            lines!(ax2, data[:, 1], data[:, 3], # This is the Fitted curve
                linewidth = 2, 
                linestyle = :dash, 
                color = :black
            )
        end

            axislegend("Wavelength", 
                position = :rt, 
                nbanks = 3, 
                framevisible = false, 
                labelsize = 15
                )

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

            ### Add Inset axis into plot for kinetics
            # ax2_1 = Axis(f[1,2], 
            #     width=Relative(0.5), 
            #     height=Relative(0.5), 
            #     halign=0.8, valign=0.5, backgroundcolor=:white, 
            #     xlabel = "Time (ns)", 
            #     ylabel = "ΔA (OD)", 
            # )

            #             CairoMakie.xlims!(-10, 20)
            #             # CairoMakie.ylims!(-0.003,0.02)

            # for (data, label) in TA_Kinetics
            #     Makie.lines!(ax2_1, data[:, 1], data[:, 2], # This is the raw kinetic Trace
            #         linewidth = 3, 
            #         label = label
            #     )
            #     lines!(ax2_1, data[:, 1], data[:, 3], # This is the Fitted curve
            #         linewidth = 2, 
            #         linestyle = :dash, 
            #         color = :black
            #     )
            # end


            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
TA_Fig()
save("JULIA/Outputs/TA_Fig.pdf",TA_Fig())






