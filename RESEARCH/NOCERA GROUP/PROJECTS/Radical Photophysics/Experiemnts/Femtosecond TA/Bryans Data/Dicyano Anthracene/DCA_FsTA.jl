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

process_Surface_csv("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical_Lifetimes/Experiemnts/Femtosecond TA/Bryans Data/Dicyano Anthracene/DCA_FsTA_Spectra.csv", make_file = true)

## IMPORT DATA
DCA_fsTA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical_Lifetimes/Experiemnts/Femtosecond TA/Bryans Data/Dicyano Anthracene/DCA_FsTA_Spectra_Surface_processed.csv", header = true, DataFrame)
DCA_fsTA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical_Lifetimes/Experiemnts/Femtosecond TA/Bryans Data/Dicyano Anthracene/DCA_FsTA_Kinetics.csv", header = true, DataFrame)

# Data Workup
####################################################################################################
DCA_fsTA_374nm_Kinetics = BiFitIRF_new(DCA_fsTA_Kinetics,2,0.0, 100,[0.0009, 0.09, 0.1548, 43.9892, -0.0056, 4.314, -0.0002],10,100,false,100)
DaDCA_fsTA_474nm_Kineticsta = MonoFitIRF_New(DCA_fsTA_Kinetics,3,0.0, 100,[0.008, 0.001, 0.022, 4.314, -0.0], 100,Inf,false,100)


### Make Figures 
####################################################################################################


function DCA_fsTA_Fig(f = Figure())

    TA_Spectra = DCA_fsTA_Spectra
    TA_Kinetics = [
        (DCA_fsTA_374nm_Kinetics[1], " 374 nm"),
        (DaDCA_fsTA_474nm_Kineticsta[1], "474 nm")
    ]
    
    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                title = "FS TA on DCA radican anion 710 nm pump",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(350,600)
            CairoMakie.ylims!(-0.001,0.011)

            for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
            
            end

            axislegend("Delay Time", position = :rt,nbanks = 4,framevisible = false,labelsize = 15)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    #### Second Panel TA KINETICS ##############################################################################

        ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "FS TA on DCA radical anion 710 nm pump",
                xlabel = "Time (ps)", 
                ylabel = "ΔA (OD)", 
                # xscale = Makie.Symlog10(1.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                xminorticks = IntervalsBetween(2)
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(-10, 2000)
            # CairoMakie.ylims!(-0.001,0.02)

        for (data, label) in TA_Kinetics
            lines!(ax2, data[:, 1], data[:, 2],linewidth = 3,label = label) # Kinetic Trace
            lines!(ax2, data[:, 1], data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
        end

            axislegend("Wavelength", position = :rt, nbanks = 4, framevisible = false, labelsize = 15)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                ### Add Inset axis into plot for kinetics
                # ax2_1 = Axis(f[1,2], width=Relative(0.5), height=Relative(0.5),halign=0.8, valign=0.5, backgroundcolor=:white, 
                #xlabel = "Time (ns)", 
                #ylabel = "ΔA (OD)", 
                # )

                #CairoMakie.xlims!(-10, 20)
                # CairoMakie.ylims!(-0.003,0.02)

                # for (data, label) in TA_Kinetics
                #     Makie.lines!(ax2_1, data[:, 1], data[:, 2],linewidth = 3,label = label)
                #     lines!(ax2_1, data[:, 1], data[:, 3],linewidth = 2, linestyle = :dash, color = :black)
                # end

            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
DCA_fsTA_Fig()
save("JULIA/Outputs/DCA_fsTA_Fig.pdf",DCA_fsTA_Fig())






