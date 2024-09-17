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
## IMPORT DATA
file1 = CSV.read("path", header = true, DataFrame)
Data = MonoFitIRF_New(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017, 0.0], 10,Inf,false,10000)
Data = BiFitIRF_new(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017,10.0, 100, 0.0], 10,Inf,false,10000)



#### Pd_TFP ################################################################################################
    Pd_TFP_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/PdTFP_TA_Spectra.csv", header = true, DataFrame)
    Pd_TFP_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/PdTFP_TA_Kinetics.csv", header = true, DataFrame)

    Pd_TFP_Kinetics_415nm = MonoFitIRF_New(Pd_TFP_Kinetics,2,-100.0, 4700, [0.001, 0.981, -69.073, 603.598, -0.001], 10,Inf,false,200)
    Pd_TFP_Kinetics_456nm = MonoFitIRF_New(Pd_TFP_Kinetics,3,-100.0, 4700, [0.001, 0.981, 17.159, 608.947, -0.0001], 10,Inf,false,200)
    Pd_TFP_Kinetics_524nm = MonoFitIRF_New(Pd_TFP_Kinetics,4,-100.0, 4700, [0.001, 0.981, -5.148, 564.724, 0.0], 10,Inf,false,200)


    τ_mean_Pd_TFP = mean([603.598,608.947,564.724])
    τ_stdev_Pd_TFP = std([603.598,608.947,564.724])

    function Pd_TFP_TA_Fig(f = Figure())

        TA_Spectra = Pd_TFP_Spectra
        TA_Kinetics = [
            (Pd_TFP_Kinetics_415nm[1], " 415 nm"),
            (Pd_TFP_Kinetics_456nm[1], "456 nm"),
            (Pd_TFP_Kinetics_524nm[1], "524 nm")
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
                # CairoMakie.xlims!(375,425)
                # CairoMakie.ylims!(nothing,0.04)
    
            ### Plot Data ##############################################################################
    
                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end
    
                axislegend("Delay Time", 
                    position = :rb, 
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

                # ### Add Inset axis into plot for Spectra
                # ax1_inset = Axis(f[1,1],palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),),
                #     width=Relative(0.5), 
                #     height=Relative(0.5), 
                #     halign=0.8, valign=0.5, backgroundcolor=:white, 
                #     xlabel = "Time (ns)", 
                #     ylabel = "ΔA (OD)", 
                # )
    
                # CairoMakie.xlims!(427, 850)
                # CairoMakie.ylims!(-0.015,nothing)

                # for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                #     lines!(ax1_inset, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                # end
  

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
                # CairoMakie.ylims!(-0.001,0.02)
    
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
    Pd_TFP_TA_Fig()
    save("JULIA/Outputs/Pd_TFP_TA_Fig.pdf",Pd_TFP_TA_Fig())


#### Pd_TFP_CH3 ################################################################################################
    Pd_TFP_CH3_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/PdTFPCH3_TA_Spectra.csv", header = true, DataFrame)
    Pd_TFP_CH3_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/PdTFPCH3_TA_Kinetics.csv", header = true, DataFrame)

    Pd_TFP_CH3_Kinetics_417nm = MonoFitIRF_New(Pd_TFP_CH3_Kinetics,2,-100.0, 4700, [0.001, 0.981, -62.275, 657.05, -0.0], 10,Inf,false,200)
    Pd_TFP_CH3_Kinetics_454nm = MonoFitIRF_New(Pd_TFP_CH3_Kinetics,3,-100.0, 4700, [0.001, 0.981, 17.115, 658.451, -0.0], 10,Inf,false,200)
    Pd_TFP_CH3_Kinetics_523nm = MonoFitIRF_New(Pd_TFP_CH3_Kinetics,4,-100.0, 4700, [0.001, 0.981, -4.933, 654.094, -0.0], 10,Inf,false,200)

    τ_mean_Pd_TFP_CH3 = mean([657.05,658.451,654.094])
    τ_stdev_Pd_TFP_CH3 = std([657.05,658.451,654.094])

    function Pd_TFP_CH3_TA_Fig(f = Figure())

        TA_Spectra = Pd_TFP_CH3_Spectra
        TA_Kinetics = [
            (Pd_TFP_CH3_Kinetics_417nm[1], " 417 nm"),
            (Pd_TFP_CH3_Kinetics_454nm[1], "454 nm"),
            (Pd_TFP_CH3_Kinetics_523nm[1], "523 nm")
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
                    position = :rb, 
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
                # CairoMakie.ylims!(-0.001,0.02)
    
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
    Pd_TFP_CH3_TA_Fig()
    save("JULIA/Outputs/Pd_TFP_CH3_TA_Fig.pdf",Pd_TFP_CH3_TA_Fig())

#### Zn_Chlorin Band 1 ################################################################################################
    Zn_Chlorin_b1_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band1_TA_Spectra.csv", header = true, DataFrame)
    Zn_Chlorin_b1_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band1_TA_Kinetics.csv", header = true, DataFrame)

    Zn_Chlorin_b1_398nm = MonoFitIRF_New(Zn_Chlorin_b1_Kinetics,3,-100.0, 4700, [0.001, 0.981, -5.406, 635.619, 0.0], 10,Inf,false,200)
    Zn_Chlorin_b1_417nm = MonoFitIRF_New(Zn_Chlorin_b1_Kinetics,4,-100.0, 4700, [0.001, 0.981, -72.844, 648.133, -0.0], 10,Inf,false,200)
    Zn_Chlorin_b1_444nm = MonoFitIRF_New(Zn_Chlorin_b1_Kinetics,5,-100.0, 4700, [0.001, 0.981, 11.748, 653.121, 0.0], 10,Inf,false,200)
    Zn_Chlorin_b1_580nm = MonoFitIRF_New(Zn_Chlorin_b1_Kinetics,6,-100.0, 4700, [0.001, 0.981, -0.854, 609.212, -0.0], 10,Inf,false,200)
    Zn_Chlorin_b1_613nm = MonoFitIRF_New(Zn_Chlorin_b1_Kinetics,7,-100.0, 4700, [0.001, 0.981, -14.809, 664.592, 0.0], 10,Inf,false,200)

    τ_mean_Zn_Chlorin_b1 = mean([648.133,653.121, 664.592])
    τ_stdev_Zn_Chlorin_b1 = std([648.133,653.121, 664.592])


    function Zn_Chlorin_b1_TA(f = Figure())

        TA_Spectra = Zn_Chlorin_b1_Spectra
        TA_Kinetics = [
            (Zn_Chlorin_b1_417nm[1], " 417 nm"),
            (Zn_Chlorin_b1_444nm[1], "444 nm"),
            (Zn_Chlorin_b1_613nm[1], "613 nm")
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
                    position = :rb, 
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
                # CairoMakie.ylims!(-0.001,0.02)
    
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
    Zn_Chlorin_b1_TA()
    save("JULIA/Outputs/TA_Fig.pdf",Zn_Chlorin_b1_TA())

#### Zn_Chlorin Band 2 ################################################################################################
    Zn_Chlorin_b2_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_TA_Spectra.csv", header = true, DataFrame)
    Zn_Chlorin_b2_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_TA_Kinetics.csv", header = true, DataFrame)

    Zn_Chlorin_b2_399nm = MonoFitIRF_New(Zn_Chlorin_b2_Kinetics,2,-100.0, 4700, [0.001, 0.981, -6.815, 754.606, 0.0], 10,Inf,false,200)
    Zn_Chlorin_b2_417nm = MonoFitIRF_New(Zn_Chlorin_b2_Kinetics,3,-100.0, 4700, [0.001, 0.981, -85.302, 753.715, -0.0], 10,Inf,false,200)
    Zn_Chlorin_b2_445nm = MonoFitIRF_New(Zn_Chlorin_b2_Kinetics,4,-100.0, 4700, [0.001, 0.981, 16.713, 740.315, 0.0], 10,Inf,false,200)
    Zn_Chlorin_b2_585nm = MonoFitIRF_New(Zn_Chlorin_b2_Kinetics,5,-100.0, 4700, [0.001, 0.981, -0.967, 628.478, 0.0], 10,Inf,false,200)
    Zn_Chlorin_b2_611nm = MonoFitIRF_New(Zn_Chlorin_b2_Kinetics,6,-100.0, 4700, [0.001, 0.981, -16.907, 779.807, 0.0], 10,Inf,false,200)

    τ_mean_Zn_Chlorin_b2 = mean([753.715,740.315, 779.807])
    τ_stdev_Zn_Chlorin_b2 = std([753.715,740.315, 779.807])


    function Zn_Chlorin_b2_TA(f = Figure())

        TA_Spectra = Zn_Chlorin_b2_Spectra
        TA_Kinetics = [
            (Zn_Chlorin_b2_417nm[1], " 417 nm"),
            (Zn_Chlorin_b2_445nm[1], "445 nm"),
            (Zn_Chlorin_b2_611nm[1], "611 nm")
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
                    position = :rb, 
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
                # CairoMakie.ylims!(-0.001,0.02)
    
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
    Zn_Chlorin_b2_TA()
    save("JULIA/Outputs/Zn_Chlorin_b2_TA.pdf",Zn_Chlorin_b2_TA())

### Make Figures 
####################################################################################################



