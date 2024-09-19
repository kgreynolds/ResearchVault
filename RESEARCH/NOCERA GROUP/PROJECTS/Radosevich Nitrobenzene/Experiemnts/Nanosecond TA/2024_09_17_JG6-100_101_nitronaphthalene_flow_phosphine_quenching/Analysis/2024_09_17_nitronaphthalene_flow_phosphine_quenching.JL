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

# Data Workup
####################################################################################################
Data = MonoFitIRF_New(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017, 0.0], 10,Inf,false,10000)
Data = BiFitIRF_new(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017,10.0, 100, 0.0], 10,Inf,false,10000)


### Nitronaphthalene quenching by Phosphine in flow (3mL/min at 405 nm ~90 µJ) 
####################################################################################################

##### NNaph 0ul Phosphine ####################################################################################################
    NNaph_0ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-1_nitronaphthalene_0uL-phosphine_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_0ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-1_nitronaphthalene_0uL-phosphine_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_0ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_0ul_Phos_Kinetics,2,-10.0, 30000, [3.068, 2.007, 62.082, 2024.024, 0.001], 10,Inf,false,200)
    NNaph_0ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_0ul_Phos_Kinetics,3,-10.0, 30000, [3.068, 2.007, 125.69, 2063.693, 0.0], 10,Inf,false,200)
    # NNaph_0ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_0ul_Phos_Kinetics,3,-10.0, 30000, [3.1466, 1.5492, 115.4766, 2582.2144, 16.8124, 767.6591, 0.0001], 10,Inf,false,200)

    function NNaph_0ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_0ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_0ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_0ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 0 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 0 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_0ul_TA_Fig()
    save("JULIA/Outputs/NNaph_0ul_TA_Fig.png",NNaph_0ul_TA_Fig())

##### NNaph 50ul Phosphine ####################################################################################################
    NNaph_50ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-2_nitronaphthalene_50uL-phosphine_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_50ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-2_nitronaphthalene_50uL-phosphine_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_50ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_50ul_Phos_Kinetics,2,-10.0, 5000, [3.708, 1.362, 8.095, 315.776, 0.01], 10,10,false,200)
    NNaph_50ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_50ul_Phos_Kinetics,3,-10.0, 5000, [3.708, 1.362, 20.046, 312.896, 0.001], 10,10,false,200)
    #NNaph_0ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_0ul_Phos_Kinetics,3,-10.0, 30000, [0.008, 0.001, 115.6171, 2582.2166, 16.8814, 767.6586, 0.0001], 10,10,false,200)

    function NNaph_50ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_50ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_50ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_50ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 50 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 50 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_50ul_TA_Fig()
    save("JULIA/Outputs/NNaph_50ul_TA_Fig.png",NNaph_50ul_TA_Fig())

##### NNaph 100ul Phosphine ####################################################################################################
    NNaph_100ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-3_nitronaphthalene_100uL-phosphine_TS_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_100ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-3_nitronaphthalene_100uL-phosphine_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_100ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_100ul_Phos_Kinetics,2,-10.0, 5000, [4.54, 0.725, 4.652, 193.32, 0.011], 10,10,false,200)
    NNaph_100ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_100ul_Phos_Kinetics,3,-10.0, 5000, [4.731, 0.766, 12.273, 193.972, 0.0], 10,10,false,200)
    #NNaph_100ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_100ul_Phos_Kinetics,3,-10.0, 30000, [0.008, 0.001, 115.6171, 2582.2166, 16.8814, 767.6586, 0.0001], 10,10,false,200)

    function NNaph_100ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_100ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_100ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_100ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 100 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 100 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_100ul_TA_Fig()
    save("JULIA/Outputs/NNaph_100ul_TA_Fig.png",NNaph_100ul_TA_Fig())

##### NNaph 200ul Phosphine ####################################################################################################
    NNaph_200ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-4_nitronaphthalene_200uL_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_200ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-4_nitronaphthalene_200uL_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_200ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_200ul_Phos_Kinetics,2,-10.0, 5000, [4.139, 0.756, 2.354, 109.336, 0.013], 10,10,false,200)
    NNaph_200ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_200ul_Phos_Kinetics,3,-10.0, 5000, [4.57, 0.773, 6.385, 101.638, -0.0], 10,10,false,200)
    #NNaph_200ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_200ul_Phos_Kinetics,3,-10.0, 30000, [0.008, 0.001, 115.6171, 2582.2166, 16.8814, 767.6586, 0.0001], 10,10,false,200)

    function NNaph_200ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_200ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_200ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_200ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 200 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 200 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_200ul_TA_Fig()
    save("JULIA/Outputs/NNaph_200ul_TA_Fig.png",NNaph_200ul_TA_Fig())

##### NNaph 300ul Phosphine ####################################################################################################
    NNaph_300ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-5_nitronaphthalene_300uL_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_300ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-5_nitronaphthalene_300uL_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_300ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_300ul_Phos_Kinetics,2,-10.0, 5000, [4.211, 0.754, 1.694, 78.186, 0.013], 10,10,false,200)
    NNaph_300ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_300ul_Phos_Kinetics,3,-10.0, 5000, [4.657, 0.769, 4.491, 71.255, -0.0001], 10,10,false,200)
    #NNaph_300ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_300ul_Phos_Kinetics,3,-10.0, 30000, [0.008, 0.001, 115.6171, 2582.2166, 16.8814, 767.6586, 0.0001], 10,10,false,200)

    function NNaph_300ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_300ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_300ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_300ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 300 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 300 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_300ul_TA_Fig()
    save("JULIA/Outputs/NNaph_300ul_TA_Fig.png",NNaph_300ul_TA_Fig())

##### NNaph 400ul Phosphine ####################################################################################################
    NNaph_400ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-6_nitronaphthalene_400uL_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_400ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-6_nitronaphthalene_400uL_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_400ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_400ul_Phos_Kinetics,2,-10.0, 5000, [4.271, 0.752, 1.266, 59.791, 0.014], 10,10,false,200)
    NNaph_400ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_400ul_Phos_Kinetics,3,-10.0, 5000, [4.638, 0.77, 3.42, 53.849, 0.0], 10,10,false,200)
    #NNaph_400ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_400ul_Phos_Kinetics,3,-10.0, 30000, [0.008, 0.001, 115.6171, 2582.2166, 16.8814, 767.6586, 0.0001], 10,10,false,200)

    function NNaph_400ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_400ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_400ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_400ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 400 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 400 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_400ul_TA_Fig()
    save("JULIA/Outputs/NNaph_400ul_TA_Fig.png",NNaph_400ul_TA_Fig())

##### NNaph 500ul Phosphine ####################################################################################################
    NNaph_500ul_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-7_nitronaphthalene_500uL_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    NNaph_500ul_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024_09_17_JG6-100_101_nitronaphthalene_flow_phosphine_quenching/Analysis/JG6-100-7_nitronaphthalene_500uL_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    NNaph_500ul_Phos_Kinetics_420nm = MonoFitIRF_New(NNaph_500ul_Phos_Kinetics,2,-10.0, 5000, [4.149, 0.755, 0.958, 48.643, 0.014], 10,10,false,200)
    NNaph_500ul_Phos_Kinetics_570nm = MonoFitIRF_New(NNaph_500ul_Phos_Kinetics,3,-10.0, 5000, [4.512, 0.709, 2.678, 43.779, -0.0], 10,10,false,200)
    #NNaph_500ul_Phos_Kinetics_570nm = BiFitIRF_new(NNaph_500ul_Phos_Kinetics,3,-10.0, 30000, [0.008, 0.001, 115.6171, 2582.2166, 16.8814, 767.6586, 0.0001], 10,10,false,200)

    function NNaph_500ul_TA_Fig(f = Figure())

        TA_Spectra = NNaph_500ul_Phos_Spectra
        TA_Kinetics = [
            (NNaph_500ul_Phos_Kinetics_420nm[1], " 420 nm"),
            (NNaph_500ul_Phos_Kinetics_570nm[1], "570 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "Nitronaphthalene + 500 µL Phosphine TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(400, 410; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 3),),
                    title = "Nitronaphthalene + 500 µL Phosphine TA Kinetcs",
                    xlabel = "Time (µs)", 
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
                lines!(ax2, data[:, 1] ./1000, data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1] ./1000, data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 2, framevisible = false, labelsize = 15)
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
    NNaph_500ul_TA_Fig()
    save("JULIA/Outputs/NNaph_500ul_TA_Fig.png",NNaph_500ul_TA_Fig())


### COmparisons 
function Kinetics_Comparison(f = Figure())

    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.winter, 6),), ## or replace ncol with Integer 
                title = "Nitronaphthalene quenching by Phosphine",
                xlabel = "Time (ns)", 
                ylabel = "Normalized ΔA", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            ax.xticks= 0:250:1200
            CairoMakie.xlims!(nothing, 1200)
            # CairoMakie.ylims!(-0.02,0.7)

        ### PlotUvVis  Data ##############################################################################

         
            # Makie.lines!(ax, NNaph_0ul_Phos_Kinetics_570nm[1][:, 1], NNaph_0ul_Phos_Kinetics_570nm[1][:, 2],linewidth = width, label = "0 µL")
            Makie.lines!(ax, NNaph_50ul_Phos_Kinetics_570nm[1][:, 1], MaxNorm(NNaph_50ul_Phos_Kinetics_570nm[1][:, 2]),linewidth = width, label = "50 µL")
            Makie.lines!(ax, NNaph_100ul_Phos_Kinetics_570nm[1][:, 1], MaxNorm(NNaph_100ul_Phos_Kinetics_570nm[1][:, 2]),linewidth = width, label = "100 µL")
            Makie.lines!(ax, NNaph_200ul_Phos_Kinetics_570nm[1][:, 1], MaxNorm(NNaph_200ul_Phos_Kinetics_570nm[1][:, 2]),linewidth = width, label = "200 µL")
            Makie.lines!(ax, NNaph_300ul_Phos_Kinetics_570nm[1][:, 1], MaxNorm(NNaph_300ul_Phos_Kinetics_570nm[1][:, 2]),linewidth = width, label = "300 µL")
            Makie.lines!(ax, NNaph_400ul_Phos_Kinetics_570nm[1][:, 1], MaxNorm(NNaph_400ul_Phos_Kinetics_570nm[1][:, 2]),linewidth = width, label = "400 µL")
            Makie.lines!(ax, NNaph_500ul_Phos_Kinetics_570nm[1][:, 1], MaxNorm(NNaph_500ul_Phos_Kinetics_570nm[1][:, 2]),linewidth = width, label = "500 µL")


            axislegend("Phosphine Concentration", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)

            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        
    resize_to_layout!(f)
    f
end
Kinetics_Comparison()
save("JULIA/Outputs/Kinetics_Comparison.png",Kinetics_Comparison())


