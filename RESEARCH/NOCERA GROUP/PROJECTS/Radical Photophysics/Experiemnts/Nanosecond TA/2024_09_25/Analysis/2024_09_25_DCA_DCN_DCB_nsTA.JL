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


#### DCB Ns TA 1mM, 291 nm Pump ################################################################################################
    DCB_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical Photophysics/Experiemnts/Nanosecond TA/2024_09_25/Analysis/DCB_1mM_MeCN_2mm_291nmpump_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    DCB_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical Photophysics/Experiemnts/Nanosecond TA/2024_09_25/Analysis/DCB_1mM_MeCN_2mm_291nmpump_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)


    DCB_TA_Kinetics_440nm = BiFitIRF_new(DCB_TA_Kinetics,2,-10.0, 5000, [4.4806, 0.404, 0.0653, 11.6751, 6.5572, 3688.3081, -0.0007], 10,300,false,200)
    DCB_TA_Kinetics_519nm = BiFitIRF_new(DCB_TA_Kinetics,4,-10.0, 5000, [4.4806, 0.404, 0.0261, 11.7496, 10.0272, 3535.3549, -0.0007], 10,100,false,200)
    DCB_TA_Kinetics_612nm = BiFitIRF_new(DCB_TA_Kinetics,3,-10.0, 5000, [4.4805, 0.404, 0.0551, 18.7118, 0.3715, 3535.4441, -0.0003], 10,100,false,200)



    function DCB_ns_TA_Fig(f = Figure())

        TA_Spectra = DCB_TA_Spectra
        TA_Kinetics = [
            (DCB_TA_Kinetics_440nm[1], " 440 nm"),
            (DCB_TA_Kinetics_519nm[1], "519 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "DCB TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                CairoMakie.ylims!(-0.002,0.01)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "DCB TA Kinetics",
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

            for (data, label) in TA_Kinetics
                lines!(ax2, data[:, 1], data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1], data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)
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
    DCB_ns_TA_Fig()
    save("JULIA/Outputs/DCB_ns_TA_Fig.png",DCB_ns_TA_Fig())


#### DCN Ns TA 500µM, 333 nm Pump ################################################################################################
    DCN_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical Photophysics/Experiemnts/Nanosecond TA/2024_09_25/Analysis/DCN_500uM_MeCN_2mm_333nmpump_TA_Spectra_Surface_processed.csv", header = true, DataFrame)
    DCN_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical Photophysics/Experiemnts/Nanosecond TA/2024_09_25/Analysis/DCN_500uM_MeCN_2mm_333nmpump_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

    DCN_TA_Kinetics_447nm = BiFitIRF_new(DCN_TA_Kinetics,3,-10.0, 5000, [0.1524, 0.9264, 0.1457, 15.9869, 19.7026, 3322.9794, -0.0001], 10,100,false,200)
    DCN_TA_Kinetics_484nm = BiFitIRF_new(DCN_TA_Kinetics,2,-10.0, 5000,[0.1524, 0.9264, 0.3506, 12.309, 15.7783, 3510.3211, 0.0], 10,1000,false,200)
    DCN_TA_Kinetics_852nm = BiFitIRF_new(DCN_TA_Kinetics,4,-10.0, 5000, [0.1372, 0.8338, 0.0529, 17.4227, 2.0174, 3514.1113, -0.0004], 10,100,false,200)

    DCN_global = GlobalIRF(DCN_TA_Kinetics, [2, 3], -10.0, 5000.0, :bi,[0.13, 0.9, 17.296, 3863.465],[1.771, 11.858, 0.0, 0.982, 20.885, -0.0], 50, 100, 2,0.1, 200)


    function DCN_ns_TA_Fig(f = Figure())

        TA_Spectra = DCN_TA_Spectra
        TA_Kinetics = [
            (DCN_TA_Kinetics_447nm[1], " 447 nm"),
            (DCN_TA_Kinetics_484nm[1], "484 nm"), 
            (DCN_TA_Kinetics_852nm[1], "852 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "DCN Ns TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(-0.002,0.01)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "DCN Ns TA Kinetics",
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

            for (data, label) in TA_Kinetics
                lines!(ax2, data[:, 1], data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1], data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)
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
    DCN_ns_TA_Fig()
    save("JULIA/Outputs/DCN_ns_TA_Fig.png",DCN_ns_TA_Fig())

#### DCA Ns TA 220 µM, 425 nm Pump ################################################################################################
    DCA_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical Photophysics/Experiemnts/Nanosecond TA/2024_09_25/Analysis/DCA_220uM_MeCN_2mm_425nmpump_YA_Spectra_Surface_processed.csv", header = true, DataFrame)
    DCA_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radical Photophysics/Experiemnts/Nanosecond TA/2024_09_25/Analysis/DCA_220uM_MeCN_2mm_425nmpump_YA_Kinetics_Surface_processed.csv", header = true, DataFrame)


    DCA_TA_Kinetics_353nm = BiFitIRF_new(DCA_TA_Kinetics,5,-10.0, 5000, [0.1841, 1.0725, 90.6849, 15.6736, -88.2289, 15.6738, -0.0011], Inf,Inf,false,200)
    DCA_TA_Kinetics_414nm = BiFitIRF_new(DCA_TA_Kinetics,4,-10.0, 5000, [0.1092, 1.0734, 89.6959, 17.2196, -89.2176, 17.2462, 0.001], Inf,Inf,false,200)
    DCA_TA_Kinetics_447nm = BiFitIRF_new(DCA_TA_Kinetics,3,-10.0, 5000, [0.1331, 1.0329, 90.1214, 18.1088, -88.7924, 18.2078, 0.0022], Inf,Inf,false,200)
    DCA_TA_Kinetics_618nm = BiFitIRF_new(DCA_TA_Kinetics,2,-10.0, 5000,[0.1464, 0.9702, 89.8437, 18.8483, -89.07, 18.9096, 0.0002], 10,100,false,200)


    DCN_global = GlobalIRF(DCA_TA_Kinetics, [2, 3, 4], -10.0, 100.0, :mono,[0.178, 1.052, 13.844],[0.809, -0.0, 1.415, 0.001, 0.501, 0.001], 100, 100, 20,0.1, 200)


    function DCA_ns_TA_Fig(f = Figure())

        TA_Spectra = DCA_TA_Spectra
        TA_Kinetics = [
            
            (DCA_TA_Kinetics_414nm[1], "414 nm"), 
            (DCA_TA_Kinetics_447nm[1], "447 nm"), 
            (DCA_TA_Kinetics_618nm[1], "618 nm")
        ]
        
        ### First Panel TA Spectra  ##############################################################################

            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 
                    title = "DCA Ns TA Spectra",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "ΔA (OD)", 
                    xminorticks = IntervalsBetween(2)
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(-0.002,0.01)

                for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                    lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
                
                end

                axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false,labelsize = 15)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        #### Second Panel TA KINETICS ##############################################################################

            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "DCA Ns TA Kinetics",
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

            for (data, label) in TA_Kinetics
                lines!(ax2, data[:, 1], data[:, 2],linewidth = 3,label = label) # Kinetic Trace
                lines!(ax2, data[:, 1], data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
            end

                axislegend("Wavelength", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)
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
    DCA_ns_TA_Fig()
    save("JULIA/Outputs/DCA_ns_TA_Fig.png",DCA_ns_TA_Fig())





