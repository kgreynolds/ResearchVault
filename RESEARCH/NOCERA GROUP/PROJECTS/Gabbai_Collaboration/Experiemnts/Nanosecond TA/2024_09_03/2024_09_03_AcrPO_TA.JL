## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Kris_Figure_Theme = merge(Theme(
    Axis = (
        alignmode = Outside(),xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)

##  MesAcr PhosphineOxide
####################################################################################################


    MesAcr_PhosphineOxide_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/UVvis/2024_09_03/Acr_PhosOx_MeCN_sample2_processed.csv", header = true, DataFrame)


    MesAcr_PhosphineOxide_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPO_TASpectra.csv", header = true, DataFrame)
    MesAcr_PhosphineOxide_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPO_TAKinetics.csv", header = true, DataFrame)


    MesAcr_PhosphineOxide_TAKinetics_391nm = MonoFitIRF_New(MesAcr_PhosphineOxide_TAKinetics,3,-10.0, 2000,[0.552, 1.085, 0.951, 2.323, -0.0], 1,10,false,400)
    MesAcr_PhosphineOxide_TAKinetics_560nm = BiFitIRF_new(MesAcr_PhosphineOxide_TAKinetics,4,-10.0, 2000, [0.4663, 0.9416, 5.6179, 457.4754, 0.2655, 1.8474, 0.0009], 10,10,false,400)
    MesAcr_PhosphineOxide_TAKinetics_595nm = BiFitIRF_new(MesAcr_PhosphineOxide_TAKinetics,5,-10.0, 2000,[-0.008, 0.9957, 6.7437, 474.0959, -0.0462, 2.7784, -0.0], 10,10,false,400)
    MesAcr_PhosphineOxide_TAKinetics_700nm = BiFitIRF_new(MesAcr_PhosphineOxide_TAKinetics,2,-10.0, 2000,[-0.0088, 0.8961, 1.4656, 426.6863, 0.1474, 2.5854, 0.0001], 10,10,false,400)

function MesAcr_PhosphineOxideUvVis_Fig(f = Figure())

    data = MesAcr_PhosphineOxide_UvVis
    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.hsv, 3),), ## or replace ncol with Integer 

                title = "Acr Phosphine Oxide UvVis",
                xlabel = "Wavelength (nm)", 
                ylabel = "A (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(350,600)
            CairoMakie.ylims!(-0.02,1.5)

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
MesAcr_PhosphineOxideUvVis_Fig()
save("JULIA/Outputs/AcrPhos_UvVis_Fig.png",MesAcr_PhosphineOxideUvVis_Fig())


function MesAcr_PhosphineOxideUvVis_FigTA_Fig(f = Figure())

    TA_Spectra = MesAcr_PhosphineOxide_TASpectra
    TA_Kinetics = [
        (MesAcr_PhosphineOxide_TAKinetics_391nm[1], " 391 nm"),
        (MesAcr_PhosphineOxide_TAKinetics_560nm[1], "560 nm"),
        (MesAcr_PhosphineOxide_TAKinetics_595nm[1], "595 nm"), 
        (MesAcr_PhosphineOxide_TAKinetics_700nm[1], "700 nm")
    ]
    
    ### First Panel TA Spectra  ##############################################################################
        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(TA_Spectra)+1),), ## or replace ncol with Integer 

                title = "Acr Phosphine Oxide TA Spectra",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            CairoMakie.ylims!(-0.02,nothing)

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
            CairoMakie.vspan!(431, 441; ymin = 0.0, ymax = 1.0, color = :lightblue)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    # ### Second Panel TA KINETICS ##############################################################################

        ax2 =    CairoMakie.Axis(
                f[1,2], palette = (color = palette(ColorSchemes.hsv, 5),),
                title = "Acr Phosphine Oxide TA Kinetics",
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
MesAcr_PhosphineOxideUvVis_FigTA_Fig()
save("JULIA/Outputs/ MesAcr_PhosphineOxide_TA_Fig.png",MesAcr_PhosphineOxideUvVis_FigTA_Fig())


##  MesAcr PhosphineOxide + Anthracene
####################################################################################################

    MesAcr_PhosphineOxide_Antrhacene_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/UVvis/2024_09_03/Acr_PhosOx_MeCN_Anthracene_sample2_processed.csv", header = true, DataFrame)


    MesAcr_PhosphineOxide_Antrhacene_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPO_Anthracene_TASpectra.csv", header = true, DataFrame)
    MesAcr_PhosphineOxide_Antrhacene_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPO_Anthracene_TAKinetics.csv", header = true, DataFrame)


    MesAcr_PhosphineOxide_Antrhacene_TAKinetics_390nm = MonoFitIRF_New(MesAcr_PhosphineOxide_Antrhacene_TAKinetics,3,-10.0, 2000,[0.47, 1.017, 0.942, 2.446, 0.002], 1,10,false,400)
    MesAcr_PhosphineOxide_Antrhacene_TAKinetics_420nm = BiFitIRF_new(MesAcr_PhosphineOxide_Antrhacene_TAKinetics,2,-10.0, 2000, [-0.0001, 0.7356, 22.3432, 295.0017, -12.2158, 146.903, 0.0005], 100,1000,false,400)
    MesAcr_PhosphineOxide_Antrhacene_TAKinetics_510nm = BiFitIRF_new(MesAcr_PhosphineOxide_Antrhacene_TAKinetics,4,-10.0, 2000,[-0.01, 1.3678, 6.969, 317.0656, 0.2459, 2.2028, 0.0006], 100,100,false,400)
    MesAcr_PhosphineOxide_Antrhacene_TAKinetics_710nm = BiFitIRF_new(MesAcr_PhosphineOxide_Antrhacene_TAKinetics,6,-10.0, 2000,[-0.001, 1.056, 0.6029, 202.7728, 0.151, 2.3662, 0.0007], 100,100,false,400)
    MesAcr_PhosphineOxide_Antrhacene_TAKinetics_895nm = BiFitIRF_new(MesAcr_PhosphineOxide_Antrhacene_TAKinetics,5,-10.0, 2000,[-0.0, 1.613, 4.2712, 311.9637, -0.0362, 1.6686, -0.0004], 10,10,false,400)


    function MesAcr_PhosphineOxide_AntrhaceneUvVis_Fig(f = Figure())

            data = MesAcr_PhosphineOxide_Antrhacene_UvVis
            ### First Panel TA Spectra  ##############################################################################

                ax =    CairoMakie.Axis(
                        f[1,1], palette = (color = palette(ColorSchemes.hsv, 3),), ## or replace ncol with Integer 

                        title = "Acr Phosphine Oxide + Anthracene UvVis",
                        xlabel = "Wavelength (nm)", 
                        ylabel = "A (OD)", 
                        xminorticks = IntervalsBetween(2)
                    )
                
                width = 3

                ### Custom Axis limits ##############################################################################
                    # ax.xticks= 1:2:17
                    CairoMakie.xlims!(350,600)
                    CairoMakie.ylims!(-0.02,1.5)

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
    MesAcr_PhosphineOxide_AntrhaceneUvVis_Fig()
    save("JULIA/Outputs/MesAcr_PhosphineOxide_AntrhaceneUvVis_Fig.png",MesAcr_PhosphineOxide_AntrhaceneUvVis_Fig())


function MesAcr_PhosphineOxide_AntrhaceneUvVis_FigTA_Fig(f = Figure())

            TA_Spectra = MesAcr_PhosphineOxide_Antrhacene_TASpectra
            TA_Kinetics = [
                (MesAcr_PhosphineOxide_Antrhacene_TAKinetics_420nm[1], "420 nm"),
                (MesAcr_PhosphineOxide_Antrhacene_TAKinetics_510nm[1], "510 nm"), 
                (MesAcr_PhosphineOxide_Antrhacene_TAKinetics_710nm[1], "710 nm"), 
                (MesAcr_PhosphineOxide_Antrhacene_TAKinetics_895nm[1], "895 nm")
            ]

            ### First Panel TA Spectra  ##############################################################################
                ax =    CairoMakie.Axis(
                        f[1,1], palette = (color = palette(ColorSchemes.inferno, 16),), ## or replace ncol with Integer 

                        title = "Acr Phosphine Oxide + Anthracene TA Spectra",
                        xlabel = "Wavelength (nm)", 
                        ylabel = "ΔA (OD)", 
                        xminorticks = IntervalsBetween(2)
                    )
                width = 3

                ### Custom Axis limits ##############################################################################
                    # ax.xticks= 1:2:17
                    # CairoMakie.xlims!(375,910)
                    CairoMakie.ylims!(-0.02,nothing)

                ### Plot Data ##############################################################################

                    for i in [2,3,4,5,6,7,8,10,15, 20, 25, 28, 30, 32, 33] #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
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
                    CairoMakie.vspan!(431, 441; ymin = 0.0, ymax = 1.0, color = :lightblue)
                    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                    rowsize!(f.layout,1, 400)
                
            # ### Second Panel TA KINETICS ##############################################################################

                ax2 =    CairoMakie.Axis(
                        f[1,2], palette = (color = palette(ColorSchemes.hsv, 5),),
                        title = "Acr Phosphine Oxide + Anthracene TA Kinetics",
                        xlabel = "Time (ns)", 
                        ylabel = "ΔA (OD)", 
                        # xscale = Makie.Symlog10(10.0), 
                        # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                        # xminorticks = IntervalsBetween(2)
                )
                width = 3

                ### Custom Axis limits ##############################################################################
                    # ax.xticks= 1:2:17
                    CairoMakie.xlims!(-10, 1500)
                    CairoMakie.ylims!(-0.01,0.03)

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
MesAcr_PhosphineOxide_AntrhaceneUvVis_FigTA_Fig()
save("JULIA/Outputs/ MesAcr_PhosphineOxide_AntrhaceneUvVis_FigTA_Fig.png",MesAcr_PhosphineOxide_AntrhaceneUvVis_FigTA_Fig())


function MesAcr_PhosphineOxide_Antrhacene_TA_Comparison(f = Figure())

    ### First Panel TA Spectra  ##############################################################################
        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.hsv, 4),), ## or replace ncol with Integer 

                title = "Acr Phosphine Oxide + Anthracene Comparison",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(-0.02,nothing)

        ### Plot Data ##############################################################################


            lines!(ax, MesAcr_PhosphineOxide_TASpectra[:,1],MesAcr_PhosphineOxide_TASpectra[:,12],linewidth = width,label = "AcrPO $(names(MesAcr_PhosphineOxide_TASpectra)[14])")
            # lines!(ax, MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,1],MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,15],linewidth = width,label = " $(names(MesAcr_PhosphineOxide_Antrhacene_TASpectra)[15])")
            lines!(ax, MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,1],MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,17],linewidth = width,label = "AcrPO + Anthracene $(names(MesAcr_PhosphineOxide_Antrhacene_TASpectra)[17])")
            # lines!(ax, MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,1],MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,25],linewidth = width,label = " $(names(MesAcr_PhosphineOxide_Antrhacene_TASpectra)[25])")
            # lines!(ax, MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,1],MesAcr_PhosphineOxide_Antrhacene_TASpectra[:,30],linewidth = width,label = " $(names(MesAcr_PhosphineOxide_Antrhacene_TASpectra)[30])")

            lines!(ax, Anthracene_TA_Spectra[:,1],Anthracene_TA_Spectra[:,22],linewidth = width,label = "Anthracene Triplet")
            # lines!(ax, MesAcr_PhosphineOxide_Antrhacene_TAKinetics_510nm[1][:,1],MesAcr_PhosphineOxide_Antrhacene_TAKinetics_510nm[1][:,2],linewidth = width,label = "AcrPO + Anthracene 510 nm Kinetics")



            axislegend(
                position = :rt, 
                nbanks = 1, 
                framevisible = false, 
                labelsize = 15
            )

            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(431, 441; ymin = 0.0, ymax = 1.0, color = :lightblue)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)

    resize_to_layout!(f)
    f
end
MesAcr_PhosphineOxide_Antrhacene_TA_Comparison()
save("JULIA/Outputs/ MesAcr_PhosphineOxide_Antrhacene_TA_Comparison.png",MesAcr_PhosphineOxide_Antrhacene_TA_Comparison())

## Need to run 2024_09_18/2024_09_18_Anthracene_TA.JL before this section
function MesAcr_PhosphineOxide_Antrhacene_PCA_Comparison(f = Figure())

    ### First Panel TA Spectra  ##############################################################################
        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.hsv, 3),), ## or replace ncol with Integer 

                title = "Acr Phosphine Oxide + Anthracene PCA Comparison",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(-0.02,nothing)

        ### Plot Data ##############################################################################


            lines!(ax, PCA_Spectra[:,1],PCA_Spectra[:,4],linewidth = width,label = "PCA 3 of AcrPO + Anthracene")
            lines!(ax, Anthracene_PCA_Spectra[:,1],Anthracene_PCA_Spectra[:,2],linewidth = width,label = "PCA 2 of Anthracene")


            axislegend(
                position = :rt, 
                nbanks = 1, 
                framevisible = false, 
                labelsize = 15
            )

            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(431, 441; ymin = 0.0, ymax = 1.0, color = :lightblue)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)

    resize_to_layout!(f)
    f
end
MesAcr_PhosphineOxide_Antrhacene_PCA_Comparison()
save("JULIA/Outputs/ MesAcr_PhosphineOxide_Antrhacene_PCA_Comparison.png",MesAcr_PhosphineOxide_Antrhacene_PCA_Comparison())


####### PCA AcrPO + Anthracene ##########################################################################################
    TA_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPO_Anthracene_Carpet.csv", header = false, DataFrame) # TA carpet data with wavelength in rows and time points as columns

    # Extract the absorbance data (excluding the first row and column)
    absorbance_data = Matrix(TA_data[2:end, 2:end])

    Wavelengths = TA_data[2:end,1]
    Delay_times = transpose(Matrix(TA_data))[2:end,1]

    function pca_ta(data::Matrix{Float64})
        # Assume data is a matrix where each column is a spectrum at a different time point
        # and each row corresponds to a wavelength

        # Center the data
        data_centered = data .- mean(data, dims=2)

        # Perform SVD
        U, S, V = svd(data_centered)
        ## U = principle components
        ## S = Singular Values
        ## V = time dependnet scores

        # Calculate explained variance ratio
        explained_variance_ratio = S.^2 / sum(S.^2)

        # Return principal components, singular values, and explained variance ratio
        return U, S, V, explained_variance_ratio
    end

    U, S, V, explained_variance_ratio = pca_ta(absorbance_data)
    PCA_Spectra = DataFrame(Wavelength_nm = Wavelengths, PC_1 = U[:,1], PC_2 = U[:,2], PC_3 = U[:,3])
    PCA_Kinetics = DataFrame(time_µs = Delay_times,Kinetic_PC1 = V[:,1],Kinetic_PC2 = V[:,2],Kinetic_PC3 = V[:,3])

    CSV.write("JULIA/Outputs/PCA_Spectra.csv",PCA_Spectra)
    PCA_1_Kinetics = BiFitIRF_new(PCA_Kinetics,2,-10.0, 200,[0.0017, 0.001, 0.0188, 0.3503, 0.0012, 0.0018, -0.0265], 10,10,false,400)
    PCA_2_Kinetics = BiFitIRF_new(PCA_Kinetics,3,-10.0, 200,[0.0095, 0.0011, 0.0314, 0.3347, 0.0071, 0.1444, -0.043], 10,10,false,400)
    PCA_3_Kinetics = BiFitIRF_new(PCA_Kinetics,4,-10.0, 200,[0.0009, 0.0152, 0.1144, 0.3534, -0.046, 0.1058, -0.0404], 10,10,false,400)


    function PCA_Fig(f = Figure())

        U, S, V, explained_variance_ratio = pca_ta(absorbance_data)
        Wavelengths = TA_data[2:end,1]
        Delay_times = transpose(Matrix(TA_data))[2:end,1]
        components = collect(1:3)
        ### 1st Panel Principal Components  ##############################################################################
            ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.hsv, 4),), ## or replace ncol with Integer
                    title = "Principal Components",
                    xlabel = "Wavelength (nm)", 
                    ylabel = "Amplitude", 
                    xminorticks = IntervalsBetween(2)
                )
            width = 3 

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(-0.02,nothing)

            ### Plot Data ##############################################################################

            for i in components #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, Wavelengths,U[:,i],linewidth = width,label = "PCA $i")
            end

                axislegend(position = :rt,nbanks = 3, framevisible = false, labelsize = 15)
                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                vlines!(ax,436, linewidth = 1,color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            
        # ### 2nd Panel SVD Kinetics ##############################################################################
            ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "Time-dependent Scores",
                    xlabel = "Time µs)", 
                    ylabel = "Score", 
                    # xscale = Makie.Symlog10(10.0), 
                    # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                    # xminorticks = IntervalsBetween(2)
            )
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(-10, 1500)
                # CairoMakie.ylims!(-0.01,0.03)

            ### Plot Kinetic Data ##############################################################################

            for i in components #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax2, Delay_times,V[:,i],linewidth = width,label = "PCA $i")
            end
                axislegend( position = :rt, nbanks = 3, framevisible = false,labelsize = 15)
                hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

        # ### 3rd Panel: Explained Varience ##############################################################################

            ax3 =    CairoMakie.Axis(f[2,1], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "Explained varience",
                xlabel = "Principal Component", 
                ylabel = "Cumulative Explained Varience Ratio", 
                # xscale = Makie.Symlog10(10.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                # xminorticks = IntervalsBetween(2)
            )
            width = 3

            ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(0, 10)
            # CairoMakie.ylims!(-0.01,0.03)

            ### Plot Kinetic Data ##############################################################################

            CairoMakie.scatter!(ax3,1:length(explained_variance_ratio), cumsum(explained_variance_ratio),markersize = 10, color = :black)
        

        # ### 4th Panel: Singular Values ##############################################################################

            ax3 =    CairoMakie.Axis(f[2,2], palette = (color = palette(ColorSchemes.hsv, 4),),
            title = "Singular Values",
            xlabel = "Principal Component", 
            ylabel = "Singular Value", 
            # xscale = Makie.Symlog10(10.0), 
            # xticks = [-5, 0, 5, 10,50, 100, 1000], 
            # xminorticks = IntervalsBetween(2)
                )
            width = 3

            ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(0, 10)
            # CairoMakie.ylims!(-0.01,0.03)

            ### Plot Kinetic Data ##############################################################################

            CairoMakie.scatter!(ax3,1:length(S), S,markersize = 10, color = :black)




        Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
        Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
        Label(f[2,1,TopLeft()],"(c)", fontsize = 20)
        Label(f[2,2,TopLeft()],"(d)", fontsize = 20)
        colgap!(f.layout,1,20)
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        rowsize!(f.layout,2, 400)

        resize_to_layout!(f)
        f
    end


PCA_Fig()
save("JULIA/Outputs/ PCA_AcrPO_Anthracene_TA_Fig.png",PCA_Fig())


####### PCA AcrPO ##########################################################################################
Acr_PO_TA_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPO_Carpet.csv", header = false, DataFrame) # TA carpet data with wavelength in rows and time points as columns

# Extract the absorbance data (excluding the first row and column)
Acr_PO_absorbance_data = Matrix(Acr_PO_TA_data[2:end, 2:end])

Wavelengths = Acr_PO_TA_data[2:end,1]
Delay_times = transpose(Matrix(Acr_PO_TA_data))[2:end,1]

function pca_ta(data::Matrix{Float64})
    # Assume data is a matrix where each column is a spectrum at a different time point
    # and each row corresponds to a wavelength

    # Center the data
    data_centered = data .- mean(data, dims=2)

    # Perform SVD
    U, S, V = svd(data_centered)
    ## U = principle components
    ## S = Singular Values
    ## V = time dependnet scores

    # Calculate explained variance ratio
    explained_variance_ratio = S.^2 / sum(S.^2)

    # Return principal components, singular values, and explained variance ratio
    return U, S, V, explained_variance_ratio
end

U, S, V, explained_variance_ratio = pca_ta(Acr_PO_absorbance_data) .* (-1)
PCA_Spectra = DataFrame(Wavelength_nm = Wavelengths, PC_1 = U[:,1], PC_2 = U[:,2], PC_3 = U[:,3])
PCA_Kinetics = DataFrame(time_µs = Delay_times,Kinetic_PC1 = V[:,1],Kinetic_PC2 = V[:,2],Kinetic_PC3 = V[:,2])

PCA_1_Kinetics = BiFitIRF_new(PCA_Kinetics,2,-10.0, 200,[0.002, 0.001, 0.0223, 0.4238, 0.0012, 0.0018, -0.0279], 10,10,false,400)
PCA_2_Kinetics = BiFitIRF_new(PCA_Kinetics,3,-10.0, 200,[0.0092, 0.0009, 0.0663, 0.405, -0.0057, 0.1747, -0.0545], 10,10,false,400)

PCA_1_Kinetics = MonoFitIRF_New(PCA_Kinetics,2,-10.0, 200,[-0.001, 0.001, 0.02, 0.254, -0.026], Inf,Inf,false,400)
PCA_2_Kinetics = MonoFitIRF_New(PCA_Kinetics,3,-20.0, 0.54,[0.01, 0.001, 0.061, 0.445, -0.055], Inf,Inf,false,400)


function Acr_PO_PCA_Fig(f = Figure())

    U, S, V, explained_variance_ratio = pca_ta(Acr_PO_absorbance_data) .* (-1)
    Wavelengths = Acr_PO_TA_data[2:end,1]
    Delay_times = transpose(Matrix(Acr_PO_TA_data))[2:end,1]
    components = collect(1:2)
    ### 1st Panel Principal Components  ##############################################################################
        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.hsv, 4),), ## or replace ncol with Integer
                title = "Principal Components",
                xlabel = "Wavelength (nm)", 
                ylabel = "Amplitude", 
                xminorticks = IntervalsBetween(2)
            )
        width = 3 

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(-0.02,nothing)

        ### Plot Data ##############################################################################

        for i in components #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
            lines!(ax, Wavelengths,U[:,i],linewidth = width,label = "PCA $i")
        end

            axislegend(position = :rt,nbanks = 3, framevisible = false, labelsize = 15)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            vlines!(ax,436, linewidth = 1,color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    # ### 2nd Panel SVD Kinetics ##############################################################################
        ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "Time-dependent Scores",
                xlabel = "Time (ns)", 
                ylabel = "Score", 
                # xscale = Makie.Symlog10(10.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                # xminorticks = IntervalsBetween(2)
        )
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(-10, 1500)
            # CairoMakie.ylims!(-0.01,0.03)

        ### Plot Kinetic Data ##############################################################################

        for i in components #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
            lines!(ax2, Delay_times,V[:,i],linewidth = width,label = "PCA $i")
        end
            axislegend( position = :rt, nbanks = 3, framevisible = false,labelsize = 15)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

    # ### 3rd Panel: Explained Varience ##############################################################################

        ax3 =    CairoMakie.Axis(f[2,1], palette = (color = palette(ColorSchemes.hsv, 4),),
            title = "Explained varience",
            xlabel = "Principal Component", 
            ylabel = "Cumulative Explained Varience Ratio", 
            # xscale = Makie.Symlog10(10.0), 
            # xticks = [-5, 0, 5, 10,50, 100, 1000], 
            # xminorticks = IntervalsBetween(2)
        )
        width = 3

        ### Custom Axis limits ##############################################################################
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(0, 10)
        # CairoMakie.ylims!(-0.01,0.03)

        ### Plot Kinetic Data ##############################################################################

        CairoMakie.scatter!(ax3,1:length(explained_variance_ratio), cumsum(explained_variance_ratio),markersize = 10, color = :black)
    

    # ### 4th Panel: Singular Values ##############################################################################

        ax3 =    CairoMakie.Axis(f[2,2], palette = (color = palette(ColorSchemes.hsv, 4),),
        title = "Singular Values",
        xlabel = "Principal Component", 
        ylabel = "Singular Value", 
        # xscale = Makie.Symlog10(10.0), 
        # xticks = [-5, 0, 5, 10,50, 100, 1000], 
        # xminorticks = IntervalsBetween(2)
            )
        width = 3

        ### Custom Axis limits ##############################################################################
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(0, 10)
        # CairoMakie.ylims!(-0.01,0.03)

        ### Plot Kinetic Data ##############################################################################

        CairoMakie.scatter!(ax3,1:length(S), S,markersize = 10, color = :black)




    Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
    Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
    Label(f[2,1,TopLeft()],"(c)", fontsize = 20)
    Label(f[2,2,TopLeft()],"(d)", fontsize = 20)
    colgap!(f.layout,1,20)
    colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    rowsize!(f.layout,2, 400)

    resize_to_layout!(f)
    f
end


Acr_PO_PCA_Fig()
save("JULIA/Outputs/ Acr_PO_PCA_Fig.png",Acr_PO_PCA_Fig())

AcrPO_PCA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPhosphineOxide_Anthracene_MeCN_TA_436nmPump_Sample2_400points -bg -chirp PrincipalSpectra.txt", header = true, DataFrame)
AcrPO_PCA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/AcrPhosphineOxide_Anthracene_MeCN_TA_436nmPump_Sample2_400points -bg -chirp PrincipalKinetics.csv", header = true, DataFrame)

function AcrPO_PCA_Fig_2(f = Figure())

    TA_Spectra = AcrPO_PCA_Spectra

    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.hsv, ncol(TA_Spectra)),), ## or replace ncol with Integer 
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

            for i in 2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i] .* (-1),linewidth = width,label = " PCA $(i)")
            
            end

            axislegend("Principal Component", position = :rt,nbanks = 1,framevisible = false,labelsize = 15)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    #### Second Panel TA KINETICS ##############################################################################

        ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
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


            lines!(ax2, AcrPO_PCA_Kinetics[:, 1], AcrPO_PCA_Kinetics[:, 2],linewidth = 3,label = "PCA 1") # Kinetic Trace
            lines!(ax2, AcrPO_PCA_Kinetics[:, 1], AcrPO_PCA_Kinetics[:, 3], linewidth = 2, label = "PCA 2") # Fit 
            lines!(ax2, AcrPO_PCA_Kinetics[:, 1], AcrPO_PCA_Kinetics[:, 4], linewidth = 2, label = "PCA 3") # Fit 


            axislegend("Wavelength", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
AcrPO_PCA_Fig_2()
save("JULIA/Outputs/AcrPO_PCA_Fig_2.pdf",AcrPO_PCA_Fig_2())
