## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")


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

    function  Pd_TFP_TA(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Pd TFP",
                palette = (color = palette(ColorSchemes.inferno, ncol(Pd_TFP_Spectra)+1),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(nothing,700)
            # CairoMakie.ylims!(nothing,0.04)
            
            for i in 2:ncol(Pd_TFP_Spectra)
                lines!(ax, Pd_TFP_Spectra[:,1],Pd_TFP_Spectra[:,i],linewidth = width,label = " $(names(Pd_TFP_Spectra)[i])")
            
            end



            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Delay Time", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)


            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
            rowsize!(f.layout,1, 400)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        # Second Panel of figure 
            ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of Pd TFP", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA (mOD)")

                lines!(Pd_TFP_Kinetics_415nm[:,1],Pd_TFP_Kinetics_415nm[:,2],linewidth = width,label = "415 nm")
                lines!(Pd_TFP_Kinetics_415nm[:,1],Pd_TFP_Kinetics_415nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(Pd_TFP_Kinetics_456nm[:,1],Pd_TFP_Kinetics_456nm[:,2],linewidth = width,label = "456 nm")
                lines!(Pd_TFP_Kinetics_456nm[:,1],Pd_TFP_Kinetics_456nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(Pd_TFP_Kinetics_524nm[:,1],Pd_TFP_Kinetics_524nm[:,2],linewidth = width,label = "524 nm")
                lines!(Pd_TFP_Kinetics_524nm[:,1],Pd_TFP_Kinetics_524nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
                text!(2000,-0.05,text = "τ = 592.42 ± 24.14 ns", color = :black, fontsize = 20)


            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)



            # rowsize!(f.layout,2, 400)

            colgap!(f.layout,1,20)

        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    Pd_TFP_TA()
    save("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/Pd_TFP_TA.pdf",Pd_TFP_TA())


#### Pd_TFP_CH3 ################################################################################################
    Pd_TFP_CH3_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/PdTFPCH3_TA_Spectra.csv", header = true, DataFrame)
    Pd_TFP_CH3_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/PdTFPCH3_TA_Kinetics.csv", header = true, DataFrame)

    Pd_TFP_CH3_Kinetics_417nm = MonoFitIRF_New(Pd_TFP_CH3_Kinetics,2,-100.0, 4700, [0.001, 0.981, -62.275, 657.05, -0.0], 10,Inf,false,200)
    Pd_TFP_CH3_Kinetics_454nm = MonoFitIRF_New(Pd_TFP_CH3_Kinetics,3,-100.0, 4700, [0.001, 0.981, 17.115, 658.451, -0.0], 10,Inf,false,200)
    Pd_TFP_CH3_Kinetics_523nm = MonoFitIRF_New(Pd_TFP_CH3_Kinetics,4,-100.0, 4700, [0.001, 0.981, -4.933, 654.094, -0.0], 10,Inf,false,200)

    τ_mean_Pd_TFP_CH3 = mean([657.05,658.451,654.094])
    τ_stdev_Pd_TFP_CH3 = std([657.05,658.451,654.094])

    function  Pd_TFP_CH3_TA(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Pd TFP-CH₃",
                palette = (color = palette(ColorSchemes.inferno, ncol(Pd_TFP_CH3_Spectra)+1),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)
            
            for i in 2:ncol(Pd_TFP_CH3_Spectra)
                lines!(ax, Pd_TFP_CH3_Spectra[:,1],Pd_TFP_CH3_Spectra[:,i],linewidth = width,label = " $(names(Pd_TFP_CH3_Spectra)[i])")
            
            end



            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Delay Time", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)


            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


        # Second Panel of figure 
            ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of of Pd TFP-CH₃", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA (mOD)")

                lines!(Pd_TFP_CH3_Kinetics_417nm[:,1],Pd_TFP_CH3_Kinetics_417nm[:,2],linewidth = width,label = "417 nm")
                lines!(Pd_TFP_CH3_Kinetics_417nm[:,1],Pd_TFP_CH3_Kinetics_417nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(Pd_TFP_CH3_Kinetics_454nm[:,1],Pd_TFP_CH3_Kinetics_454nm[:,2],linewidth = width,label = "454 nm")
                lines!(Pd_TFP_CH3_Kinetics_454nm[:,1],Pd_TFP_CH3_Kinetics_454nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(Pd_TFP_CH3_Kinetics_523nm[:,1],Pd_TFP_CH3_Kinetics_523nm[:,2],linewidth = width,label = "523 nm")
                lines!(Pd_TFP_CH3_Kinetics_523nm[:,1],Pd_TFP_CH3_Kinetics_523nm[:,3],linewidth = 2, linestyle = :dash,color = :black)
                axislegend("Wavelength", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
                text!(2000,-0.05,text = "τ = 656.53 ± 2.22 ns", color = :black, fontsize = 20)


            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            rowsize!(f.layout,1, 400)


            # colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    Pd_TFP_CH3_TA()
    save("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/Pd_TFP_CH3_TA.pdf",Pd_TFP_CH3_TA())


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

    function  Zn_Chlorin_b1_TA(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Zn Chlorin band 1",
                palette = (color = palette(ColorSchemes.inferno, ncol(Zn_Chlorin_b1_Spectra)+1),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)
            
            for i in 2:ncol(Zn_Chlorin_b1_Spectra)
                lines!(ax, Zn_Chlorin_b1_Spectra[:,1],Zn_Chlorin_b1_Spectra[:,i],linewidth = width,label = " $(names(Zn_Chlorin_b1_Spectra)[i])")
            
            end



            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Delay Time", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)


            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


        # Second Panel of figure 
            ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of of Zn Chlorin band 1", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA (mOD)")

                lines!(Zn_Chlorin_b1_417nm[:,1],Zn_Chlorin_b1_417nm[:,2],linewidth = width,label = "417 nm")
                lines!(Zn_Chlorin_b1_417nm[:,1],Zn_Chlorin_b1_417nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(Zn_Chlorin_b1_444nm[:,1],Zn_Chlorin_b1_444nm[:,2],linewidth = width,label = "444 nm")
                lines!(Zn_Chlorin_b1_444nm[:,1],Zn_Chlorin_b1_444nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(Zn_Chlorin_b1_613nm[:,1],Zn_Chlorin_b1_613nm[:,2],linewidth = width,label = "613 nm")
                lines!(Zn_Chlorin_b1_613nm[:,1],Zn_Chlorin_b1_613nm[:,3],linewidth = 2, linestyle = :dash,color = :black)
                text!(2000,-0.05,text = "τ = 655.28 ± 8.44 ns", color = :black, fontsize = 20)


            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Wavelength", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

            rowsize!(f.layout,1, 400)


            # colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    Zn_Chlorin_b1_TA()
    save("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/Zn_Chlorin_b1_TA.pdf",Zn_Chlorin_b1_TA())


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

    function  Zn_Chlorin_b2_TA(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Zn Chlorin band 2",
                palette = (color = palette(ColorSchemes.inferno, ncol(Zn_Chlorin_b2_Spectra)+1),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)
            
            for i in 2:ncol(Zn_Chlorin_b2_Spectra)
                lines!(ax, Zn_Chlorin_b2_Spectra[:,1],Zn_Chlorin_b2_Spectra[:,i],linewidth = width,label = " $(names(Zn_Chlorin_b2_Spectra)[i])")
            
            end



            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Delay Time", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)


            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


        # Second Panel of figure 
            ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of Zn Chlorin band 2", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA (mOD)")

                lines!(Zn_Chlorin_b2_417nm[:,1],Zn_Chlorin_b2_417nm[:,2],linewidth = width,label = "417 nm")
                lines!(Zn_Chlorin_b2_417nm[:,1],Zn_Chlorin_b2_417nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(Zn_Chlorin_b2_445nm[:,1],Zn_Chlorin_b2_445nm[:,2],linewidth = width,label = "445 nm")
                lines!(Zn_Chlorin_b2_445nm[:,1],Zn_Chlorin_b2_445nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(Zn_Chlorin_b2_611nm[:,1],Zn_Chlorin_b2_611nm[:,2],linewidth = width,label = "611 nm")
                lines!(Zn_Chlorin_b2_611nm[:,1],Zn_Chlorin_b2_611nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rb,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
                text!(2000,-0.05,text = "τ = 757.94 ± 20.08 ns", color = :black, fontsize = 20)

            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            rowsize!(f.layout,1, 400)

            # colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    Zn_Chlorin_b2_TA()
    save("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/Zn_Chlorin_b2_TA.pdf",Zn_Chlorin_b2_TA())


### Make Figures 
####################################################################################################



