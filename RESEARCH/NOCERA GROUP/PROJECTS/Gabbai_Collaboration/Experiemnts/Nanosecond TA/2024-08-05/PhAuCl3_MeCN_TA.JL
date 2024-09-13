## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")



Plots.plot(PhAuCl3_Photolyzed_TA_Spectra[:,1],MaxNorm(PhAuCl3_Photolyzed_TA_Spectra[:,8]))
Plots.plot!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],MaxNorm(iPrAuCl3_Photolyzed_TA_Spectra[:,3]))

###### PhAuCl3 Sample 1  ##############################################################################################
    PhAuCl3_2mm_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/PhAuCl3_2mm_TA_Spectra.csv", header = true, DataFrame)
    PhAuCl3_2mm_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/PhAuCl3_2mm_TA_Kinetics.csv", header = true, DataFrame)

    PhAuCl3_2mm_TA_Kinetics_390nm = MonoFitIRF_New(PhAuCl3_2mm_TA_Kinetics,2,0.0,1000,[0.1, 0.65, 0.249, 16.791, -0.001],100,Inf,false,200)

    # PhAuCl3_2mm_TA_Kinetics_418nm = MonoFitIRF_New(PhAuCl3_2mm_TA_Kinetics,3,0.0,1000,[0.1, 0.65, 0.249, 16.791, -0.001],100,Inf,false,200)
    PhAuCl3_2mm_TA_Kinetics_418nm = BiFitIRF_new(PhAuCl3_2mm_TA_Kinetics,3,0.0,1000,[0.1, 0.65, 0.0975, 13.3662, -2.4516, 401.6066, -0.0016],100,Inf,false,200)

    # PhAuCl3_2mm_TA_Kinetics_436nm = MonoFitIRF_New(PhAuCl3_2mm_TA_Kinetics,4,0.0,1000,[0.1, 0.65, -4.115, 656.558, -0.0],100,Inf,false,200)
    PhAuCl3_2mm_TA_Kinetics_436nm = BiFitIRF_new(PhAuCl3_2mm_TA_Kinetics,4,0.0,1000,[0.1, 0.65, 0.0716, 23.9207, -2.5448, 417.5258, -0.0012],100,Inf,false,200)

    # PhAuCl3_2mm_TA_Kinetics_496nm = MonoFitIRF_New(PhAuCl3_2mm_TA_Kinetics,5,0.0,1000,[0.1, 0.65, -4.115, 656.558, -0.0],100,Inf,false,200)
    
    PhAuCl3_2mm_TA_Kinetics_496nm = BiFitIRF_new(PhAuCl3_2mm_TA_Kinetics,5,0.0,1000,[0.1, 0.65, 0.0321, 14.3195, 5.8572, 357.267, 0.001],100,Inf,false,200)

    # PhAuCl3_2mm_TA_Kinetics_890nm = MonoFitIRF_New(PhAuCl3_2mm_TA_Kinetics,7,0.0,1000,[0.1, 0.65, -4.115, 656.558, -0.0],100,Inf,false,200)
    PhAuCl3_2mm_TA_Kinetics_890nm = BiFitIRF_new(PhAuCl3_2mm_TA_Kinetics,7,0.0,1000,[0.1, 0.65, 0.0341, 3.9793, 3.3033, 365.2296, 0.0006],100,Inf,false,200)

    function PhAuCl3_2mm_Spectra(f = Figure())
        width = 2.5
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "NS TA of PhAu(III)Cl₃ in MeCN",
                palette = (color = palette(ColorSchemes.inferno, 14),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            CairoMakie.ylims!(nothing,0.026)
            
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,2],linewidth = width,label = "2 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,3],linewidth = width,label = "5 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,4],linewidth = width,label = "10 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,5],linewidth = width,label = "20 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,6],linewidth = width,label = "30 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,7],linewidth = width,label = "50 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,8],linewidth = width,label = "100 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,9],linewidth = width,label = "150 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,10],linewidth = width,label = "200 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,11],linewidth = width,label = "300 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,12],linewidth = width,label = "450 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,13],linewidth = width,label = "700 ns")
            lines!(PhAuCl3_2mm_TA_Spectra[:,1],PhAuCl3_2mm_TA_Spectra[:,14],linewidth = width,label = "1000 ns")


            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)

                # Second Panel of figure 
                ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of PhAu(III)Cl₃ in MeCN", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA (mOD)")

                lines!(PhAuCl3_2mm_TA_Kinetics_418nm[:,1],PhAuCl3_2mm_TA_Kinetics_418nm[:,2],linewidth = width,label = "418 nm")
                lines!(PhAuCl3_2mm_TA_Kinetics_418nm[:,1],PhAuCl3_2mm_TA_Kinetics_418nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(PhAuCl3_2mm_TA_Kinetics_496nm[:,1],PhAuCl3_2mm_TA_Kinetics_496nm[:,2],linewidth = width,label = "496 nm")
                lines!(PhAuCl3_2mm_TA_Kinetics_496nm[:,1],PhAuCl3_2mm_TA_Kinetics_496nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(PhAuCl3_2mm_TA_Kinetics_890nm[:,1],PhAuCl3_2mm_TA_Kinetics_890nm[:,2],linewidth = width,label = "890 nm")
                lines!(PhAuCl3_2mm_TA_Kinetics_890nm[:,1],PhAuCl3_2mm_TA_Kinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)


            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            rowsize!(f.layout,1, 400)
            # rowsize!(f.layout,2, 400)

            colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    PhAuCl3_2mm_Spectra()
    save("JULIA/Outputs/PhAuCl3_2mm_Spectra.pdf",PhAuCl3_2mm_Spectra())


###### PhAuCl3 Sample 2 Photolyzed  ##############################################################################################

    PhAuCl3_Photolyzed_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/PhAuCl3 in MeCN (2mm)_sample2_PostPhotolysis_TA_Spectra.csv", header = true, DataFrame)
    PhAuCl3_Photolyzed_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/PhAuCl3 in MeCN (2mm)_sample2_PostPhotolysis_TA_Kinetics.csv", header = true, DataFrame)

    # PhAuCl3_2mm_TA_Kinetics_390nm = BiFitIRF_new(PhAuCl3_Photolyzed_TA_Kinetics,2,0.0,2000,[0.1, 0.65, 0.2205, 10.0743, -44.9313, 23790.6048, 0.0009],100,Inf,false,200)
    PhAuCl3_Photolyzed_TA_Kinetics_435nm = BiFitIRF_new(PhAuCl3_Photolyzed_TA_Kinetics,3,-10.0,2000,[0.1164, 0.4375, 0.0587, 18.8678, -5.9085, 448.0558, -0.0007],100,Inf,false,200)
    PhAuCl3_Photlyzed_TA_Kinetics_490nm = BiFitIRF_new(PhAuCl3_Photolyzed_TA_Kinetics,4,-10.0,2000,[0.1, 0.65, 0.0856, 20.0389, 13.7137, 454.2208, -0.0003],100,Inf,false,200)
    PhAuCl3_Photlyzed_TA_Kinetics_890nm = BiFitIRF_new(PhAuCl3_Photolyzed_TA_Kinetics,5,-10.0,2000,[0.1, 0.65, 0.07, 8.8991, 7.8039, 447.3663, 0.0002],100,Inf,false,200)

    function PhAuCl3_Photolyzed(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Photolyzed PhAu(III)Cl₃ in MeCN",
                palette = (color = palette(ColorSchemes.inferno, 17),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            CairoMakie.ylims!(nothing,0.04)
            
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,2],linewidth = width,label = "5 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,3],linewidth = width,label = "10 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,4],linewidth = width,label = "15 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,5],linewidth = width,label = "20 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,6],linewidth = width,label = "40 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,7],linewidth = width,label = "75 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,8],linewidth = width,label = "100 ns")
            # lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,9],linewidth = width,label = "150 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,10],linewidth = width,label = "200 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,11],linewidth = width,label = "300 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,12],linewidth = width,label = "400 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,13],linewidth = width,label = "600 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,14],linewidth = width,label = "800 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,15],linewidth = width,label = "1000 ns")
            # lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,16],linewidth = width,label = "1500 ns")
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],PhAuCl3_Photolyzed_TA_Spectra[:,17],linewidth = width,label = "2000 ns")


            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)

        # Second Panel of figure 
            ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of Photolyzed PhAu(III)Cl₃ in MeCN", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA (mOD)")

                lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,2],linewidth = width,label = "435 nm")
                lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,3],linewidth = 2, linestyle = :dash, color = :black)
 
                lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,2],linewidth = width,label = "490 nm")
                lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,2],linewidth = width,label = "890 nm")
                lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)


            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            rowsize!(f.layout,1, 400)
            # rowsize!(f.layout,2, 400)

            colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    PhAuCl3_Photolyzed()
    save("JULIA/Outputs/PhAuCl3_Photolyzed.pdf",PhAuCl3_Photolyzed())

###### iPrAuCl3  ##############################################################################################


    iPrAuCl3_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/iPrAuCl3 in MeCN_TA_Spectra.csv", header = true, DataFrame)

    function iPrAuCl3_Spectra(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of iPrAu(III)Cl₃ in MeCN",
                palette = (color = palette(ColorSchemes.inferno, 6),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)
            
            lines!(iPrAuCl3_TA_Spectra[:,1],iPrAuCl3_TA_Spectra[:,2],linewidth = width,label = "5 ns")
            lines!(iPrAuCl3_TA_Spectra[:,1],iPrAuCl3_TA_Spectra[:,3],linewidth = width,label = "10 ns")
            lines!(iPrAuCl3_TA_Spectra[:,1],iPrAuCl3_TA_Spectra[:,4],linewidth = width,label = "15 ns")
            lines!(iPrAuCl3_TA_Spectra[:,1],iPrAuCl3_TA_Spectra[:,5],linewidth = width,label = "20 ns")
            lines!(iPrAuCl3_TA_Spectra[:,1],iPrAuCl3_TA_Spectra[:,6],linewidth = width,label = "40 ns")



            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)

        # Second Panel of figure 
            # ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of Photolyzed PhAu(III)Cl₃ in MeCN", 
            #     palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
            #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
            #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
            #     xlabel = "Time (ns)", 
            #     ylabel = "ΔA (mOD)")

            #     lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,2],linewidth = width,label = "435 nm")
            #     lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

            #     lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,2],linewidth = width,label = "490 nm")
            #     lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

            #     lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,2],linewidth = width,label = "890 nm")
            #     lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

            #     axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)


            # Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            # Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            # hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            # rowsize!(f.layout,1, 400)
            # # rowsize!(f.layout,2, 400)

            # colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        # colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    iPrAuCl3_Spectra()
    save("JULIA/Outputs/iPrAuCl3_Spectra.pdf",iPrAuCl3_Spectra())

###### iPrAuCl3 Sample Photolyzed  ##############################################################################################

    iPrAuCl3_Photolyzed_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/iPrAuCl3 in MeCN (2mm) photolyzed_TA_Spectra.csv", header = true, DataFrame)
    iprAuCl3_Photolyzed_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024-08-05/iPrAuCl3 in MeCN (2mm) photolyzed_TA_Kinetics.csv", header = true, DataFrame)

    # PhAuCl3_2mm_TA_Kinetics_390nm = BiFitIRF_new(PhAuCl3_Photolyzed_TA_Kinetics,2,0.0,2000,[0.1, 0.65, 0.2205, 10.0743, -44.9313, 23790.6048, 0.0009],100,Inf,false,200)
    iPrAuCl3_Photolyzed_TA_Kinetics_415nm = BiFitIRF_new(iprAuCl3_Photolyzed_TA_Kinetics,2,-10.0,50000,[0.1, 0.65, -23.9382, 1886.0605, -50.1771, 5388.9329, -0.0001],10,10,false,200)
    iPrAuCl3_Photlyzed_TA_Kinetics_490nm = BiFitIRF_new(iprAuCl3_Photolyzed_TA_Kinetics,3,-10.0,50000,[0.1, 0.65, 96.9478, 2594.1662, 116.3791, 8427.9195, -0.0001],10,10,false,200)
    iPrAuCl3_Photlyzed_TA_Kinetics_890nm = BiFitIRF_new(iprAuCl3_Photolyzed_TA_Kinetics,5,-10.0,50000,[0.1, 0.65, 32.9046, 2241.3455, 143.7315, 9459.9812, -0.0002],10,10,false,200)

    function iPrAuCl3_Photolyzed(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Photolyzed iPrAu(III)Cl₃ in MeCN",
                palette = (color = palette(ColorSchemes.inferno, 13),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            CairoMakie.ylims!(nothing,0.06)
            
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,2],linewidth = width,label = "10 ns")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,3],linewidth = width,label = "100 ns")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,4],linewidth = width,label = "400 ns")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,5],linewidth = width,label = "800 ns")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,6],linewidth = width,label = "1.5 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,7],linewidth = width,label = "2.5 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,8],linewidth = width,label = "4 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,9],linewidth = width,label = "6 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,10],linewidth = width,label = "10 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,11],linewidth = width,label = "20 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,12],linewidth = width,label = "30 µs")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],iPrAuCl3_Photolyzed_TA_Spectra[:,13],linewidth = width,label = "50 µs")


            CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)

        # Second Panel of figure 
            ax2 = CairoMakie.Axis(f[1,2], title = "TA Kinetics of Photolyzed iPrAu(III)Cl₃ in MeCN", 
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (µs)", 
                ylabel = "ΔA (mOD)")

                lines!(iPrAuCl3_Photolyzed_TA_Kinetics_415nm[:,1] ./ 1000,iPrAuCl3_Photolyzed_TA_Kinetics_415nm[:,2],linewidth = width,label = "415 nm")
                lines!(iPrAuCl3_Photolyzed_TA_Kinetics_415nm[:,1]./ 1000,iPrAuCl3_Photolyzed_TA_Kinetics_415nm[:,3],linewidth = 2, linestyle = :dash, color = :black)
 
                lines!(iPrAuCl3_Photlyzed_TA_Kinetics_490nm[:,1]./ 1000,iPrAuCl3_Photlyzed_TA_Kinetics_490nm[:,2],linewidth = width,label = "490 nm")
                lines!(iPrAuCl3_Photlyzed_TA_Kinetics_490nm[:,1]./ 1000,iPrAuCl3_Photlyzed_TA_Kinetics_490nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(iPrAuCl3_Photlyzed_TA_Kinetics_890nm[:,1]./ 1000,iPrAuCl3_Photlyzed_TA_Kinetics_890nm[:,2],linewidth = width,label = "890 nm")
                lines!(iPrAuCl3_Photlyzed_TA_Kinetics_890nm[:,1]./ 1000,iPrAuCl3_Photlyzed_TA_Kinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)

                # text!(20,0.025,text = "τ₁ = 1757 ns, τ₂ = 5124", color = :red, fontsize = 20)


            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)

            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


            rowsize!(f.layout,1, 400)
            # rowsize!(f.layout,2, 400)

            colgap!(f.layout,1,20)

        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    iPrAuCl3_Photolyzed()
    save("JULIA/Outputs/iPrAuCl3_Photolyzed.pdf",iPrAuCl3_Photolyzed())


    function Ph_iPr_Comparisons(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "Normalized TA Spectra of Photolyzed PhAuCl₃ and iPrAu(III)Cl₃",
                palette = (color = palette(ColorSchemes.inferno, 3),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (A.U.)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.06)
            
            lines!(PhAuCl3_Photolyzed_TA_Spectra[:,1],MaxNorm(PhAuCl3_Photolyzed_TA_Spectra[:,8]),linewidth = width,label = "PhAuCl₃")
            lines!(iPrAuCl3_Photolyzed_TA_Spectra[:,1],MaxNorm(iPrAuCl3_Photolyzed_TA_Spectra[:,3]),linewidth = width,label = "iPrAuCl₃")



            CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend( position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)



        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        # colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
    Ph_iPr_Comparisons()
    save("JULIA/Outputs/Ph_iPr_Comparisons.pdf",Ph_iPr_Comparisons())

