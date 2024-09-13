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

# Open and read MatLab Files (For Magnitude)
File = EnVisionKinetic(2,"path")

## IMPORT DATA
file1 = CSV.read("path", header = true, DataFrame)

# Data Workup
####################################################################################################
Data = MonoFitIRF_New(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017, 0.0], 10,Inf,false,10000)
Data = BiFitIRF_new(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017,10.0, 100, 0.0], 10,Inf,false,10000)

MesAcr_PPh3_And_iPrAuCl_Absorption = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_iPrAuCl and PPh3 UvVis.csv", header = true, DataFrame)

function MesAcr_PPh3_And_iPrAuCl_UVVis(f = Figure())



    ### First Panel ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, 3),), ## or replace ncol with Integer 
                xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.01, 0.01),
                title = "MesAcriPrAu(I)Cl Uv-Vis",
                xlabel = "Wavelength (nm)", 
                ylabel = "A (OD)", 
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(300,600)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

        lines!(MesAcr_PPh3_And_iPrAuCl_Absorption[:,1] ,MesAcr_PPh3_And_iPrAuCl_Absorption[:,2],linewidth = width,label = "Pre TA")
        lines!(MesAcr_PPh3_And_iPrAuCl_Absorption[:,1],MesAcr_PPh3_And_iPrAuCl_Absorption[:,3],linewidth = width, label = "Post TA")



            axislegend( position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        

    # # ### Second Panel ##############################################################################

        ax2 =    CairoMakie.Axis(
                f[1,2], palette = (color = palette(ColorSchemes.inferno, 3),),
                title = "MesAcrP(Ph)₃ Uv-Vis",
                xlabel = "Wavelength (nm)", 
                ylabel = "A (OD)", 
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(350,600)
            CairoMakie.ylims!(nothing,1.6)

        ### Plot Data ##############################################################################

            lines!(MesAcr_PPh3_And_iPrAuCl_Absorption[:,1],MesAcr_PPh3_And_iPrAuCl_Absorption[:,4],linewidth = width,label = "Pre TA ")
            lines!(MesAcr_PPh3_And_iPrAuCl_Absorption[:,1],MesAcr_PPh3_And_iPrAuCl_Absorption[:,5],linewidth = width,label = "Post TA")

            axislegend( position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        linkaxes!(ax,ax2)
    resize_to_layout!(f)
    f
end
MesAcr_PPh3_And_iPrAuCl_UVVis()
save("JULIA/Outputs/MesAcr_PPh3_And_iPrAuCl_UVVis.pdf",MesAcr_PPh3_And_iPrAuCl_UVVis())

MesAcr_Phosphineoxide_DPA_Absorption = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_phosphineOxide and DPA UvVis.csv", header = true, DataFrame)

function MesAcr_Phosphineoxide_DPA_UVVis(f = Figure())



    ### First Panel ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, 3),), ## or replace ncol with Integer 
                xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.01, 0.01),
                title = "MesAcr Phosphine Oxide Uv-Vis",
                xlabel = "Wavelength (nm)", 
                ylabel = "A (OD)", 
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(300,600)
            CairoMakie.ylims!(nothing,0.6)

        ### Plot Data ##############################################################################

        lines!(MesAcr_Phosphineoxide_DPA_Absorption[:,1] ,MesAcr_Phosphineoxide_DPA_Absorption[:,2],linewidth = width,label = "Pre TA")
        lines!(MesAcr_Phosphineoxide_DPA_Absorption[:,1],MesAcr_Phosphineoxide_DPA_Absorption[:,3],linewidth = width, label = "Post TA")



            axislegend( position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        

    # # ### Second Panel ##############################################################################

        ax2 =    CairoMakie.Axis(
                f[1,2], palette = (color = palette(ColorSchemes.inferno, 3),),
                title = "MesAcr Phosphine Oxide + DPA Uv-Vis",
                xlabel = "Wavelength (nm)", 
                ylabel = "A (OD)", 
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(350,600)
            CairoMakie.ylims!(nothing,2.0)

        ### Plot Data ##############################################################################

            lines!(MesAcr_Phosphineoxide_DPA_Absorption[:,1],MesAcr_Phosphineoxide_DPA_Absorption[:,4],linewidth = width,label = "Pre TA ")
            lines!(MesAcr_Phosphineoxide_DPA_Absorption[:,1],MesAcr_Phosphineoxide_DPA_Absorption[:,5],linewidth = width,label = "Post TA")

            axislegend( position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        linkxaxes!(ax,ax2)

    resize_to_layout!(f)
    f
end
MesAcr_Phosphineoxide_DPA_UVVis()
save("JULIA/Outputs/MesAcr_Phosphineoxide_DPA_UVVis.pdf",MesAcr_Phosphineoxide_DPA_UVVis())



##  MesAcr PPh3 
####################################################################################################

    MesAcr_PPh3_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_Phosphine_Ph3_2mm_200ns_nostirr_425nmpump_TA_Spectra.csv", header = true, DataFrame)
    MesAcr_PPh3_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_Phosphine_Ph3_2mm_200ns_nostirr_425nmpump_TA_Kinetics.csv", header = true, DataFrame)

    MesAcr_PPh3_TAKinetics_387nm = MonoFitIRF_New(MesAcr_PPh3_TAKinetics,2,-10.0, 200, [-1.128, 1.049, 0.231, 16.489, 0.001], 10,100,false,200)
    # MesAcr_PPh3_TAKinetics_387nm = BiFitIRF_new(MesAcr_PPh3_TAKinetics,2,-10.0, 200, [-0.024, 0.6, -19.042, 14.1147, 19.256, 14.1416, 0.0012], 10,Inf,false,200)


    MesAcr_PPh3_TAKinetics_472nm = MonoFitIRF_New(MesAcr_PPh3_TAKinetics,3,-10.0, 200, [-0.052, 1.223, 0.068, 9.243, 0.001], 10,10,false,200)
    # MesAcr_PPh3_TAKinetics_490nm = BiFitIRF_new(MesAcr_PPh3_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, 407.538, 29353.8024, 142.6518, 7492.1869, 0.0009], 10,25,false,390)


    MesAcr_iPrAuCl_TAKinetics_549nm = MonoFitIRF_New(MesAcr_PPh3_TAKinetics,4,-10.0,  200, [-0.048, 0.208, -0.115, 30.825, 0.001], 10,10,false,200)
    MesAcr_iPrAuCl_TAKinetics_892nm = MonoFitIRF_New(MesAcr_PPh3_TAKinetics,5,-10.0,  200, [-0.048, 0.206, 0.076, 18.819, 0.0], 10,10,false,200)


    function MesAcr_PPh3_Spectra(f = Figure())

        Data = MesAcr_PPh3_TASpectra

        ### First Panel ##############################################################################

            ax =    CairoMakie.Axis(
                    f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                    xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.01, 0.01),
                    title = "MesAcrPPh₃ TA Spectra",
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
                CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            

        # # ### Second Panel ##############################################################################

            ax2 =    CairoMakie.Axis(
                    f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "MesAcrPPh₃ TA Kinetics",
                    xlabel = "Time (ns)", 
                    ylabel = "ΔA (OD)", 
                )

            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

            ### Plot Data ##############################################################################

                lines!(MesAcr_PPh3_TAKinetics_387nm[:,1] ,MesAcr_PPh3_TAKinetics_387nm[:,2],linewidth = width,label = "387 nm")
                lines!(MesAcr_PPh3_TAKinetics_387nm[:,1],MesAcr_PPh3_TAKinetics_387nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(MesAcr_PPh3_TAKinetics_472nm[:,1],MesAcr_PPh3_TAKinetics_472nm[:,2],linewidth = width,label = "472 nm")
                lines!(MesAcr_PPh3_TAKinetics_472nm[:,1],MesAcr_PPh3_TAKinetics_472nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(MesAcr_iPrAuCl_TAKinetics_549nm[:,1],MesAcr_iPrAuCl_TAKinetics_549nm[:,2],linewidth = width,label = "549 nm")
                lines!(MesAcr_iPrAuCl_TAKinetics_549nm[:,1],MesAcr_iPrAuCl_TAKinetics_549nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
                Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
                colgap!(f.layout,1,20)
                colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
MesAcr_PPh3_Spectra()
save("JULIA/Outputs/MesAcr_PPh3_Spectra.pdf",MesAcr_PPh3_Spectra())

##  MesAcr iPrAu(I)Cl 425 nm Pump
####################################################################################################

    MesAcr_iPrAuCl_425nm_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_iPrAuCl_2mm_425nmPump_TA_Spectra.csv", header = true, DataFrame)
    MesAcr_iPrAuCl_425nm_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_iPrAuCl_2mm_425nmPump_TA_Kinetics.csv", header = true, DataFrame)


    # MesAcr_iPrAuCl_425nm_TAKinetics_382nm = MonoFitIRF_New(MesAcr_iPrAuCl_425nm_TAKinetics,2,-10.0, 2000, [-0.048, 1.2, 0.12, 11.225, 0.002], 10,10,false,1000)
    MesAcr_iPrAuCl_425nm_TAKinetics_382nm = BiFitIRF_new(MesAcr_iPrAuCl_425nm_TAKinetics,2,-10.0, 200000, [-0.048, 0.9619, 0.1259, 7.1711, 95.222, 42721.8561, 0.0002], 10,10,false,10000)


    # MesAcr_iPrAuCl_425nm_TAKinetics_415nm = MonoFitIRF_New(MesAcr_iPrAuCl_425nm_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, -183.999, 13502.95, -0.0], 100,Inf,false,390)
    MesAcr_iPrAuCl_425nm_TAKinetics_415nm = BiFitIRF_new(MesAcr_iPrAuCl_425nm_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, -218.2452, 15833.4779, -63.8465, 3824.9275, 0.0], 1,1,false,1000)


    # MesAcr_iPrAuCl_TAKinetics_490nm = MonoFitIRF_New(MesAcr_iPrAuCl_425nm_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, 470.679, 16170.317, 0.001], Inf,Inf,false,390)
    MesAcr_iPrAuCl_425nm_TAKinetics_485nm = BiFitIRF_new(MesAcr_iPrAuCl_425nm_TAKinetics,4,-10.0, 200000, [-0.024, 0.6, 495.3847, 17564.2006, 155.2832, 4245.3874, 0.0014], 1,1,false,1000)


    MesAcr_iPrAuCl_425nm_TAKinetics_880nm = BiFitIRF_new(MesAcr_iPrAuCl_425nm_TAKinetics,5,-10.0, 200000, [-0.024, 0.6, 228.1282, 19489.8409, 99.6803, 4344.6876, 0.0001], 1,1,false,1000)


    function MesAcr_iPrAuCl_425nm_Spectra(f = Figure())

        Data = MesAcr_iPrAuCl_425nm_TASpectra

        ### First Panel ##############################################################################

            ax =    CairoMakie.Axis(
                    f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                    xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.4, 0.05),
                    title = "MesAcr-iPrAu(I)Cl TA Spectra",
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

                axislegend("Delay Time", position = :rb,nbanks = 4,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(420, 430; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            

        # ### Second Panel ##############################################################################

            ax2 =    CairoMakie.Axis(
                    f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "MesAcr-iPrAu(I)Cl TA Kinetics",
                    xlabel = "Time (µs)", 
                    ylabel = "ΔA (OD)", 
                )
            
            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(-0.02,0.2)
                # CairoMakie.ylims!(nothing,0.04)

            ### Plot Data ##############################################################################

                lines!(MesAcr_iPrAuCl_425nm_TAKinetics_415nm[:,1]./ 1000,MesAcr_iPrAuCl_425nm_TAKinetics_415nm[:,2],linewidth = width,label = "415 nm")
                lines!(MesAcr_iPrAuCl_425nm_TAKinetics_415nm[:,1]./ 1000,MesAcr_iPrAuCl_425nm_TAKinetics_415nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(MesAcr_iPrAuCl_425nm_TAKinetics_485nm[:,1]./ 1000,MesAcr_iPrAuCl_425nm_TAKinetics_485nm[:,2],linewidth = width,label = "485 nm")
                lines!(MesAcr_iPrAuCl_425nm_TAKinetics_485nm[:,1]./ 1000,MesAcr_iPrAuCl_425nm_TAKinetics_485nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(MesAcr_iPrAuCl_425nm_TAKinetics_880nm[:,1]./ 1000,MesAcr_iPrAuCl_425nm_TAKinetics_880nm[:,2],linewidth = width,label = "890 nm")
                lines!(MesAcr_iPrAuCl_425nm_TAKinetics_880nm[:,1]./ 1000,MesAcr_iPrAuCl_425nm_TAKinetics_880nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rb,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                ax2_inset = Axis(f[1,2], width=Relative(0.5), height=Relative(0.4), halign=0.9, valign=0.8, palette = (color = palette(ColorSchemes.cool, 4),),
                xlabel = "Time (ns)", 
                ylabel = "ΔA (OD)", )
                CairoMakie.xlims!(-10,200)

                lines!(ax2_inset, MesAcr_iPrAuCl_425nm_TAKinetics_382nm[:,1],MesAcr_iPrAuCl_425nm_TAKinetics_382nm[:,2],linewidth = width,label = "382 nm")
                lines!(ax2_inset,MesAcr_iPrAuCl_425nm_TAKinetics_382nm[:,1],MesAcr_iPrAuCl_425nm_TAKinetics_382nm[:,3],linewidth = 2, linestyle = :dash, color = :black)
                axislegend(position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)

                Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
                Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
                colgap!(f.layout,1,20)
                colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
MesAcr_iPrAuCl_425nm_Spectra()
save("JULIA/Outputs/MesAcr_iPrAuCl_425nm_Spectra.pdf",MesAcr_iPrAuCl_425nm_Spectra())


##  MesAcr iPrAu(I)Cl 
####################################################################################################

    MesAcr_iPrAuCl_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_iPrAuCl_2mm_375nmPump_TA_Spectra.csv", header = true, DataFrame)
    MesAcr_iPrAuCl_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/MesAcr_iPrAuCl_2mm_375nmPump_TA_Kinetics.csv", header = true, DataFrame)

    MesAcr_iPrAuCl_TAKinetics_420nm = MonoFitIRF_New(MesAcr_iPrAuCl_TAKinetics,4,-10.0, 200000, [-0.024, 0.6, -183.999, 13502.95, -0.0], 100,Inf,false,390)
    MesAcr_iPrAuCl_TAKinetics_420nm = BiFitIRF_new(MesAcr_iPrAuCl_TAKinetics,4,-10.0, 200000, [-0.024, 0.6, -127.8041, 36078.4206, -98.649, 8773.7352, 0.0002], 10,100,false,398)


    # MesAcr_iPrAuCl_TAKinetics_490nm = MonoFitIRF_New(MesAcr_iPrAuCl_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, 470.679, 16170.317, 0.001], Inf,Inf,false,390)
    MesAcr_iPrAuCl_TAKinetics_490nm = BiFitIRF_new(MesAcr_iPrAuCl_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, 399.8178, 29993.1574, 150.7399, 7950.9206, 0.0009], 10,10,false,398)


    MesAcr_iPrAuCl_TAKinetics_890nm = BiFitIRF_new(MesAcr_iPrAuCl_TAKinetics,2,-10.0, 200000, [-0.024, 0.6, 136.2565, 53423.7825, 129.0576, 10226.2776, 0.0005], 10,10,false,398)


    function MesAcr_iPrAuCl_Spectra(f = Figure())

        Data = MesAcr_iPrAuCl_TASpectra

        ### First Panel ##############################################################################

            ax =    CairoMakie.Axis(
                    f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                    xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.4, 0.05),
                    title = "MesAcr-iPrAu(I)Cl TA Spectra",
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

                axislegend("Delay Time", position = :rb,nbanks = 4,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            

        # ### Second Panel ##############################################################################

            ax2 =    CairoMakie.Axis(
                    f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "MesAcr-iPrAu(I)Cl TA Kinetics",
                    xlabel = "Time (µs)", 
                    ylabel = "ΔA (OD)", 
                )

            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

            ### Plot Data ##############################################################################

                lines!(MesAcr_iPrAuCl_TAKinetics_420nm[:,1] ./ 1000,MesAcr_iPrAuCl_TAKinetics_420nm[:,2],linewidth = width,label = "420 nm")
                lines!(MesAcr_iPrAuCl_TAKinetics_420nm[:,1]./ 1000,MesAcr_iPrAuCl_TAKinetics_420nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(MesAcr_iPrAuCl_TAKinetics_490nm[:,1]./ 1000,MesAcr_iPrAuCl_TAKinetics_490nm[:,2],linewidth = width,label = "490 nm")
                lines!(MesAcr_iPrAuCl_TAKinetics_490nm[:,1]./ 1000,MesAcr_iPrAuCl_TAKinetics_490nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(MesAcr_iPrAuCl_TAKinetics_890nm[:,1]./ 1000,MesAcr_iPrAuCl_TAKinetics_890nm[:,2],linewidth = width,label = "890 nm")
                lines!(MesAcr_iPrAuCl_TAKinetics_890nm[:,1]./ 1000,MesAcr_iPrAuCl_TAKinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
                Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
                colgap!(f.layout,1,20)
                colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
MesAcr_iPrAuCl_Spectra()
save("JULIA/Outputs/MesAcr_iPrAuCl_Spectra.pdf",MesAcr_iPrAuCl_Spectra())

##  MesAcr PhosphineOxide
####################################################################################################

    MesAcr_PhosphineOxide_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/Mesacr_PhosphineOxide_2mm_MeCN_20us_TA_Spectra.csv", header = true, DataFrame)
    MesAcr_PhosphineOxide_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/Mesacr_PhosphineOxide_2mm_MeCN_20us_TA_Kinetics.csv", header = true, DataFrame)

    # MesAcr_PhosphineOxide_TAKinetics_523nm = MonoFitIRF_New(MesAcr_PhosphineOxide_TAKinetics,2,-10.0, 2000, [3.286, 0.248, 0.839, 312.471, 0.0], 100,Inf,false,1000)
    MesAcr_PhosphineOxide_TAKinetics_523nm = BiFitIRF_new(MesAcr_PhosphineOxide_TAKinetics,2,-10.0, 2000, [4.5742, 0.3891, 0.887, 378.0788, 0.0305, 1.8661, 0.0003], 10,10,false,1000)

    MesAcr_PhosphineOxide_TAKinetics_389nm = MonoFitIRF_New(MesAcr_PhosphineOxide_TAKinetics,3,-10.0, 2000, [3.637, 1.114, 0.124, 1.997, 0.0], 10,10,false,1000)
    # MesAcr_PhosphineOxide_TAKinetics_490nm = BiFitIRF_new(MesAcr_PhosphineOxide_TAKinetics,3,-10.0, 200000, [-0.024, 0.6, 407.538, 29353.8024, 142.6518, 7492.1869, 0.0009], 10,25,false,390)


    # MesAcr_PhosphineOxide_TAKinetics_701nm = MonoFitIRF_New(MesAcr_PhosphineOxide_TAKinetics,4,-10.0,  2000,[3.637, 1.114, 0.124, 1.997, 0.0], 100,100,false,200)


    function MesAcr_PhosphineOxide_Spectra(f = Figure())

        Data = MesAcr_PhosphineOxide_TASpectra

        ### First Panel ##############################################################################

            ax =    CairoMakie.Axis(
                    f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                    xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.01, 0.01),
                    title = "MesAcr Phosphine Oxide TA Spectra",
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
                CairoMakie.vspan!(431, 441; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            

        # # ### Second Panel ##############################################################################

            ax2 =    CairoMakie.Axis(
                    f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                    title = "MesAcr Phosphine Oxide TA Kinetics",
                    xlabel = "Time (ns)", 
                    ylabel = "ΔA (OD)", 
                )

            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

            ### Plot Data ##############################################################################

                lines!(MesAcr_PhosphineOxide_TAKinetics_389nm[:,1] ,MesAcr_PhosphineOxide_TAKinetics_389nm[:,2],linewidth = width,label = "389 nm")
                lines!(MesAcr_PhosphineOxide_TAKinetics_389nm[:,1],MesAcr_PhosphineOxide_TAKinetics_389nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(MesAcr_PhosphineOxide_TAKinetics_523nm[:,1],MesAcr_PhosphineOxide_TAKinetics_523nm[:,2],linewidth = width,label = "523 nm")
                lines!(MesAcr_PhosphineOxide_TAKinetics_523nm[:,1],MesAcr_PhosphineOxide_TAKinetics_523nm[:,3],linewidth = 2,linestyle = :dash,color = :black)



                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
                Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
                colgap!(f.layout,1,20)
                colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
MesAcr_PhosphineOxide_Spectra()
save("JULIA/Outputs/MesAcr_PhosphineOxide_Spectra.pdf",MesAcr_PhosphineOxide_Spectra())

##  MesAcr PhosphineOxide + Diphenyl Anthracene (DPA)
####################################################################################################

    MesAcr_PhosphineOxide_DPA_TASpectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/Mesacr_PhosphineOxide_DPA_2mm_MeCN_1us_TA_Spectra.csv", header = true, DataFrame)
    MesAcr_PhosphineOxide_DPA_TAKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_08_21/NSTA/Analysis/Mesacr_PhosphineOxide_DPA_2mm_MeCN_1us_TA_Kinetics.csv", header = true, DataFrame)

    MesAcr_PhosphineOxide_DPA_TAKinetics_386nm = TriFitIRF_new(MesAcr_PhosphineOxide_DPA_TAKinetics,2,-10.0, 1000,[3.919, 0.62, -0.208, 39.144, 0.281, 2.384, 1.444, 375.357, -0.001], 1,10,false,200)

    MesAcr_PhosphineOxide_DPA_TAKinetics_442nm = BiFitIRF_new(MesAcr_PhosphineOxide_DPA_TAKinetics,3,-20.0, 1000,[2.5825, 0.3319, -18.9062, 156.7413, 23.0934, 198.193, 0.0], 10,1,false,200)



    MesAcr_PhosphineOxide_DPA_TAKinetics_501nm = BiFitIRF_new(MesAcr_PhosphineOxide_DPA_TAKinetics,4,-20.0, 1000,[2.7392, 0.3179, 2.1021, 198.3503, 0.2503, 1.9692, 0.0009], 10,10,false,200)
    MesAcr_PhosphineOxide_DPA_TAKinetics_501nm = TriFitIRF_new(MesAcr_PhosphineOxide_DPA_TAKinetics,4,-20.0, 1000,[3.563, 0.656, -1.377, 141.067, 0.159, 1.914, 1.749, 218.012, 0.001], 10,100,false,200)


    MesAcr_PhosphineOxide_DPA_TAKinetics_700nm = BiFitIRF_new(MesAcr_PhosphineOxide_DPA_TAKinetics,5,-20.0, 1000,[3.6887, 0.6112, 0.1165, 192.2196, 0.1128, 2.5425, 0.0012], 10,10,false,200)
    MesAcr_PhosphineOxide_DPA_TAKinetics_891nm = BiFitIRF_new(MesAcr_PhosphineOxide_DPA_TAKinetics,6,-20.0, 1000,[3.3198, 0.6723, 1.1048, 172.9976, -0.0374, 2.2883, 0.0003], 100,100,false,200)



    function MesAcr_PhosphineOxide_DPA_Spectra(f = Figure())

        Data = MesAcr_PhosphineOxide_DPA_TASpectra

        ### First Panel ##############################################################################

            ax =    CairoMakie.Axis(
                    f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                    xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.01, 0.16),
                    title = "MesAcr Phosphine Oxide + DPA TA Spectra",
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

                axislegend("Delay Time", position = :rt,nbanks = 4,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

                hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
                # text!(x,y,text = "text", color = :red, textsize = 20)
                CairoMakie.vspan!(431, 441; ymin = 0.0, ymax = 1.0, color = :lightblue)

                colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
                rowsize!(f.layout,1, 400)
            

        # ### Second Panel ##############################################################################

            ax2 =    CairoMakie.Axis(
                    f[1,2], palette = (color = palette(ColorSchemes.hsv, 5),),
                    title = "MesAcr Phosphine Oxide + DPA TA Kinetics",
                    xlabel = "Time (ns)", 
                    ylabel = "ΔA (OD)", 
                )

            width = 3

            ### Custom Axis limits ##############################################################################
                # ax.xticks= 1:2:17
                # CairoMakie.xlims!(375,910)
                # CairoMakie.ylims!(nothing,0.04)

            ### Plot Data ##############################################################################

                lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_386nm[:,1] ,MesAcr_PhosphineOxide_DPA_TAKinetics_386nm[:,2],linewidth = width,label = "386 nm")
                lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_386nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_386nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

                lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_442nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_442nm[:,2],linewidth = width,label = "442 nm")
                lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_442nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_442nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

                lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_501nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_501nm[:,2],linewidth = width,label = "501 nm")
                lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_501nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_501nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

                # lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_700nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_700nm[:,2],linewidth = width,label = "700 nm")
                # lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_700nm[:,1],MesAcr_PhosphineOxide_DPA_TAKinetics_700nm[:,3],linewidth = 2, linestyle = :dash,color = :black)


                axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
                hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
                Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
                colgap!(f.layout,1,20)
                colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

        resize_to_layout!(f)
        f
    end
MesAcr_PhosphineOxide_DPA_Spectra()
save("JULIA/Outputs/MesAcr_PhosphineOxide_DPA_Spectra.pdf",MesAcr_PhosphineOxide_DPA_Spectra())













function Comparisons(f = Figure())


    ### First Panel ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.hsv, 5),), ## or replace ncol with Integer 
                title = "Title",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
            )
        
        width = 3

        # lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_386nm[:,1] ,MesAcr_PhosphineOxide_DPA_TAKinetics_386nm[:,2],linewidth = width,label = "336 nm")
        # lines!(MesAcr_PhosphineOxide_TAKinetics_389nm[:,1] ,MesAcr_PhosphineOxide_TAKinetics_389nm[:,2],linewidth = width,label = "389 nm")

        lines!(MesAcr_PhosphineOxide_DPA_TAKinetics_501nm[:,1] ,MaxNorm(MesAcr_PhosphineOxide_DPA_TAKinetics_501nm[:,2]),linewidth = width,label = "501 DPA nm")
        lines!(MesAcr_PhosphineOxide_TAKinetics_523nm[:,1] ,MaxNorm(MesAcr_PhosphineOxide_TAKinetics_523nm[:,2]),linewidth = width,label = "523 PO nm")

        # lines!(MesAcr_PhosphineOxide_DPA_TASpectra[:,1] ,MaxNorm(MesAcr_PhosphineOxide_DPA_TASpectra[:,15]),linewidth = width,label = "+ DPA")
        # lines!(MesAcr_PhosphineOxide_TASpectra[:,1] ,MaxNorm(MesAcr_PhosphineOxide_TASpectra[:,10]),linewidth = width,label = "Phosphine Oxide")

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,500)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################


            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        

    # ### Second Panel ##############################################################################

    #     ax2 =    CairoMakie.Axis(
    #             f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
    #             title = "Title",
    #             xlabel = "Time (ns)", 
    #             ylabel = "ΔA (OD)", 
    #         )

    #     width = 3

    #     ### Custom Axis limits ##############################################################################
    #         # ax.xticks= 1:2:17
    #         # CairoMakie.xlims!(375,910)
    #         # CairoMakie.ylims!(nothing,0.04)

    #     ### Plot Data ##############################################################################

    #         lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,2],linewidth = width,label = "nm")
    #         lines!(PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,1],PhAuCl3_Photolyzed_TA_Kinetics_435nm[:,3],linewidth = 2, linestyle = :dash, color = :black)

    #         lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,2],linewidth = width,label = "nm")
    #         lines!(PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_490nm[:,3],linewidth = 2,linestyle = :dash,color = :black)

    #         lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,2],linewidth = width,label = "nm")
    #         lines!(PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,1],PhAuCl3_Photlyzed_TA_Kinetics_890nm[:,3],linewidth = 2, linestyle = :dash,color = :black)

    #         axislegend("Wavelength", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    #         hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

    #         Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
    #         Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
    #         colgap!(f.layout,1,20)
    #         colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
Comparisons()
save("JULIA/Outputs/Comparisons.pdf",Comparisons())


