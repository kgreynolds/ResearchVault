
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

###########################################################################

#### ReBr2 in DMF
    ReBr2_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/RNR/2024_08_02/Analysis/ReBr_inDMF_needleptamount_run2_spectra.csv", header = true, DataFrame)
    ReBr2_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/RNR/2024_08_02/Analysis/ReBr_inDMF_needleptamount_run2_Kinetics.csv", header = true, DataFrame)

    ReBr2_Kinetics_430nm = MonoFitIRF_New(ReBr2_Kinetics,3,0.0,1500, [0.008, 1.5, 4.544, 157.644, 0.001], 5,Inf,false, 199)
    #ReBr2_Kinetics_430nm = BiFitIRF_new(ReBr2_Kinetics,3,0.0,1500, [0.008, 1.4998, 4.5565, 161.929, 0.0374, 8.2608, 0.0005], 10,Inf,false, 200)

    ReBr2_Kinetics_470nm = MonoFitIRF_New(ReBr2_Kinetics,4,0.0,1500, [0.008, 1.5, 4.655, 162.207, 0.0], 10,Inf,false, 199)
    ReBr2_Kinetics_700nm = MonoFitIRF_New(ReBr2_Kinetics,5,0.0,1500, [0.008, 1.5, 2.267, 167.807, 0.0], 100,Inf,false, 199)
    ReBr2_Kinetics_700nm = BiFitIRF_new(ReBr2_Kinetics,5,0.0,1500, [0.008, 1.499, 2.2256, 177.678, 0.0791, 54.5673, 0.0004], 10,Inf,false, 199)

    mean([162.2,167.8,177.7])
    std([162.2,167.8,177.7])

    function ReBr2_TA_Spectra(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of ReBr₂ in DMF",
                palette = (color = palette(ColorSchemes.inferno, 12),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,2],linewidth = width, label = "5 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,3],linewidth = width, label = "10 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,4],linewidth = width, label = "25 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,5],linewidth = width, label = "50 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,6],linewidth = width, label = "75 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,7],linewidth = width, label = "100 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,8],linewidth = width, label = "150 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,9],linewidth = width, label = "200 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,10],linewidth = width, label = "300 ns")
            # lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,11],linewidth = width, label = "400 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,12],linewidth = width, label = "500 ns")
            lines!(ReBr2_Spectra[:,1],ReBr2_Spectra[:,13],linewidth = width, label = "1000 ns")


    
    
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
    
    
            # vlines!(ax, 436, linewidth = 2,color = :blue,linestyle = :dash)
            hlines!(ax,0, linewidth = 2,color = :black,linestyle = :dash)
    
            axislegend("Delay Time", position = :rt,nbanks = 3, framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    
            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    ReBr2_TA_Spectra()
    save("JULIA/Outputs/ReBr2_TA_Spectra.pdf",ReBr2_TA_Spectra())

    function ReBr2_TA_Kinetics(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Kinetics of ReBr₂ in DMF",
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            lines!(ReBr2_Kinetics_430nm[:,1],MaxNorm(ReBr2_Kinetics_430nm[:,2]),linewidth = width, label = "430 nm")
            # lines!(ReBr2_Kinetics_430nm[:,1],MaxNorm(ReBr2_Kinetics_430nm[:,3]),linewidth = 2, color = :lightblue)

            lines!(ReBr2_Kinetics_470nm[:,1],MaxNorm(ReBr2_Kinetics_470nm[:,2]),linewidth = width, label = "470 nm")
            # lines!(ReBr2_Kinetics_470nm[:,1],MaxNorm(ReBr2_Kinetics_470nm[:,3]),linewidth = 2, color = :lightblue)

            lines!(ReBr2_Kinetics_700nm[:,1],MaxNorm(ReBr2_Kinetics_700nm[:,2]),linewidth = width, label = "70 nm")


    
    
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
    
    
            # vlines!(ax, 436, linewidth = 2,color = :blue,linestyle = :dash)
            hlines!(ax,0, linewidth = 2,color = :black,linestyle = :dash)
    
            axislegend("Wavelength", position = :rt,nbanks = 1, framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    
            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    ReBr2_TA_Kinetics()
    save("JULIA/Outputs/ReBr2_TA_Kinetics.pdf",ReBr2_TA_Kinetics())

#### 50 µM photo beta 2 in pH 8.3 buffer
    Photobeta_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/RNR/2024_08_02/Analysis/Photob2_50uM_8p3_Spectra.csv", header = true, DataFrame)
    Photobeta_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/RNR/2024_08_02/Analysis/Photob2_50uM_8p3_Kinetics.csv", header = true, DataFrame)

    Photobeta_Kinetics_423nm = BiFitIRF_new(Photobeta_Kinetics,3,0.0,1500, [0.008, 1.4998, 5.6708, 326.1197, 0.6039, 37.1014, 0.0014], 10,Inf,false, 200)
    Photobeta_Kinetics_500nm = BiFitIRF_new(Photobeta_Kinetics,4,0.0,1500,[0.008, 1.4998, 4.8537, 348.6673, 0.5491, 43.4274, 0.0009], 10,Inf,false, 200)
    Photobeta_Kinetics_550nm = BiFitIRF_new(Photobeta_Kinetics,5,0.0,1500,[0.008, 1.4998, 4.0704, 359.3464, 0.722, 41.9728, 0.0011], 10,Inf,false, 200)

    mean([326.1, 348.7, 359.3])
    std([326.1, 348.7, 359.3])

    mean([37.1, 43.4, 41.9])
    std([37.1, 43.4, 41.9])


    function Photobeta_TA_Spectra(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Photo-β (50 µM)",
                palette = (color = palette(ColorSchemes.inferno, 10),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            CairoMakie.ylims!(-0.001,0.037)
            
            
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,3],linewidth = width, label = "5 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,4],linewidth = width, label = "10 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,5],linewidth = width, label = "20 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,6],linewidth = width, label = "40 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,7],linewidth = width, label = "80 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,8],linewidth = width, label = "150 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,9],linewidth = width, label = "400 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,10],linewidth = width, label = "1000 ns")
            lines!(Photobeta_Spectra[:,1],Photobeta_Spectra[:,11],linewidth = width, label = "4000 ns")



    
    
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
    
    
            # vlines!(ax, 423, linewidth = 2,color = :blue,linestyle = :dash)
            # vlines!(ax, 500, linewidth = 2,color = :blue,linestyle = :dash)
            # vlines!(ax, 550, linewidth = 2,color = :blue,linestyle = :dash)

            hlines!(ax,0, linewidth = 2,color = :black,linestyle = :dash)
    
            axislegend("Delay Time", position = :rt,nbanks = 2, framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    
            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    Photobeta_TA_Spectra()
    save("JULIA/Outputs/Photobeta_TA_Spectra.pdf",Photobeta_TA_Spectra())

    function Photobeta_TA_Kinetics(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Kinetics of Photo-β (50 µM)",
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            lines!(Photobeta_Kinetics_423nm[:,1],MaxNorm(Photobeta_Kinetics_423nm[:,2]),linewidth = width, label = "423 nm")
            # lines!(ReBr2_Kinetics_430nm[:,1],MaxNorm(ReBr2_Kinetics_430nm[:,3]),linewidth = 2, color = :lightblue)

            lines!(Photobeta_Kinetics_500nm[:,1],MaxNorm(Photobeta_Kinetics_500nm[:,2]),linewidth = width, label = "500 nm")
            # lines!(ReBr2_Kinetics_470nm[:,1],MaxNorm(ReBr2_Kinetics_470nm[:,3]),linewidth = 2, color = :lightblue)

            lines!(Photobeta_Kinetics_550nm[:,1],MaxNorm(Photobeta_Kinetics_550nm[:,2]),linewidth = width, label = "550 nm")


    
    
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
    
    
            # vlines!(ax, 436, linewidth = 2,color = :blue,linestyle = :dash)
            hlines!(ax,0, linewidth = 2,color = :black,linestyle = :dash)
    
            axislegend("Wavelength", position = :rt,nbanks = 1, framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    
            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    Photobeta_TA_Kinetics()
    save("JULIA/Outputs/Photobeta_TA_Kinetics.pdf",Photobeta_TA_Kinetics())



#### 50 µM photo beta 2 + 10 mM RuHex in pH 8.3 buffer
    Photobeta_RuHex_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/RNR/2024_08_02/Analysis/Photob2_50uM_10mM_RuHex_8p3_3sInt_Spectra.csv", header = true, DataFrame)
    Photobeta_RuHex_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/RNR/2024_08_02/Analysis/Photob2_50uM_10mM_RuHex_8p3_3sInt_Kinetics.csv", header = true, DataFrame)

    Photobeta_RuHex_Kinetics_423nm = MonoFitIRF_New(Photobeta_RuHex_Kinetics,3,0.0,500, [0.008, 1.5, 1.164, 55.092, 0.002], 5,Inf,false, 199)
    Photobeta_RuHex_Kinetics_411nm = MonoFitIRF_New(Photobeta_RuHex_Kinetics,4,0.0,500,[0.008, 1.5, 0.998, 53.859, 0.003], 5,Inf,false, 199)
    Photobeta_RuHex_Kinetics_530nm = MonoFitIRF_New(Photobeta_RuHex_Kinetics,5,0.0,500, [0.008, 1.5, 0.946, 71.357, 0.001], 5,Inf,false, 199)

    Photobeta_RuHex_Kinetics_423nm = BiFitIRF_new(Photobeta_RuHex_Kinetics,3,0.0,500, [0.008, 1.4998, 5.6708, 326.1197, 0.6039, 37.1014, 0.0014], 10,Inf,false, 200)
    Photobeta_RuHex_Kinetics_411nm = BiFitIRF_new(Photobeta_RuHex_Kinetics,4,0.0,500,[0.008, 1.4998, 4.8537, 348.6673, 0.5491, 43.4274, 0.0009], 10,Inf,false, 200)
    Photobeta_RuHex_Kinetics_530nm = BiFitIRF_new(Photobeta_RuHex_Kinetics,5,0.0,500,[0.008, 1.4998, 4.0704, 359.3464, 0.722, 41.9728, 0.0011], 10,Inf,false, 200)


    function Photobeta_RuHex_TA_Spectra(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Spectra of Photo-β (50 µM) + RuHex (10 mM)",
                palette = (color = palette(ColorSchemes.inferno, 10),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (mOD)"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            CairoMakie.ylims!(-0.001,0.037)
            
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,3],linewidth = width, label = "5 ns")
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,4],linewidth = width, label = "15 ns")
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,5],linewidth = width, label = "35 ns")
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,6],linewidth = width, label = "75 ns")
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,7],linewidth = width, label = "150 ns")
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,8],linewidth = width, label = "12.5 µs")
            lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,9],linewidth = width, label = "30 µs")
            # lines!(Photobeta_RuHex_Spectra[:,1],Photobeta_RuHex_Spectra[:,10],linewidth = width, label = "50000 ns")

            CairoMakie.vspan!(705, 720; ymin = 0.0, ymax = 1.0, color = :lightblue)
    
            vlines!(ax, 423, linewidth = 2,color = :blue,linestyle = :dash)
            # vlines!(ax, 500, linewidth = 2,color = :blue,linestyle = :dash)
            # vlines!(ax, 550, linewidth = 2,color = :blue,linestyle = :dash)

            hlines!(ax,0, linewidth = 2,color = :black,linestyle = :dash)
    
            axislegend("Delay Time", position = :rt,nbanks = 1, framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    
            text!(450,0.03,text = "τ₁ ≈ 60 ns", color = :blue, fontsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    Photobeta_RuHex_TA_Spectra()
    save("JULIA/Outputs/Photobeta_RuHex_TA_Spectra.pdf",Photobeta_RuHex_TA_Spectra())

    function Photobeta_RuHex_TA_Kinetics(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "TA Kinetics of Photo-β",
                palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Time (ns)", 
                ylabel = "ΔA"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            lines!(Photobeta_Kinetics_423nm[:,1],MaxNorm(Photobeta_Kinetics_423nm[:,2]),linewidth = width, label = "423 nm")
            # lines!(ReBr2_Kinetics_430nm[:,1],MaxNorm(ReBr2_Kinetics_430nm[:,3]),linewidth = 2, color = :lightblue)

            lines!(Photobeta_Kinetics_500nm[:,1],MaxNorm(Photobeta_Kinetics_500nm[:,2]),linewidth = width, label = "500 nm")
            # lines!(ReBr2_Kinetics_470nm[:,1],MaxNorm(ReBr2_Kinetics_470nm[:,3]),linewidth = 2, color = :lightblue)

            lines!(Photobeta_Kinetics_550nm[:,1],MaxNorm(Photobeta_Kinetics_550nm[:,2]),linewidth = width, label = "550 nm")


    
    
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)
    
    
            # vlines!(ax, 436, linewidth = 2,color = :blue,linestyle = :dash)
            hlines!(ax,0, linewidth = 2,color = :black,linestyle = :dash)
    
            axislegend("Delay Time", position = :rt,nbanks = 3, framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
    
            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    Photobeta_RuHex_TA_Kinetics()
    save("JULIA/Outputs/Photobeta_RuHex_TA_Kinetics.pdf",PhotobetaPhotobeta_RuHex_TA_Kinetics_TA_Kinetics())