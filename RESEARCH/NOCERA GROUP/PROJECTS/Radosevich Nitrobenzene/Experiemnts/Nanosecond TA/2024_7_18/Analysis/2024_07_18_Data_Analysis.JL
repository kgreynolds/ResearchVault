## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

# Open and read MatLab Files (For Magnitude)
#File = EnVisionKinetic(2,"path")
function smoothMat(Mat,SmoothPts)
    """This function takes and input matrix Mar and generates a new matrix sMat with the moving average over
        some number of points defined by the SmoothPts argument"""
        sMat = zeros(nrow(Mat)-SmoothPts,ncol(Mat))

    for i in 1:nrow(Mat)-SmoothPts
        for j in 1:ncol(Mat)
            sMat[i,j] = mean(skipmissing(Mat[i:i+SmoothPts,j])) # This works! calculates the mean of the i to i+2th row in the jth columns and replaces the E[i,j]th value with that mean
        end
    end
    return sMat #, display(Plots.plot(sMat[:,1], sMat[:,2],xlabel = "Time", ylabel = " ΔA"))
end

## IMPORT DATA

Nb_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-3_nitrobenzene 300pt -bg -chirp_TA_Spectra.csv",header = true, DataFrame)
Nb_Spectra_Smoothed = smoothMat(Nb_Spectra,3)

Nb_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-3_nitrobenzene 300pt -bg -chirp_TA_Kinetics.csv",header = true, DataFrame)
Nb_Kinetics = hcat(Nb_Kinetics[:,1]*1000,Nb_Kinetics[:,2:end])


Nb_Si_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG-45-1_Nitrobenzene+Silane inCPME_(2mm)_390_TA_Spectra.csv", header = true, DataFrame)
Nb_Si_Spectra = smoothMat(Nb_Si_Spectra,3)


Nb_Si_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG-45-1_Nitrobenzene+Silane inCPME_(2mm)_390_TA_Kinetics.csv", header = true, DataFrame)
Nb_Si_Kinetics = hcat(Nb_Si_Kinetics[:,1]*1000,Nb_Si_Kinetics[:,3:end])



Nb_Si_Spectra_100ns = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-1_nitrobenzene +silane in CPME (2mm)_390nm 100ns -bg -bg -chirp_TA_Spectra.csv", header = true, DataFrame)
Nb_Si_Kinetics_100ns = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-1_nitrobenzene +silane in CPME (2mm)_390nm 100ns -bg -bg -chirp_TA_Kinetics.csv", header = true, DataFrame)
Nb_Si_Kinetics_100ns = hcat(Nb_Si_Kinetics_100ns[:,1]*1000,Nb_Si_Kinetics_100ns[:,2:end])


Nb_Phos_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-4_nitrobenzene_phosphine 200pt_100ns -bg -chirp_TA_Spectra.csv", header = true, DataFrame)
Nb_Phos_Spectra_Smoothed = smoothMat(Nb_Phos_Spectra,3)

Nb_Phos_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-4_nitrobenzene_phosphine 200pt_100ns -bg -chirp_TA_Kinetics.csv", header = true, DataFrame)
Nb_Phos_Kinetics = hcat(Nb_Phos_Kinetics[:,1]*1000,Nb_Phos_Kinetics[:,2:end])

Nb_DIPEA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-07-17/Sample4_Nitrobenzene_DIPEA inCPME (1cm)_390nm_1ms -bg -t0 TA Kinetics.csv",header = true, DataFrame)
Nb_DIPEA_Kinetics_Smoothed = smoothMat(Nb_DIPEA_Kinetics,3)


Nb_DIPEA_Spectra = CSV.read("RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-07-17/Sample4_Nitrobenzene_DIPEA inCPME (1cm)_390nm_1ms -bg -t0 TA Spectra.csv",header = true, DataFrame)
Nb_DIPEA_Spectra_smooth = smoothMat(Nb_DIPEA_Spectra,3)


NNaph_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-5_nitro-naphthalene_1cm_5us -bg -chirp_TA_Spectra.csv", header = true, DataFrame)
NNaph_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Radosevich Nitrobenzene/Experiemnts/Nanosecond TA/2024-7-18/Analysis/Data to analyze/JG6-45-5_nitro-naphthalene_1cm_5us -bg -chirp_TA_Kinetics.csv", header = true, DataFrame)
NNaph_Kinetics = hcat(NNaph_Kinetics[:,1]*1000,NNaph_Kinetics[:,2:end])


# Data Workup
####################################################################################################

Nb_Kinetics_433nm = MonoFitIRF_New(Nb_Kinetics,3,-2,18.0,[0.499, 0.636, 0.043, 0.895, 0.001], 0,0,false,120)
Nb_Kinetics_700nm = MonoFitIRF_New(Nb_Kinetics,5,-2,18.0,[0.499, 0.636, 0.033, 0.911, -0.0], 0,0,false,120)

Nb_Si_Kinetics_413nm = MonoFitIRF_New(Nb_Si_Kinetics,2,-2,10.0,[0.499, 0.635, 0.041, 0.898, 0.002], 0,0,false,120)
Nb_Si_Kinetics_432nm = MonoFitIRF_New(Nb_Si_Kinetics,3,-2,10.0,[0.499, 0.635, 0.044, 0.893, 0.002], 0,0,false,120)
Nb_Si_Kinetics_449nm = MonoFitIRF_New(Nb_Si_Kinetics,4,-2,10.0,[0.499, 0.635, 0.045, 0.883, 0.001], 0,0,false,120)
Nb_Si_Kinetics_700nm = MonoFitIRF_New(Nb_Si_Kinetics,5,-2,10.0,[0.499, 0.635, 0.033, 0.875, 0.0], 0,0,false,120)

Nb_Si_Kinetics_432nm = BiFitIRF_new(Nb_Si_Kinetics,3,-5,18.0,[0.2894, 0.5164, -0.4223, 0.617, 0.4584, 0.6737, 0.0019], 0,0,false,120)

Nb_DIPEA_Kinetics_427nm = MonoFitIRF_New(Nb_DIPEA_Kinetics,3,-2,1000000.0,[0.499, 0.636, 3134.024, 304575.279, -0.0], 0,0,false,178)
Nb_DIPEA_Kinetics_468nm = MonoFitIRF_New(Nb_DIPEA_Kinetics,4,-2,1000000.0,[0.499, 0.636, 789.339, 211647.609, -0.0], 0,0,false,178)


NNaph_Kinetics_400nm = MonoFitIRF_New(NNaph_Kinetics,2,-2,4800.0,[0.486, 0.636, 41.227, 1917.454, 0.001], 0,0,false,300)
NNaph_Kinetics_570nm = MonoFitIRF_New(NNaph_Kinetics,3,-2,4800.0,[0.486, 0.636, 69.335, 1944.009, 0.002], 0,0,false,300)
NNaph_Kinetics_815nm = MonoFitIRF_New(NNaph_Kinetics,4,-2,4800.0,[0.486, 0.636, 12.89, 1896.758, -0.0], 0,0,false,300)


##Nitrobenzene + Silane
####################################################################################################

function Nb_Si_TA_Spectra(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Silane TA Spectra in CPME",
            palette = (color = palette(ColorSchemes.inferno, 8),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(395,850)
        CairoMakie.ylims!(-0.001,nothing)
        
        # lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,3],linewidth = width,label = "-1 ns")
        # lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,4],linewidth = width,label = "0 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,5],linewidth = width,label = "1 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,6],linewidth = width,label = "2 ns")
        # lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,7],linewidth = width,label = "2.1 ns")
        # lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,8],linewidth = width,label = "2.2 ns")
        # lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,9],linewidth = width,label = "2.3 ns")
        # lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,10],linewidth = width,label = "2.4 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,11],linewidth = width,label = "2.5 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,16],linewidth = width,label = "3 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,17],linewidth = width,label = "3.5 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,18],linewidth = width,label = "4 ns")
        lines!(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,19],linewidth = width,label = "5 ns")

        ## SMOOTHED SPECTRA

        # # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,3],linewidth = width,label = "-1 ns")
        # # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,4],linewidth = width,label = "0 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,5],linewidth = width,label = "1 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,6],linewidth = width,label = "2 ns")
        # # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,7],linewidth = width,label = "2.1 ns")
        # # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,8],linewidth = width,label = "2.2 ns")
        # # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,9],linewidth = width,label = "2.3 ns")
        # # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,10],linewidth = width,label = "2.4 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,11],linewidth = width,label = "2.5 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,16],linewidth = width,label = "3 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,17],linewidth = width,label = "3.5 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,18],linewidth = width,label = "4 ns")
        # lines!(Nb_Si_Spectra_smooth[:,1],Nb_Si_Spectra_smooth[:,19],linewidth = width,label = "5 ns")


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Si_TA_Spectra()
save("JULIA/Outputs/2024_07_18_Nb_Si_TA_Spectra.png",Nb_Si_TA_Spectra())

function Nb_Si_TA_Kinetics(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Silane TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        # lines!(Nb_Si_Kinetics_413nm[:,1],Nb_Si_Kinetics_413nm[:,2],linewidth = width,label = "413 nm , τ = 0.898 ns")
        # lines!(Nb_Si_Kinetics_413nm[:,1],Nb_Si_Kinetics_413nm[:,3],linewidth = width)

        # lines!(Nb_Si_Kinetics_432nm[:,1],Nb_Si_Kinetics_432nm[:,2],linewidth = width,label = "432 nm, τ = 0.893 ns")
        # lines!(Nb_Si_Kinetics_432nm[:,1],Nb_Si_Kinetics_432nm[:,3],linewidth = width)

        lines!(Nb_Si_Kinetics_449nm[:,1],Nb_Si_Kinetics_449nm[:,2],linewidth = width,label = "449 nm, τ = 0.883 ns")
        lines!(Nb_Si_Kinetics_449nm[:,1],Nb_Si_Kinetics_449nm[:,3],linewidth = 2)

        lines!(Nb_Si_Kinetics_700nm[:,1],Nb_Si_Kinetics_700nm[:,2],linewidth = width,label = "700 nm, τ = 0.875ns")
        lines!(Nb_Si_Kinetics_700nm[:,1],Nb_Si_Kinetics_700nm[:,3],linewidth = 2)


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Si_TA_Kinetics()
save("JULIA/Outputs/2024_07_18_Nb_Si_TA_Kinetics.png",Nb_Si_TA_Kinetics())

# Nitrobenzene + Silane 100ns delay
####################################################################################################

function Nb_Si_TA_Spectra_100ns(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Silane TA Spectra in CPME",
            palette = (color = palette(ColorSchemes.inferno, 8),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(395,850)
        CairoMakie.ylims!(-0.005,nothing)
        
        # lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,3],linewidth = width,label = "-1 ns")
        # lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,4],linewidth = width,label = "0 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,5],linewidth = width,label = "1 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,6],linewidth = width,label = "2 ns")
        # lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,7],linewidth = width,label = "2.1 ns")
        # lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,8],linewidth = width,label = "2.2 ns")
        # lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,9],linewidth = width,label = "2.3 ns")
        # lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,10],linewidth = width,label = "2.4 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,11],linewidth = width,label = "2.5 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,16],linewidth = width,label = "3 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,17],linewidth = width,label = "3.5 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,18],linewidth = width,label = "4 ns")
        lines!(Nb_Si_Spectra_100ns[:,1],Nb_Si_Spectra_100ns[:,19],linewidth = width,label = "5 ns")


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Si_TA_Spectra_100ns()
save("JULIA/Outputs/2024_07_18_Nb_Si_TA_Spectra_100ns.png",Nb_Si_TA_Spectra_100ns())

function Nb_Si_TA_Kinetics_100ns(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Silane TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 4),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        lines!(Nb_Si_Kinetics_100ns[:,1],Nb_Si_Kinetics_100ns[:,2],linewidth = width,label = "417 nm")
        lines!(Nb_Si_Kinetics_100ns[:,1],Nb_Si_Kinetics_100ns[:,3],linewidth = width,label = "432 nm")
        lines!(Nb_Si_Kinetics_100ns[:,1],Nb_Si_Kinetics_100ns[:,4],linewidth = width,label = "450 nm")
        # lines!(Nb_Si_Kinetics_100ns[:,1],Nb_Si_Kinetics_100ns[:,5],linewidth = width,label = "700 nm")




        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Si_TA_Kinetics_100ns()
save("JULIA/Outputs/2024_07_18_Nb_Si_TA_Kinetics_100ns.png",Nb_Si_TA_Kinetics_100ns())

# Nitrobenzene alone
####################################################################################################

function Nb_TA_Spectra(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene TA Spectra in CPME",
            palette = (color = palette(ColorSchemes.inferno, 6),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(395,850)
        CairoMakie.ylims!(-0.001,nothing)
        
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,2],linewidth = width,label = "-1 ns")
        lines!(Nb_Spectra[:,1],Nb_Spectra[:,3],linewidth = width,label = "0 ns")
        lines!(Nb_Spectra[:,1],Nb_Spectra[:,4],linewidth = width,label = "1 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,5],linewidth = width,label = "1.5 ns")
        lines!(Nb_Spectra[:,1],Nb_Spectra[:,6],linewidth = width,label = "2 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,7],linewidth = width,label = "2.1 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,8],linewidth = width,label = "2.2 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,9],linewidth = width,label = "2.5 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,10],linewidth = width,label = "2.7 ns")
        lines!(Nb_Spectra[:,1],Nb_Spectra[:,11],linewidth = width,label = "3 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,12],linewidth = width,label = "3.5 ns")
        # lines!(Nb_Spectra[:,1],Nb_Spectra[:,13],linewidth = width,label = "4 ns")
        lines!(Nb_Spectra[:,1],Nb_Spectra[:,14],linewidth = width,label = "5 ns")

        ## SMOOTHED SPECTRA

                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,2],linewidth = width,label = "-1 ns")
                # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,3],linewidth = width,label = "0 ns")
                # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,4],linewidth = width,label = "1 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,5],linewidth = width,label = "1.5 ns")
                # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,6],linewidth = width,label = "2 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,7],linewidth = width,label = "2.1 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,8],linewidth = width,label = "2.2 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,9],linewidth = width,label = "2.5 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,10],linewidth = width,label = "2.7 ns")
                # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,11],linewidth = width,label = "3 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,12],linewidth = width,label = "3.5 ns")
                # # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,13],linewidth = width,label = "4 ns")
                # lines!(Nb_Spectra_Smoothed[:,1],Nb_Spectra_Smoothed[:,14],linewidth = width,label = "5 ns")


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_TA_Spectra()
save("JULIA/Outputs/2024_07_18_Nb_TA_Spectra.png",Nb_TA_Spectra())

function Nb_TA_Kinetics(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        lines!(Nb_Kinetics_433nm[:,1],Nb_Kinetics_433nm[:,2],linewidth = width,label = "417 nm τ = 0.895 ns")
        lines!(Nb_Kinetics_433nm[:,1],Nb_Kinetics_433nm[:,3],linewidth = 3)

        lines!(Nb_Kinetics_700nm[:,1],Nb_Kinetics_700nm[:,2],linewidth = width,label = "450 nm, τ = 0.903 ns")
        lines!(Nb_Kinetics_700nm[:,1],Nb_Kinetics_700nm[:,3],linewidth = 3)




        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_TA_Kinetics()
save("JULIA/Outputs/2024_07_18_Nb_TA_Kinetics.png",Nb_TA_Kinetics())

# Nitrobenzene + Phosphine 
####################################################################################################

function Nb_Phos_TA_Spectra(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Phosphine TA Spectra in CPME",
            palette = (color = palette(ColorSchemes.inferno, 7),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,2],linewidth = width,label = "-1 ns")
        # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,3],linewidth = width,label = "0 ns")
        # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,4],linewidth = width,label = "1 ns")
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,5],linewidth = width,label = "1.5 ns")
        # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,6],linewidth = width,label = "2 ns")
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,7],linewidth = width,label = "2.1 ns")
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,8],linewidth = width,label = "2.2 ns")
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,9],linewidth = width,label = "2.5 ns")
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,10],linewidth = width,label = "2.7 ns")
        # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,11],linewidth = width,label = "3 ns")
        # # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,12],linewidth = width,label = "3.5 ns")
        # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,13],linewidth = width,label = "4 ns")
        # lines!(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,14],linewidth = width,label = "5 ns")

        ## SMOOTHED SPECTRA

                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,2],linewidth = width,label = "-1 ns")
                lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,3],linewidth = width,label = "0 ns")
                lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,4],linewidth = width,label = "1 ns")
                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,5],linewidth = width,label = "1.5 ns")
                lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,6],linewidth = width,label = "2 ns")
                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,7],linewidth = width,label = "2.1 ns")
                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,8],linewidth = width,label = "2.2 ns")
                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,9],linewidth = width,label = "2.5 ns")
                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,10],linewidth = width,label = "2.7 ns")
                lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,11],linewidth = width,label = "3 ns")
                # lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,12],linewidth = width,label = "3.5 ns")
                lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,13],linewidth = width,label = "4 ns")
                lines!(Nb_Phos_Spectra_Smoothed[:,1],Nb_Phos_Spectra_Smoothed[:,14],linewidth = width,label = "5 ns")



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Phos_TA_Spectra()
save("JULIA/Outputs/2024_07_18_Nb_Phos_TA_Spectra.png",Nb_Phos_TA_Spectra())

function Nb_Phos_TA_Kinetics(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Phosphine TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 6),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        lines!(Nb_Phos_Kinetics[:,1],Nb_Phos_Kinetics[:,3],linewidth = width,label = "402 nm")
        lines!(Nb_Phos_Kinetics[:,1],Nb_Phos_Kinetics[:,4],linewidth = width,label = "416 nm")
        lines!(Nb_Phos_Kinetics[:,1],Nb_Phos_Kinetics[:,5],linewidth = width,label = "430 nm")
        lines!(Nb_Phos_Kinetics[:,1],Nb_Phos_Kinetics[:,6],linewidth = width,label = "535 nm")
        lines!(Nb_Phos_Kinetics[:,1],Nb_Phos_Kinetics[:,7],linewidth = width,label = "701 nm")





        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Phos_TA_Kinetics()
save("JULIA/Outputs/2024_07_18_Nb_Phos_TA_Kinetics.png",Nb_Phos_TA_Kinetics())

# Nitrobenzene queched by DIPEA
####################################################################################################

function Nb_DIPEA_TA_Spectra(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitrobenzene DIPEA Flash Quench TA Spectra in CPME",
            palette = (color = palette(ColorSchemes.inferno, 10),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        # lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,3],linewidth = width,label = "0 ns")
        # lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,4],linewidth = width,label = "10 ns")
        # lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,5],linewidth = width,label = "50 ns")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,6],linewidth = width,label = "10 μs")
        # lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,7],linewidth = width,label = "20 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,8],linewidth = width,label = "50 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,9],linewidth = width,label = "100 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,10],linewidth = width,label = "200 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,11],linewidth = width,label = "300 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,12],linewidth = width,label = "400 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,13],linewidth = width,label = "500 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,14],linewidth = width,label = "750 μs")
        lines!(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,15],linewidth = width,label = "1000 μs")





        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_DIPEA_TA_Spectra()
save("JULIA/Outputs/2024_07_18_Nb_DIPEA_TA_Spectra.png",Nb_DIPEA_TA_Spectra())

function Nb_DIPEA_TA_Kinetics(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronaphthalene TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (μs)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        

        lines!(Nb_DIPEA_Kinetics_427nm[:,1] ./ 1000,Nb_DIPEA_Kinetics_427nm[:,2],linewidth = width,label = "427 nm, τ = 305 μs")
        lines!(Nb_DIPEA_Kinetics_427nm[:,1]./ 1000,Nb_DIPEA_Kinetics_427nm[:,3],linewidth = 3)
        

        lines!(Nb_DIPEA_Kinetics_468nm[:,1]./ 1000,Nb_DIPEA_Kinetics_468nm[:,2],linewidth = width,label = "468 nm, τ = 211 μs")
        lines!(Nb_DIPEA_Kinetics_468nm[:,1]./ 1000,Nb_DIPEA_Kinetics_468nm[:,3],linewidth = 3)

        # lines!(Nb_DIPEA_Kinetics_427nm[:,1]./ 1000,MaxNorm(Nb_DIPEA_Kinetics_427nm[:,3]),linewidth = 3)
        # lines!(Nb_DIPEA_Kinetics_468nm[:,1]./ 1000,MaxNorm(Nb_DIPEA_Kinetics_468nm[:,3]),linewidth = 3)



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_DIPEA_TA_Kinetics()
save("JULIA/Outputs/2024_07_18_Nb_DIPEA_TA_Kinetics.png",Nb_DIPEA_TA_Kinetics())

####################################################################################################

# NitroNaphthalene 
####################################################################################################

function NNaph_TA_Spectra(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronaphthalene TA Spectra in CPME",
            palette = (color = palette(ColorSchemes.inferno, 9),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,2],linewidth = width,label = "-2 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,3],linewidth = width,label = "0 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,4],linewidth = width,label = "1 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,5],linewidth = width,label = "2 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,6],linewidth = width,label = "3 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,7],linewidth = width,label = "5 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,8],linewidth = width,label = "10 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,9],linewidth = width,label = "50 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,10],linewidth = width,label = "100 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,11],linewidth = width,label = "200 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,12],linewidth = width,label = "300 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,13],linewidth = width,label = "400 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,14],linewidth = width,label = "500 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,15],linewidth = width,label = "750 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,16],linewidth = width,label = "1000 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,17],linewidth = width,label = "2000 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,18],linewidth = width,label = "2500 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,19],linewidth = width,label = "3000 ns")
        # lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,20],linewidth = width,label = "4000 ns")
        lines!(NNaph_Spectra[:,1],NNaph_Spectra[:,21],linewidth = width,label = "5000 ns")



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
NNaph_TA_Spectra()
save("JULIA/Outputs/2024_07_18_NNaph_TA_Spectra.png",NNaph_TA_Spectra())

function NNaph_TA_Kinetics(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronaphthalene TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 7),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        # lines!(NNaph_Kinetics[:,1],NNaph_Kinetics[:,2],linewidth = width,label = "400 nm")
        # lines!(NNaph_Kinetics[:,1],NNaph_Kinetics[:,3],linewidth = width,label = "570 nm")
        # lines!(NNaph_Kinetics[:,1],NNaph_Kinetics[:,4],linewidth = width,label = "815 nm")

        lines!(NNaph_Kinetics_400nm[:,1],NNaph_Kinetics_400nm[:,2],linewidth = width,label = "400 nm, τ = 1917 ns")
        lines!(NNaph_Kinetics_400nm[:,1],NNaph_Kinetics_400nm[:,3],linewidth = 3)
        

        lines!(NNaph_Kinetics_570nm[:,1],NNaph_Kinetics_570nm[:,2],linewidth = width,label = "570 nm, τ = 1944 ns")
        lines!(NNaph_Kinetics_570nm[:,1],NNaph_Kinetics_570nm[:,3],linewidth = 3)

        lines!(NNaph_Kinetics_815nm[:,1],NNaph_Kinetics_815nm[:,2],linewidth = width,label = "815 nm, τ = 1896 ns")
        lines!(NNaph_Kinetics_815nm[:,1],NNaph_Kinetics_815nm[:,3],linewidth = 3)




        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
NNaph_TA_Kinetics()
save("JULIA/Outputs/2024_07_18_NNaph_TA_Kinetics.png",NNaph_TA_Kinetics())

####################################################################################################

function Overlay(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronaphthalene TA Kinetics in CPME",
            palette = (color = palette(ColorSchemes.rainbow, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(400,850)
        CairoMakie.ylims!(-0.15,1.2)

        zNorm = 415

        lines!(Nb_Si_Spectra[:,1],Norm(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,5],zNorm),linewidth = width,label = " Si 4 ns")
        lines!(Nb_Phos_Spectra[:,1],Norm(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,6],zNorm),linewidth = width,label = "³NB + NB•⁻")  #NB + Phos 2 ns
        lines!(Nb_DIPEA_Spectra[:,1],Norm(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,6],zNorm),linewidth = width,label = "NB•⁻")

        lines!(Nb_Spectra[:,1],Norm(Nb_Spectra[:,1],Nb_Spectra[:,4],zNorm),linewidth = width,label = "³NB")

        # lines!(Nb_Si_Spectra[:,1],Norm(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,5],zNorm),linewidth = width,label = " Nitrobenzene + Si")



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend( position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Overlay()
save("JULIA/Outputs/2024_07_18_Overlay.png",Overlay())

function Nb_Si_Overlay(f = Figure())
    width = 4
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitrobenzene with and without Silane at 1 ns delay",
            palette = (color = palette(ColorSchemes.inferno, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "Norm. ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(400,850)
        CairoMakie.ylims!(-0.1,1.2)

        zNorm = 415

        lines!(Nb_Spectra[:,1],Norm(Nb_Spectra[:,1],Nb_Spectra[:,4],zNorm),linewidth = width,label = "Nitrobenzne")
        lines!(Nb_Si_Spectra[:,1],Norm(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,5],zNorm),linewidth = width,label = "Nitrobenzene + Si")



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend( position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Nb_Si_Overlay()
save("JULIA/Outputs/2024_07_18_Nb_Si_Overlay.png",Nb_Si_Overlay())

function NB_Quench_Overlay(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Quenching of Nitrobenzene by Phosphine",
            palette = (color = palette(ColorSchemes.inferno, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "Norm. ΔA"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(400,850)
        CairoMakie.ylims!(-0.15,1.2)

        zNorm = 415

        lines!(Nb_Si_Spectra[:,1],Norm(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,19],420),linewidth = width,label = " NB + Si 5 ns")
        lines!(Nb_Phos_Spectra[:,1],Norm(Nb_Phos_Spectra[:,1],Nb_Phos_Spectra[:,6],415),linewidth = width,label = "NB + Phos")  #NB + Phos 2 ns
        lines!(Nb_DIPEA_Spectra[:,1],Norm(Nb_DIPEA_Spectra[:,1],Nb_DIPEA_Spectra[:,6],430),linewidth = width,label = "NB•⁻")

        # lines!(Nb_Spectra[:,1],Norm(Nb_Spectra[:,1],Nb_Spectra[:,4],zNorm),linewidth = width,label = "³NB")



        # lines!(Nb_Si_Spectra[:,1],Norm(Nb_Si_Spectra[:,1],Nb_Si_Spectra[:,5],zNorm),linewidth = width,label = " Nitrobenzene + Si")



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend( position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
NB_Quench_Overlay()
save("JULIA/Outputs/2024_07_18_NB_Quench_Overlay.png",NB_Quench_Overlay())