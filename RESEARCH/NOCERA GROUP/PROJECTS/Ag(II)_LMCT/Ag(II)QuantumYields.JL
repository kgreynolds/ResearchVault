include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

function Flux(vector,Wavelength)
    # Calculates the mean flux in moles photons per second of a 50 Hz pulsed laser at some wavelenght
    h = 6.62607015e-34
    c = 2.99792458e+8
    NA = 6.0221408e23 
    
    E_photon = h*c/(Wavelength*10^(-9))

    Avg_energy = mean(vector) # J per pulse
    Avg_Power = Avg_energy*50 # power in Watts at 50 Hz rep rate
    Avg_Flux = Avg_Power/E_photon/NA # mol photons /s
    return Avg_Flux
end

function QY_method_2(k, ϵ, volume, init_flux, init_A)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value of the sample
    I_abs = init_flux*(1-10^(-init_A))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end

function QY_method_3(k, ϵ_450, volume, init_flux, init_A_pump, init_A_450)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value of the sample
    I_abs = init_flux*(1-10^(-init_A_pump))/volume

    Quantum_yield = k * init_A_450 / (I_abs * ϵ_450) * 100

    return Quantum_yield
end

#### Sample #########################################################################################################

    ## Sample UvVis
    Sample_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/Ferrioxalate_Absorbance__0__10-31-59-317.csv",header = false, DataFrame))
    Sample_UvVis[1,120]
    Sample_UvVis[1,139]
    Sample_UvVis[1,197]

    Sample_UvVis_Kinetics = hcat(Sample_UvVis[15:end,1], Sample_UvVis[15:end,120])

    Sample_Absorbances = transpose(Sample_UvVis)
    Sample_initial_Abs_450 = Sample_Absorbances[197,2]
    Sample_initial_Abs_pump = Sample_Absorbances[139,2]

    Plots.plot(Sample_Absorbances[2:end,1],Sample_Absorbances[2:end,2:end])
    Plots.plot(Sample_UvVis[15:end,1],Sample_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    Sample_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/afterstage-blank.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Sample_Blank[:,1],Sample_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Sample_Blank_avg = mean(Sample_Blank[:,2])

    ## Fitting for rate constant
    Sample_Fit = MonoFitIRF_New(Sample_UvVis_Kinetics,2,0.0, 400, [0.08, 0.0001, 132.659, 210.211, 0.031], 10,Inf,false,10000)
    k_Sample = 1/210.211

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Sample = Flux(Sample_Blank_avg,450)

    QY_Sample= QY_method_3(k_Sample, 1297, 0.003,FLux_Blank_Sample,Sample_initial_Abs_pump, Sample_initial_Abs_450)


    function Sample_Fig(f = Figure())

        Data = PhAuCl3_Photolyzed_TA_Spectra
    
        ### First Panel ##############################################################################
    
            ax =    CairoMakie.Axis(
                    f[1,1], palette = (color = palette(ColorSchemes.inferno, ncol(Data)+1),), ## or replace ncol with Integer 
                    title = "Title",
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
    Sample_Fig()
    save("JULIA/Outputs/Sample_Fig.pdf",Sample_Fig())
    