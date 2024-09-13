include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Outside(),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)

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
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value
    I_abs = init_flux*(1-10^(-init_A))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end

function QY_method_3(k, ϵ, volume, init_flux, init_A_pump, init_A)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value 
    I_abs = init_flux*(1-10^(-init_A_pump))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end

### Summary ###### RUN LINES BELOW First

## Mean QY for Ag_Chloro at 525nm 
mean_QY_Ag_Chloro_525 = mean([QY_Ag_Chloro_525nm_s1,QY_Ag_Chloro_525nm_s2,QY_Ag_Chloro_525nm_s3])
stdev_QY_Ag_Chloro_525 = std([QY_Ag_Chloro_525nm_s1,QY_Ag_Chloro_525nm_s2,QY_Ag_Chloro_525nm_s3])

## Mean QY for Ag_Chloro at 450nm 
mean_QY_Ag_Chloro_450 = mean([QY_Ag_Chloro_450nm_s1,QY_Ag_Chloro_450nm_s2,QY_Ag_Chloro_450nm_s3])
stdev_QY_Ag_Chloro_450 = std([QY_Ag_Chloro_450nm_s1,QY_Ag_Chloro_450nm_s2,QY_Ag_Chloro_450nm_s3])

## Mean QY for Ag_Chloro at 450nm 
mean_QY_Ag_Chloro_405 = mean([QY_Ag_Chloro_405nm_s1,QY_Ag_Chloro_405nm_s2,QY_Ag_Chloro_405nm_s3])
stdev_QY_Ag_Chloro_405 = std([QY_Ag_Chloro_405nm_s1,QY_Ag_Chloro_405nm_s2,QY_Ag_Chloro_405nm_s3])

function Ag_Chloro_QY_Action_Spectrum(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax =    CairoMakie.Axis(
            f[1,1], palette = (color = palette(ColorSchemes.inferno, 3),), ## or replace ncol with Integer 
            title = "Ag(II) Chloro Action Spectrum",
            xlabel = "Wavelength (nm)", 
            ylabel = "QY", 
        )

        width = 3
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(nothing,600)
        # CairoMakie.ylims!(-0.001,nothing)
        
        
        CairoMakie.errorbars!([405,450,525],[mean_QY_Ag_Chloro_405,mean_QY_Ag_Chloro_450,mean_QY_Ag_Chloro_525], [stdev_QY_Ag_Chloro_405,stdev_QY_Ag_Chloro_450,stdev_QY_Ag_Chloro_525], color = :red, whiskerwidth = 10)
        CairoMakie.scatter!([405,450,525],[mean_QY_Ag_Chloro_405,mean_QY_Ag_Chloro_450,mean_QY_Ag_Chloro_525], markersize = 10, color = :black)



        # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


        # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        # axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)
        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        resize_to_layout!(f)
   f
        
end
Ag_Chloro_QY_Action_Spectrum()
save("JULIA/Outputs/Ag_Chloro_QY_Action_Spectrum.pdf",Ag_Chloro_QY_Action_Spectrum())

#### FerriOxalate #########################################################################################################

    ## FerriOxalate UvVis
    FerriOxalate_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/Ferrioxalate_Absorbance__0__10-31-59-317.csv",header = false, DataFrame))
    FerriOxalate_UvVis[1,134]

    FerriOxalate_UvVis_Kinetics = hcat(FerriOxalate_UvVis[5:end,1], FerriOxalate_UvVis[5:end,134])

    FerriOxalate_Absorbances = transpose(FerriOxalate_UvVis)

    Plots.plot(FerriOxalate_Absorbances[2:end,1],FerriOxalate_Absorbances[2:end,2:end])
    Plots.plot(FerriOxalate_UvVis[15:end,1],FerriOxalate_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    FerriOxalate_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/afterstage-blank.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(FerriOxalate_Blank[:,1],FerriOxalate_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    FerriOxalate_Blank_avg = mean(FerriOxalate_Blank[:,2])

    ## Fitting for rate constant
    FerriOxalate_Fit = MonoFitIRF_New(FerriOxalate_UvVis_Kinetics,2,0.0, 400, [0.08, 0.0001, 132.659, 210.211, 0.031], 10,Inf,false,10000)
    k = 1/210.211

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_FerriOxalate = Flux(FerriOxalate_Blank_avg,436)

    QY_FerriOxalate= QY_method_2(k, 50.26, 0.002,FLux_Blank_FerriOxalate,0.15078)


    function FerriOxalate_Abs_decay(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "Ferrioxalate, 3mM QY",
                palette = (color = palette(ColorSchemes.inferno, 7),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA"
            )
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            lines!(FerriOxalate_Absorbances[2:end,1],FerriOxalate_Absorbances[2:end,2],linewidth = width, label = "1 s")




            CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


            # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    FerriOxalate_Abs_decay()
    save("JULIA/Outputs/FerriOxalate_Abs_decay.pdf",FerriOxalate_Abs_decay())

#### Ag_Chloro_450nm_s1 #########################################################################################################

    ## ϵ_450nm = 1305 for Ag(II)_Chloro
    ## Ag_Chloro_450nm_s1 UvVis
    Ag_Chloro_450nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv",header = false, DataFrame))
 
    Ag_Chloro_450nm_s1_UvVis[1,134]
    Ag_Chloro_450nm_s1_UvVis[1,197]

    Ag_Chloro_450nm_s1_UvVis_Kinetics = hcat(Ag_Chloro_450nm_s1_UvVis[2:end,1], Ag_Chloro_450nm_s1_UvVis[2:end,134])

    Ag_Chloro_450nm_s1_Absorbances = transpose(Ag_Chloro_450nm_s1_UvVis)


    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header = true, DataFrame)
    Ag_Chloro_450nm_SteadyState_Uvvis[151,1]
    Ag_Chloro_450nm_s1_initial_Abs_450 = Ag_Chloro_450nm_SteadyState_Uvvis[151,3] - Ag_Chloro_450nm_SteadyState_Uvvis[151,2]
    Ag_Chloro_450nm_s1_initial_Abs_pump = Ag_Chloro_450nm_SteadyState_Uvvis[151,3] - Ag_Chloro_450nm_SteadyState_Uvvis[151,2]

    # Plots.plot(Ag_Chloro_450nm_s1_Absorbances[2:end,1],Ag_Chloro_450nm_s1_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_450nm_s1_UvVis[2:end,1],Ag_Chloro_450nm_s1_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    Ag_Chloro_450nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_450nm_s1_Blank[:,1],Ag_Chloro_450nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_450nm_s1_Blank_avg = mean(Ag_Chloro_450nm_s1_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_450nm_s1_Fit = MonoFitIRF_New(Ag_Chloro_450nm_s1_UvVis_Kinetics,2,5, 600, [0.08, 0.001, 20.335, 129.17, 0.031], 10,Inf,false,10000)
    k_Ag_Chloro_450nm_s1 = 1/129.17

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_450nm_s1 = Flux(Ag_Chloro_450nm_s1_Blank_avg,450)

    QY_Ag_Chloro_450nm_s1= QY_method_3(k_Ag_Chloro_450nm_s1, 1305, 0.003,FLux_Blank_Ag_Chloro_450nm_s1,Ag_Chloro_450nm_s1_initial_Abs_pump, Ag_Chloro_450nm_s1_initial_Abs_450)



#### Ag_Chloro_450nm_s2 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## Ag_Chloro_450nm_s2 UvVis
    Ag_Chloro_450nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s2/Absorbance_AgTGBP_Cl_450nm_s2.csv",header = false, DataFrame))
    
    Ag_Chloro_450nm_s2_UvVis[1,134]
    Ag_Chloro_450nm_s2_UvVis[1,197]

    Ag_Chloro_450nm_s2_UvVis_Kinetics = hcat(Ag_Chloro_450nm_s2_UvVis[2:end,1], Ag_Chloro_450nm_s2_UvVis[2:end,134])

    Ag_Chloro_450nm_s2_Absorbances = transpose(Ag_Chloro_450nm_s2_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header = true, DataFrame)
    Ag_Chloro_450nm_s2_initial_Abs_450 = Ag_Chloro_450nm_SteadyState_Uvvis[151,4] - Ag_Chloro_450nm_SteadyState_Uvvis[151,2]
    Ag_Chloro_450nm_s2_initial_Abs_pump = Ag_Chloro_450nm_SteadyState_Uvvis[151,4] - Ag_Chloro_450nm_SteadyState_Uvvis[151,2]

    # Plots.plot(Ag_Chloro_450nm_s2_Absorbances[2:end,1],Ag_Chloro_450nm_s2_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_450nm_s2_UvVis[2:end,1],Ag_Chloro_450nm_s2_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    Ag_Chloro_450nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_450nm_s2_Blank[:,1],Ag_Chloro_450nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_450nm_s2_Blank_avg = mean(Ag_Chloro_450nm_s2_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_450nm_s2_Fit = MonoFitIRF_New(Ag_Chloro_450nm_s2_UvVis_Kinetics,2,5.0, 600, [0.08, 0.001, 27.396, 135.575, -0.002], 10,Inf,false,10000)
    k_Ag_Chloro_450nm_s2 = 1/135.575

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_450nm_s2 = Flux(Ag_Chloro_450nm_s2_Blank_avg,450)

    QY_Ag_Chloro_450nm_s2= QY_method_3(k_Ag_Chloro_450nm_s2, 1305, 0.003,FLux_Blank_Ag_Chloro_450nm_s2,Ag_Chloro_450nm_s2_initial_Abs_pump, Ag_Chloro_450nm_s2_initial_Abs_450)



#### Ag_Chloro_450nm_s3 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## Ag_Chloro_450nm_s3 UvVis
    Ag_Chloro_450nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s3/Absorbance_AgTGBP_Cl_450nm_s3.csv",header = false, DataFrame))

    Ag_Chloro_450nm_s3_UvVis[1,134]
    Ag_Chloro_450nm_s3_UvVis[1,197]

    Ag_Chloro_450nm_s3_UvVis_Kinetics = hcat(Ag_Chloro_450nm_s3_UvVis[2:end,1], Ag_Chloro_450nm_s3_UvVis[2:end,134])

    Ag_Chloro_450nm_s3_Absorbances = transpose(Ag_Chloro_450nm_s3_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header = true, DataFrame)
    Ag_Chloro_450nm_s3_initial_Abs_450 = Ag_Chloro_450nm_SteadyState_Uvvis[151,5] - Ag_Chloro_450nm_SteadyState_Uvvis[151,2]
    Ag_Chloro_450nm_s3_initial_Abs_pump = Ag_Chloro_450nm_SteadyState_Uvvis[151,5] - Ag_Chloro_450nm_SteadyState_Uvvis[151,2]

    Plots.plot(Ag_Chloro_450nm_s3_Absorbances[2:end,1],Ag_Chloro_450nm_s3_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_450nm_s3_UvVis[2:end,1],Ag_Chloro_450nm_s3_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    Ag_Chloro_450nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_450nm_s3_Blank[:,1],Ag_Chloro_450nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_450nm_s3_Blank_avg = mean(Ag_Chloro_450nm_s3_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_450nm_s3_Fit = MonoFitIRF_New(Ag_Chloro_450nm_s3_UvVis_Kinetics,2,10.0, 800, [0.08, 0.001, 21.63, 130.265, 0.02], 10,Inf,false,10000)
    k_Ag_Chloro_450nm_s3 = 1/ 130.265

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_450nm_s3 = Flux(Ag_Chloro_450nm_s3_Blank_avg,450)

    QY_Ag_Chloro_450nm_s3= QY_method_3(k_Ag_Chloro_450nm_s3, 1305, 0.003,FLux_Blank_Ag_Chloro_450nm_s3,Ag_Chloro_450nm_s3_initial_Abs_pump, Ag_Chloro_450nm_s3_initial_Abs_450)



#### Ag_Chloro_405nm_s1 #########################################################################################################

    # ϵ_405nm = 2203.7 for Ag(II)TFA
    ## Ag_Chloro_405nm_s1 UvVis
    Ag_Chloro_405nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/BC-4-031/s1/Ag_Chloro_s1_Absorbances_405nm.csv",header = false, DataFrame))
    Ag_Chloro_405nm_s1_UvVis[1,108]
    Ag_Chloro_405nm_s1_UvVis[1,139]


    Ag_Chloro_405nm_s1_UvVis_Kinetics = hcat(Ag_Chloro_405nm_s1_UvVis[2:end,1], Ag_Chloro_405nm_s1_UvVis[2:end,108])

    Ag_Chloro_405nm_s1_Absorbances = transpose(Ag_Chloro_405nm_s1_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_405nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/SteadyState_UvVis_Ag_Chloro_405nm.csv", header = true, DataFrame)
    Ag_Chloro_405nm_SteadyState_Uvvis[251,1]
    Ag_Chloro_405nm_s1_initial_Abs_450 = Ag_Chloro_405nm_SteadyState_Uvvis[251,3] - Ag_Chloro_405nm_SteadyState_Uvvis[251,2]
    Ag_Chloro_405nm_s1_initial_Abs_pump = Ag_Chloro_405nm_SteadyState_Uvvis[206,3] - Ag_Chloro_405nm_SteadyState_Uvvis[206,2]

    Plots.plot(Ag_Chloro_405nm_s1_Absorbances[2:end,1],Ag_Chloro_405nm_s1_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_405nm_s1_UvVis[2:end,1],Ag_Chloro_405nm_s1_UvVis[2:end,108])

    ## Per-pulse energy after the Blank
    Ag_Chloro_405nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/Ag_Chloro_Blank_Energy_405nm.csv", skipto = 6, header = true, DataFrame)
    
    Plots.plot(Ag_Chloro_405nm_s1_Blank[:,1],Ag_Chloro_405nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_405nm_s1_Blank_avg = mean(Ag_Chloro_405nm_s1_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_405nm_s1_Fit = MonoFitIRF_New(Ag_Chloro_405nm_s1_UvVis_Kinetics,2,10.0, 600, [0.08, 0.001, 9.801, 71.984, -0.007], 10,Inf,false,10000)
    k_Ag_Chloro_405nm_s1 = 1/  71.984

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_405nm_s1 = Flux(Ag_Chloro_405nm_s1_Blank_avg,405)

    QY_Ag_Chloro_405nm_s1= QY_method_3(k_Ag_Chloro_405nm_s1, 1305, 0.003,FLux_Blank_Ag_Chloro_405nm_s1,Ag_Chloro_405nm_s1_initial_Abs_pump, Ag_Chloro_405nm_s1_initial_Abs_450)


 

#### Ag_Chloro_405nm_s2 #########################################################################################################

    # ϵ_405nm = 2203.7 for Ag(II)TFA
    ## Ag_Chloro_405nm_s2 UvVis
    Ag_Chloro_405nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/BC-4-031/s2/Ag_Chloro_s2_Absorbances_405nm.csv",header = false, DataFrame))

    Ag_Chloro_405nm_s2_UvVis[1,139]
    Ag_Chloro_405nm_s2_UvVis[1,108]

    Ag_Chloro_405nm_s2_UvVis_Kinetics = hcat(Ag_Chloro_405nm_s2_UvVis[2:end,1], Ag_Chloro_405nm_s2_UvVis[2:end,108])

    Ag_Chloro_405nm_s2_Absorbances = transpose(Ag_Chloro_405nm_s2_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_405nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/SteadyState_UvVis_Ag_Chloro_405nm.csv", header = true, DataFrame)
    Ag_Chloro_405nm_SteadyState_Uvvis[206,1]
    Ag_Chloro_405nm_s2_initial_Abs_450 = Ag_Chloro_405nm_SteadyState_Uvvis[251,4] - Ag_Chloro_405nm_SteadyState_Uvvis[251,2]
    Ag_Chloro_405nm_s2_initial_Abs_pump = Ag_Chloro_405nm_SteadyState_Uvvis[206,4] - Ag_Chloro_405nm_SteadyState_Uvvis[206,2]

    Plots.plot(Ag_Chloro_405nm_s2_Absorbances[2:end,1],Ag_Chloro_405nm_s2_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_405nm_s2_UvVis[2:end,1],Ag_Chloro_405nm_s2_UvVis[2:end,108])

    ## Per-pulse energy after the Blank
    Ag_Chloro_405nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/Ag_Chloro_Blank_Energy_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_405nm_s2_Blank[:,1],Ag_Chloro_405nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_405nm_s2_Blank_avg = mean(Ag_Chloro_405nm_s2_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_405nm_s2_Fit = MonoFitIRF_New(Ag_Chloro_405nm_s2_UvVis_Kinetics,2,10.0, 400, [0.08, 0.001, 9.606, 67.124, -0.018], 10,Inf,false,10000)
    k_Ag_Chloro_405nm_s2 = 1/67.124

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_405nm_s2 = Flux(Ag_Chloro_405nm_s2_Blank_avg,405)

    QY_Ag_Chloro_405nm_s2= QY_method_3(k_Ag_Chloro_405nm_s2, 1305, 0.003,FLux_Blank_Ag_Chloro_405nm_s2,Ag_Chloro_405nm_s2_initial_Abs_pump, Ag_Chloro_405nm_s2_initial_Abs_450)



#### Ag_Chloro_405nm_s3 #########################################################################################################

    ## Ag_Chloro_405nm_s3 UvVis
    Ag_Chloro_405nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/BC-4-031/s3/Ag_Chloro_s3_Absorbances_405nm.csv",header = false, DataFrame))

    Ag_Chloro_405nm_s3_UvVis[1,108]
    Ag_Chloro_405nm_s3_UvVis[1,139]

    Ag_Chloro_405nm_s3_UvVis_Kinetics = hcat(Ag_Chloro_405nm_s3_UvVis[2:end,1], Ag_Chloro_405nm_s3_UvVis[2:end,108])

    Ag_Chloro_405nm_s3_Absorbances = transpose(Ag_Chloro_405nm_s3_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_405nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/SteadyState_UvVis_Ag_Chloro_405nm.csv", header = true, DataFrame)
    Ag_Chloro_405nm_s3_initial_Abs_450 = Ag_Chloro_405nm_SteadyState_Uvvis[251,5] - Ag_Chloro_405nm_SteadyState_Uvvis[251,2]
    Ag_Chloro_405nm_s3_initial_Abs_pump = Ag_Chloro_405nm_SteadyState_Uvvis[206,5] - Ag_Chloro_405nm_SteadyState_Uvvis[206,2]

    Plots.plot(Ag_Chloro_405nm_s3_Absorbances[2:end,1],Ag_Chloro_405nm_s3_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_405nm_s3_UvVis[2:end,1],Ag_Chloro_405nm_s3_UvVis[2:end,108])

    ## Per-pulse energy after the Blank
    Ag_Chloro_405nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_22/Ag_Chloro_405nm_QY/Ag_Chloro_Blank_Energy_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_405nm_s3_Blank[:,1],Ag_Chloro_405nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_405nm_s3_Blank_avg =  mean(Ag_Chloro_405nm_s3_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_405nm_s3_Fit = MonoFitIRF_New(Ag_Chloro_405nm_s3_UvVis_Kinetics,2,0.0, 400, [0.08, 0.001, 6.008, 61.176, -0.006], 10,Inf,false,10000)
    k_Ag_Chloro_405nm_s3 = 1/61.176

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_405nm_s3 = Flux(Ag_Chloro_405nm_s3_Blank_avg,405)

    QY_Ag_Chloro_405nm_s3= QY_method_3(k_Ag_Chloro_405nm_s3, 1305, 0.003,FLux_Blank_Ag_Chloro_405nm_s3,Ag_Chloro_405nm_s3_initial_Abs_pump, Ag_Chloro_405nm_s3_initial_Abs_450)


#### Ag_Chloro_525nm_s1 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## Tracked Kinetics at 400nm, calculated concentration at 450nm
    ## Ag_Chloro_525nm_s1 UvVis
    Ag_Chloro_525nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/s1/AgTGBP_Cl_525nm_s1_Absorbances.csv",header = false, DataFrame))
    Ag_Chloro_525nm_s1_UvVis[1,293]
    Ag_Chloro_525nm_s1_UvVis[1,197]
    Ag_Chloro_525nm_s1_UvVis[1,134]

    Ag_Chloro_525nm_s1_UvVis_Kinetics = hcat(Ag_Chloro_525nm_s1_UvVis[2:end,1], Ag_Chloro_525nm_s1_UvVis[2:end,134])

    Ag_Chloro_525nm_s1_Absorbances = transpose(Ag_Chloro_525nm_s1_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_525nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/AgTGBP_Cl_Steadystate_UvVis_525nmex.csv", header = true, DataFrame)
    Ag_Chloro_525nm_s1_initial_Abs_450 = Ag_Chloro_525nm_SteadyState_Uvvis[151,3] - Ag_Chloro_525nm_SteadyState_Uvvis[151,2]
    Ag_Chloro_525nm_s1_initial_Abs_pump = Ag_Chloro_525nm_SteadyState_Uvvis[226,3] - Ag_Chloro_525nm_SteadyState_Uvvis[226,2]

    # Plots.plot(Ag_Chloro_525nm_s1_Absorbances[2:end,1],Ag_Chloro_525nm_s1_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_525nm_s1_UvVis[2:end,1],Ag_Chloro_525nm_s1_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    Ag_Chloro_525nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/Poststage_blank_power_525nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_525nm_s1_Blank[:,1],Ag_Chloro_525nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_525nm_s1_Blank_avg = mean(Ag_Chloro_525nm_s1_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_525nm_s1_Fit = MonoFitIRF_New(Ag_Chloro_525nm_s1_UvVis_Kinetics,2,100.0, 1300, [0.08, 0.001, 303.679, 632.312, 0.017], 10,Inf,false,10000)
    k_Ag_Chloro_525nm_s1 = 1/632.312

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_525nm_s1 = Flux(Ag_Chloro_525nm_s1_Blank_avg,525)

    QY_Ag_Chloro_525nm_s1= QY_method_3(k_Ag_Chloro_525nm_s1, 1305, 0.003,FLux_Blank_Ag_Chloro_525nm_s1,Ag_Chloro_525nm_s1_initial_Abs_pump, Ag_Chloro_525nm_s1_initial_Abs_450)


#### Ag_Chloro_525nm_s2 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## Ag_Chloro_525nm_s2 UvVis
    Ag_Chloro_525nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/s2/AgTGBP_Cl_525nm_s2_Absorbances.csv",header = false, DataFrame))
    Ag_Chloro_525nm_s2_UvVis[1,293]
    Ag_Chloro_525nm_s2_UvVis[1,134]
    Ag_Chloro_525nm_s2_UvVis[1,197]

    Ag_Chloro_525nm_s2_UvVis_Kinetics = hcat(Ag_Chloro_525nm_s2_UvVis[2:end,1], Ag_Chloro_525nm_s2_UvVis[2:end,134])

    Ag_Chloro_525nm_s2_Absorbances = transpose(Ag_Chloro_525nm_s2_UvVis)

    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_525nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/AgTGBP_Cl_Steadystate_UvVis_525nmex.csv", header = true, DataFrame)
    Ag_Chloro_525nm_s2_initial_Abs_450 = Ag_Chloro_525nm_SteadyState_Uvvis[151,4] - Ag_Chloro_525nm_SteadyState_Uvvis[151,2]
    Ag_Chloro_525nm_s2_initial_Abs_pump = Ag_Chloro_525nm_SteadyState_Uvvis[226,4] - Ag_Chloro_525nm_SteadyState_Uvvis[226,2]

    # Plots.plot(Ag_Chloro_525nm_s2_Absorbances[2:end,1],Ag_Chloro_525nm_s2_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_525nm_s2_UvVis[2:end,1],Ag_Chloro_525nm_s2_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    Ag_Chloro_525nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/Poststage_blank_power_525nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_525nm_s2_Blank[:,1],Ag_Chloro_525nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_525nm_s2_Blank_avg = mean(Ag_Chloro_525nm_s2_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_525nm_s2_Fit = MonoFitIRF_New(Ag_Chloro_525nm_s2_UvVis_Kinetics,2,100.0, 3000, [0.08, 0.001, 320.319, 636.647, -0.005], 10,Inf,false,10000)
    k_Ag_Chloro_525nm_s2 = 1/636.647

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_525nm_s2 = Flux(Ag_Chloro_525nm_s2_Blank_avg,525)

    QY_Ag_Chloro_525nm_s2= QY_method_3(k_Ag_Chloro_525nm_s2, 1305, 0.003,FLux_Blank_Ag_Chloro_525nm_s2,Ag_Chloro_525nm_s2_initial_Abs_pump, Ag_Chloro_525nm_s2_initial_Abs_450)



#### Ag_Chloro_525nm_s3 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## Ag_Chloro_525nm_s3 UvVis
    Ag_Chloro_525nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/s3/AgTGBP_Cl_525nm_s3_Absorbances.csv",header = false, DataFrame))
    Ag_Chloro_525nm_s3_UvVis[1,293]
    Ag_Chloro_525nm_s3_UvVis[1,134]
    Ag_Chloro_525nm_s3_UvVis[1,197]

    Ag_Chloro_525nm_s3_UvVis_Kinetics = hcat(Ag_Chloro_525nm_s3_UvVis[2:end,1], Ag_Chloro_525nm_s3_UvVis[2:end,134])

    Ag_Chloro_525nm_s3_Absorbances = transpose(Ag_Chloro_525nm_s3_UvVis)
    
    #### SteadyState UvVis #########################################################################################################
    Ag_Chloro_525nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/AgTGBP_Cl_Steadystate_UvVis_525nmex.csv", header = true, DataFrame)
    Ag_Chloro_525nm_SteadyState_Uvvis[226,1]
    Ag_Chloro_525nm_s3_initial_Abs_450 = Ag_Chloro_525nm_SteadyState_Uvvis[151,5] - Ag_Chloro_525nm_SteadyState_Uvvis[151,2]
    Ag_Chloro_525nm_s3_initial_Abs_pump = Ag_Chloro_525nm_SteadyState_Uvvis[226,5] - Ag_Chloro_525nm_SteadyState_Uvvis[226,2]
    # Plots.plot(Ag_Chloro_525nm_s3_Absorbances[2:end,1],Ag_Chloro_525nm_s3_Absorbances[2:end,2:end])
    Plots.plot(Ag_Chloro_525nm_s3_UvVis[2:end,1],Ag_Chloro_525nm_s3_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    Ag_Chloro_525nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/Poststage_blank_power_525nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Ag_Chloro_525nm_s3_Blank[:,1],Ag_Chloro_525nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    Ag_Chloro_525nm_s3_Blank_avg = mean(Ag_Chloro_525nm_s3_Blank[:,2])

    ## Fitting for rate constant
    Ag_Chloro_525nm_s3_Fit = MonoFitIRF_New(Ag_Chloro_525nm_s3_UvVis_Kinetics,2,100.0, 3000, [0.08, 0.001, 300.092, 627.175, -0.005], 10,Inf,false,10000)
    k_Ag_Chloro_525nm_s3 = 1/ 627.175

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_Ag_Chloro_525nm_s3 = Flux(Ag_Chloro_525nm_s3_Blank_avg,525)

    QY_Ag_Chloro_525nm_s3= QY_method_3(k_Ag_Chloro_525nm_s3, 1305, 0.003,FLux_Blank_Ag_Chloro_525nm_s3,Ag_Chloro_525nm_s3_initial_Abs_pump, Ag_Chloro_525nm_s3_initial_Abs_450)

