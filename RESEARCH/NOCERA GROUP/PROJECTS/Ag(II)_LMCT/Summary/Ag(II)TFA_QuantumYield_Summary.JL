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

## Mean QY for AgTFA at 525nm 
mean_QY_AgTFA_525 = mean([QY_AgTFA_525nm_s1,QY_AgTFA_525nm_s2,QY_AgTFA_525nm_s3])
stdev_QY_AgTFA_525 = std([QY_AgTFA_525nm_s1,QY_AgTFA_525nm_s2,QY_AgTFA_525nm_s3])

## Mean QY for AgTFA at 450nm 
mean_QY_AgTFA_450 = mean([QY_AgTFA_450nm_s1,QY_AgTFA_450nm_s2,QY_AgTFA_450nm_s3])
stdev_QY_AgTFA_450 = std([QY_AgTFA_450nm_s1,QY_AgTFA_450nm_s2,QY_AgTFA_450nm_s3])

## Mean QY for AgTFA at 450nm 
mean_QY_AgTFA_405 = mean([QY_AgTFA_405nm_s1,QY_AgTFA_405nm_s2,QY_AgTFA_405nm_s3])
stdev_QY_AgTFA_405 = std([QY_AgTFA_405nm_s1,QY_AgTFA_405nm_s2,QY_AgTFA_405nm_s3])

function AgTFA_QY_Action_Spectrum(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax =    CairoMakie.Axis(
            f[1,1], palette = (color = palette(ColorSchemes.inferno, 3),), ## or replace ncol with Integer 
            title = "Ag(II)TFA Action Spectrum",
            xlabel = "Wavelength (nm)", 
            ylabel = "QY", 
        )

        width = 3
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(nothing,600)
        # CairoMakie.ylims!(-0.001,nothing)
        
        
        CairoMakie.errorbars!([405,450,525],[mean_QY_AgTFA_405,mean_QY_AgTFA_450,mean_QY_AgTFA_525], [stdev_QY_AgTFA_405,stdev_QY_AgTFA_450,stdev_QY_AgTFA_525], color = :red, whiskerwidth = 10)
        CairoMakie.scatter!([405,450,525],[mean_QY_AgTFA_405,mean_QY_AgTFA_450,mean_QY_AgTFA_525], markersize = 10, color = :black)



        # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


        # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        # axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)
        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        resize_to_layout!(f)
   f
        
end
AgTFA_QY_Action_Spectrum()
save("JULIA/Outputs/AgTFA_QY_Action_Spectrum.pdf",AgTFA_QY_Action_Spectrum())

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

#### AgTFA_450nm_s1 #########################################################################################################

    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## AgTFA_450nm_s1 UvVis
    AgTFA_450nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s3-4sided/sample_Absorbance__0__001.csv",header = false, DataFrame))
 
    AgTFA_450nm_s1_UvVis[1,134]
    AgTFA_450nm_s1_UvVis[1,197]

    AgTFA_450nm_s1_UvVis_Kinetics = hcat(AgTFA_450nm_s1_UvVis[2:end,1], AgTFA_450nm_s1_UvVis[2:end,134])

    AgTFA_450nm_s1_Absorbances = transpose(AgTFA_450nm_s1_UvVis)
    AgTFA_450nm_s1_initial_Abs_450 = AgTFA_450nm_s1_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_450nm_s1_initial_Abs_pump = AgTFA_450nm_s1_Absorbances[197,2]

    Plots.plot(AgTFA_450nm_s1_Absorbances[2:end,1],AgTFA_450nm_s1_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_450nm_s1_UvVis[15:end,1],AgTFA_450nm_s1_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_450nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s3-4sided/poststage_Blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_450nm_s1_Blank[:,1],AgTFA_450nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_450nm_s1_Blank_avg = mean(AgTFA_450nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_450nm_s1_Fit = MonoFitIRF_New(AgTFA_450nm_s1_UvVis_Kinetics,2,5, 600, [0.08, 0.001, 49.45, 308.689, -0.007], 10,Inf,false,10000)
    k_AgTFA_450nm_s1 = 1/308.689

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_450nm_s1 = Flux(AgTFA_450nm_s1_Blank_avg,450)

    QY_AgTFA_450nm_s1= QY_method_3(k_AgTFA_450nm_s1, 1297, 0.003,FLux_Blank_AgTFA_450nm_s1,AgTFA_450nm_s1_initial_Abs_pump, AgTFA_450nm_s1_initial_Abs_450)



#### AgTFA_450nm_s2 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## AgTFA_450nm_s2 UvVis
    AgTFA_450nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s4/real_Absorbance__0__049.csv",header = false, DataFrame))
    
    AgTFA_450nm_s2_UvVis[1,134]
    AgTFA_450nm_s2_UvVis[1,197]

    AgTFA_450nm_s2_UvVis_Kinetics = hcat(AgTFA_450nm_s2_UvVis[2:end,1], AgTFA_450nm_s2_UvVis[2:end,134])

    AgTFA_450nm_s2_Absorbances = transpose(AgTFA_450nm_s2_UvVis)
    AgTFA_450nm_s2_initial_Abs_450 = AgTFA_450nm_s2_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_450nm_s2_initial_Abs_pump = AgTFA_450nm_s2_Absorbances[197,2]

    Plots.plot(AgTFA_450nm_s2_Absorbances[2:end,1],AgTFA_450nm_s2_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_450nm_s2_UvVis[2:end,1],AgTFA_450nm_s2_UvVis[2:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_450nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s4/afterstage-blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_450nm_s2_Blank[:,1],AgTFA_450nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_450nm_s2_Blank_avg = mean(AgTFA_450nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_450nm_s2_Fit = MonoFitIRF_New(AgTFA_450nm_s2_UvVis_Kinetics,2,5.0, 600, [0.08, 0.001, 41.288, 288.675, 0.026], 10,Inf,false,10000)
    k_AgTFA_450nm_s2 = 1/288.675

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_450nm_s2 = Flux(AgTFA_450nm_s2_Blank_avg,450)

    QY_AgTFA_450nm_s2= QY_method_3(k_AgTFA_450nm_s2, 1297, 0.003,FLux_Blank_AgTFA_450nm_s2,AgTFA_450nm_s2_initial_Abs_pump, AgTFA_450nm_s2_initial_Abs_450)



#### AgTFA_450nm_s3 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## AgTFA_450nm_s3 UvVis
    AgTFA_450nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s5/real_Absorbance__0__001.csv",header = false, DataFrame))

    AgTFA_450nm_s3_UvVis[1,134]
    AgTFA_450nm_s3_UvVis[1,197]

    AgTFA_450nm_s3_UvVis_Kinetics = hcat(AgTFA_450nm_s3_UvVis[2:end,1], AgTFA_450nm_s3_UvVis[2:end,134])

    AgTFA_450nm_s3_Absorbances = transpose(AgTFA_450nm_s3_UvVis)
    AgTFA_450nm_s3_initial_Abs_450 = AgTFA_450nm_s3_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_450nm_s3_initial_Abs_pump = AgTFA_450nm_s3_Absorbances[197,2]

    Plots.plot(AgTFA_450nm_s3_Absorbances[2:end,1],AgTFA_450nm_s3_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_450nm_s3_UvVis[2:end,1],AgTFA_450nm_s3_UvVis[2:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_450nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s5/afterstage-blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_450nm_s3_Blank[:,1],AgTFA_450nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_450nm_s3_Blank_avg = mean(AgTFA_450nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_450nm_s3_Fit = MonoFitIRF_New(AgTFA_450nm_s3_UvVis_Kinetics,2,10.0, 800, [0.08, 0.001, 50.593, 307.115, 0.0], 10,Inf,false,10000)
    k_AgTFA_450nm_s3 = 1/ 307.115

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_450nm_s3 = Flux(AgTFA_450nm_s3_Blank_avg,450)

    QY_AgTFA_450nm_s3= QY_method_3(k_AgTFA_450nm_s3, 1297, 0.003,FLux_Blank_AgTFA_450nm_s3,AgTFA_450nm_s3_initial_Abs_pump, AgTFA_450nm_s3_initial_Abs_450)



#### AgTFA_405nm_s1 #########################################################################################################

    # ϵ_405nm = 2203.7 for Ag(II)TFA
    ## AgTFA_405nm_s1 UvVis
    AgTFA_405nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/s1_real/Ag(II)TFA_Absorbances_405nm_s1.csv",header = false, DataFrame))
    AgTFA_405nm_s1_UvVis[1,108]
    AgTFA_405nm_s1_UvVis[1,139]


    AgTFA_405nm_s1_UvVis_Kinetics = hcat(AgTFA_405nm_s1_UvVis[15:end,1], AgTFA_405nm_s1_UvVis[15:end,108])

    AgTFA_405nm_s1_Absorbances = transpose(AgTFA_405nm_s1_UvVis)
    AgTFA_405nm_s1_initial_Abs_380 = AgTFA_405nm_s1_Absorbances[108,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_405nm_s1_initial_Abs_pump = AgTFA_405nm_s1_Absorbances[139,2]

    Plots.plot(AgTFA_405nm_s1_Absorbances[2:end,1],AgTFA_405nm_s1_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_405nm_s1_UvVis[15:end,1],AgTFA_405nm_s1_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/Ag(II)TFA_poststage_blank_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_405nm_s1_Blank[:,1],AgTFA_405nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s1_Blank_avg = mean(AgTFA_405nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s1_Fit = MonoFitIRF_New(AgTFA_405nm_s1_UvVis_Kinetics,2,10.0, 350, [0.08, 0.001, 3.156, 82.671, 0.045], 10,Inf,false,10000)
    k_AgTFA_405nm_s1 = 1/ 82.671

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s1 = Flux(AgTFA_405nm_s1_Blank_avg,405)

    QY_AgTFA_405nm_s1= QY_method_3(k_AgTFA_405nm_s1, 2203.7, 0.003,FLux_Blank_AgTFA_405nm_s1,AgTFA_405nm_s1_initial_Abs_pump, AgTFA_405nm_s1_initial_Abs_380)


 

#### AgTFA_405nm_s2 #########################################################################################################

    # ϵ_405nm = 2203.7 for Ag(II)TFA
    ## AgTFA_405nm_s2 UvVis
    AgTFA_405nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/s2/Ag(II)TFA_Absorbances_405nm_s2.csv",header = false, DataFrame))

    AgTFA_405nm_s2_UvVis[1,139]
    AgTFA_405nm_s2_UvVis[1,108]

    AgTFA_405nm_s2_UvVis_Kinetics = hcat(AgTFA_405nm_s2_UvVis[2:end,1], AgTFA_405nm_s2_UvVis[2:end,108])

    AgTFA_405nm_s2_Absorbances = transpose(AgTFA_405nm_s2_UvVis)
    AgTFA_405nm_s2_initial_Abs_380 = AgTFA_405nm_s2_Absorbances[108,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_405nm_s2_initial_Abs_pump = AgTFA_405nm_s2_Absorbances[139,2]

    Plots.plot(AgTFA_405nm_s2_Absorbances[2:end,1],AgTFA_405nm_s2_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_405nm_s2_UvVis[2:end,1],AgTFA_405nm_s2_UvVis[2:end,108])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/Ag(II)TFA_poststage_blank_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_405nm_s2_Blank[:,1],AgTFA_405nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s2_Blank_avg = mean(AgTFA_405nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s2_Fit = MonoFitIRF_New(AgTFA_405nm_s2_UvVis_Kinetics,2,5.0, 400, [0.08, 0.001, 10.223, 98.068, 0.023], 10,Inf,false,10000)
    k_AgTFA_405nm_s2 = 1/98.068

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s2 = Flux(AgTFA_405nm_s2_Blank_avg,405)

    QY_AgTFA_405nm_s2= QY_method_3(k_AgTFA_405nm_s2, 2203.7, 0.003,FLux_Blank_AgTFA_405nm_s2,AgTFA_405nm_s2_initial_Abs_pump, AgTFA_405nm_s2_initial_Abs_380)



#### AgTFA_405nm_s3 #########################################################################################################

    ## AgTFA_405nm_s3 UvVis
    AgTFA_405nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/s3/Ag(II)TFA_Absorbances_405nm_s3.csv",header = false, DataFrame))

    AgTFA_405nm_s3_UvVis[1,108]
    AgTFA_405nm_s3_UvVis[1,139]

    AgTFA_405nm_s3_UvVis_Kinetics = hcat(AgTFA_405nm_s3_UvVis[2:end,1], AgTFA_405nm_s3_UvVis[2:end,108])

    AgTFA_405nm_s3_Absorbances = transpose(AgTFA_405nm_s3_UvVis)
    AgTFA_405nm_s3_initial_Abs_380 = AgTFA_405nm_s3_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_405nm_s3_initial_Abs_pump = AgTFA_405nm_s3_Absorbances[139,2]

    Plots.plot(AgTFA_405nm_s3_Absorbances[2:end,1],AgTFA_405nm_s3_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_405nm_s3_UvVis[2:end,1],AgTFA_405nm_s3_UvVis[2:end,108])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/Ag(II)TFA_poststage_blank_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_405nm_s3_Blank[:,1],AgTFA_405nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s3_Blank_avg = mean(AgTFA_405nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s3_Fit = MonoFitIRF_New(AgTFA_405nm_s3_UvVis_Kinetics,2,0.0, 400, [0.08, 0.001, 13.436, 119.074, -0.008], 10,Inf,false,10000)
    k_AgTFA_405nm_s3 = 1/119.074

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s3 = Flux(AgTFA_405nm_s3_Blank_avg,405)

    QY_AgTFA_405nm_s3= QY_method_3(k_AgTFA_405nm_s3, 1297, 0.003,FLux_Blank_AgTFA_405nm_s3,AgTFA_405nm_s3_initial_Abs_pump, AgTFA_405nm_s3_initial_Abs_380)


#### AgTFA_525nm_s1 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## Tracked Kinetics at 400nm, calculated concentration at 450nm
    ## AgTFA_525nm_s1 UvVis
    AgTFA_525nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S1/Ag(II)TFA_Absorbances_525nm_s1.csv",header = false, DataFrame))
    AgTFA_525nm_s1_UvVis[1,293]
    AgTFA_525nm_s1_UvVis[1,197]
    AgTFA_525nm_s1_UvVis[1,134]

    AgTFA_525nm_s1_UvVis_Kinetics = hcat(AgTFA_525nm_s1_UvVis[2:end,1], AgTFA_525nm_s1_UvVis[2:end,134])

    AgTFA_525nm_s1_Absorbances = transpose(AgTFA_525nm_s1_UvVis)
    AgTFA_525nm_s1_initial_Abs_450 = AgTFA_525nm_s1_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_525nm_s1_initial_Abs_pump = AgTFA_525nm_s1_Absorbances[293,2]

    Plots.plot(AgTFA_525nm_s1_Absorbances[2:end,1],AgTFA_525nm_s1_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_525nm_s1_UvVis[2:end,1],AgTFA_525nm_s1_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_525nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S1/Ag(II)TFA_poststage_blank_525nm_s1.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_525nm_s1_Blank[:,1],AgTFA_525nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_525nm_s1_Blank_avg = mean(AgTFA_525nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_525nm_s1_Fit = MonoFitIRF_New(AgTFA_525nm_s1_UvVis_Kinetics,2,0.0, 1000, [0.08, 0.001, 256.294, 1394.101, 0.026], 10,Inf,false,10000)
    k_AgTFA_525nm_s1 = 1/1394.101

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_525nm_s1 = Flux(AgTFA_525nm_s1_Blank_avg,525)

    QY_AgTFA_525nm_s1= QY_method_3(k_AgTFA_525nm_s1, 1297, 0.003,FLux_Blank_AgTFA_525nm_s1,AgTFA_525nm_s1_initial_Abs_pump, AgTFA_525nm_s1_initial_Abs_450)


#### AgTFA_525nm_s2 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## AgTFA_525nm_s2 UvVis
    AgTFA_525nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S2/Ag(II)TFA_Absorbances_525nm_s2.csv",header = false, DataFrame))
    AgTFA_525nm_s2_UvVis[1,293]
    AgTFA_525nm_s2_UvVis[1,134]
    AgTFA_525nm_s2_UvVis[1,197]

    AgTFA_525nm_s2_UvVis_Kinetics = hcat(AgTFA_525nm_s2_UvVis[2:end,1], AgTFA_525nm_s2_UvVis[2:end,134])

    AgTFA_525nm_s2_Absorbances = transpose(AgTFA_525nm_s2_UvVis)
    AgTFA_525nm_s2_initial_Abs_450 = AgTFA_525nm_s2_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_525nm_s2_initial_Abs_pump = AgTFA_525nm_s2_Absorbances[293,2]

    # Plots.plot(AgTFA_525nm_s2_Absorbances[2:end,1],AgTFA_525nm_s2_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_525nm_s2_UvVis[2:end,1],AgTFA_525nm_s2_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_525nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S2/Ag(II)TFA_poststage_blank_525nm_s2.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_525nm_s2_Blank[:,1],AgTFA_525nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_525nm_s2_Blank_avg = mean(AgTFA_525nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_525nm_s2_Fit = MonoFitIRF_New(AgTFA_525nm_s2_UvVis_Kinetics,2,100.0, 3000, [0.08, 0.001, 378.858, 1539.248, 0.01], 10,Inf,false,10000)
    k_AgTFA_525nm_s2 = 1/1539.248

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_525nm_s2 = Flux(AgTFA_525nm_s2_Blank_avg,525)

    QY_AgTFA_525nm_s2= QY_method_3(k_AgTFA_525nm_s2, 1297, 0.003,FLux_Blank_AgTFA_525nm_s2,AgTFA_525nm_s2_initial_Abs_pump, AgTFA_525nm_s2_initial_Abs_450)



#### AgTFA_525nm_s3 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
    ## AgTFA_525nm_s3 UvVis
    AgTFA_525nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S3/Ag(II)TFA_Absorbances_525nm_s3.csv",header = false, DataFrame))
    AgTFA_525nm_s3_UvVis[1,293]
    AgTFA_525nm_s3_UvVis[1,134]
    AgTFA_525nm_s3_UvVis[1,197]

    AgTFA_525nm_s3_UvVis_Kinetics = hcat(AgTFA_525nm_s3_UvVis[2:end,1], AgTFA_525nm_s3_UvVis[2:end,134])

    AgTFA_525nm_s3_Absorbances = transpose(AgTFA_525nm_s3_UvVis)
    AgTFA_525nm_s3_initial_Abs_450 = AgTFA_525nm_s3_Absorbances[197,2] ## Absorbance value at the wavelength for which you have ϵ
    AgTFA_525nm_s3_initial_Abs_pump = AgTFA_525nm_s3_Absorbances[293,2]

    # Plots.plot(AgTFA_525nm_s3_Absorbances[2:end,1],AgTFA_525nm_s3_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_525nm_s3_UvVis[2:end,1],AgTFA_525nm_s3_UvVis[2:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_525nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S3/Ag(II)TFA_poststage_blank_525nm_s3.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_525nm_s3_Blank[:,1],AgTFA_525nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_525nm_s3_Blank_avg = mean(AgTFA_525nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_525nm_s3_Fit = MonoFitIRF_New(AgTFA_525nm_s3_UvVis_Kinetics,2,100.0, 3000, [0.08, 0.001, 749.428, 1854.956, -0.034], 10,Inf,false,10000)
    k_AgTFA_525nm_s3 = 1/1854.956

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_525nm_s3 = Flux(AgTFA_525nm_s3_Blank_avg,525)

    QY_AgTFA_525nm_s3= QY_method_3(k_AgTFA_525nm_s3, 1297, 0.003,FLux_Blank_AgTFA_525nm_s3,AgTFA_525nm_s3_initial_Abs_pump, AgTFA_525nm_s3_initial_Abs_450)

