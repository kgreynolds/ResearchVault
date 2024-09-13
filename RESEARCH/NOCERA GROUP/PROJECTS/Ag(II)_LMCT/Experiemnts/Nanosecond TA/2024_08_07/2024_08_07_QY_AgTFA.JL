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

#### AgTFA_1p5mJ_450nm #########################################################################################################
    # 40000 shots of 1.5 mJ of 450 nm at 50 Hz rep rate
    ## AgTFA_1p5mJ_450nm UvVis
    AgTFA_1p5mJ_450nm_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-017/Ag(II)TFA_Absorbances_1p5mJ_450nm.csv",header = false, DataFrame))
    AgTFA_1p5mJ_450nm_UvVis[1,134]
 
    AgTFA_1p5mJ_450nm_UvVis[1,197]
    AgTFA_1p5mJ_450nm_UvVis_Kinetics = hcat(AgTFA_1p5mJ_450nm_UvVis[15:end,1], AgTFA_1p5mJ_450nm_UvVis[15:end,134])

    AgTFA_1p5mJ_450nm_Absorbances = transpose(AgTFA_1p5mJ_450nm_UvVis)
    AgTFA_1p5mJ_450nm_initial_Abs_450 = AgTFA_1p5mJ_450nm_Absorbances[197,2]
    AgTFA_1p5mJ_450nm_initial_Abs_pump = AgTFA_1p5mJ_450nm_Absorbances[197,2]

    Plots.plot(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_1p5mJ_450nm_UvVis[15:end,1],AgTFA_1p5mJ_450nm_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_1p5mJ_450nm_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-017/450nm_Postblank_power.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_1p5mJ_450nm_Blank[:,1],AgTFA_1p5mJ_450nm_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_1p5mJ_450nm_Blank_avg = mean(AgTFA_1p5mJ_450nm_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_1p5mJ_450nm_Fit = MonoFitIRF_New(AgTFA_1p5mJ_450nm_UvVis_Kinetics,2,0.0, 800, [0.08, 0.001, 10.57, 155.194, 0.035], 10,Inf,false,807)
    k_AgTFA_1p5mJ_450nm = 1/155.194

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_1p5mJ_450nm = Flux(AgTFA_1p5mJ_450nm_Blank_avg,450)

    QY_AgTFA_1p5mJ_450nm= QY_method_3(k_AgTFA_1p5mJ_450nm, 1297, 0.003,FLux_Blank_AgTFA_1p5mJ_450nm,AgTFA_1p5mJ_450nm_initial_Abs_pump,AgTFA_1p5mJ_450nm_initial_Abs_450)


    function AgTFA_QY_Action(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "Ag(II)TFA QY Action Spectrum",
                palette = (color = palette(ColorSchemes.inferno, 7),), xautolimitmargin = (0.1, 0.1),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "QY"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            # lines!(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2],linewidth = width, label = "1 s")
            CairoMakie.errorbars!([405,450,525],[4.375,2.267,0.659 ], [0.104,0.087,0.011], color = :red, whiskerwidth = 10)
            CairoMakie.scatter!([405,450,525],[4.375,2.267,0.659], markersize = 10, color = :black)



            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


            # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            # axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    AgTFA_QY_Action()
    save("JULIA/Outputs/AgTFA_QY_Action.pdf",AgTFA_QY_Action())

#### AgTFA_525nm_s1 #########################################################################################################
    # 60000 shots of 1 mJ of 525 nm at 50 Hz 
    ## AgTFA_525nm_s1 UvVis
    AgTFA_525nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S1/Ag(II)TFA_Absorbances_525nm_s1.csv",header = false, DataFrame))
    AgTFA_525nm_s1_UvVis[1,134]

    AgTFA_525nm_s1_UvVis[1,197]
    AgTFA_525nm_s1_UvVis[1,293]

    AgTFA_525nm_s1_UvVis_Kinetics = hcat(AgTFA_525nm_s1_UvVis[15:end,1], AgTFA_525nm_s1_UvVis[15:end,134])

    AgTFA_525nm_s1_Absorbances = transpose(AgTFA_525nm_s1_UvVis)
    AgTFA_525nm_s1_initial_Abs_450nm = AgTFA_525nm_s1_Absorbances[197,2]
    AgTFA_525nm_s1_initial_Abs_pump = AgTFA_525nm_s1_Absorbances[293,2]

    Plots.plot(AgTFA_525nm_s1_Absorbances[2:end,1],AgTFA_525nm_s1_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_525nm_s1_UvVis[15:end,1],AgTFA_525nm_s1_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_525nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S1/Ag(II)TFA_poststage_blank_525nm_s1.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_525nm_s1_Blank[:,1],AgTFA_525nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_525nm_s1_Blank_avg = mean(AgTFA_525nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_525nm_s1_Fit = MonoFitIRF_New(AgTFA_525nm_s1_UvVis_Kinetics,2,20.0, 1000, [20.012, 0.001, 255.941, 1405.288, 0.025], 10,Inf,false,1008)
    k_AgTFA_525nm_s1 = 1/1405.288

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_525nm_s1 = Flux(AgTFA_525nm_s1_Blank_avg,525)

    QY_AgTFA_525nm_s1= QY_method_3(k_AgTFA_525nm_s1, 1297, 0.003,FLux_Blank_AgTFA_525nm_s1,AgTFA_525nm_s1_initial_Abs_pump,AgTFA_525nm_s1_initial_Abs_450nm)


#### AgTFA_525nm_s2 #########################################################################################################
    #180000 shots of 1 mJ of 525 nm at 50 Hz
    ## AgTFA_525nm_s2 UvVis
    AgTFA_525nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S2/Ag(II)TFA_Absorbances_525nm_s2.csv",header = false, DataFrame))
    AgTFA_525nm_s2_UvVis[1,134]
    AgTFA_525nm_s2_UvVis[1,293]
    AgTFA_525nm_s2_UvVis[1,197]

    AgTFA_525nm_s2_UvVis_Kinetics = hcat(AgTFA_525nm_s2_UvVis[15:end,1], AgTFA_525nm_s2_UvVis[15:end,134])

    AgTFA_525nm_s2_Absorbances = transpose(AgTFA_525nm_s2_UvVis)
    AgTFA_525nm_s2_initial_Abs_450 = AgTFA_525nm_s2_Absorbances[197,2]
    AgTFA_525nm_s2_initial_Abs_pump = AgTFA_525nm_s2_Absorbances[293,2]

    # Plots.plot(AgTFA_525nm_s2_Absorbances[2:end,1],AgTFA_525nm_s2_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_525nm_s2_UvVis[15:end,1],AgTFA_525nm_s2_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_525nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S2/Ag(II)TFA_poststage_blank_525nm_s2.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_525nm_s2_Blank[:,1],AgTFA_525nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_525nm_s2_Blank_avg = mean(AgTFA_525nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_525nm_s2_Fit = MonoFitIRF_New(AgTFA_525nm_s2_UvVis_Kinetics,2,100.0, 3000, [100.011, 0.001, 355.073, 1539.35, 0.01], 10,Inf,false,3000)
    k_AgTFA_525nm_s2 = 1/1539.35

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_525nm_s2 = Flux(AgTFA_525nm_s2_Blank_avg,525)

    QY_AgTFA_525nm_s2= QY_method_3(k_AgTFA_525nm_s2, 1297, 0.003,FLux_Blank_AgTFA_525nm_s2,AgTFA_525nm_s2_initial_Abs_pump, AgTFA_525nm_s2_initial_Abs_450)

#### AgTFA_525nm_s3 #########################################################################################################

    ## AgTFA_525nm_s3 UvVis
    AgTFA_525nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S3/Ag(II)TFA_Absorbances_525nm_s3.csv",header = false, DataFrame))
    AgTFA_525nm_s3_UvVis[1,134]
    AgTFA_525nm_s3_UvVis[1,293]
    AgTFA_525nm_s3_UvVis[1,197]

    AgTFA_525nm_s3_UvVis_Kinetics = hcat(AgTFA_525nm_s3_UvVis[15:end,1], AgTFA_525nm_s3_UvVis[15:end,134])

    AgTFA_525nm_s3_Absorbances = transpose(AgTFA_525nm_s3_UvVis)
    AgTFA_525nm_s3_initial_Abs_450 = AgTFA_525nm_s3_Absorbances[197,2]
    AgTFA_525nm_s3_initial_Abs_pump = AgTFA_525nm_s3_Absorbances[293,2]

    # Plots.plot(AgTFA_525nm_s3_Absorbances[2:end,1],AgTFA_525nm_s3_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_525nm_s3_UvVis[15:end,1],AgTFA_525nm_s3_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    AgTFA_525nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-18_525nm/S3/Ag(II)TFA_poststage_blank_525nm_s3.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_525nm_s3_Blank[:,1],AgTFA_525nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_525nm_s3_Blank_avg = mean(AgTFA_525nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_525nm_s3_Fit = MonoFitIRF_New(AgTFA_525nm_s3_UvVis_Kinetics,2,120.0, 3000, [120.0, 0.001, 699.981, 1849.473, -0.033], 10,Inf,false,3000)
    k_AgTFA_525nm_s3 = 1/1849.473

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_525nm_s3 = Flux(AgTFA_525nm_s3_Blank_avg,525)

    QY_AgTFA_525nm_s3 = QY_method_3(k_AgTFA_525nm_s3, 1297, 0.003,FLux_Blank_AgTFA_525nm_s3,AgTFA_525nm_s3_initial_Abs_pump, AgTFA_525nm_s3_initial_Abs_450)

    ## mean QY for 525nm excitation
    QY_mean_525 = mean([QY_AgTFA_525nm_s1,QY_AgTFA_525nm_s2,QY_AgTFA_525nm_s3])
    QY_stdev_525 = std([QY_AgTFA_525nm_s1,QY_AgTFA_525nm_s2,QY_AgTFA_525nm_s3])

#### AgTFA_405nm_s1 #########################################################################################################

    ## AgTFA_405nm_s1 UvVis
    AgTFA_405nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/s1_real/Ag(II)TFA_Absorbances_405nm_s1.csv",header = false, DataFrame))
    AgTFA_405nm_s1_UvVis[1,120]
    AgTFA_405nm_s1_UvVis[1,139]
    AgTFA_405nm_s1_UvVis[1,108]

    AgTFA_405nm_s1_UvVis_Kinetics = hcat(AgTFA_405nm_s1_UvVis[15:end,1], AgTFA_405nm_s1_UvVis[15:end,120])

    AgTFA_405nm_s1_Absorbances = transpose(AgTFA_405nm_s1_UvVis)
    AgTFA_405nm_s1_initial_Abs_380 = AgTFA_405nm_s1_Absorbances[108,2]
    AgTFA_405nm_s1_initial_Abs_pump = AgTFA_405nm_s1_Absorbances[139,2]

    Plots.plot(AgTFA_405nm_s1_Absorbances[2:end,1],AgTFA_405nm_s1_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_405nm_s1_UvVis[15:end,1],AgTFA_405nm_s1_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/Ag(II)TFA_poststage_blank_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_405nm_s1_Blank[:,1],AgTFA_405nm_s1_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s1_Blank_avg = mean(AgTFA_405nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s1_Fit = MonoFitIRF_New(AgTFA_405nm_s1_UvVis_Kinetics,2,0.0, 400,[0.08, 0.001, 3.069, 80.328, 0.045] , 10,Inf,false,10000)
    k_AgTFA_405nm_s1 = 1/80.328

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s1 = Flux(AgTFA_405nm_s1_Blank_avg,405)

    # calculate at 380nm with ϵ_380 = 2203.7
    QY_AgTFA_405nm_s1= QY_method_3(k_AgTFA_405nm_s1, 2203.7, 0.003,FLux_Blank_AgTFA_405nm_s1,AgTFA_405nm_s1_initial_Abs_pump, AgTFA_405nm_s1_initial_Abs_380)

#### AgTFA_405nm_s2 #########################################################################################################

    ## AgTFA_405nm_s2 UvVis
    AgTFA_405nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/s2/Ag(II)TFA_Absorbances_405nm_s2.csv",header = false, DataFrame))
    AgTFA_405nm_s2_UvVis[1,120]
    AgTFA_405nm_s2_UvVis[1,139]
    AgTFA_405nm_s2_UvVis[1,108]

    AgTFA_405nm_s2_UvVis_Kinetics = hcat(AgTFA_405nm_s2_UvVis[15:end,1], AgTFA_405nm_s2_UvVis[15:end,120])

    AgTFA_405nm_s2_Absorbances = transpose(AgTFA_405nm_s2_UvVis)
    AgTFA_405nm_s2_initial_Abs_380 = AgTFA_405nm_s2_Absorbances[108,2]
    AgTFA_405nm_s2_initial_Abs_pump = AgTFA_405nm_s2_Absorbances[139,2]

    Plots.plot(AgTFA_405nm_s2_Absorbances[2:end,1],AgTFA_405nm_s2_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_405nm_s2_UvVis[15:end,1],AgTFA_405nm_s2_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/Ag(II)TFA_poststage_blank_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_405nm_s2_Blank[:,1],AgTFA_405nm_s2_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s2_Blank_avg = mean(AgTFA_405nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s2_Fit = MonoFitIRF_New(AgTFA_405nm_s2_UvVis_Kinetics,2,20.0, 400, [20.0, 0.001, 7.458, 99.7, 0.021], 1,Inf,false,408)
    k_AgTFA_405nm_s2 = 1/99.7

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s2 = Flux(AgTFA_405nm_s2_Blank_avg,405)
 # calculate at 380nm with ϵ_380 = 2203.7
    QY_AgTFA_405nm_s2= QY_method_3(k_AgTFA_405nm_s2, 2203.7, 0.003,  FLux_Blank_AgTFA_405nm_s2,  AgTFA_405nm_s2_initial_Abs_pump,  AgTFA_405nm_s2_initial_Abs_380)

#### AgTFA_405nm_s3 #########################################################################################################

    ## AgTFA_405nm_s3 UvVis
    AgTFA_405nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/s3/Ag(II)TFA_Absorbances_405nm_s3.csv",header = false, DataFrame))
    AgTFA_405nm_s3_UvVis[1,120]
    AgTFA_405nm_s3_UvVis[1,139]
    AgTFA_405nm_s3_UvVis[1,108]

    AgTFA_405nm_s3_UvVis_Kinetics = hcat(AgTFA_405nm_s3_UvVis[15:end,1], AgTFA_405nm_s3_UvVis[15:end,120])

    AgTFA_405nm_s3_Absorbances = transpose(AgTFA_405nm_s3_UvVis)
    AgTFA_405nm_s3_initial_Abs_380 = AgTFA_405nm_s3_Absorbances[108,2]
    AgTFA_405nm_s3_initial_Abs_pump = AgTFA_405nm_s3_Absorbances[139,2]

    Plots.plot(AgTFA_405nm_s3_Absorbances[2:end,1],AgTFA_405nm_s3_Absorbances[2:end,2:end])
    Plots.plot(AgTFA_405nm_s3_UvVis[15:end,1],AgTFA_405nm_s3_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_07/BC-4-19_405nm/Ag(II)TFA_poststage_blank_405nm.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(AgTFA_405nm_s3_Blank[:,1],AgTFA_405nm_s3_Blank[:,2])

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s3_Blank_avg = mean(AgTFA_405nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s3_Fit = MonoFitIRF_New(AgTFA_405nm_s3_UvVis_Kinetics,2,20.0, 300,[19.972, 0.001, 8.331, 103.882, -0.005], 1,Inf,false,408)
    k_AgTFA_405nm_s3 = 1/103.882

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s3 = Flux(AgTFA_405nm_s3_Blank_avg,405)

    # calculate at 380nm with ϵ_380 = 2203.7
    QY_AgTFA_405nm_s3= QY_method_3(k_AgTFA_405nm_s3, 2203.7, 0.003,FLux_Blank_AgTFA_405nm_s3,AgTFA_405nm_s3_initial_Abs_pump, AgTFA_405nm_s3_initial_Abs_380)

    AgTFA_405nm_s2_Absorbances[800,3]


    ## mean QY for 405nm excitation
    QY_mean_405 = mean([QY_AgTFA_405nm_s3,QY_AgTFA_405nm_s2,QY_AgTFA_405nm_s3])
    QY_stdev_405 = std([QY_AgTFA_405nm_s3,QY_AgTFA_405nm_s2,QY_AgTFA_405nm_s3])

    ## Normalized absorbances baseline corrected at different wavelngths. at low concentrations (405 nm exitation) the error in inital concentration estimaition is large
    Plots.plot(AgTFA_405nm_s2_Absorbances[2:end,1],Norm(AgTFA_405nm_s2_Absorbances[2:end,1],AgTFA_405nm_s2_Absorbances[2:end,2].- AgTFA_405nm_s2_Absorbances[800,3], 450))
    Plots.plot!(AgTFA_525nm_s3_Absorbances[2:end,1],Norm(AgTFA_525nm_s3_Absorbances[2:end,1],AgTFA_525nm_s3_Absorbances[2:end,2] .- AgTFA_525nm_s3_Absorbances[800,3],450))
    Plots.plot!(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],Norm(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2] .-AgTFA_1p5mJ_450nm_Absorbances[800,3],450))

    Plots.plot(AgTFA_405nm_s2_Absorbances[2:end,1],Norm(AgTFA_405nm_s2_Absorbances[2:end,1],AgTFA_405nm_s2_Absorbances[2:end,2].- AgTFA_405nm_s2_Absorbances[800,3], 400))
    Plots.plot!(AgTFA_525nm_s3_Absorbances[2:end,1],Norm(AgTFA_525nm_s3_Absorbances[2:end,1],AgTFA_525nm_s3_Absorbances[2:end,2].- AgTFA_525nm_s3_Absorbances[800,3],400))
    Plots.plot!(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],Norm(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2].- AgTFA_1p5mJ_450nm_Absorbances[800,3],400))

    Plots.plot(AgTFA_405nm_s2_Absorbances[2:end,1],Norm(AgTFA_405nm_s2_Absorbances[2:end,1],AgTFA_405nm_s2_Absorbances[2:end,2].- AgTFA_405nm_s2_Absorbances[800,3], 380))
    Plots.plot!(AgTFA_525nm_s3_Absorbances[2:end,1],Norm(AgTFA_525nm_s3_Absorbances[2:end,1],AgTFA_525nm_s3_Absorbances[2:end,2].- AgTFA_525nm_s3_Absorbances[800,3],380))
    Plots.plot!(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],Norm(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2].- AgTFA_1p5mJ_450nm_Absorbances[800,3],380))

    Plots.plot(AgTFA_405nm_s2_Absorbances[2:end,1],Norm(AgTFA_405nm_s2_Absorbances[2:end,1],AgTFA_405nm_s2_Absorbances[2:end,2], 800))
    Plots.plot!(AgTFA_525nm_s3_Absorbances[2:end,1],Norm(AgTFA_525nm_s3_Absorbances[2:end,1],AgTFA_525nm_s3_Absorbances[2:end,2],800))
    Plots.plot!(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],Norm(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2],800))