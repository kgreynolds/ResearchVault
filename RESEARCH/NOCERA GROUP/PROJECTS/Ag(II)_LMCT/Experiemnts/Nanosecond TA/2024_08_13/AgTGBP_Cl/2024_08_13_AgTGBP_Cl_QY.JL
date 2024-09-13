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

function QY_method_3(k, ϵ, volume, init_flux, init_A_pump, init_A)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value of the sample
    I_abs = init_flux*(1-10^(-init_A_pump))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end


#### AgTGBP_Cl_525nm_s1 #########################################################################################################

    ## AgTGBP_Cl_525nm_s1 UvVis
    AgTGBP_Cl_525nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/s1/AgTGBP_Cl_525nm_s1_Absorbances.csv",header = false, DataFrame))
    AgTGBP_Cl_525nm_s1_UvVis[1,108]

    AgTGBP_Cl_525nm_s1_UvVis_Kinetics = hcat(AgTGBP_Cl_525nm_s1_UvVis[15:end,1], AgTGBP_Cl_525nm_s1_UvVis[15:end,108])

    AgTGBP_Cl_525nm_s1_Absorbances = transpose(AgTGBP_Cl_525nm_s1_UvVis)

 #### SteadyState UvVis #########################################################################################################
    AgTGBP_Cl_525nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/AgTGBP_Cl_Steadystate_UvVis_525nmex.csv", header = true, DataFrame)
    AgTGBP_Cl_525nm_SteadyState_Uvvis[151,1]
    AgTGBP_Cl_525nm_s1_initial_Abs_450 = AgTGBP_Cl_525nm_SteadyState_Uvvis[151,3] - AgTGBP_Cl_525nm_SteadyState_Uvvis[151,2]
    AgTGBP_Cl_525nm_s1_initial_Abs_pump = AgTGBP_Cl_525nm_SteadyState_Uvvis[226,3] - AgTGBP_Cl_525nm_SteadyState_Uvvis[226,2]

    Plots.plot(AgTGBP_Cl_525nm_SteadyState_Uvvis[:,1],AgTGBP_Cl_525nm_SteadyState_Uvvis[:,3] - AgTGBP_Cl_525nm_SteadyState_Uvvis[:,2])

    Plots.plot(AgTGBP_Cl_525nm_s1_UvVis[15:end,1],AgTGBP_Cl_525nm_s1_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_Cl_525nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/Poststage_blank_power_525nm.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_Cl_525nm_s1_Blank_avg = mean(AgTGBP_Cl_525nm_s1_Blank[:,2])
    AgTGBP_Cl_525nm_s1_Blank_stdev = std(AgTGBP_Cl_525nm_s1_Blank[:,2])


    ## Fitting for rate constant
    AgTGBP_Cl_525nm_s1_Fit = MonoFitIRF_New(AgTGBP_Cl_525nm_s1_UvVis_Kinetics,2,100.0, 1200,[0.001, 0.001, 370.631, 627.204, 0.015], 10,Inf,false,1200)
    k_AgTGBP_Cl_525nm_s1 = 1/627.204

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_Cl_525nm_s1 = Flux(AgTGBP_Cl_525nm_s1_Blank_avg,525)


    QY_AgTGBP_Cl_525nm_s1= QY_method_3(k_AgTGBP_Cl_525nm_s1,1305 , 0.003,FLux_Blank_AgTGBP_Cl_525nm_s1,AgTGBP_Cl_525nm_s1_initial_Abs_pump, AgTGBP_Cl_525nm_s1_initial_Abs_450)

#### AgTGBP_Cl_525nm_s2 #########################################################################################################

    ## AgTGBP_Cl_525nm_s2 UvVis
    AgTGBP_Cl_525nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/s2/AgTGBP_Cl_525nm_s2_Absorbances.csv",header = false, DataFrame))
    AgTGBP_Cl_525nm_s2_UvVis[1,108]

    AgTGBP_Cl_525nm_s2_UvVis_Kinetics = hcat(AgTGBP_Cl_525nm_s2_UvVis[15:end,1], AgTGBP_Cl_525nm_s2_UvVis[15:end,108])

    AgTGBP_Cl_525nm_s2_Absorbances = transpose(AgTGBP_Cl_525nm_s2_UvVis)

    #### SteadyState UvVis #########################################################################################################
    AgTGBP_Cl_525nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/AgTGBP_Cl_Steadystate_UvVis_525nmex.csv", header = true, DataFrame)
    AgTGBP_Cl_525nm_s2_initial_Abs_450 = AgTGBP_Cl_525nm_SteadyState_Uvvis[151,4] - AgTGBP_Cl_525nm_SteadyState_Uvvis[151,2]
    AgTGBP_Cl_525nm_s2_initial_Abs_pump = AgTGBP_Cl_525nm_SteadyState_Uvvis[226,4] - AgTGBP_Cl_525nm_SteadyState_Uvvis[226,2]

    Plots.plot(AgTGBP_Cl_525nm_SteadyState_Uvvis[:,1],AgTGBP_Cl_525nm_SteadyState_Uvvis[:,4] - AgTGBP_Cl_525nm_SteadyState_Uvvis[:,2])


    Plots.plot(AgTGBP_Cl_525nm_s2_UvVis[15:end,1],AgTGBP_Cl_525nm_s2_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_Cl_525nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/Poststage_blank_power_525nm.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_Cl_525nm_s2_Blank_avg = mean(AgTGBP_Cl_525nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_Cl_525nm_s2_Fit = MonoFitIRF_New(AgTGBP_Cl_525nm_s2_UvVis_Kinetics,2,100.0, 1200, [0.01, 0.001, 391.206, 632.802, -0.005], 10,Inf,false,1200)
    k_AgTGBP_Cl_525nm_s2 = 1/632.802

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_Cl_525nm_s2 = Flux(AgTGBP_Cl_525nm_s2_Blank_avg,525)

    QY_AgTGBP_Cl_525nm_s2= QY_method_3(k_AgTGBP_Cl_525nm_s2,1305 , 0.003,  FLux_Blank_AgTGBP_Cl_525nm_s2,  AgTGBP_Cl_525nm_s2_initial_Abs_pump,  AgTGBP_Cl_525nm_s2_initial_Abs_450)

#### AgTGBP_Cl_525nm_s3 #########################################################################################################

    ## AgTGBP_Cl_525nm_s3 UvVis
    AgTGBP_Cl_525nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/s3/AgTGBP_Cl_525nm_s3_Absorbances.csv",header = false, DataFrame))
    AgTGBP_Cl_525nm_s3_UvVis[1,108]

    AgTGBP_Cl_525nm_s3_UvVis_Kinetics = hcat(AgTGBP_Cl_525nm_s3_UvVis[15:end,1], AgTGBP_Cl_525nm_s3_UvVis[15:end,108])

    AgTGBP_Cl_525nm_s3_Absorbances = transpose(AgTGBP_Cl_525nm_s3_UvVis)

    #### SteadyState UvVis #########################################################################################################
    AgTGBP_Cl_525nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/AgTGBP_Cl_Steadystate_UvVis_525nmex.csv", header = true, DataFrame)
    AgTGBP_Cl_525nm_s3_initial_Abs_450 = AgTGBP_Cl_525nm_SteadyState_Uvvis[151,5] - AgTGBP_Cl_525nm_SteadyState_Uvvis[151,2]
    AgTGBP_Cl_525nm_s3_initial_Abs_pump = AgTGBP_Cl_525nm_SteadyState_Uvvis[226,5] - AgTGBP_Cl_525nm_SteadyState_Uvvis[226,2]

    Plots.plot(AgTGBP_Cl_525nm_SteadyState_Uvvis[:,1],AgTGBP_Cl_525nm_SteadyState_Uvvis[:,5] - AgTGBP_Cl_525nm_SteadyState_Uvvis[:,2])

    Plots.plot(AgTGBP_Cl_525nm_s3_UvVis[15:end,1],AgTGBP_Cl_525nm_s3_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_Cl_525nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/525nm/Poststage_blank_power_525nm.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_Cl_525nm_s3_Blank_avg = mean(AgTGBP_Cl_525nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_Cl_525nm_s3_Fit = MonoFitIRF_New(AgTGBP_Cl_525nm_s3_UvVis_Kinetics,2,100.0, 1200,[0.01, 0.001, 367.562, 622.949, -0.006], 1,Inf,false,1200)
    k_AgTGBP_Cl_525nm_s3 = 1/622.949

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_Cl_525nm_s3 = Flux(AgTGBP_Cl_525nm_s3_Blank_avg,525)

    
    QY_AgTGBP_Cl_525nm_s3= QY_method_3(k_AgTGBP_Cl_525nm_s3,1305 , 0.003,FLux_Blank_AgTGBP_Cl_525nm_s3,AgTGBP_Cl_525nm_s3_initial_Abs_pump, AgTGBP_Cl_525nm_s3_initial_Abs_450)

## mean QY for 525nm excitation
QY_mean_525 = mean([QY_AgTGBP_Cl_525nm_s3,QY_AgTGBP_Cl_525nm_s2,QY_AgTGBP_Cl_525nm_s3])
QY_stdev_525 = std([QY_AgTGBP_Cl_525nm_s3,QY_AgTGBP_Cl_525nm_s2,QY_AgTGBP_Cl_525nm_s3])











#### AgTGBP_Cl_10C_450nm_s1 #########################################################################################################

    ## AgTGBP_Cl_10C_450nm_s1 UvVis
    AgTGBP_Cl_10C_450nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/10C/AgTGBP_Cl_45nm_10C_Absorbances.csv",header = false, DataFrame))
    AgTGBP_Cl_10C_450nm_s1_UvVis[1,108]

    AgTGBP_Cl_10C_450nm_s1_UvVis_Kinetics = hcat(AgTGBP_Cl_10C_450nm_s1_UvVis[15:end,1], AgTGBP_Cl_10C_450nm_s1_UvVis[15:end,108])

    AgTGBP_Cl_10C_450nm_s1_Absorbances = transpose(AgTGBP_Cl_10C_450nm_s1_UvVis)

    #### SteadyState UvVis #########################################################################################################
    AgTGBP_Cl_10C_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/AgTGBP_Cl_450nm_Steadystate_UvVis.csv", header = true, DataFrame)
    AgTGBP_Cl_10C_450nm_s1_initial_Abs_450 = AgTGBP_Cl_10C_450nm_SteadyState_Uvvis[151,3] - AgTGBP_Cl_10C_450nm_SteadyState_Uvvis[151,2]
    AgTGBP_Cl_10C_450nm_s1_initial_Abs_pump = AgTGBP_Cl_10C_450nm_SteadyState_Uvvis[151,3] - AgTGBP_Cl_10C_450nm_SteadyState_Uvvis[151,2]


    Plots.plot(AgTGBP_Cl_10C_450nm_s1_UvVis[15:end,1],AgTGBP_Cl_10C_450nm_s1_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_Cl_10C_450nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/Poststage_blank_450nm_power.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_Cl_10C_450nm_s1_Blank_avg = mean(AgTGBP_Cl_10C_450nm_s1_Blank[:,2])
    AgTGBP_Cl_10C_450nm_s1_Blank_stdev = std(AgTGBP_Cl_10C_450nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_Cl_10C_450nm_s1_Fit = MonoFitIRF_New(AgTGBP_Cl_10C_450nm_s1_UvVis_Kinetics,2,0.0, 620,[0.001, 0.001, 45.288, 155.206, 0.018], 10,Inf,false,600)
    k_AgTGBP_Cl_10C_450nm_s1 = 1/155.206

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_Cl_10C_450nm_s1 = Flux(AgTGBP_Cl_10C_450nm_s1_Blank_avg,450)

    
    QY_AgTGBP_Cl_10C_450nm_s1= QY_method_3(k_AgTGBP_Cl_10C_450nm_s1, 1305, 0.003,FLux_Blank_AgTGBP_Cl_10C_450nm_s1,AgTGBP_Cl_10C_450nm_s1_initial_Abs_pump, AgTGBP_Cl_10C_450nm_s1_initial_Abs_450)




#### AgTGBP_Cl_15C_450nm_s2 #########################################################################################################

    ## AgTGBP_Cl_15C_450nm_s1 UvVis
    AgTGBP_Cl_15C_450nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/15C/AgTGBP_Cl_45nm_15C_Absorbances.csv",header = false, DataFrame))
    AgTGBP_Cl_15C_450nm_s1_UvVis[1,108]

    AgTGBP_Cl_15C_450nm_s1_UvVis_Kinetics = hcat(AgTGBP_Cl_15C_450nm_s1_UvVis[15:end,1], AgTGBP_Cl_15C_450nm_s1_UvVis[15:end,108])

    AgTGBP_Cl_15C_450nm_s1_Absorbances = transpose(AgTGBP_Cl_15C_450nm_s1_UvVis)

     #### SteadyState UvVis #########################################################################################################
    AgTGBP_Cl_15C_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/AgTGBP_Cl_450nm_Steadystate_UvVis.csv", header = true, DataFrame)
    AgTGBP_Cl_15C_450nm_s1_initial_Abs_450 = AgTGBP_Cl_15C_450nm_SteadyState_Uvvis[151,4] - AgTGBP_Cl_15C_450nm_SteadyState_Uvvis[151,2]
    AgTGBP_Cl_15C_450nm_s1_initial_Abs_pump = AgTGBP_Cl_15C_450nm_SteadyState_Uvvis[151,4] - AgTGBP_Cl_15C_450nm_SteadyState_Uvvis[151,2]


    Plots.plot(AgTGBP_Cl_15C_450nm_s1_UvVis[15:end,1],AgTGBP_Cl_15C_450nm_s1_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_Cl_15C_450nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/Poststage_blank_450nm_power.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_Cl_15C_450nm_s1_Blank_avg = mean(AgTGBP_Cl_15C_450nm_s1_Blank[:,2])
    AgTGBP_Cl_15C_450nm_s1_Blank_stdev = std(AgTGBP_Cl_15C_450nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_Cl_15C_450nm_s1_Fit = MonoFitIRF_New(AgTGBP_Cl_15C_450nm_s1_UvVis_Kinetics,2,0.0, 620,[0.001, 0.001, 29.647, 136.117, 0.033], 10,Inf,false,600)
    k_AgTGBP_Cl_15C_450nm_s1 = 1/136.117

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_Cl_15C_450nm_s1 = Flux(AgTGBP_Cl_15C_450nm_s1_Blank_avg,450)

    
    QY_AgTGBP_Cl_15C_450nm_s1= QY_method_3(k_AgTGBP_Cl_15C_450nm_s1, 1305, 0.003,FLux_Blank_AgTGBP_Cl_15C_450nm_s1,AgTGBP_Cl_15C_450nm_s1_initial_Abs_pump, AgTGBP_Cl_15C_450nm_s1_initial_Abs_450)

#### AgTGBP_Cl_5C_450nm_s3 #########################################################################################################

    ## AgTGBP_Cl_5C_450nm_s3 UvVis
    AgTGBP_Cl_5C_450nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/5C/AgTGBP_Cl_45nm_5C_Absorbances.csv",header = false, DataFrame))
    AgTGBP_Cl_5C_450nm_s3_UvVis[1,108]

    AgTGBP_Cl_5C_450nm_s3_UvVis_Kinetics = hcat(AgTGBP_Cl_5C_450nm_s3_UvVis[15:end,1], AgTGBP_Cl_5C_450nm_s3_UvVis[15:end,108])

    AgTGBP_Cl_5C_450nm_s3_Absorbances = transpose(AgTGBP_Cl_5C_450nm_s3_UvVis)


    #### SteadyState UvVis #########################################################################################################
    AgTGBP_Cl_5C_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/AgTGBP_Cl_450nm_Steadystate_UvVis.csv", header = true, DataFrame)
    AgTGBP_Cl_5C_450nm_s3_initial_Abs_450 = AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[151,5] - AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[151,2]
    AgTGBP_Cl_5C_450nm_s3_initial_Abs_pump = AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[151,5] - AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[151,2]

    Plots.plot(AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[:,1],AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[:,5] - AgTGBP_Cl_5C_450nm_SteadyState_Uvvis[:,2])

    Plots.plot(AgTGBP_Cl_5C_450nm_s3_UvVis[15:end,1],AgTGBP_Cl_5C_450nm_s3_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_Cl_5C_450nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_13/AgTGBP_Cl/450nm/Poststage_blank_450nm_power.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_Cl_5C_450nm_s3_Blank_avg = mean(AgTGBP_Cl_5C_450nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_Cl_5C_450nm_s3_Fit = MonoFitIRF_New(AgTGBP_Cl_5C_450nm_s3_UvVis_Kinetics,2,00.0, 620,[0.001, 0.001, 37.456, 149.025, 0.008], 1,Inf,false,600)
    k_AgTGBP_Cl_5C_450nm_s3 = 1/149.025

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_Cl_5C_450nm_s3 = Flux(AgTGBP_Cl_5C_450nm_s3_Blank_avg,450) 

    
    QY_AgTGBP_Cl_5C_450nm_s3= QY_method_3(k_AgTGBP_Cl_5C_450nm_s3, 1305, 0.003,FLux_Blank_AgTGBP_Cl_5C_450nm_s3,AgTGBP_Cl_5C_450nm_s3_initial_Abs_pump, AgTGBP_Cl_5C_450nm_s3_initial_Abs_450)

    function AgTGBP_Cl_QY_Temp(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "Ag(II) Tgbp Cl Temperature Depedance of QY",
                palette = (color = palette(ColorSchemes.inferno, 7),), xautolimitmargin = (0.1, 0.1),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Temperature (˚C)", 
                ylabel = "QY"
            )
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            # lines!(AgTFA_1p5mJ_450nm_Absorbances[2:end,1],AgTFA_1p5mJ_450nm_Absorbances[2:end,2],linewidth = width, label = "1 s")
            # CairoMakie.errorbars!([405,450,525],[4.375,2.267,0.659 ], [0.104,0.087,0.011], color = :red, whiskerwidth = 10)
            CairoMakie.scatter!([5,10,15,23],[3.74,3.69,4.01,4.54], markersize = 15, color = :black)



            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


            # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            # axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    AgTGBP_Cl_QY_Temp()
    save("JULIA/Outputs/AgTGBP_Cl_QY_Temp.pdf",AgTGBP_Cl_QY_Temp())
