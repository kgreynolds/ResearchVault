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


#### AgTGBP_450nm_s1 #########################################################################################################

    ## AgTGBP_450nm_s1 UvVis
    AgTGBP_450nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv",header = false, DataFrame))
    AgTGBP_450nm_s1_UvVis[1,108]

    AgTGBP_450nm_s1_UvVis_Kinetics = hcat(AgTGBP_450nm_s1_UvVis[15:end,1], AgTGBP_450nm_s1_UvVis[15:end,108])

    AgTGBP_450nm_s1_Absorbances = transpose(AgTGBP_450nm_s1_UvVis)

 #### SteadyState UvVis #########################################################################################################
    AgTGBP_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header = true, DataFrame)
    AgTGBP_450nm_s1_initial_Abs_450 = AgTGBP_450nm_SteadyState_Uvvis[151,3] - AgTGBP_450nm_SteadyState_Uvvis[151,2]
    AgTGBP_450nm_s1_initial_Abs_pump = AgTGBP_450nm_SteadyState_Uvvis[151,3] - AgTGBP_450nm_SteadyState_Uvvis[151,2]


    Plots.plot(AgTGBP_450nm_s1_UvVis[15:end,1],AgTGBP_450nm_s1_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTGBP_450nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_450nm_s1_Blank_avg = mean(AgTGBP_450nm_s1_Blank[:,2])
    AgTGBP_450nm_s1_Blank_stdev = std(AgTGBP_450nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_450nm_s1_Fit = MonoFitIRF_New(AgTGBP_450nm_s1_UvVis_Kinetics,2,0.0, 620,[0.001, 0.001, 24.956, 127.473, 0.032], 10,Inf,false,600)
    k_AgTGBP_450nm_s1 = 1/127.473

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_450nm_s1 = Flux(AgTGBP_450nm_s1_Blank_avg,450)

    
    QY_AgTGBP_450nm_s1= QY_method_3(k_AgTGBP_450nm_s1, 1305, 0.003,FLux_Blank_AgTGBP_450nm_s1,AgTGBP_450nm_s1_initial_Abs_pump, AgTGBP_450nm_s1_initial_Abs_450)

#### AgTGBP_450nm_s2 #########################################################################################################

    ## AgTGBP_450nm_s2 UvVis
    AgTGBP_450nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s2/Absorbance_AgTGBP_Cl_450nm_s2.csv",header = false, DataFrame))
    AgTGBP_450nm_s2_UvVis[1,108]


    AgTGBP_450nm_s2_UvVis_Kinetics = hcat(AgTGBP_450nm_s2_UvVis[15:end,1], AgTGBP_450nm_s2_UvVis[15:end,108])

    AgTGBP_450nm_s2_Absorbances = transpose(AgTGBP_450nm_s2_UvVis)

    #### SteadyState UvVis #########################################################################################################
    AgTGBP_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header = true, DataFrame)
    AgTGBP_450nm_s2_initial_Abs_450 = AgTGBP_450nm_SteadyState_Uvvis[151,4] - AgTGBP_450nm_SteadyState_Uvvis[151,2]
    AgTGBP_450nm_s2_initial_Abs_pump = AgTGBP_450nm_SteadyState_Uvvis[151,4] - AgTGBP_450nm_SteadyState_Uvvis[151,2]

    Plots.plot(AgTGBP_450nm_SteadyState_Uvvis[:,1],AgTGBP_450nm_SteadyState_Uvvis[:,4] - AgTGBP_450nm_SteadyState_Uvvis[:,2])


    Plots.plot(AgTGBP_450nm_s2_UvVis[15:end,1],AgTGBP_450nm_s2_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTGBP_450nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_450nm_s2_Blank_avg = mean(AgTGBP_450nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_450nm_s2_Fit = MonoFitIRF_New(AgTGBP_450nm_s2_UvVis_Kinetics,2,00.0, 620, [0.01, 0.001, 33.573, 133.981, -0.001], 10,Inf,false,600)
    k_AgTGBP_450nm_s2 = 1/133.981

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_450nm_s2 = Flux(AgTGBP_450nm_s2_Blank_avg,450)
 
    QY_AgTGBP_450nm_s2= QY_method_3(k_AgTGBP_450nm_s2, 1305, 0.003,  FLux_Blank_AgTGBP_450nm_s2,  AgTGBP_450nm_s2_initial_Abs_pump,  AgTGBP_450nm_s2_initial_Abs_450)

#### AgTGBP_450nm_s3 #########################################################################################################

    ## AgTGBP_450nm_s3 UvVis
    AgTGBP_450nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s3/Absorbance_AgTGBP_Cl_450nm_s3.csv",header = false, DataFrame))
    AgTGBP_450nm_s3_UvVis[1,108]

    AgTGBP_450nm_s3_UvVis_Kinetics = hcat(AgTGBP_450nm_s3_UvVis[15:end,1], AgTGBP_450nm_s3_UvVis[15:end,108])

    AgTGBP_450nm_s3_Absorbances = transpose(AgTGBP_450nm_s3_UvVis)


    #### SteadyState UvVis #########################################################################################################
    AgTGBP_450nm_SteadyState_Uvvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header = true, DataFrame)
    AgTGBP_450nm_s3_initial_Abs_450 = AgTGBP_450nm_SteadyState_Uvvis[151,5] - AgTGBP_450nm_SteadyState_Uvvis[151,2]
    AgTGBP_450nm_s3_initial_Abs_pump = AgTGBP_450nm_SteadyState_Uvvis[151,5] - AgTGBP_450nm_SteadyState_Uvvis[151,2]

    Plots.plot(AgTGBP_450nm_SteadyState_Uvvis[:,1],AgTGBP_450nm_SteadyState_Uvvis[:,5] - AgTGBP_450nm_SteadyState_Uvvis[:,2])

    Plots.plot(AgTGBP_450nm_s3_UvVis[15:end,1],AgTGBP_450nm_s3_UvVis[15:end,108])

    ## Per-pulse energy after the Blank
    AgTGBP_450nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTGBP_450nm_s3_Blank_avg = mean(AgTGBP_450nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTGBP_450nm_s3_Fit = MonoFitIRF_New(AgTGBP_450nm_s3_UvVis_Kinetics,2,00.0, 620,[0.01, 0.001, 26.542, 129.256, 0.021], 1,Inf,false,600)
    k_AgTGBP_450nm_s3 = 1/129.256

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTGBP_450nm_s3 = Flux(AgTGBP_450nm_s3_Blank_avg,450)

    
    QY_AgTGBP_450nm_s3= QY_method_3(k_AgTGBP_450nm_s3, 1305, 0.003,FLux_Blank_AgTGBP_450nm_s3,AgTGBP_450nm_s3_initial_Abs_pump, AgTGBP_450nm_s3_initial_Abs_450)

## mean QY for 450nm excitation
QY_mean_450 = mean([QY_AgTGBP_450nm_s3,QY_AgTGBP_450nm_s2,QY_AgTGBP_450nm_s3])
QY_stdev_450 = std([QY_AgTGBP_450nm_s3,QY_AgTGBP_450nm_s2,QY_AgTGBP_450nm_s3])

