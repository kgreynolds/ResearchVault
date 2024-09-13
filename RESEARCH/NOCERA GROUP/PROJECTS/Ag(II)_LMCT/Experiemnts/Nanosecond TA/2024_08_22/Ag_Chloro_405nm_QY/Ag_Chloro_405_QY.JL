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

## Mean QY for Ag_Chloro at 450nm 
mean_QY_Ag_Chloro_405 = mean([QY_Ag_Chloro_405nm_s1,QY_Ag_Chloro_405nm_s2,QY_Ag_Chloro_405nm_s3])
stdev_QY_Ag_Chloro_405 = std([QY_Ag_Chloro_405nm_s1,QY_Ag_Chloro_405nm_s2,QY_Ag_Chloro_405nm_s3])

#### Experimental Notes #########################################################################################################
"""
All samples prepared in glovebox
excitation energy around 1 mJ at 405 nm with a 50 Hz Ns pulsed laser (EKSPLA NT230-50)
ND filter started to get damaged so more light got transmitted that wasnt accounted for in samples 2 and 3 leading to artificially high QY
Initial concentrations were determined by steady state UVVis
"""

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
