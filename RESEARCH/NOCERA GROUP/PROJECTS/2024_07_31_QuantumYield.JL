## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

function Flux(vector,Wavelength)

    h = 6.62607015e-34
    c = 2.99e+8
    NA = 6.023e23 
    
    E_photon = h*c/(Wavelength*10^(-9))

    Avg_energy = mean(vector) # J per pulse
    Avg_Power = Avg_energy*50 # power in Watts at 50 Hz rep rate
    Avg_Flux = Avg_Power/E_photon/NA # mol photons /s

end

function QY(ΔA, ϵ, volume, flux, time)

    Δn = ΔA/(ϵ)*volume
    Quantum_yield = Δn / (flux * time) * 100

    return Quantum_yield
end

function QY_method_2(k, ϵ, volume, init_flux, init_A)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value
    I_abs = init_flux*(1-10^(-init_A))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end


#### SAMPLE 1 #########################################################################################################
    ## Sample 1 UvVis
    Sample1_UvVis_pre = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/Absorbance__0__10-16-28-098.csv", header = true, DataFrame)
    Sample1_UvVis_post = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/Absorbance__1__10-38-26-189.csv", header = true, DataFrame)
    Plots.plot(Sample1_UvVis_pre[:,1],Sample1_UvVis_pre[:,2])
    Plots.plot!(Sample1_UvVis_post[:,1],Sample1_UvVis_post[:,2])

    ## Per-pulse energy after the Blank
    Post_stage_wBlank_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/Afterstage-blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Post_stage_wBlank_avg = mean(Post_stage_wBlank_2024_07_30[:,2])
    Plots.plot(Post_stage_wBlank_2024_07_30[:,1],Post_stage_wBlank_2024_07_30[:,2])

    ## Per-pulse energy after the Sample for 20000 shots
    Post_stage_Sample_20000_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/Afterstage-sample-20000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot!(Post_stage_Sample_20000_2024_07_30[:,1],Post_stage_Sample_20000_2024_07_30[:,2])

    ## calculated absorbance from the ΔP 
    Post_stage_Sample_20000_abs = hcat(Post_stage_Sample_20000_2024_07_30[:,1], - log.(Post_stage_Sample_20000_2024_07_30[:,2] ./Post_stage_wBlank_avg))
    mean(Post_stage_Sample_20000_abs[1:10,2])
    Plots.plot(Post_stage_Sample_20000_abs[:,1],Post_stage_Sample_20000_abs[:,2])


    Sample1_Fit = MonoFitIRF_New(Post_stage_Sample_20000_abs,2,0.0, 400, [0.008, 0.0001, 54.335, 242.318, 0.02], 10,50,false,10000)
    k1 = 1/242.318

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_1 = Flux(Post_stage_wBlank_avg,450)
    Flux_Sample_1 = Flux(mean(Post_stage_Sample_20000_2024_07_30[1:10,2]),450)
    ΔFlux_1 = FLux_Blank_1 - Flux_Sample_1

    # Quantum Yield based on the measured inital absorbance
    QY_linear_Sample_Flux_20000_2 = QY_method_2(k1, 1297, 0.003,FLux_Blank_1,0.087)


#### SAMPLE 2 #########################################################################################################

    ## Sample 2 UvVis
    Sample2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s2-4sided/Sample/Sample_real_Absorbance__0__106.csv",header = false, DataFrame))
    UvVis_kinetics_400nm = hcat(Sample2_UvVis[15:end,1], Sample2_UvVis[15:end,135])
    UvVis_kinetics_350nm = hcat(Sample2_UvVis[15:end,1], Sample2_UvVis[15:end,70])


    Sample2_UvVis[1,135]
    Sample2_UvVis[1,70]

    Sample2_Absorbances = transpose(Sample2_UvVis)

    Plots.plot(Sample2_Absorbances[2:end,1],Sample2_Absorbances[2:end,2:end])
    
    Plots.plot(Sample2_UvVis[15:end,1],Sample2_UvVis[15:end,135])

    Plots.plot(Sample2_UvVis[15:end,1],Sample2_UvVis[15:end,70])


    ## Per-pulse energy after the Blank
    Post_stage_wBlank1_2024_07_30_sample2 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s2-4sided/Blank1-1000shots.csv", skipto = 6, header = true, DataFrame)
    Post_stage_wBlank2_2024_07_30_sample2 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s2-4sided/Blank2-1000shots.csv", skipto = 6, header = true, DataFrame)
    Post_stage_wBlank3_2024_07_30_sample2 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s2-4sided/blank3-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_wBlank1_2024_07_30_sample2[:,1],Post_stage_wBlank1_2024_07_30_sample2[:,2])

    Post_stage_wBlank_Sample2_avg = mean(Post_stage_wBlank1_2024_07_30_sample2[:,2])
    Post_stage_wBlank_Sample3_avg = mean(Post_stage_wBlank2_2024_07_30_sample2[:,2])
    Post_stage_wBlank_Sample3_avg = mean(Post_stage_wBlank3_2024_07_30_sample2[:,2])


    ## Per-pulse energy after the Sample for 20000 shots, Threshold energy sensor value at 1 mJ some pulses went over causing spikes
    Post_stage_Sample_30000_sample2_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s2-4sided/Sample-30000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_30000_sample2_2024_07_30[:,1],Post_stage_Sample_30000_sample2_2024_07_30[:,2])

    ## calculated absorbance from the ΔP 
    Post_stage_Sample_30000_sample2_2024_07_30_abs = hcat(Post_stage_Sample_30000_sample2_2024_07_30[:,1], - log.(Post_stage_Sample_30000_sample2_2024_07_30[:,2] ./Post_stage_wBlank_Sample2_avg))
    mean(Post_stage_Sample_20000_abs[1:10,2])
    Plots.plot(Post_stage_Sample_30000_sample2_2024_07_30_abs[:,1],Post_stage_Sample_30000_sample2_2024_07_30_abs[:,2])



    ## Fitting for rate constant
    Sample2_Fit = MonoFitIRF_New(UvVis_kinetics_400nm,2,0.0, 500, [0.008, 0.0001, 35.431, 260.486, 0.012], 100,Inf,false,10000)
    k2 = 1/260.486

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_2_1 = Flux(Post_stage_wBlank_Sample2_avg,450)
    FLux_Blank_2_2 = Flux(Post_stage_wBlank_Sample3_avg,450)
    FLux_Blank_2_3 = Flux(Post_stage_wBlank_Sample3_avg,450)

    QY_Sample_2 = QY_method_2(k2, 1297, 0.003,FLux_Blank_2_1,0.116)



#### SAMPLE 3 #########################################################################################################

    ## Sample 3 UvVis
    Sample3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s3-4sided/sample_Absorbance__0__001.csv",header = false, DataFrame))
    Sample3_UvVis_kinetics_400nm = hcat(Sample3_UvVis[15:end,1], Sample3_UvVis[15:end,135])
    Sample3_UvVis_kinetics_350nm = hcat(Sample3_UvVis[15:end,1], Sample3_UvVis[15:end,75])

    Sample3_UvVis[1,135]

    Sample3_Absorbances = transpose(Sample3_UvVis)

    Plots.plot(Sample3_Absorbances[2:end,1],Sample3_Absorbances[2:end,2:end])
    Plots.plot!(Sample3_UvVis[15:end,1],Sample3_UvVis[15:end,135])

    ## Per-pulse energy after the Blank
    Post_stage_wBlank_sample3 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s3-4sided/poststage_Blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_wBlank_sample3[:,1],Post_stage_wBlank_sample3[:,2])

    Post_stage_wBlank_Sample3_avg = mean(Post_stage_wBlank_sample3[:,2])



    ## Per-pulse energy after the Sample for 20000 shots, Threshold energy sensor value at 1 mJ some pulses went over causing spikes
    Post_stage_Sample_30000_sample3 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s3-4sided/afterstage-sample-30000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_30000_sample3[:,1],Post_stage_Sample_30000_sample3[:,2])

    ## calculated absorbance from the ΔP 
    Post_stage_Sample_30000_sample3_abs = hcat(Post_stage_Sample_30000_sample3[:,1], - log.(Post_stage_Sample_30000_sample3[:,2] ./Post_stage_wBlank_Sample3_avg))
    mean(Post_stage_Sample_30000_sample3_abs[1:10,2])
    Plots.plot(Post_stage_Sample_30000_sample3_abs[:,1],Post_stage_Sample_30000_sample3_abs[:,2])



    ## Fitting for rate constant
    Sample3_Fit = MonoFitIRF_New(Sample3_UvVis_kinetics_400nm,2,0.0, 500, [0.008, 0.0001, 50.16, 314.061, -0.006], 10,Inf,false,10000)
    k3 = 1/314.061


    ## Calculation of the Flux and Quantum yields
    FLux_Blank_3 = Flux(Post_stage_wBlank_Sample3_avg,450)

    QY_Sample_3 = QY_method_2(k3, 1297, 0.003,FLux_Blank_3,0.116)

#### SAMPLE 4 #########################################################################################################

    ## Sample 4 UvVis
    Sample4_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s4/real_Absorbance__0__049.csv",header = false, DataFrame))
    Sample4_UvVis_kinetics_400nm = hcat(Sample4_UvVis[15:end,1], Sample4_UvVis[15:end,135])

    Sample4_UvVis[1,135]

    Sample4_Absorbances = transpose(Sample4_UvVis)

    Plots.plot(Sample4_Absorbances[2:end,1],Sample4_Absorbances[2:end,2:end])
    Plots.plot(Sample4_UvVis[15:end,1],Sample4_UvVis[15:end,135])

    ## Per-pulse energy after the Blank
    Post_stage_wBlank_sample4 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s4/afterstage-blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_wBlank_sample4[:,1],Post_stage_wBlank_sample4[:,2])

    Post_stage_wBlank_Sample4_avg = mean(Post_stage_wBlank_sample4[:,2])



    ## Per-pulse energy after the Sample for 40000 shots, Threshold energy sensor value at 1 mJ some pulses went over causing spikes
    Post_stage_Sample_40000_sample4 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s4/afterstage-sample-40000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_40000_sample4[:,1],Post_stage_Sample_30000_sample4[:,2])

    ## calculated absorbance from the ΔP 
    Post_stage_Sample_40000_sample4_abs = hcat(Post_stage_Sample_40000_sample4[:,1], - log.(Post_stage_Sample_40000_sample4[:,2] ./Post_stage_wBlank_Sample4_avg))
    mean(Post_stage_Sample_40000_sample4_abs[1:10,2])
    Plots.plot(Post_stage_Sample_40000_sample4_abs[:,1],Post_stage_Sample_40000_sample4_abs[:,2])



    ## Fitting for rate constant
    Sample4_Fit = MonoFitIRF_New(Sample4_UvVis_kinetics_400nm,2,0.0, 800, [0.008, 0.0, 42.184, 295.011, 0.026], 10,Inf,false,10000)
    k4 = 1/295.011


    ## Calculation of the Flux and Quantum yields
    FLux_Blank_4 = Flux(Post_stage_wBlank_Sample4_avg,450)

    QY_Sample_4 = QY_method_2(k4, 1297, 0.003,FLux_Blank_4,0.137)

#### SAMPLE 5 #########################################################################################################

    ## Sample 5 UvVis
    Sample5_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s5/real_Absorbance__0__001.csv",header = false, DataFrame))
    Sample5_UvVis_kinetics_400nm = hcat(Sample5_UvVis[15:end,1], Sample5_UvVis[15:end,135])

    Sample5_UvVis[1,135]

    Sample5_Absorbances = transpose(Sample5_UvVis)

    Plots.plot(Sample5_Absorbances[2:end,1],Sample5_Absorbances[2:end,2:end])
    Plots.plot(Sample5_UvVis[15:end,1],Sample5_UvVis[15:end,135])

    ## Per-pulse energy after the Blank
    Post_stage_wBlank_sample5 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s5/afterstage-blank-1000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_wBlank_sample5[:,1],Post_stage_wBlank_sample5[:,2])

    Post_stage_wBlank_Sample5_avg = mean(Post_stage_wBlank_sample5[:,2])



    ## Per-pulse energy after the Sample for 40000 shots, Threshold energy sensor value at 1 mJ some pulses went over causing spikes
    Post_stage_Sample_40000_sample5 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_23/QY/BC-4-008/s5/afterstage-sample-40000shots.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_40000_sample5[:,1],Post_stage_Sample_40000_sample5[:,2])

    ## calculated absorbance from the ΔP 
    Post_stage_Sample_40000_sample5_abs = hcat(Post_stage_Sample_40000_sample5[:,1], - log.(Post_stage_Sample_40000_sample5[:,2] ./Post_stage_wBlank_Sample5_avg))
    mean(Post_stage_Sample_40000_sample5_abs[1:10,2])
    Plots.plot(Post_stage_Sample_40000_sample5_abs[:,1],Post_stage_Sample_40000_sample5_abs[:,2])



    ## Fitting for rate constant
    Sample5_Fit = MonoFitIRF_New(Sample5_UvVis_kinetics_400nm,2,0.0, 800, [0.008, 0.0001, 51.63, 313.736, 0.0], 10,Inf,false,10000)
    k5 = 1/313.736

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_5 = Flux(Post_stage_wBlank_Sample5_avg,450)

    QY_Sample_5 = QY_method_2(k5, 1297, 0.003,FLux_Blank_5,0.128)


    1/k1
    1/k2
    1/k3
    1/k4
    1/k5

    QY_Sample_2
    QY_Sample_3
    QY_Sample_4
    QY_Sample_5

mean([QY_Sample_2,QY_Sample_3,QY_Sample_4,QY_Sample_5])
std([QY_Sample_2,QY_Sample_3,QY_Sample_4,QY_Sample_5])

mean([QY_Sample_3,QY_Sample_4,QY_Sample_5])
std([QY_Sample_3,QY_Sample_4,QY_Sample_5])
