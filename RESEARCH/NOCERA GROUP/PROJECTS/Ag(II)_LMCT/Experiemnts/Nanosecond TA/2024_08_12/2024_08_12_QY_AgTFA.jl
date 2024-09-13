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


#### AgTFA_405nm_s1 #########################################################################################################

    ## AgTFA_405nm_s1 UvVis
    AgTFA_405nm_s1_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/s1/Absorbances_AgTFA_405_s1.csv",header = false, DataFrame))
    AgTFA_405nm_s1_UvVis[1,120]

    AgTFA_405nm_s1_UvVis_Kinetics = hcat(AgTFA_405nm_s1_UvVis[15:end,1], AgTFA_405nm_s1_UvVis[15:end,120])

    AgTFA_405nm_s1_Absorbances = transpose(AgTFA_405nm_s1_UvVis)

    ## SteadyState UvVis #########################################################################################################
    AgTFA_SteadyState_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/AgTFA_SteadyState_UvVis.csv", header = true, DataFrame)
    AgTFA_SteadyState_UvVis[81,1]
    AgTFA_405nm_s1_initial_Abs_380 = AgTFA_SteadyState_UvVis[81,3] - AgTFA_SteadyState_UvVis[81,2] 
    AgTFA_405nm_s1_initial_Abs_pump = AgTFA_SteadyState_UvVis[106,3] - AgTFA_SteadyState_UvVis[106,2] 

    Plots.plot(AgTFA_405nm_s1_UvVis[15:end,1],AgTFA_405nm_s1_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s1_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/Poststage_Spectrocell_Blank.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s1_Blank_avg = mean(AgTFA_405nm_s1_Blank[:,2])
    AgTFA_405nm_s1_Blank_stdev = std(AgTFA_405nm_s1_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s1_Fit = MonoFitIRF_New(AgTFA_405nm_s1_UvVis_Kinetics,2,75, 620,[75.022, 0.001, 8.235, 118.295, 0.015], 10,Inf,false,600)
    k_AgTFA_405nm_s1 = 1/118.295

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s1 = Flux(AgTFA_405nm_s1_Blank_avg,405)

    # calculate at 380nm with ϵ_380 = 2203.7
    QY_AgTFA_405nm_s1= QY_method_3(k_AgTFA_405nm_s1, 2182, 0.003,FLux_Blank_AgTFA_405nm_s1,AgTFA_405nm_s1_initial_Abs_pump, AgTFA_405nm_s1_initial_Abs_380)

#### AgTFA_405nm_s2 #########################################################################################################

    ## AgTFA_405nm_s2 UvVis
    AgTFA_405nm_s2_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/s2/Absorbances_AgTFA_405_s2.csv",header = false, DataFrame))
    AgTFA_405nm_s2_UvVis[1,120]


    AgTFA_405nm_s2_UvVis_Kinetics = hcat(AgTFA_405nm_s2_UvVis[15:end,1], AgTFA_405nm_s2_UvVis[15:end,120])

    AgTFA_405nm_s2_Absorbances = transpose(AgTFA_405nm_s2_UvVis)

    ## SteadyState UvVis #########################################################################################################
    AgTFA_SteadyState_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/AgTFA_SteadyState_UvVis.csv", header = true, DataFrame)
    AgTFA_405nm_s2_initial_Abs_380 = AgTFA_SteadyState_UvVis[81,4] - AgTFA_SteadyState_UvVis[81,2]
    AgTFA_405nm_s2_initial_Abs_pump = AgTFA_SteadyState_UvVis[106,4] - AgTFA_SteadyState_UvVis[106,2] 

    Plots.plot(AgTFA_405nm_s2_UvVis[15:end,1],AgTFA_405nm_s2_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s2_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/Poststage_Spectrocell_Blank.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s2_Blank_avg = mean(AgTFA_405nm_s2_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s2_Fit = MonoFitIRF_New(AgTFA_405nm_s2_UvVis_Kinetics,2,00.0, 620, [0.01, 0.001, 11.136, 106.945, 0.02], 10,Inf,false,600)
    k_AgTFA_405nm_s2 = 1/106.945

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s2 = Flux(AgTFA_405nm_s2_Blank_avg,405)

    # calculate at 380nm with ϵ_380 = 2203.7
    QY_AgTFA_405nm_s2= QY_method_3(k_AgTFA_405nm_s2, 2182, 0.003,  FLux_Blank_AgTFA_405nm_s2,  AgTFA_405nm_s2_initial_Abs_pump,  AgTFA_405nm_s2_initial_Abs_380)

#### AgTFA_405nm_s3 #########################################################################################################

    ## AgTFA_405nm_s3 UvVis
    AgTFA_405nm_s3_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/s3/Absorbances_AgTFA_405_s3.csv",header = false, DataFrame))
    AgTFA_405nm_s3_UvVis[1,120]

    AgTFA_405nm_s3_UvVis_Kinetics = hcat(AgTFA_405nm_s3_UvVis[15:end,1], AgTFA_405nm_s3_UvVis[15:end,120])

    AgTFA_405nm_s3_Absorbances = transpose(AgTFA_405nm_s3_UvVis)

    ## SteadyState UvVis #########################################################################################################
    AgTFA_SteadyState_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/AgTFA_SteadyState_UvVis.csv", header = true, DataFrame)
    AgTFA_405nm_s3_initial_Abs_380 = AgTFA_SteadyState_UvVis[81,5] - AgTFA_SteadyState_UvVis[81,2]
    AgTFA_405nm_s3_initial_Abs_pump = AgTFA_SteadyState_UvVis[106,5] - AgTFA_SteadyState_UvVis[106,2] 

    Plots.plot(AgTFA_405nm_s3_UvVis[15:end,1],AgTFA_405nm_s3_UvVis[15:end,120])

    ## Per-pulse energy after the Blank
    AgTFA_405nm_s3_Blank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTFA/Poststage_Spectrocell_Blank.csv", skipto = 6, header = true, DataFrame)

    ## Average Per-pulse energy after the Blank
    AgTFA_405nm_s3_Blank_avg = mean(AgTFA_405nm_s3_Blank[:,2])

    ## Fitting for rate constant
    AgTFA_405nm_s3_Fit = MonoFitIRF_New(AgTFA_405nm_s3_UvVis_Kinetics,2,20.0, 620,[19.972, 0.001, 10.286, 107.629, 0.0], 1,Inf,false,600)
    k_AgTFA_405nm_s3 = 1/107.629

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_AgTFA_405nm_s3 = Flux(AgTFA_405nm_s3_Blank_avg,405)

    # calculate at 380nm with ϵ_380 = 2203.7
    QY_AgTFA_405nm_s3= QY_method_3(k_AgTFA_405nm_s3, 2182, 0.003,FLux_Blank_AgTFA_405nm_s3,AgTFA_405nm_s3_initial_Abs_pump, AgTFA_405nm_s3_initial_Abs_380)

## mean QY for 405nm excitation
QY_mean_405 = mean([QY_AgTFA_405nm_s3,QY_AgTFA_405nm_s2,QY_AgTFA_405nm_s3])
QY_stdev_405 = std([QY_AgTFA_405nm_s3,QY_AgTFA_405nm_s2,QY_AgTFA_405nm_s3])




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
        CairoMakie.errorbars!([405,450,525],[4.53,2.267,0.659 ], [0.04,0.087,0.011], color = :red, whiskerwidth = 10)
        CairoMakie.scatter!([405,450,525],[4.53,2.267,0.659], markersize = 10, color = :black)



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

with_theme(AgTFA_QY_Action, theme_latexfonts())