include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Outside(),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)


# Constants
const H = 6.62607015e-34  # Define Planck's constant
const C = 2.99792458e+8   # Define speed of light
const NA = 6.0221408e23   # Define Avogadro's number

# Function to calculate mean flux of pulsed laser
function calculate_flux(vector, wavelength)
    E_photon = H * C / (wavelength * 1e-9)  # Calculate energy of a single photon
    avg_power = mean(vector) * 50  # Calculate average power (50 Hz rep rate)
    return avg_power / (E_photon * NA)  # Return calculated flux
end

# Function to calculate quantum yield (method 3)
function quantum_yield(k, ϵ, volume, init_flux, init_A_pump, init_A)
    # k = rate constant, ϵ = molar absorbtivity at wavelength of kinetic, 
    # volume = volume of sample, init_flux = init photon flux of pump, init_A_pump = absorbance at pump WL, init_A = Absorbance at kinetic WL
    I_abs = init_flux * (1 - 10^(-init_A_pump)) / volume  # Calculate absorbed photon flux
    return k * init_A / (I_abs * ϵ) * 100  # Return calculated quantum yield as a percentage
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


#### Ag_Chloro_450nm_s1 #########################################################################################################

    ## ϵ_450nm = 1305 for Ag(II)_Chloro
        function process_ag_Chloro_450nm_s1(ϵ, volume, kinetic_wavelength, pump_wavelength)
            # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
            
            # Read UV-Vis data from CSV file
            # Kinetics obtained from Ocean optic series of UvVi's obtained every second
            uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
            # Extract kinetics data from UV-Vis data
            kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, 134])
            
            # Read steady-state data from CSV file
            steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
            
            # Find row index for the desired Wavelength for kinetic
            target_kinetic_wavelength = kinetic_wavelength
            row_index_kinetic = argmin(abs.(steady_state_data[:, 1] .- target_kinetic_wavelength))
        
            # Find row index for the desired Wavelength for pump
            target_pump_wavelength = pump_wavelength
            row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
        
            # Calculate initial absorbance at kinetic trace
            initial_abs_kinetic = steady_state_data[row_index_kinetic, 3] - steady_state_data[row_index_kinetic, 2]
            # Calculate initial absorbance at pump wavelength
            initial_abs_pump = steady_state_data[row_index_pump, 3] - steady_state_data[row_index_pump, 2]
        
            # Read blank data from CSV file
            blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto=6, header=true, DataFrame)
            # Calculate average of blank data
            blank_avg = mean(blank_data[:, 2])
        
            # Fit kinetics data using MonoFitIRF_New function
            fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 600, [0.08, 0.001, 20.335, 129.17, 0.031], 10, Inf, false, 10000)
            # Calculate rate constant from fit result
            k = 1 / fit_result[2][4]
            # Calculate flux using blank average
            flux = calculate_flux(blank_avg, pump_wavelength)
            # Calculate quantum yield
            qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)
        
            return qy  # Return calculated quantum yield
        end
        
        # Run the analysis
        QY_chloro_450_s1 = process_ag_Chloro_450nm_s1(1305, 0.003,450, 450)  # Call the processing function and store the result
        

#### Ag_Chloro_450nm_s2 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA
        function process_ag_Chloro_450nm_s2(ϵ, volume, kinetic_wavelength, pump_wavelength)
            # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
            
            # Read UV-Vis data from CSV file
            # Kinetics obtained from Ocean optic series of UvVi's obtained every second
            uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s2/Absorbance_AgTGBP_Cl_450nm_s2.csv", header=false, DataFrame))
            # Extract kinetics data from UV-Vis data
            kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, 134])
            
            # Read steady-state data from CSV file
            steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
            
            # Find row index for the desired Wavelength for kinetic
            target_kinetic_wavelength = kinetic_wavelength
            row_index_kinetic = argmin(abs.(steady_state_data[:, 1] .- target_kinetic_wavelength))
        
            # Find row index for the desired Wavelength for pump
            target_pump_wavelength = pump_wavelength
            row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
        
            # Calculate initial absorbance at kinetic trace
            initial_abs_kinetic = steady_state_data[row_index_kinetic, 4] - steady_state_data[row_index_kinetic, 2]
            # Calculate initial absorbance at pump wavelength
            initial_abs_pump = steady_state_data[row_index_pump, 4] - steady_state_data[row_index_pump, 2]
        
            # Read blank data from CSV file
            blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/Spectrocell_poststage_blank_power_450nm.csv", skipto=6, header=true, DataFrame)
            # Calculate average of blank data
            blank_avg = mean(blank_data[:, 2])
        
            # Fit kinetics data using MonoFitIRF_New function
            fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 600,[0.08, 0.001, 27.396, 135.575, -0.002], 10, Inf, false, 10000)
            # Calculate rate constant from fit result
            k = 1 / fit_result[2][4]
            # Calculate flux using blank average
            flux = calculate_flux(blank_avg, pump_wavelength)
            # Calculate quantum yield
            qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)
        
            return qy  # Return calculated quantum yield
        end
        
        # Run the analysis
        QY_chloro_450_s2 = process_ag_Chloro_450nm_s2(1305, 0.003,450, 450)  # Call the processing function and store the result
        


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

