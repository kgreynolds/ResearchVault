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

# Function for data processing of Ag_Bromo_450nm_s1 experiment
function process_ag_bromo_450nm_s1(ϵ, volume, kinetic_wavelength, pump_wavelength)
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
QY = process_ag_bromo_450nm_s1(1305, 0.003,450, 450)  # Call the processing function and store the result


### Summary ###### RUN LINES BELOW First

### STILLE NEED 405nm, 450nm, and 525nm

## Mean QY for Ag_Bromo at 525nm 
mean_QY_Ag_Bromo_525 = mean([QY_Ag_Bromo_525nm_s1,QY_Ag_Bromo_525nm_s2,QY_Ag_Bromo_525nm_s3])
stdev_QY_Ag_Bromo_525 = std([QY_Ag_Bromo_525nm_s1,QY_Ag_Bromo_525nm_s2,QY_Ag_Bromo_525nm_s3])

## Mean QY for Ag_Bromo at 450nm 
mean_QY_Ag_Bromo_450 = mean([QY_Ag_Bromo_450nm_s1,QY_Ag_Bromo_450nm_s2,QY_Ag_Bromo_450nm_s3])
stdev_QY_Ag_Bromo_450 = std([QY_Ag_Bromo_450nm_s1,QY_Ag_Bromo_450nm_s2,QY_Ag_Bromo_450nm_s3])

## Mean QY for Ag_Bromo at 450nm 
mean_QY_Ag_Bromo_405 = mean([QY_Ag_Bromo_405nm_s1,QY_Ag_Bromo_405nm_s2,QY_Ag_Bromo_405nm_s3])
stdev_QY_Ag_Bromo_405 = std([QY_Ag_Bromo_405nm_s1,QY_Ag_Bromo_405nm_s2,QY_Ag_Bromo_405nm_s3])

function Ag_Bromo_QY_Action_Spectrum(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax =    CairoMakie.Axis(
            f[1,1], palette = (color = palette(ColorSchemes.inferno, 3),), ## or replace ncol with Integer 
            title = "Ag(II) Bromo Action Spectrum",
            xlabel = "Wavelength (nm)", 
            ylabel = "QY", 
        )

        width = 3
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(nothing,600)
        # CairoMakie.ylims!(-0.001,nothing)
        
        
        CairoMakie.errorbars!([405,450,525],[mean_QY_Ag_Bromo_405,mean_QY_Ag_Bromo_450,mean_QY_Ag_Bromo_525], [stdev_QY_Ag_Bromo_405,stdev_QY_Ag_Bromo_450,stdev_QY_Ag_Bromo_525], color = :red, whiskerwidth = 10)
        CairoMakie.scatter!([405,450,525],[mean_QY_Ag_Bromo_405,mean_QY_Ag_Bromo_450,mean_QY_Ag_Bromo_525], markersize = 10, color = :black)



        # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


        # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        # axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)
        colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        resize_to_layout!(f)
   f
        
end
Ag_Bromo_QY_Action_Spectrum()
save("JULIA/Outputs/Ag_Bromo_QY_Action_Spectrum.pdf",Ag_Bromo_QY_Action_Spectrum())



#### Ag_Bromo_450nm #########################################################################################################

    ## ϵ_450nm = 1305 for Ag(II)_Bromo
        function process_ag_bromo_450nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
            # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
        
            # Read UV-Vis data from CSV file
            # Kinetics obtained from Ocean optic series of UvVi's obtained every second
            uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
            # Extract kinetics data from UV-Vis data
            target_kinetic_wavelength = kinetic_wavelength
            row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
            kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
                
            # Read steady-state data from CSV file
            steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
            
            # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
            target_C0_wavelength = C0_wavelength
            row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
        
            # Find row index for the desired Wavelength for pump
            target_pump_wavelength = pump_wavelength
            row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
        
            # Calculate initial absorbance at kinetic trace
            initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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

        QY_ag_bromo_450nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 450, 400, 450)
        function process_ag_bromo_450nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
            # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
        
            # Read UV-Vis data from CSV file
            # Kinetics obtained from Ocean optic series of UvVi's obtained every second
            uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
            # Extract kinetics data from UV-Vis data
            target_kinetic_wavelength = kinetic_wavelength
            row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
            kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
                
            # Read steady-state data from CSV file
            steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
            
            # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
            target_C0_wavelength = C0_wavelength
            row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
        
            # Find row index for the desired Wavelength for pump
            target_pump_wavelength = pump_wavelength
            row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
        
            # Calculate initial absorbance at kinetic trace
            initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
        QY_ag_bromo_450nm_s2 = process_ag_bromo_450nm_s1(1305, 0.003, 450, 400, 450)

        function process_ag_bromo_450nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
            # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
        
            # Read UV-Vis data from CSV file
            # Kinetics obtained from Ocean optic series of UvVi's obtained every second
            uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
            # Extract kinetics data from UV-Vis data
            target_kinetic_wavelength = kinetic_wavelength
            row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
            kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
                
            # Read steady-state data from CSV file
            steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
            
            # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
            target_C0_wavelength = C0_wavelength
            row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
        
            # Find row index for the desired Wavelength for pump
            target_pump_wavelength = pump_wavelength
            row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
        
            # Calculate initial absorbance at kinetic trace
            initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
        QY_ag_bromo_450nm_s3 = process_ag_bromo_450nm_s1(1305, 0.003, 450, 400, 450)


#### Ag_Bromo_405nm_s1 #########################################################################################################

    # ϵ_405nm = 2203.7 for Ag(II)TFA
    function process_ag_bromo_405nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
    
        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
    
        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
    
        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
    QY_ag_bromo_405nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 405, 400, 405)

    function process_ag_bromo_405nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
    
        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
    
        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
    
        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
    QY_ag_bromo_405nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 405, 400, 405)

    function process_ag_bromo_405nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
    
        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
    
        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
    
        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
    QY_ag_bromo_405nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 405, 400, 405)


#### Ag_Bromo_525nm_s1 #########################################################################################################
    ## ϵ_450nm = 1297 for Ag(II)TFA

    function process_ag_bromo_525nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
    
        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
    
        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
    
        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
    QY_ag_bromo_525nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 450, 400, 525)

    function process_ag_bromo_525nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
    
        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
    
        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
    
        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
    QY_ag_bromo_525nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 450, 400, 525)

    function process_ag_bromo_525nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample
    
        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/s1/Absorbance_AgTGBP_Cl_450nm_s1.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_12/AgTGBP_Cl/AgTGBP_Cl_SteadyState_UvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))
    
        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))
    
        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 3] - steady_state_data[row_index_C0, 2]
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
    QY_ag_bromo_525nm_s1 = process_ag_bromo_450nm_s1(1305, 0.003, 450, 400, 525)

