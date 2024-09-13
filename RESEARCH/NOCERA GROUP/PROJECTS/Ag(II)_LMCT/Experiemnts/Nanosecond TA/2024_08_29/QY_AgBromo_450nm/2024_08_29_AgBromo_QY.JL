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


### Ag(II) Bromo 450 nm 10 deg QY ################################################################################################

#  ϵ_450 = 1358

    function process_ag_bromo_450nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/s1/AgBromo_450nm_s1_Absorbances.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/AgBromo_450_10deg_SteadystateUvVis.csv", header=true, DataFrame)
        
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
        blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/AgBromo_450nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
        # Calculate average of blank data
        blank_avg = mean(blank_data[:, 2])

        # Fit kinetics data using MonoFitIRF_New function
        fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 450, [0.08, 0.001, 20.418, 103.607, 0.033], 10, Inf, false, 10000)
        # Calculate rate constant from fit result
        k = 1 / fit_result[2][4]
        # Calculate flux using blank average
        flux = calculate_flux(blank_avg, pump_wavelength)
        # Calculate quantum yield
        qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

        return qy  # Return calculated quantum yield
    end
    QY_ag_bromo_450nm_s1 = process_ag_bromo_450nm_s1(1358, 0.003,450, 400, 450)  # Call the processing function and store the result

    function process_ag_bromo_450nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/s2/AgBromo_450nm_s2_Absorbances.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/AgBromo_450_10deg_SteadystateUvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))

        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))

        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 4] - steady_state_data[row_index_C0, 2]
        # Calculate initial absorbance at pump wavelength
        initial_abs_pump = steady_state_data[row_index_pump, 4] - steady_state_data[row_index_pump, 2]

        # Read blank data from CSV file
        blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/AgBromo_450nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
        # Calculate average of blank data
        blank_avg = mean(blank_data[:, 2])

        # Fit kinetics data using MonoFitIRF_New function
        fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 600, [0.08, 0.001, 21.005, 99.377, 0.038], 10, Inf, false, 10000)
        # Calculate rate constant from fit result
        k = 1 / fit_result[2][4]
        # Calculate flux using blank average
        flux = calculate_flux(blank_avg, pump_wavelength)
        # Calculate quantum yield
        qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

        return qy  # Return calculated quantum yield
    end
    QY_ag_bromo_450nm_s2 = process_ag_bromo_450nm_s2(1358, 0.003,450, 400, 450)  # Call the processing function and store the result

    function process_ag_bromo_450nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/s3/AgBromo_450nm_s3_Absorbances.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/AgBromo_450_10deg_SteadystateUvVis.csv", header=true, DataFrame)
        
        # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
        target_C0_wavelength = C0_wavelength
        row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))

        # Find row index for the desired Wavelength for pump
        target_pump_wavelength = pump_wavelength
        row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))

        # Calculate initial absorbance at kinetic trace
        initial_abs_kinetic = steady_state_data[row_index_C0, 5] - steady_state_data[row_index_C0, 2]
        # Calculate initial absorbance at pump wavelength
        initial_abs_pump = steady_state_data[row_index_pump, 5] - steady_state_data[row_index_pump, 2]

        # Read blank data from CSV file
        blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_450nm/AgBromo_450nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
        # Calculate average of blank data
        blank_avg = mean(blank_data[:, 2])

        # Fit kinetics data using MonoFitIRF_New function
        fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 600,[0.08, 0.001, 21.588, 101.735, 0.035], 10, Inf, false, 10000)
        # Calculate rate constant from fit result
        k = 1 / fit_result[2][4]
        # Calculate flux using blank average
        flux = calculate_flux(blank_avg, pump_wavelength)
        # Calculate quantum yield
        qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

        return qy  # Return calculated quantum yield
    end
    QY_ag_bromo_450nm_s3 = process_ag_bromo_450nm_s3(1358, 0.003,450, 400, 450)  # Call the processing function and store the result

    mean([QY_ag_bromo_450nm_s1,QY_ag_bromo_450nm_s2,QY_ag_bromo_450nm_s3])
    std([QY_ag_bromo_450nm_s1,QY_ag_bromo_450nm_s2,QY_ag_bromo_450nm_s3])


### Ag(II) Bromo 525 nm 10 deg QY ################################################################################################
#  ϵ_450 = 1358

function process_ag_bromo_525nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
    # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

    # Read UV-Vis data from CSV file
    # Kinetics obtained from Ocean optic series of UvVi's obtained every second
    uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/s1/AgBromo_525nm_s1_Absorbances.csv", header=false, DataFrame))
    # Extract kinetics data from UV-Vis data
    target_kinetic_wavelength = kinetic_wavelength
    row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
    kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
        
    # Read steady-state data from CSV file
    steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/AgBromo_525nm_SteadyStateAbsorbancescsv.csv", header=true, DataFrame)
    
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
    blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/AgBromo_525nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
    # Calculate average of blank data
    blank_avg = mean(blank_data[:, 2])

    # Fit kinetics data using MonoFitIRF_New function
    fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 600, [0.08, 0.001, 20.418, 103.607, 0.033], 10, Inf, false, 10000)
    # Calculate rate constant from fit result
    k = 1 / fit_result[2][4]
    # Calculate flux using blank average
    flux = calculate_flux(blank_avg, pump_wavelength)
    # Calculate quantum yield
    qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

    return qy  # Return calculated quantum yield
end
QY_ag_bromo_525nm_s1 = process_ag_bromo_525nm_s1(1358, 0.003,450, 400, 525)  # Call the processing function and store the result

function process_ag_bromo_525nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
    # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

    # Read UV-Vis data from CSV file
    # Kinetics obtained from Ocean optic series of UvVi's obtained every second
    uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/s2/QY_AgBromo_525nm_s2_Absorbances.csv", header=false, DataFrame))
    # Extract kinetics data from UV-Vis data
    target_kinetic_wavelength = kinetic_wavelength
    row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
    kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
        
    # Read steady-state data from CSV file
    steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/AgBromo_525nm_SteadyStateAbsorbancescsv.csv", header=true, DataFrame)
    
    # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
    target_C0_wavelength = C0_wavelength
    row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))

    # Find row index for the desired Wavelength for pump
    target_pump_wavelength = pump_wavelength
    row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))

    # Calculate initial absorbance at kinetic trace
    initial_abs_kinetic = steady_state_data[row_index_C0, 4] - steady_state_data[row_index_C0, 2]
    # Calculate initial absorbance at pump wavelength
    initial_abs_pump = steady_state_data[row_index_pump, 4] - steady_state_data[row_index_pump, 2]

    # Read blank data from CSV file
    blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/AgBromo_525nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
    # Calculate average of blank data
    blank_avg = mean(blank_data[:, 2])

    # Fit kinetics data using MonoFitIRF_New function
    fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 1000, [0.08, 0.001, 234.315, 549.492, 0.051], 10, Inf, false, 10000)
    # Calculate rate constant from fit result
    k = 1 / fit_result[2][4]
    # Calculate flux using blank average
    flux = calculate_flux(blank_avg, pump_wavelength)
    # Calculate quantum yield
    qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

    return qy  # Return calculated quantum yield
end
QY_ag_bromo_525nm_s2 = process_ag_bromo_525nm_s2(1358, 0.003,450, 400, 525)  # Call the processing function and store the result

function process_ag_bromo_525nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
    # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

    # Read UV-Vis data from CSV file
    # Kinetics obtained from Ocean optic series of UvVi's obtained every second
    uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/s3/QY_AgBromo_525nm_s3_Absorbances.csv", header=false, DataFrame))
    # Extract kinetics data from UV-Vis data
    target_kinetic_wavelength = kinetic_wavelength
    row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
    kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
        
    # Read steady-state data from CSV file
    steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/AgBromo_525nm_SteadyStateAbsorbancescsv.csv", header=true, DataFrame)
    
    # Find row index for the desired Wavelength for calculating inital concentration, this is the wavelengh for which ϵ is known
    target_C0_wavelength = C0_wavelength
    row_index_C0 = argmin(abs.(steady_state_data[:, 1] .- target_C0_wavelength))

    # Find row index for the desired Wavelength for pump
    target_pump_wavelength = pump_wavelength
    row_index_pump = argmin(abs.(steady_state_data[:, 1] .- target_pump_wavelength))

    # Calculate initial absorbance at kinetic trace
    initial_abs_kinetic = steady_state_data[row_index_C0, 5] - steady_state_data[row_index_C0, 2]
    # Calculate initial absorbance at pump wavelength
    initial_abs_pump = steady_state_data[row_index_pump, 5] - steady_state_data[row_index_pump, 2]

    # Read blank data from CSV file
    blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_08_29/QY_AgBromo_525nm/AgBromo_525nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
    # Calculate average of blank data
    blank_avg = mean(blank_data[:, 2])

    # Fit kinetics data using MonoFitIRF_New function
    fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 1500,[0.08, 0.001, 241.354, 567.394, 0.051], 10, Inf, false, 10000)
    # Calculate rate constant from fit result
    k = 1 / fit_result[2][4]
    # Calculate flux using blank average
    flux = calculate_flux(blank_avg, pump_wavelength)
    # Calculate quantum yield
    qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

    return qy  # Return calculated quantum yield
end
QY_ag_bromo_525nm_s3 = process_ag_bromo_525nm_s3(1358, 0.003,450, 400, 525)  # Call the processing function and store the result

mean([QY_ag_bromo_525nm_s1,QY_ag_bromo_525nm_s2,QY_ag_bromo_525nm_s3])
std([QY_ag_bromo_525nm_s1,QY_ag_bromo_525nm_s2,QY_ag_bromo_525nm_s3])