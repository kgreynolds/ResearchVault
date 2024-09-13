include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Inside(),
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

uv_vis_data_Ag_DFA_450nm =  Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/s1/QY_ag_CF2_450nm_s1_Absorbances.csv", header=false, DataFrame))
uv_vis_data_Ag_DFA_525nm = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/s1/Ag_DFA_525nm_s1_Absorbances.csv", header=false, DataFrame))

x = uv_vis_data_Ag_DFA_450nm[1,2:end]


function UvVis_Fig(f = Figure())

    data = uv_vis_data_Ag_DFA_525nm
    x = uv_vis_data_Ag_DFA_525nm[1,2:end]

    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.winter, 12),), ## or replace ncol with Integer 

                title = "Ag(II)CF₂H Absorbance Decay 525 nm pump",
                xlabel = "Wavelength (nm)", 
                ylabel = "A (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(300,600)
            # CairoMakie.ylims!(-0.02,0.7)

        ### PlotUvVis  Data ##############################################################################


            Makie.lines!(ax, x, data[2, 2:end],linewidth = 3, label = "0 s")
            Makie.lines!(ax, x, data[21, 2:end],linewidth = 3, label = "20 s")
            Makie.lines!(ax, x, data[41, 2:end],linewidth = 3, label = "40 s")
            Makie.lines!(ax, x, data[61, 2:end],linewidth = 3, label = "60 s")
            Makie.lines!(ax, x, data[81, 2:end],linewidth = 3, label = "80 s")
            Makie.lines!(ax, x, data[101, 2:end],linewidth = 3, label = "100 s")
            Makie.lines!(ax, x, data[126, 2:end],linewidth = 3, label = "125 s")
            Makie.lines!(ax, x, data[151, 2:end],linewidth = 3, label = "150 s")
            Makie.lines!(ax, x, data[201, 2:end],linewidth = 3, label = "200 s")
            Makie.lines!(ax, x, data[251, 2:end],linewidth = 3, label = "250 s")
            Makie.lines!(ax, x, data[351, 2:end],linewidth = 3, label = "350 s")

  

            axislegend("Delay time",
                position = :rt, 
                nbanks = 3, 
                framevisible = false, 
                labelsize = 15
            )

            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            # CairoMakie.vspan!(520, 530; ymin = 0.0, ymax = 0.5, color = :lightblue)

            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
        
    resize_to_layout!(f)
    f
end
UvVis_Fig()
save("JULIA/Outputs/Ag(II)_DFA Absorbance Decay 525 nm.png",UvVis_Fig())


### Ag(II) CF2 450 nm 10 deg QY ################################################################################################
 # ϵ_450  = 1214
    function process_ag_CF2_450nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/s1/QY_ag_CF2_450nm_s1_Absorbances.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/QY_ag_CF2_450nm_Steadystate.csv", header=true, DataFrame)
        
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
        blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/Ag_DFA_450nm_blank_energy.csv", skipto=6, header=true, DataFrame)
        # Calculate average of blank data
        blank_avg = mean(blank_data[:, 2])

        # Fit kinetics data using MonoFitIRF_New function
        fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 200, [0.08, 0.001, 6.332, 36.919, 0.048], 10, Inf, false, 10000)
        # Calculate rate constant from fit result
        k = 1 / fit_result[2][4]
        # Calculate flux using blank average
        flux = calculate_flux(blank_avg, pump_wavelength)
        # Calculate quantum yield
        qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

        return qy  # Return calculated quantum yield
    end
    QY_ag_CF2_450nm_s1 = process_ag_CF2_450nm_s1(1214, 0.003,450, 410, 450)  # Call the processing function and store the result

    function process_ag_CF2_450nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/s2/QY_ag_CF2_450nm_s2_absorbances.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/QY_ag_CF2_450nm_Steadystate.csv", header=true, DataFrame)
        
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
        blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/Ag_DFA_450nm_blank_energy.csv", skipto=6, header=true, DataFrame)
        # Calculate average of blank data
        blank_avg = mean(blank_data[:, 2])

        # Fit kinetics data using MonoFitIRF_New function
        fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 200,[0.08, 0.001, 6.583, 36.494, 0.055], 10, Inf, false, 10000)
        # Calculate rate constant from fit result
        k = 1 / fit_result[2][4]
        # Calculate flux using blank average
        flux = calculate_flux(blank_avg, pump_wavelength)
        # Calculate quantum yield
        qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

        return qy  # Return calculated quantum yield
    end
    QY_ag_CF2_450nm_s2 = process_ag_CF2_450nm_s2(1214, 0.003,450, 410, 450)  # Call the processing function and store the result

    function process_ag_CF2_450nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
        # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

        # Read UV-Vis data from CSV file
        # Kinetics obtained from Ocean optic series of UvVi's obtained every second
        uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/s3/QY_ag_CF2_450nm_s3_absorbances.csv", header=false, DataFrame))
        # Extract kinetics data from UV-Vis data
        target_kinetic_wavelength = kinetic_wavelength
        row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
        kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
            
        # Read steady-state data from CSV file
        steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/QY_ag_CF2_450nm_Steadystate.csv", header=true, DataFrame)
        
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
        blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_450nm_QY/Ag_DFA_450nm_blank_energy.csv", skipto=6, header=true, DataFrame)
        # Calculate average of blank data
        blank_avg = mean(blank_data[:, 2])

        # Fit kinetics data using MonoFitIRF_New function
        fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 200,[0.08, 0.001, 6.085, 36.762, 0.059], 10, Inf, false, 10000)
        # Calculate rate constant from fit result
        k = 1 / fit_result[2][4]
        # Calculate flux using blank average
        flux = calculate_flux(blank_avg, pump_wavelength)
        # Calculate quantum yield
        qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

        return qy  # Return calculated quantum yield
    end
    QY_ag_CF2_450nm_s3 = process_ag_CF2_450nm_s3(1214, 0.003,450, 410, 450)  # Call the processing function and store the result

    mean([QY_ag_CF2_450nm_s1,QY_ag_CF2_450nm_s2,QY_ag_CF2_450nm_s3])
    std([QY_ag_CF2_450nm_s1,QY_ag_CF2_450nm_s2,QY_ag_CF2_450nm_s3])


### Ag(II) CF2 525 nm 10 deg QY ################################################################################################
 # ϵ_450  = 1214
 function process_ag_CF2_525nm_s1(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
    # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

    # Read UV-Vis data from CSV file
    # Kinetics obtained from Ocean optic series of UvVi's obtained every second
    uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/s1/Ag_DFA_525nm_s1_Absorbances.csv", header=false, DataFrame))
    # Extract kinetics data from UV-Vis data
    target_kinetic_wavelength = kinetic_wavelength
    row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
    kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
        
    # Read steady-state data from CSV file
    steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/Ag_DFA_525nm_SteadyState_UvVis.csv", header=true, DataFrame)
    
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
    blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/Ag_DFA_525nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
    # Calculate average of blank data
    blank_avg = mean(blank_data[:, 2])

    # Fit kinetics data using MonoFitIRF_New function
    fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 200,[0.08, 0.001, 35.117, 121.222, 0.089], 10, Inf, false, 10000)
    # Calculate rate constant from fit result
    k = 1 / fit_result[2][4]
    # Calculate flux using blank average
    flux = calculate_flux(blank_avg, pump_wavelength)
    # Calculate quantum yield
    qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

    return qy  # Return calculated quantum yield
end
QY_ag_CF2_525nm_s1 = process_ag_CF2_525nm_s1(1214, 0.003,450, 450, 525)  # Call the processing function and store the result

function process_ag_CF2_525nm_s2(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
    # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

    # Read UV-Vis data from CSV file
    # Kinetics obtained from Ocean optic series of UvVi's obtained every second
    uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/s2/Ag_DFA_525nm_s2_Absorbances.csv", header=false, DataFrame))
    # Extract kinetics data from UV-Vis data
    target_kinetic_wavelength = kinetic_wavelength
    row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
    kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
        
    # Read steady-state data from CSV file
    steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/Ag_DFA_525nm_SteadyState_UvVis.csv", header=true, DataFrame)
    
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
    blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/Ag_DFA_525nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
    # Calculate average of blank data
    blank_avg = mean(blank_data[:, 2])

    # Fit kinetics data using MonoFitIRF_New function
    fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 400,[0.08, 0.001, 38.23, 122.618, 0.1], 10, Inf, false, 10000)
    # Calculate rate constant from fit result
    k = 1 / fit_result[2][4]
    # Calculate flux using blank average
    flux = calculate_flux(blank_avg, pump_wavelength)
    # Calculate quantum yield
    qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

    return qy  # Return calculated quantum yield
end
QY_ag_CF2_525nm_s2 = process_ag_CF2_525nm_s2(1214, 0.003,450, 450, 525)  # Call the processing function and store the result

function process_ag_CF2_525nm_s3(ϵ, volume, C0_wavelength,kinetic_wavelength, pump_wavelength)
    # ϵ = molar absorbtivity at wavelength of kinetic, volume = volume of sample

    # Read UV-Vis data from CSV file
    # Kinetics obtained from Ocean optic series of UvVi's obtained every second
    uv_vis_data = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/s3/Ag_DFA_525nm_s3_Absorbances.csv", header=false, DataFrame))
    # Extract kinetics data from UV-Vis data
    target_kinetic_wavelength = kinetic_wavelength
    row_index_kinetic = argmin(abs.(uv_vis_data[1, :] .- target_kinetic_wavelength))
    kinetics_data = hcat(uv_vis_data[2:end, 1], uv_vis_data[2:end, row_index_kinetic])
        
    # Read steady-state data from CSV file
    steady_state_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/Ag_DFA_525nm_SteadyState_UvVis.csv", header=true, DataFrame)
    
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
    blank_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_09_09/Ag_CF2H_525nm_QY/Ag_DFA_525nm_Blank_energy.csv", skipto=6, header=true, DataFrame)
    # Calculate average of blank data
    blank_avg = mean(blank_data[:, 2])

    # Fit kinetics data using MonoFitIRF_New function
    fit_result = MonoFitIRF_New(kinetics_data, 2, 5, 400,[0.08, 0.001, 34.116, 116.913, 0.1], 10, Inf, false, 10000)
    # Calculate rate constant from fit result
    k = 1 / fit_result[2][4]
    # Calculate flux using blank average
    flux = calculate_flux(blank_avg, pump_wavelength)
    # Calculate quantum yield
    qy = quantum_yield(k, ϵ, volume, flux, initial_abs_pump, initial_abs_kinetic)

    return qy  # Return calculated quantum yield
end
QY_ag_CF2_525nm_s3 = process_ag_CF2_525nm_s3(1214, 0.003,450, 450, 525)  # Call the processing function and store the result

mean([QY_ag_CF2_525nm_s1,QY_ag_CF2_525nm_s2,QY_ag_CF2_525nm_s3])
std([QY_ag_CF2_525nm_s1,QY_ag_CF2_525nm_s2,QY_ag_CF2_525nm_s3])