using Interpolations,Plots,Statistics, DataFrames, CSV, Optim,LinearAlgebra, ColorSchemes,MAT,LineSearches,SpecialFunctions, EasyFit,CairoMakie; CairoMakie.activate!()
plotly()


function MonoFitIRF_Global(file, columns::Vector{Int}, t0, uppertimebound, shared_params, initial_amplitudes, bound1, bound2, guess,num_optimizations,noise_level, size)
    
    # Resize function to interpolate and resize the kinetic trace
    function Resize(file, column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1], file[end,1], size))
        Interp  = LinearInterpolation(file[:,1], file[:,column])
        y_out = Interp(x_out)
        ResizedDF = DataFrame(x = x_out, y = y_out)
        return ResizedDF
    end

    # Model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁)
        @. ((p[4] + A₁ / p[3] * exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])) * 
            (erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2)) + 1) / 2)
    end

    # Cut the data within time bounds and resize
    cut = t0 .< file[:,1] .< uppertimebound
    t = file[:,1][cut]
    
    # Store resized data for all selected columns
    KineticTraces = []
    for col in columns
        trace_data = file[:,col][cut]
        Cut_Data = hcat(t, trace_data)
        Resized_Cut_Data = Resize(Cut_Data, 2, size)
        push!(KineticTraces, Resized_Cut_Data)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = params[1:4] # Shared lifetimes (xc, w, τ₁, A_inf)
        amplitudes = params[5:end] # Independent amplitudes for each trace
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            model_values = Y(shared_coeffs, x, amplitudes[i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess == true
        # Plot each trace and initial guess
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            p_raw = Plots.plot!(x, KineticTrace, label="Data Trace $i")
            p_guess = Plots.plot!(x, Y(shared_params, x, initial_amplitudes[i]), label="Initial Guess")
            display(p_raw)
            display(p_guess)
        end
    else

        # Set optimization bounds for shared lifetimes and amplitudes
        lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1/100)),(shared_params[2] - abs(shared_params[2]) * (bound1/100)), (shared_params[3] - abs(shared_params[3] * (bound2/100))), -Inf]
        upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1/100)), (shared_params[2] + abs(shared_params[2]) * (bound1/100)),  (shared_params[3] + abs(shared_params[3] * (bound2/100))), Inf]
        
        lower_amplitudes = fill(-Inf, length(columns)) # Allow amplitudes to vary
        upper_amplitudes = fill(Inf, length(columns))

        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        # num_optimizations = 10
        # noise_level = 0.25
        results = []
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        for _ in 1:num_optimizations
            # Perform global fit using the optimizer
            # Add noise to the initial guess
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        
            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds,noisy_coeffs,Fminbox(inner_optimizer),Optim.Options(time_limit = 60,g_tol = 1e-12))
            best_fit_params = Optim.minimizer(result)
            push!(results, Optim.minimizer(result))
            print(";Coeffs: ", noisy_coeffs[3] )
        end

        avg_values = round.(mean(results), digits = 3)
        std_values = round.(std(results), digits  =3)

        # Extract the best-fit parameters
        best_lifetimes = round.(avg_values[1:4], digits = 3)
        best_amplitudes = round.(avg_values[5:end], digits = 3)

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = round(1 - total_ss_residual / total_ss, digits = 3)

        # Output the fit coefficients and R²
        println("Best Fit Lifetimes: ", best_lifetimes)
        println("Best Fit Amplitudes: ", best_amplitudes)
        println("Global R²: ", R²)

        # Plot the data, initial guess, and final fit for each trace
        Fit_Fig = Figure(font = "", figure_padding = 25, fontsize = 20)
        width = 3
        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])" , subtitle = "τ₁ = $(best_lifetimes[3]) ± $(std_values[3])",
        palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),),
        xlabel = "Time")

        ax2 = CairoMakie.Axis(Fit_Fig[2,1],palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),), title = "Residuals, R² = $R²")

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values
        
            lines!(ax1, x, KineticTrace, linewidth = width, label = "Data $i")
            lines!(ax1, x, fit_values, linewidth = width/2, linestyle = :dash,color = :black)

            # Calculate and plot Bollinger Bands
            upper_band = fit_values .+ 2000 .* std_values[3] .* abs.(fit_values)
            lower_band = fit_values .- 2000 .* std_values[3] .* abs.(fit_values)
            CairoMakie.band!(ax1,x, lower_band, upper_band, color = (:blue, 0.2))

            
            lines!(ax2, x, residuals, linewidth = width)

            linkxaxes!(ax1, ax2)

        end
        axislegend(ax1, position = :rb,nbanks =3, framevisible = false)

            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,130)
            resize_to_layout!(Fit_Fig)
        display(Fit_Fig)
        # save("JULIA/Outputs/GlobalFit_Fig.png", Fit_Fig)
        return best_lifetimes, best_amplitudes, R²
    end
end

function MonoFitIRF_Global_mod(file, columns::Vector{Int}, t0, uppertimebound, shared_params, initial_amplitudes, bound1, bound2, guess,num_optimizations,noise_level, size)
    
    # Resize function to interpolate and resize the kinetic trace
    function Resize(file, column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1], file[end,1], size))
        Interp  = LinearInterpolation(file[:,1], file[:,column])
        y_out = Interp(x_out)
        ResizedDF = DataFrame(x = x_out, y = y_out)
        return ResizedDF
    end

    # Model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A)
        @. ((p[4] + A[1] / p[3] * exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])) * 
            (erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2)) + 1) / 2)
    end

    # Cut the data within time bounds and resize
    cut = t0 .< file[:,1] .< uppertimebound
    t = file[:,1][cut]
    
    # Store resized data for all selected columns
    KineticTraces = []
    for col in columns
        trace_data = file[:,col][cut]
        Cut_Data = hcat(t, trace_data)
        Resized_Cut_Data = Resize(Cut_Data, 2, size)
        push!(KineticTraces, Resized_Cut_Data)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = params[1:4] # Shared lifetimes (xc, w, τ₁, A_inf)
        amplitudes = params[5:end] # Independent amplitudes for each trace
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            model_values = Y(shared_coeffs, x, amplitudes[i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess == true
        # Plot each trace and initial guess
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            p_raw = Plots.plot!(x, KineticTrace, label="Data Trace $i")
            p_guess = Plots.plot!(x, Y(shared_params, x, initial_amplitudes[i]), label="Initial Guess")
            display(p_raw)
            display(p_guess)
        end
    else

        # Set optimization bounds for shared lifetimes and amplitudes
        lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1/100)),(shared_params[2] - abs(shared_params[2]) * (bound1/100)), (shared_params[3] - abs(shared_params[3] * (bound2/100))), -Inf]
        upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1/100)), (shared_params[2] + abs(shared_params[2]) * (bound1/100)),  (shared_params[3] + abs(shared_params[3] * (bound2/100))), Inf]
        
        lower_amplitudes = fill(-Inf, length(columns)) # Allow amplitudes to vary
        upper_amplitudes = fill(Inf, length(columns))

        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        # num_optimizations = 10
        # noise_level = 0.25
        results = []
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        for _ in 1:num_optimizations
            # Perform global fit using the optimizer
            # Add noise to the initial guess
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        
            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds,noisy_coeffs,Fminbox(inner_optimizer),Optim.Options(time_limit = 60,g_tol = 1e-12))
            best_fit_params = Optim.minimizer(result)
            push!(results, Optim.minimizer(result))
            print(";Coeffs: ", noisy_coeffs[3] )
        end

        avg_values = round.(mean(results), digits = 3)
        std_values = round.(std(results), digits  =3)

        # Extract the best-fit parameters
        best_lifetimes = round.(avg_values[1:4], digits = 3)
        best_amplitudes = round.(avg_values[5:end], digits = 3)

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = round(1 - total_ss_residual / total_ss, digits = 3)

        # Output the fit coefficients and R²
        println("Best Fit Lifetimes: ", best_lifetimes)
        println("Best Fit Amplitudes: ", best_amplitudes)
        println("Global R²: ", R²)

        # Plot the data, initial guess, and final fit for each trace
        Fit_Fig = Figure(font = "", figure_padding = 25, fontsize = 20)
        width = 3
        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])" , subtitle = "τ₁ = $(best_lifetimes[3]) ± $(std_values[3])",
        palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),),
        xlabel = "Time")

        ax2 = CairoMakie.Axis(Fit_Fig[2,1],palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),), title = "Residuals, R² = $R²")

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values
        
            lines!(ax1, x, KineticTrace, linewidth = width, label = "Data $i")
            lines!(ax1, x, fit_values, linewidth = width/2, linestyle = :dash,color = :black)

            # Calculate and plot Bollinger Bands
            upper_band = fit_values .+ 2000 .* std_values[3] .* abs.(fit_values)
            lower_band = fit_values .- 2000 .* std_values[3] .* abs.(fit_values)
            CairoMakie.band!(ax1,x, lower_band, upper_band, color = (:blue, 0.2))

            
            lines!(ax2, x, residuals, linewidth = width)

            linkxaxes!(ax1, ax2)

        end
        axislegend(ax1, position = :rb,nbanks =3, framevisible = false)

            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,130)
            resize_to_layout!(Fit_Fig)
        display(Fit_Fig)
        # save("JULIA/Outputs/GlobalFit_Fig.png", Fit_Fig)
        return best_lifetimes, best_amplitudes, R²
    end
end

function biFitIRF_Global_mod(file, columns::Vector{Int}, t0, uppertimebound, shared_params, initial_amplitudes, bound1, bound2, guess,num_optimizations,noise_level, size)
    
    # Resize function to interpolate and resize the kinetic trace
    function Resize(file, column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1], file[end,1], size))
        Interp  = LinearInterpolation(file[:,1], file[:,column])
        y_out = Interp(x_out)
        ResizedDF = DataFrame(x = x_out, y = y_out)
        return ResizedDF
    end

    # Model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A)
        @. ((p[5]+ A[1] / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 
        + (A[2] /p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])*(erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))
    end

    # @. ((p[5]+ A[1] / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 
    # + (A[2] /p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])*(erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))

    # Cut the data within time bounds and resize
    cut = t0 .< file[:,1] .< uppertimebound
    t = file[:,1][cut]
    
    # Store resized data for all selected columns
    KineticTraces = []
    for col in columns
        trace_data = file[:,col][cut]
        Cut_Data = hcat(t, trace_data)
        Resized_Cut_Data = Resize(Cut_Data, 2, size)
        push!(KineticTraces, Resized_Cut_Data)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = params[1:5] # Shared lifetimes (xc, w, τ₁, A_inf)
        amplitudes = params[6:end] # Independent amplitudes for each trace
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            model_values = Y(shared_coeffs, x, amplitudes[:,i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess == true
        # Plot each trace and initial guess
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            p_raw = Plots.plot!(x, KineticTrace, label="Data Trace $i")
            p_guess = Plots.plot!(x, Y(shared_params, x, initial_amplitudes[:,i]), label="Initial Guess")
            display(p_raw)
            display(p_guess)
        end
    else

        # Set optimization bounds for shared lifetimes and amplitudes
        lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1/100)),(shared_params[2] - abs(shared_params[2]) * (bound1/100)), (shared_params[3] - abs(shared_params[3] * (bound2/100))),(shared_params[4] - abs(shared_params[4] * (bound2/100))), -Inf]
        upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1/100)), (shared_params[2] + abs(shared_params[2]) * (bound1/100)),  (shared_params[3] + abs(shared_params[3] * (bound2/100))),(shared_params[4] + abs(shared_params[4] * (bound2/100))), Inf]
        
        lower_amplitudes = fill(-Inf, 2*length(columns)) # Allow amplitudes to vary
        upper_amplitudes = fill(Inf, 2*length(columns))

        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        # num_optimizations = 10
        # noise_level = 0.25
        results = []
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        for _ in 1:num_optimizations
            # Perform global fit using the optimizer
            # Add noise to the initial guess
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        
            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds,noisy_coeffs,Fminbox(inner_optimizer),Optim.Options(time_limit = 60,g_tol = 1e-12))
            best_fit_params = Optim.minimizer(result)
            push!(results, Optim.minimizer(result))
            print(";Coeffs: ", noisy_coeffs[3] )
        end

        avg_values = round.(mean(results), digits = 3)
        std_values = round.(std(results), digits  =3)

        # Extract the best-fit parameters
        best_lifetimes = round.(avg_values[1:5], digits = 3)
        best_amplitudes = round.(avg_values[6:end], digits = 3)

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[:,i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = round(1 - total_ss_residual / total_ss, digits = 3)

        # Output the fit coefficients and R²
        println("Best Fit Lifetimes: ", best_lifetimes)
        println("Best Fit Amplitudes: ", best_amplitudes)
        println("Global R²: ", R²)

        # Plot the data, initial guess, and final fit for each trace
        Fit_Fig = Figure(font = "", figure_padding = 25, fontsize = 20)
        width = 3
        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])" , subtitle = "τ₁ = $(best_lifetimes[3]) ± $(std_values[3])",
        palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),),
        xlabel = "Time")

        ax2 = CairoMakie.Axis(Fit_Fig[2,1],palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),), title = "Residuals, R² = $R²")

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[:,i])
            residuals = KineticTrace .- fit_values
        
            lines!(ax1, x, KineticTrace, linewidth = width, label = "Data $i")
            lines!(ax1, x, fit_values, linewidth = width/2, linestyle = :dash,color = :black)

            # Calculate and plot Bollinger Bands
            upper_band = fit_values .+ 2000 .* std_values[3] .* abs.(fit_values)
            lower_band = fit_values .- 2000 .* std_values[3] .* abs.(fit_values)
            CairoMakie.band!(ax1,x, lower_band, upper_band, color = (:blue, 0.2))

            
            lines!(ax2, x, residuals, linewidth = width)

            linkxaxes!(ax1, ax2)

        end
        axislegend(ax1, position = :rb,nbanks =3, framevisible = false)

            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,130)
            resize_to_layout!(Fit_Fig)
        display(Fit_Fig)
        # save("JULIA/Outputs/GlobalFit_Fig.png", Fit_Fig)
        return best_lifetimes, best_amplitudes, R²
    end
end

function BiExpFitIRF_Global(file, columns::Vector{Int}, t0, uppertimebound, shared_params, initial_amplitudes, bound1, bound2, guess, num_optimizations, noise_level, size)
    
    # Resize function to interpolate and resize the kinetic trace
    function Resize(file, column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1], file[end,1], size))
        Interp  = LinearInterpolation(file[:,1], file[:,column])
        y_out = Interp(x_out)
        ResizedDF = DataFrame(x = x_out, y = y_out)
        return ResizedDF
    end

    # Bi-exponential model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁, A₂)
        @. ((p[5]+ A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
        (A₂ /p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])*(erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))
    end

    # @. ((p[5]+ A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
    #  (A₂ /p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])*(erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))

    # Cut the data within time bounds and resize
    cut = t0 .< file[:,1] .< uppertimebound
    t = file[:,1][cut]
    
    # Store resized data for all selected columns
    KineticTraces = []
    for col in columns
        trace_data = file[:,col][cut]
        Cut_Data = hcat(t, trace_data)
        Resized_Cut_Data = Resize(Cut_Data, 2, size)
        push!(KineticTraces, Resized_Cut_Data)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = params[1:5] # Shared lifetimes (xc, w, τ₁, τ₂, A_inf)
        amplitudes = params[6:end] # Independent amplitudes for each trace
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            model_values = Y(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess == true
        # Plot each trace and initial guess
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            p_raw = Plots.plot!(x, KineticTrace, label="Data Trace $i")
            p_guess = Plots.plot!(x, Y(shared_params, x, initial_amplitudes[2i-1], initial_amplitudes[2i]), label="Initial Guess")
            display(p_raw)
            display(p_guess)
        end
    else
        # Set optimization bounds for shared lifetimes and amplitudes
        lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1/100)),
                           (shared_params[2] - abs(shared_params[2]) * (bound1/100)),
                           (shared_params[3] - abs(shared_params[3]) * (bound2/100)),
                           (shared_params[4] - abs(shared_params[4]) * (bound2/100)),
                           -Inf]
        upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1/100)),
                           (shared_params[2] + abs(shared_params[2]) * (bound1/100)),
                           (shared_params[3] + abs(shared_params[3]) * (bound2/100)),
                           (shared_params[4] + abs(shared_params[4]) * (bound2/100)),
                           Inf]
        
        lower_amplitudes = fill(-Inf, length(initial_amplitudes))
        upper_amplitudes = fill(Inf, length(initial_amplitudes))

        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        results = []
        noisy_coeffs = []
        for _ in 1:num_optimizations
            # Add noise to the initial guess
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))

            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer), Optim.Options(time_limit = 60, g_tol = 1e-12))
            push!(results, Optim.minimizer(result))
            # print(";Coeffs: ", noisy_coeffs)
        end

        avg_values = round.(mean(results), digits = 3)
        std_values = round.(std(results), digits = 3)

        # Extract the best-fit parameters
        best_lifetimes = round.(avg_values[1:5], digits = 3)
        best_amplitudes = round.(avg_values[6:end], digits = 3)

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = round(1 - total_ss_residual / total_ss, digits = 3)

        # Output the fit coefficients and R²
        println("Best Fit Lifetimes: ", best_lifetimes)
        println("Best Fit Amplitudes: ", best_amplitudes)
        println("Stdev: ", std_values)
        println("Global R²: ", R²)

        # Plot the data, initial guess, and final fit for each trace
        Fit_Fig = Figure(font = "", figure_padding = 25, fontsize = 20)
        width = 3
        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(best_lifetimes[1]) ± $(std_values[1]); FWHM = $(best_lifetimes[2])± $(std_values[2])" , 
                              subtitle = "τ₁ = $(best_lifetimes[3]) ± $(std_values[3]); τ₂ = $(best_lifetimes[4]) ± $(std_values[4])",
                              palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),),
                              xlabel = "Time")

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], palette = (color = palette(ColorSchemes.hsv, length(KineticTraces) +1),), 
                              title = "Residuals, R² = $R²")

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
            residuals = KineticTrace .- fit_values
        
            lines!(ax1, x, KineticTrace, linewidth = width, label = "Data $i")
            lines!(ax1, x, fit_values, linewidth = width/2, linestyle = :dash, color = :black)

            # Calculate and plot Bollinger Bands
            # upper_band = fit_values .+ 2 .* (std_values[3] + std_values[4]) .* abs.(fit_values)
            # lower_band = fit_values .- 2 .* (std_values[3] + std_values[4]) .* abs.(fit_values)


            # CairoMakie.band!(ax1, x, lower_band, upper_band, color = (:black, 0.2))
            
            lines!(ax2, x, residuals, linewidth = width)

            linkxaxes!(ax1, ax2)
        end
        axislegend(ax1, position = :rb, nbanks = 3, framevisible = false)

        colsize!(Fit_Fig.layout, 1, Aspect(1, 1.5)) # Set aspect ratio of the y vs. x axis
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout, 2, 130)
        resize_to_layout!(Fit_Fig)
        display(Fit_Fig)
        # save("JULIA/Outputs/GlobalFit_Fig.png", Fit_Fig)
        # return best_lifetimes, best_amplitudes, R²
    end
end

globalfit_test  = MonoFitIRF_Global(Zn_Chlorin_b2_Kinetics, [2, 3, 4], -10.0, 3500, [4.675, 1.998, 763.324, -0.0], [-6.729, -86.094, 17.29], 100, Inf, false,10,0.1, 200)
globalfit_test  = MonoFitIRF_Global_mod(Zn_Chlorin_b2_Kinetics, [2, 3], -10.0, 3500, [4.675, 1.998, 763.324, -0.0], [-6.729, -86.094], 100, Inf, false,10,0.1, 200)
globalfit_test  = BiExpFitIRF_Global(Zn_Chlorin_b2_Kinetics, [2,3,4], -100.0, 4500, [3.612, 1.571, 522.744, 875.771, 0.0],[-2.094, -4.846, -19.742, -67.96, 4.078, 13.307], 10, 1000, false,5,0.1, 200)


Zn_Chlorin_b2_399nm = MonoFitIRF_New_Global(Zn_Chlorin_b2_Kinetics,[2,4],-100.0, 4700, [0.001, 0.981, -6.815, 754.606, 0.0], Inf,Inf,true,200)
Zn_Chlorin_b2_399nm = MonoFitIRF_New_Global(Zn_Chlorin_b2_Kinetics,[2,3,4],-100.0, 4700, [5.542, 1.236, 417.25, -0.0], [6, 100.0], Inf,Inf,false,200)

Zn_Chlorin_b2_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_TA_Kinetics.csv", header = true, DataFrame)
Zn_Chlorin_b2_AIRFREE_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_AIRFREE_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)

file = Zn_Chlorin_b2_Kinetics
t0 = -10.0
uppertimebound = 3500
columns = [2, 3, 4,6]
Size = 200
bound1 = 100
bound2 = Inf

    # Resize function to interpolate and resize the kinetic trace
    function Resize(file, column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1], file[end,1], size))
        Interp  = LinearInterpolation(file[:,1], file[:,column])
        y_out = Interp(x_out)
        ResizedDF = DataFrame(x = x_out, y = y_out)
        return ResizedDF
    end

    # Model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁)
        @. ((p[4] + A₁ / p[3] * exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])) * 
            (erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2)) + 1) / 2)
    end

    # Cut the data within time bounds and resize
    cut = t0 .< file[:,1] .< uppertimebound
    t = file[:,1][cut]
    
    # Store resized data for all selected columns
    KineticTraces = []
    for col in columns
        trace_data = file[:,col][cut]
        Cut_Data = hcat(t, trace_data)
        Resized_Cut_Data = Resize(Cut_Data, 2, Size)
        push!(KineticTraces, Resized_Cut_Data)
    end
    KineticTraces
    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        lifetimes = params[1:4] # Shared lifetimes (xc, w, τ₁)
        amplitudes = params[5:end] # Independent amplitudes for each trace
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            model_values = Y(lifetimes, x, amplitudes[i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes

    shared_lifetimes = [4.681482105411228, 2.00889735130653, 764.268248757977, 3.0959655939649345e-6]
    initial_amplitudes = [-1.865798620424938, -56.28262362041595, 7.169060293822703, -26.19677916193525]
    initial_params = vcat(shared_lifetimes, initial_amplitudes)

  
        # Plot each trace and initial guess
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            p_raw = Plots.plot!(x, KineticTrace, label="Data Trace $i")
            p_guess = Plots.plot!(x, Y(shared_lifetimes, x, initial_amplitudes[i]), label="Initial Guess")
            display(p_raw)
            display(p_guess)
        end

        # Set optimization bounds for shared lifetimes and amplitudes
        lower_lifetimes = [(shared_lifetimes[1] - abs(shared_lifetimes[1]) * (bound1/100)),(shared_lifetimes[2] - abs(shared_lifetimes[2]) * (bound1/100)), (shared_lifetimes[3] - abs(shared_lifetimes[3] * (bound2/100))), -Inf]
        upper_lifetimes = [(shared_lifetimes[1] + abs(shared_lifetimes[1]) * (bound1/100)), (shared_lifetimes[2] + abs(shared_lifetimes[2]) * (bound1/100)), (shared_lifetimes[3] + abs(shared_lifetimes[3] * (bound2/100))), Inf]
        
        lower_amplitudes = fill(-Inf, length(columns)) # Allow amplitudes to vary
        upper_amplitudes = fill(Inf, length(columns))

        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        num_optimizations = 10
        noise_level = 0.25
        results = []
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        for _ in 1:num_optimizations
            # Perform global fit using the optimizer
            # Add noise to the initial guess
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        
            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds,noisy_coeffs,Fminbox(inner_optimizer),Optim.Options(time_limit = 60,g_tol = 1e-12))
            best_fit_params = Optim.minimizer(result)
            push!(results, Optim.minimizer(result))
            print(";Coeffs: ", noisy_coeffs[3] )
        end
        results
    
    
        avg_values = mean(results)
        std_values = std(results)

        # Extract the best-fit parameters
        best_lifetimes = avg_values[1:4]
        best_amplitudes = avg_values[5:end]

        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            p_raw = Plots.plot!(x, KineticTrace, label="Data Trace $i")
            p_guess = Plots.plot!(x, Y(best_lifetimes, x, best_amplitudes[i]), label="FIt")
            display(p_raw)
            display(p_guess)
        end

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:,1]
            KineticTrace = trace[:,2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = 1 - total_ss_residual / total_ss

# Output the fit coefficients and R²
println("Best Fit Lifetimes: ", best_lifetimes)
println("Best Fit Amplitudes: ", best_amplitudes)
println("Global R²: ", R²)

KineticTraces[1]
t
    # Plot the results with Bollinger Bands
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Intensity")

    # Plot raw data
    lines!(ax, KineticTraces[1][:,1], KineticTraces[1][:,2], label = "Raw Data")

    # Plot average fit
    avg_fit = Y(best_lifetimes, KineticTraces[1][:,1], best_amplitudes[1])
    lines!(ax, KineticTraces[1][:,1], avg_fit, label = "Average Fit")

    # Calculate and plot Bollinger Bands
    upper_band = avg_fit .+ 2 .* std_values[3] .* abs.(avg_fit)
    lower_band = avg_fit .- 2 .* std_values[3] .* abs.(avg_fit)
    band!(ax, KineticTraces[1][:,1], lower_band, upper_band, color = (:blue, 0.2), label = "Bollinger Bands")

    axislegend()
    display(fig)























































# Convolved bi-exponential with Gaussian IRF (includes infinite constant A_inf)
function convolved_biexponential(t, A1, τ1, A2, τ2, A_inf, t0, sigma)
    conv_exp(A, τ) = @. ((A / 2) * exp(sigma^2 / (2 * τ^2) - (t - t0) / τ) .* erf.((t - t0) / (sqrt(2) * sigma) - sigma / (sqrt(2) * τ)))
    return conv_exp(A1, τ1) + conv_exp(A2, τ2) .+ A_inf
end



# Objective function for optimization
function residuals(p, t, data, sigma)
    A1, τ1, A2, τ2, A_inf, t0 = p
    model = convolved_biexponential(t, A1, τ1, A2, τ2, A_inf, t0, sigma)
    return sum((data .- model) .^ 2)
end

# Calculate reduced chi-squared (χ²) for goodness of fit
function reduced_chi_squared(p, t, data, sigma, num_params, uncertainty)
    model = convolved_biexponential(t, p..., sigma)
    chi_squared = sum(((data .- model) ./ uncertainty) .^ 2)
    dof = length(t) - num_params  # Degrees of freedom
    return chi_squared / dof
end

# Perform the optimization with bounds
function fit_kinetics(t, data, sigma, p0, lower_bounds, upper_bounds)
    obj = p -> residuals(p, t, data, sigma)
    result = optimize(obj, lower_bounds, upper_bounds, p0, Fminbox(BFGS()))
    return result.minimizer, result
end

# Example usage with synthetic data
t = collect(0.0:0.1:100.0)
p_true = [1.0, 20.0, 0.5, 50.0, 0.3, 0.0]
sigma, noise_level = 1.0, 0.05
data = convolved_biexponential(t, p_true..., sigma) .+ noise_level .* randn(length(t))
Plots.plot(t,data)

p0 = [0.8, 10.0, 0.4, 65.0, 0.2, 0.0]  # Initial guess
lower_bounds = [0.0, 1.0, 0.0, 1.0, 0.0, -10.0]  # Lower bounds on the parameters
upper_bounds = [2.0, 100.0, 2.0, 100.0, 1.0, 10.0]  # Upper bounds on the parameters


##### my sample data
Zn_Chlorin_b2_AIRFREE_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_AIRFREE_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)
Zn_Chlorin_b2_AIRFREE_400nm = BiFitIRF_new(Zn_Chlorin_b2_AIRFREE_Kinetics,3,-10.0, 4700, [0.0009, 0.0517, -3.0371, 525.5928, -0.5392, 155.1067, 0.0001], 10,Inf,false,200)

t = Zn_Chlorin_b2_AIRFREE_Kinetics[:,1]
data = Zn_Chlorin_b2_AIRFREE_Kinetics[:,3]
sigma = 0.0517
p0 = [-0.015, 25.5928, -0.001392, 505.1067, 0.01, 0.0009]  # Initial guess
lower_bounds = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]  # Lower bounds on the parameters
upper_bounds = [Inf, Inf, Inf, Inf, Inf, Inf]  # Upper bounds on the parameters


# Fit the data with bounds
fitted_params, result = fit_kinetics(t, data, sigma, p0, lower_bounds, upper_bounds)
fitted_params
fitted_data  = convolved_biexponential(t, fitted_params...,sigma)
# Calculate reduced chi-squared
uncertainty = noise_level  # Assuming constant uncertainty equal to noise level
num_params = length(p0)
chi_sq_reduced = reduced_chi_squared(fitted_params, t, data, sigma, num_params, uncertainty)

# Print the results
Plots.plot(t, data)
Plots.plot!(t,fitted_data)
println("Fitted parameters: ", fitted_params)
println("Reduced Chi-Squared: ", chi_sq_reduced)
