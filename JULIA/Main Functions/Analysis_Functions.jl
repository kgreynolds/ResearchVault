using Interpolations,Plots,Statistics, DataFrames, CSV, Optim,LinearAlgebra, ColorSchemes,MAT,LineSearches,SpecialFunctions, EasyFit,CairoMakie; CairoMakie.activate!()
plotly()
set_theme!(theme_latexfonts())


function GlobalIRF(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing, initial_amplitudes=nothing,bound1=10, bound2=10, num_optimizations=10, noise_level=0.1, size=1000)
    # Constants
    sqrt2_inv = 1 / sqrt(2)  # Precompute for efficiency

    # Resize function to interpolate and resize the kinetic trace to uniform size
    function Resize(x, y, size)
        x_out = collect(range(first(x), last(x), length=size))  # New x-axis with uniform size
        interp_func = LinearInterpolation(x, y)  # Linear interpolation
        y_out = interp_func(x_out)  # Interpolated y-values
        return x_out, y_out
    end

    # Define the tri-exponential, bi-exponential, and mono-exponential models
    @inline function Y_tri(p, x, A₁, A₂, A₃, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5 * (p[2]/p[3])^2 - (x - p[1]) / p[3])) *
        (erf(((x - p[1]) / p[2] - p[2] / p[3]) * sqrt2_inv) + 1) / 2 +
        (A₂ / p[4] * exp(0.5 * (p[2]/p[4])^2 - (x - p[1]) / p[4]) *
        (erf(((x - p[1]) / p[2] - p[2] / p[4]) * sqrt2_inv) + 1) / 2) +
        (A₃ / p[5] * exp(0.5 * (p[2]/p[5])^2 - (x - p[1]) / p[5]) *
        (erf(((x - p[1]) / p[2] - p[2] / p[5]) * sqrt2_inv) + 1) / 2))
    end

    @inline function Y_bi(p, x, A₁, A₂, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5 * (p[2]/p[3])^2 - (x - p[1]) / p[3])) *
        (erf(((x - p[1]) / p[2] - p[2] / p[3]) * sqrt2_inv) + 1) / 2 +
        (A₂ / p[4] * exp(0.5 * (p[2]/p[4])^2 - (x - p[1]) / p[4]) *
        (erf(((x - p[1]) / p[2] - p[2] / p[4]) * sqrt2_inv) + 1) / 2))
    end

    @inline function Y_mono(p, x, A₁, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5 * (p[2]/p[3])^2 - (x - p[1]) / p[3])) *
        (erf(((x - p[1]) / p[2] - p[2] / p[3]) * sqrt2_inv) + 1) / 2)
    end

    # Step 1: Data selection based on the time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)  # Apply the time limits
    t = file[cut, 1]  # Extract the time column

    # Preallocate kinetic traces for each selected column
    KineticTraces = [DataFrame() for _ in columns]

    # Step 2: Resize data for each selected column
    for (i, col) in enumerate(columns)
        trace_data = file[cut, col]  # Extract raw kinetic trace data for each column
        t_resized, trace_resized = Resize(t, trace_data, size)  # Resize the trace
        KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)  # Store resized data
    end

    # Step 3: Initialize guess parameters if not provided
    if shared_params === nothing
        t0_guess = 0.2  # Initial guess for t₀
        FWHM_guess = 1.5  # Initial guess for Full Width at Half Maximum (FWHM)
        τ₁_guess = 0.2 * (maximum(t) - minimum(t))  # Initial guess for first lifetime

        # Add guesses for other lifetimes based on the model type
        if model_type == :bi
            τ₂_guess = 0.8 * (maximum(t) - minimum(t))
            shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess]
        elseif model_type == :tri
            τ₂_guess = 0.6 * (maximum(t) - minimum(t))
            τ₃_guess = 1.0 * (maximum(t) - minimum(t))
            shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess, τ₃_guess]
        else
            shared_params = [t0_guess, FWHM_guess, τ₁_guess]
        end
    end

    # Step 4: Initialize amplitudes if not provided
    if initial_amplitudes === nothing
        if model_type == :bi
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces])
        elseif model_type == :tri
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces])
        else
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces])
        end
    end

    # Step 5: Define global fitting function
    function global_fit(params)
        # Extract shared lifetime parameters and independent amplitudes
        shared_coeffs = if model_type == :bi params[1:4] elseif model_type == :tri params[1:5] else params[1:3] end
        amplitudes = params[if model_type == :bi 5 elseif model_type == :tri 6 else 4 end:end]

        total_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace.x  # Time values
            kinetic_trace = trace.y  # Kinetic trace values

            # Calculate model values based on the selected model type
            model_values = if model_type == :bi
                Y_bi(shared_coeffs, x, amplitudes[3i-2], amplitudes[3i-1], amplitudes[3i])
            elseif model_type == :tri
                Y_tri(shared_coeffs, x, amplitudes[4i-3], amplitudes[4i-2], amplitudes[4i-1], amplitudes[4i])
            else
                Y_mono(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
            end

            # Calculate residuals
            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Step 6: Optimization process with bounds
    initial_params = vcat(shared_params, initial_amplitudes)

    # Define bounds for shared lifetimes and amplitudes
    lower_lifetimes = [(shared_params[i] - abs(shared_params[i]) * (bound1 / 100)) for i in 1:length(shared_params)]
    upper_lifetimes = [(shared_params[i] + abs(shared_params[i]) * (bound1 / 100)) for i in 1:length(shared_params)]

    lower_amplitudes = fill(-Inf, length(initial_amplitudes))
    upper_amplitudes = fill(Inf, length(initial_amplitudes))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    results = zeros(num_optimizations, length(initial_params))

    for i in 1:num_optimizations
        # Add small noise to initial parameters to avoid local minima
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer))
        results[i, :] = Optim.minimizer(result)
    end

    # Step 7: Calculate mean and std of the optimization results
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    # Step 8: Best-fit parameters for shared lifetimes and independent amplitudes
    best_lifetimes = round.(avg_values[1:(model_type == :bi ? 4 : model_type == :tri ? 5 : 3)], digits=3)
    best_amplitudes = round.(avg_values[(model_type == :bi ? 5 : model_type == :tri ? 6 : 4):end], digits=3)

    # Step 9: Compute R² and organize results
    total_ss = 0.0
    total_ss_residual = 0.0
    df_list = []

    final_df = DataFrame(Time = KineticTraces[1].x)  # Start final_df with the time values

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x  # Time
        kinetic_trace = trace.y  # Raw data

        # Generate fit based on the best-fit lifetimes and amplitudes
        fit_values = if model_type == :bi
            Y_bi(best_lifetimes, x, best_amplitudes[3i-2], best_amplitudes[3i-1], best_amplitudes[3i])
        elseif model_type == :tri
            Y_tri(best_lifetimes, x, best_amplitudes[4i-3], best_amplitudes[4i-2], best_amplitudes[4i-1], best_amplitudes[4i])
        else
            Y_mono(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
        end

        residuals = kinetic_trace .- fit_values  # Calculate residuals
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)

        # Add raw data, fitted data, and residuals to the output DataFrame
        final_df[!, "Raw_Trace_$i"] = kinetic_trace
        final_df[!, "Fit_Trace_$i"] = fit_values
        final_df[!, "Residuals_Trace_$i"] = residuals
    end

    # Step 9: Calculate R² and Adjusted R²
    R² = round(1 - total_ss_residual / total_ss, digits=3)
    N = nrow(final_df)  # Number of data points
    p = length(initial_params)  # Number of parameters
    Adjusted_R² = round(1 - (1 - R²) * (N - 1) / (N - p - 1), digits=3)

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x  # Time
        kinetic_trace = trace.y  # Raw data

        # Generate fit based on the best-fit lifetimes and amplitudes
        fit_values = if model_type == :bi
            Y_bi(best_lifetimes, x, best_amplitudes[3i-2], best_amplitudes[3i-1], best_amplitudes[3i])
        elseif model_type == :tri
            Y_tri(best_lifetimes, x, best_amplitudes[4i-3], best_amplitudes[4i-2], best_amplitudes[4i-1], best_amplitudes[4i])
        else
            Y_mono(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
        end

        residuals = kinetic_trace .- fit_values  # Calculate residuals
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)

        # Add raw data, fitted data, and residuals to the output DataFrame
        final_df[!, "Raw_Trace_$i"] = kinetic_trace
        final_df[!, "Fit_Trace_$i"] = fit_values
        final_df[!, "Residuals_Trace_$i"] = residuals

        # Store each trace's data for separate analysis
        trace_df = DataFrame(Time = x, Raw_Data = kinetic_trace, Fit = fit_values, Residuals = residuals)
        push!(df_list, trace_df)
    end

    # Calculate R² (goodness of fit)
    R² = round(1 - total_ss_residual / total_ss, digits=3)

    # Step 10: Prepare parameter labels and output DataFrames
    param_labels = ["t₀", "FWHM", "τ₁"]
    if model_type == :bi
        push!(param_labels, "τ₂")
    elseif model_type == :tri
        push!(param_labels, "τ₂", "τ₃")
    end

    amplitude_labels = []
    amp_per_trace = if model_type == :bi 3 else model_type == :tri ? 4 : 2 end
    for i in 1:length(columns)
        for j in 1:(amp_per_trace - 1)
            push!(amplitude_labels, "A" * string(j) * "_Trace_" * string(i))
        end
        push!(amplitude_labels, "A_inf_Trace_" * string(i))
    end

    param_labels = vcat(param_labels, amplitude_labels)

    param_amp_df = DataFrame(Parameter = param_labels, Value = vcat(best_lifetimes, best_amplitudes), Stdev = std_values)

    # Step 11: Visualization using CairoMakie
    Fit_Fig = Figure(font="", figure_padding=25, fontsize=20)
    width = 3  # Line width

    color_palette = palette(ColorSchemes.hsv, 5)  # Color scheme

    ax1 = CairoMakie.Axis(Fit_Fig[1, 1],
        title = "IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])",
        subtitle = if model_type == :tri
                        "τ₁ = $(best_lifetimes[3]); τ₂ = $(best_lifetimes[4]); τ₃ = $(best_lifetimes[5]); R² = $R²"
                    elseif model_type == :bi
                        "τ₁ = $(best_lifetimes[3]); τ₂ = $(best_lifetimes[4]); R² = $R²"
                    else
                        "τ₁ = $(best_lifetimes[3]); R² = $R²"
                    end,
        palette = (color = color_palette,),
        xlabel = "", xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false, xminorticksvisible = true,
        ylabel = "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false, yminorticksvisible = true)

    # Plot each trace and its fit
    for (i, trace) in enumerate(KineticTraces)
        x = trace.x
        kinetic_trace = trace.y
        fit_values = if model_type == :bi
            Y_bi(best_lifetimes, x, best_amplitudes[3i-2], best_amplitudes[3i-1], best_amplitudes[3i])
        elseif model_type == :tri
            Y_tri(best_lifetimes, x, best_amplitudes[4i-3], best_amplitudes[4i-2], best_amplitudes[4i-1], best_amplitudes[4i])
        else
            Y_mono(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
        end

        residuals = kinetic_trace .- fit_values  # Calculate residuals

        lines!(Fit_Fig[1, 1], x, kinetic_trace, linewidth = width, label = "Trace $i", color = color_palette[i])
        lines!(Fit_Fig[1, 1], x, fit_values, linewidth = width*0.5, linestyle = :dash, color = :black)
    end

    axislegend(ax1, position = :rb, nbanks = 3, framevisible = false, fontsize = 15)

    # Residuals plot
    ax2 = CairoMakie.Axis(Fit_Fig[2, 1], title = "Residuals",
        xlabel = "", xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false, xminorticksvisible = true,
        ylabel = "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false, yminorticksvisible = true)

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x
        fit_values = if model_type == :bi
            Y_bi(best_lifetimes, x, best_amplitudes[3i-2], best_amplitudes[3i-1], best_amplitudes[3i])
        elseif model_type == :tri
            Y_tri(best_lifetimes, x, best_amplitudes[4i-3], best_amplitudes[4i-2], best_amplitudes[4i-1], best_amplitudes[4i])
        else
            Y_mono(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
        end

        residuals = trace.y .- fit_values

        lines!(Fit_Fig[2, 1], x, residuals, linewidth = width, color = :black)
    end

    linkxaxes!(ax1, ax2)

    colsize!(Fit_Fig.layout, 1, Aspect(1, 1.5))
    rowsize!(Fit_Fig.layout, 1, 300)
    rowsize!(Fit_Fig.layout, 2, 130)
    resize_to_layout!(Fit_Fig)
    display(Fit_Fig)
    save("JULIA/Outputs/TriFit_Fig.png", Fit_Fig)

    # Return the final DataFrame with traces, parameters, and R²
    return final_df, param_amp_df, best_lifetimes, best_amplitudes, R², Adjusted_R²
end


function MonoFitIRF_New(file,column,t0,uppertimebound,coeffs, bound1,bound2,guess,size)

    function Resize(file,column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1],file[end,1],size))
        Interp  = LinearInterpolation(file[:,1],file[:,column])
        y_out = Interp(x_out)
    
        ResizedDF = DataFrame(x = x_out,y = y_out)
        return ResizedDF
    end

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< file[:,1] .< uppertimebound
    
    t = file[:,1][cut]
    I = file[:,column][cut]

    Cut_Data = hcat(t,I)

    Resized_Cut_Data = Resize(Cut_Data,2,size)

    x = Resized_Cut_Data[:,1]
    KineticTrace = Resized_Cut_Data[:,2]

    ##coeffs = [xc, w, A1, t1, y0]
    #z1 = (x-xc)/w - w/t1
    #y = y0 + A1/t1 * exp(0.5*(w/t1)^2-(x-xc)/t1)*(erf(z1/sqrt(2))+1)/2

    # Z1(p) = @. ((x-p[1])/p[2] - p[2]/p[4])

    # z1 = Z1(coeffs)
    # Y(p) = @. ((p[5]+ p[3]/p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])) * (erf(z1/sqrt(2))+1)/2)

    Y(p) = @. ((p[5]+ p[3]/p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])) * (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2)


    if guess == true

        Plots.plot(x,KineticTrace)
        Plots.plot!(x,Y(coeffs))
    
    else

        ## FITTING FUNCTION
        peakfit0(coeffs) = sum((KineticTrace .- Y(coeffs)).^2)

        lower = [(coeffs[1]-abs(coeffs[1])*(bound1/100)),(coeffs[2]-abs(coeffs[2])*(bound1/100)),-Inf,(coeffs[4]-abs(coeffs[4]*(bound2/100))),-Inf]#  ,(coeffs[6]-abs(coeffs[6])*(bound1/100)),(coeffs[7]-abs(coeffs[7])*(bound1/100))]
        upper = [(coeffs[1]+abs(coeffs[1])*(bound1/100)),(coeffs[2]+abs(coeffs[2])*(bound1/100)),Inf,(coeffs[4]+abs(coeffs[4]*(bound2/100))),Inf]# ,(coeffs[6]+abs(coeffs[6])*(bound1/100)),(coeffs[7]+abs(coeffs[7])*(bound1/100))]
        
        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        res = optimize(peakfit0, lower, upper, coeffs, Fminbox(inner_optimizer),Optim.Options(time_limit = 60,g_tol = 1e-12))
        
        print("Fit Resuluts:", res)
        values = round.(Optim.minimizer(res),digits = 3) # Kinetics results
        residual = KineticTrace .- Y(values)

        RSS =  sum(((residual)) .^ 2) # residual sum of squares
        TSS = sum((KineticTrace .- mean(KineticTrace)) .^ 2) # Total sum of squares
        R² = round(1 - RSS/TSS,digits = 4)

        print("Fit Coeff:", values)
        print("R²:", R²)

        df = DataFrame(Time = x,RawData = KineticTrace, Fit =Y(values),Residual = residual)


        #######################################################################################################
            Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20)
            width = 3

            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title ="IRF: t₀ = $(values[1]); FWHM = $(values[2])" , subtitle = "A₁ = $(values[3]); τ₁ = $(values[4]); R² = $R²",
            palette = (color = palette(ColorSchemes.hsv, 4),),
            xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)

            lines!(Fit_Fig[1,1],x, KineticTrace,linewidth = width, label = "Data")            
            lines!(Fit_Fig[1,1],x,Y(coeffs), linewidth = width, label = "Guess")
            lines!(Fit_Fig[1,1],x, Y(values), linewidth = width, label = "Fit")

            axislegend(ax1, position = :rt, framevisible = false)

            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

            lines!(Fit_Fig[2,1],x, residual, linewidth = width, color = :black)

            linkxaxes!(ax1,ax2)
            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,130)
            resize_to_layout!(Fit_Fig)

            display(Fit_Fig)
            save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)


        return  df, values
    end
end

function BiFitIRF_new(file,column,t0,uppertimebound,coeffs,bound1, bound2,guess,size)

    function Resize(file,column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1],file[end,1],size))
        Interp  = LinearInterpolation(file[:,1],file[:,column])
        y_out = Interp(x_out)
    
        ResizedDF = DataFrame(x = x_out,y = y_out)
        return ResizedDF
    end

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< file[:,1] .< uppertimebound
    
    t = file[:,1][cut]
    I = file[:,column][cut]

    Cut_Data = hcat(t,I)

    Resized_Cut_Data = Resize(Cut_Data,2,size)

    x = Resized_Cut_Data[:,1]
    KineticTrace = Resized_Cut_Data[:,2]

    ##coeffs = [xc, w, A1, t1, A2, t2, y0]
    #z1 = (x-xc)/w - w/t1
    #z2 = (x-xc)/w - w/t2
    #y = y0 + A1/t1 * exp(0.5*(w/t1)^2-(x-xc)/t1)*(erf(z1/sqrt(2))+1)/2 + A2/t2 * exp(0.5*(w/t2)^2-(x-xc)/t2)*(erf(z2/sqrt(2))+1)/2

    # Z1(p) = @. ((x-p[1])/p[2] - p[2]/p[4])
    # Z2(p) = @. ((x-p[1])/p[2] - p[2]/p[6])
    # z1 = Z1(coeffs)
    # z2 = Z2(coeffs)
    #Y(p) = @. ((p[7]+ p[3]/p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])) * (erf((z1)/sqrt(2))+1)/2 + (p[5]/p[6] * exp(0.5*(p[2]/p[6])^2-(x-p[1])/p[6])*(erf(z2/sqrt(2))+1)/2))

    Y(p) = @. ((p[7]+ p[3]/p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])) * (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2 + (p[5]/p[6] * exp(0.5*(p[2]/p[6])^2-(x-p[1])/p[6])*(erf(((x-p[1])/p[2] - p[2]/p[6])/sqrt(2))+1)/2))

   


    Y(coeffs)

    if guess == true

        Plots.plot(x,KineticTrace)
        Plots.plot!(x,Y(coeffs))
    
    else

        ## FITTING FUNCTION
        peakfit0(coeffs) = sum((KineticTrace .- Y(coeffs)).^2)

        lower = [(coeffs[1]-abs(coeffs[1])*(bound1/100)),(coeffs[2]-abs(coeffs[2])*(bound1/100)),-Inf,(coeffs[4]-abs(coeffs[4]*(bound2/100))),-Inf,(coeffs[6]-abs(coeffs[6]*(bound2/100))),-Inf]
        upper = [(coeffs[1]+abs(coeffs[1])*(bound1/100)),(coeffs[2]+abs(coeffs[2])*(bound1/100)),Inf,(coeffs[4]+abs(coeffs[4]*(bound2/100))),Inf,(coeffs[6]+abs(coeffs[6]*(bound2/100))),Inf]

        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        res = optimize(peakfit0, lower, upper, coeffs, Fminbox(inner_optimizer),Optim.Options(time_limit = 60,g_tol = 1e-12))

        print("Fit Resuluts:", res)
        values =  round.(Optim.minimizer(res),digits = 4) # Kinetics results
        residual = KineticTrace .- Y(values)
        

        RSS =  sum(((residual)) .^ 2)
        TSS = sum((KineticTrace .- mean(KineticTrace)) .^ 2)
        R² = round(1 - RSS/TSS,digits = 4)

        # Additional code for reduced chi-squared
        # n = length(KineticTrace) # number of data points
        # p = length(values) # number of parameters
        # dof = n - p # degrees of freedom
        # reduced_chi_squared = RSS / dof

        print("Fit Coeff:", values)
        print("R²:", R²)
        # println("Reduced chi-squared: ", reduced_chi_squared)

        df = DataFrame(Time = x,RawData = KineticTrace, Fit = Y(values),Residual = residual)

        #######################################################################################################################
            Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20)
            width = 3

            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(values[1]); FWHM = $(values[2])" , subtitle = "τ₁ = $(values[4]); τ₂ = $(values[6]); R² = $R²",
            palette = (color = palette(ColorSchemes.hsv, 4),),
            xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)

            lines!(Fit_Fig[1,1],x, KineticTrace,linewidth = width, label = "Data")
            lines!(Fit_Fig[1,1],x,Y(coeffs), linewidth = width, label = "Guess")
            lines!(Fit_Fig[1,1],x, Y(values), linewidth = width, label = "Fit")

            axislegend(ax1, position = :rt, framevisible = false)

            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

            lines!(Fit_Fig[2,1],x, residual, linewidth = width, color = :black)

            linkxaxes!(ax1,ax2)
            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,130)
            resize_to_layout!(Fit_Fig)

            display(Fit_Fig)
            save("JULIA/Outputs/BiFit_Fig.png",Fit_Fig)
        return df, values
    end
end

function TriFitIRF_new(file,column,t0,uppertimebound,coeffs,bound1,bound2,guess,size)

    function Resize(file,column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1],file[end,1],size))
        Interp  = LinearInterpolation(file[:,1],file[:,column])
        y_out = Interp(x_out)
    
        ResizedDF = DataFrame(x = x_out,y = y_out)
        return ResizedDF
    end

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< file[:,1] .< uppertimebound
    
    t = file[:,1][cut]
    I = file[:,column][cut]

    Cut_Data = hcat(t,I)

    Resized_Cut_Data = Resize(Cut_Data,2,size)

    x = Resized_Cut_Data[:,1]
    KineticTrace = Resized_Cut_Data[:,2]


    ##coeffs = [xc, w, A1, t1, A2, t2, A3, t3 y0]
    #z1 = (x-xc)/w - w/t1
    #z2 = (x-xc)/w - w/t2
    #z3 = (x-xc)/w - w/t3

    #y = y0 + A1/t1 * exp(0.5*(w/t1)^2-(x-xc)/t1)*(erf(z1/sqrt(2))+1)/2 + A2/t2 * exp(0.5*(w/t2)^2-(x-xc)/t2)*(erf(z2/sqrt(2))+1)/2 +A3/t3 * exp(0.5*(w/t3)^2-(x-xc)/t3)*(erf(z3/sqrt(2))+1)/2


    # Z1(p) = @. ((x-p[1])/p[2] - p[2]/p[4])
    # Z2(p) = @. ((x-p[1])/p[2] - p[2]/p[6])
    # Z3(p) = @. ((x-p[1])/p[2] - p[2]/p[8])

    # z1 = Z1(coeffs)
    # z2 = Z2(coeffs)
    # z3 = Z3(coeffs)
    # Y(p) = @. ((p[9]+ p[3]/p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])) * (erf(z1/sqrt(2))+1)/2 + (p[5]/p[6] * exp(0.5*(p[2]/p[6])^2-(x-p[1])/p[6])*(erf(z2/sqrt(2))+1)/2) + (p[7]/p[8] * exp(0.5*(p[2]/p[8])^2-(x-p[1])/p[8])*(erf(z3/sqrt(2))+1)/2))

    Y(p) = @. ((p[9]+ p[3]/p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4])) * (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2 + (p[5]/p[6] * exp(0.5*(p[2]/p[6])^2-(x-p[1])/p[6])*(erf(((x-p[1])/p[2] - p[2]/p[6])/sqrt(2))+1)/2) + (p[7]/p[8] * exp(0.5*(p[2]/p[8])^2-(x-p[1])/p[8])*(erf(((x-p[1])/p[2] - p[2]/p[8])/sqrt(2))+1)/2))
    Y(coeffs)

    if guess == true

        Plots.plot(x,KineticTrace)
        Plots.plot!(x,Y(coeffs))
    
    else

        ## FITTING FUNCTION
        peakfit0(coeffs) = sum((KineticTrace .- Y(coeffs)).^2)

        lower = [(coeffs[1]-abs(coeffs[1])*(bound1/100)),(coeffs[2]-abs(coeffs[2])*(bound1/100)),-Inf,(coeffs[4]-abs(coeffs[4]*(bound2/100))),-Inf,(coeffs[6]-abs(coeffs[6]*(bound2/100))),-Inf,(coeffs[8]-abs(coeffs[8]*(bound2/100))),-Inf ]
        upper = [(coeffs[1]+abs(coeffs[1])*(bound1/100)),(coeffs[2]+abs(coeffs[2])*(bound1/100)),Inf,(coeffs[4]+abs(coeffs[4]*(bound2/100))),Inf,(coeffs[6]+abs(coeffs[6]*(bound2/100))),Inf,(coeffs[8]+abs(coeffs[8]*(bound2/100))),Inf]

        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking(order=3))
        res = optimize(peakfit0, lower, upper, coeffs, Fminbox(inner_optimizer),Optim.Options(time_limit = 200,g_tol = 1e-12))

        print("Fit Resuluts:", res)
        values = round.(Optim.minimizer(res),digits = 3) # Kinetics results
        residual = KineticTrace .- Y(values)


        RSS =  sum(((residual)) .^ 2)
        TSS = sum((KineticTrace .- mean(KineticTrace)) .^ 2)
        R² = round(1 - RSS/TSS,digits = 4)

        print("Fit Coeff:", values)
        print("R²:", R²)

        df = DataFrame(Time = x,RawData = KineticTrace, Fit = Y(values),Residual = residual)

        #######################################################################################################################
            Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20)
            width = 3

            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(values[1]); FWHM = $(values[2])" , subtitle = "τ₁ = $(values[4]); τ₂ = $(values[6]); τ₃ = $(values[8]);R² = $R²",
        
            palette = (color = palette(ColorSchemes.hsv, 4),),
            xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)

            lines!(Fit_Fig[1,1],x, KineticTrace,linewidth = width, label = "Data")
            lines!(Fit_Fig[1,1],x,Y(coeffs), linewidth = width, label = "Guess")
            lines!(Fit_Fig[1,1],x, Y(values), linewidth = width, label = "Fit")

            axislegend(ax1, position = :rt, framevisible = false)

            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

            lines!(Fit_Fig[2,1],x, residual, linewidth = width, color = :black)

            linkxaxes!(ax1,ax2)
            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,130)
            # resize_to_layout!(Fit_Fig)

            display(Fit_Fig)
            save("JULIA/Outputs/TriFit_Fig.png",Fit_Fig)
        return df, values
    end
end

function LinearFit(x,y,residuals, ReturnFit)
    
    time = x
    KineticTrace = y

    ## FITTING FUNCTION
    fitfunction = mexp(p) = p[1] .* x .+ p[2]
    p0 = [-1.1380000000000215e-5, 0.08230333333333385]

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual = KineticTrace .- fitfunction(values)
    print("Fit Coeff:", values)

    if ReturnFit == true
        return hcat(time,KineticTrace,fitfunction(values),residual)
    

    elseif residuals == true
        Fit_Fig = Figure(figure_padding = 25,fontsize =20, resolution = (600,600))
        width = 2

        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
        xlabel = "x",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
        # xlims!(ax1,350,650)
        # ylims!(nothing,)

        lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
        lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
        lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
        axislegend(ax1, position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
        xlabel = "x",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        Fit_Fig
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
    else

        Plots.plot(time, KineticTrace, color = :blue, label = "Data")
        Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
        Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")

    end

    # return fitfunction(values)
end

function smoothMat(Mat,SmoothPts)
    """This function takes and input matrix Mar and generates a new matrix sMat with the moving average over
        some number of points defined by the SmoothPts argument"""
    sMat = zeros(size(Mat,1)-SmoothPts,size(Mat,2))

    for i in 1:size(Mat,1)-SmoothPts
        for j in 1:size(Mat,2)
            sMat[i,j] = mean(skipmissing(Mat[i:i+SmoothPts,j])) # This works! calculates the mean of the i to i+2th row in the jth columns and replaces the E[i,j]th value with that mean
        end
    end
    return sMat #, display(Plots.plot(sMat[:,1], sMat[:,2],xlabel = "Time", ylabel = " ΔA"))
end

function Resize(file,column, size)
    """Resize a column of a dataframe to desired length of Data points """
    x_out = collect(range(file[1,1],file[end,1],size))
    Interp  = LinearInterpolation(file[:,1],file[:,column])
    y_out = Interp(x_out)

    ResizedDF = DataFrame(x = x_out,y = y_out)
    return ResizedDF
end

#Normalization to x=z
function Norm(x,y,z)
    """Normalize Data to a value Z """
    x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
    y_interp = x_interp(x) # compute y values for the x vector
    y_norm = y_interp/x_interp(z) #normalize to desired value
    return y_norm
end

# Normalize to Maximum
function MaxNorm(y)
    """Normalize vector to its absolute largest value"""
    y_norm = y/maximum(abs.(y)) #normalize to desired value
    return y_norm
end

# Interpolate dataframe to different x values (interval must be inside the bounds of the data)
function InterpOnRange(Data::AbstractMatrix, start::Number, stop::Number, interval::Number)::DataFrame
    # Check if start and stop are within the bounds of the data
    if start < minimum(Data[:,1]) || stop > maximum(Data[:,1])
        throw(ArgumentError("Start and stop must be within the bounds of the data"))
    end

    # Generate the x values at which to interpolate
    x_out = collect(start:interval:stop)
    
    # Perform linear interpolation
    Interp  = LinearInterpolation(Data[:,1], Data[:,2])
    y_out = Interp(x_out)

    # Create and return the interpolated DataFrame
    df = DataFrame(x = x_out, y = y_out)
    return df
end

function ReadFilesinDir(directory::String)
    """ Read files in a directory, recurisley, of same file structure and read the data into one matrix"""
    files = []
    for (root, dirs, filenames) in walkdir(directory)
    for filename in filenames
        if isfile(joinpath(root, filename))
        push!(files, readdlm(joinpath(root, filename),skipstart = 60))
        end
    end
    end

    return files
end

function EnVisionKinetic(SmoothPts,path)
    function smoothMat(Mat,SmoothPts)
        """This function takes and input matrix Mar and generates a new matrix sMat with the moving average over
            some number of points defined by the SmoothPts argument"""
        sMat = zeros(size(Mat,1)-SmoothPts,size(Mat,2))
    
        for i in 1:size(Mat,1)-SmoothPts
            for j in 1:size(Mat,2)
                sMat[i,j] = mean(skipmissing(Mat[i:i+SmoothPts,j])) # This works! calculates the mean of the i to i+2th row in the jth columns and replaces the E[i,j]th value with that mean
            end
        end
        return sMat #, display(Plots.plot(sMat[:,1], sMat[:,2],xlabel = "Time", ylabel = " ΔA"))
    end

    matfile = matopen(path)
    newmat = read(matfile,"data")

    KineticTrace = hcat(newmat[:,1],newmat[:,2]*1000) # convert to mOD in ΔA
    Smoothed_Kinetic = smoothMat(KineticTrace,SmoothPts)

    df = DataFrame(Time_ns = Smoothed_Kinetic[:,1],ΔA_mOD = Smoothed_Kinetic[:,2])

    return df
end

function process_CaryUvVis_csv(file_path::String; make_file::Bool=false)
    # Read the CSV file
    df = CSV.read(file_path, DataFrame)

    # Drop the first row and last 42 rows
    df_trimmed = df[2:end-42, 1:end-1]

    # Sort the DataFrame by the first column in ascending order
    df_sorted = sort(df_trimmed, 1)

    # If make_file is true, write the result to a CSV
    if make_file
        output_path = replace(file_path, r"\.csv$" => "_processed.csv")
        CSV.write(output_path, df_sorted)
        println("Processed file saved as: ", output_path)
    end

    # Print the first few rows of the sorted and trimmed dataframe
    # println("First 5 rows of processed data:")
    # println(first(df_sorted, 5))

    return df_sorted
end

function process_Surface_csv(file_path::String; make_file::Bool=false)
    # Read the CSV file
    df = CSV.read(file_path, DataFrame)
    
    # Find all columns with "wavelength, nm" or "t, ns" in the title
    target_cols = findall(x -> occursin("Wavelength, nm", x) || occursin("t, ", x), names(df))
    
    # Keep only the first target column and drop the rest
    if length(target_cols) > 1
        cols_to_drop = target_cols[2:end]
        select!(df, Not(cols_to_drop))
    end
    
    # Remove "Data, " from the remaining column names and replace "ns" with " ns" and "us" with " us"
    new_names = [replace(replace(replace(replace(String(name), "Data, " => ""), "ns" => " ns"), "us" => " µs"), "nm" => " nm") for name in names(df)]
    rename!(df, new_names)

    if make_file
        output_path = replace(file_path, r"\.csv$" => "_Surface_processed.csv")
        CSV.write(output_path, df)
        println("Processed file saved as: ", output_path)
    end
    
    return df
end

