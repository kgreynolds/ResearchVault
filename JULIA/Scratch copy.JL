using Interpolations, Plots, Statistics, DataFrames,ForwardDiff, CSV, Optim, LinearAlgebra, ColorSchemes, MAT, LineSearches, SpecialFunctions,SharedArrays, EasyFit, CairoMakie
CairoMakie.activate!()

Zn_Chlorin_b2_AIRFREE_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_AIRFREE_TA_Kinetics_Surface_processed.csv", header = true, DataFrame)
Zn_Chlorin_b2_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Dilek_Porphyrins/Dilek_Porphyrins/FUELS/Analysis/ZnChlorin_band2_TA_Kinetics.csv", header = true, DataFrame)


function MonoFitIRF_Global_GPTv1(file, columns::Vector{Int}, t0, uppertimebound, shared_params, initial_amplitudes, bound1, bound2, guess, num_optimizations, noise_level, size)
    
    # Resize function to interpolate and resize the kinetic trace
    function Resize(x, y, size)
        """Resize x and y arrays to desired length of Data points"""
        x_out = collect(range(first(x), last(x), length=size))
        interp_func = LinearInterpolation(x, y)
        y_out = interp_func(x_out)
        return x_out, y_out
    end

    # Model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁)
        @. ((p[4] + A₁ / p[3] * exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])) * 
            (erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2)) + 1) / 2)
    end

    # Cut and resize data based on time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)
    t = file[cut, 1]

    # Preallocate kinetic traces
    KineticTraces = [DataFrame() for _ in columns]

    # Resize data for each column
    for (i, col) in enumerate(columns)
        trace_data = file[cut, col]
        t_resized, trace_resized = Resize(t, trace_data, size)
        KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = params[1:4]  # Shared lifetimes
        amplitudes = params[5:end]   # Independent amplitudes
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            model_values = Y(shared_coeffs, x, amplitudes[i])
            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess
        # Plot each trace and initial guess
        plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            plot!(x, kinetic_trace, label="Data Trace $i")
            plot!(x, Y(shared_params, x, initial_amplitudes[i]), label="Initial Guess $i")
        end
        return
    end

    # Set optimization bounds for shared lifetimes and amplitudes
    lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1 / 100)),
                       (shared_params[2] - abs(shared_params[2]) * (bound1 / 100)),
                       (shared_params[3] - abs(shared_params[3]) * (bound2 / 100)), -Inf]
    
    upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1 / 100)),
                       (shared_params[2] + abs(shared_params[2]) * (bound1 / 100)),
                       (shared_params[3] + abs(shared_params[3]) * (bound2 / 100)), Inf]

    lower_amplitudes = fill(-Inf, length(columns))  # Allow amplitudes to vary
    upper_amplitudes = fill(Inf, length(columns))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    # Preallocate results
    results = zeros(num_optimizations, length(initial_params))

    # Optimize the model using multiple initial guesses with noise
    for i in 1:num_optimizations
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))

        # Perform global fit using the optimizer
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(LBFGS()))
        results[i, :] = Optim.minimizer(result)
    end

    # Calculate mean and standard deviation of the optimized parameters
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    # Extract the best-fit parameters
    best_lifetimes = round.(avg_values[1:4], digits=3)
    best_amplitudes = round.(avg_values[5:end], digits=3)

    # Calculate global R² and residuals
    total_ss = 0.0
    total_ss_residual = 0.0

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x
        kinetic_trace = trace.y
        fit_values = Y(best_lifetimes, x, best_amplitudes[i])
        residuals = kinetic_trace .- fit_values
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)
    end

    R² = round(1 - total_ss_residual / total_ss, digits=3)

    # Output the fit coefficients and R²
    println("Best Fit Lifetimes: ", best_lifetimes)
    println("Best Fit Amplitudes: ", best_amplitudes)
    println("Global R²: ", R²)

    # Plot the data and final fit for each trace
    fig = Figure(fontsize=20)
    ax1 = Axis(fig[1, 1], title="IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])", xlabel="Time", ylabel="Intensity")
    ax2 = Axis(fig[2, 1], title="Residuals", xlabel="Time", ylabel="Residuals")

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x
        kinetic_trace = trace.y
        fit_values = Y(best_lifetimes, x, best_amplitudes[i])
        residuals = kinetic_trace .- fit_values

        lines!(ax1, x, kinetic_trace, label="Data $i", linewidth=2)
        lines!(ax1, x, fit_values, label="Fit $i", linewidth=1, linestyle=:dash)
        lines!(ax2, x, residuals, linewidth=2)
    end

    axislegend(ax1)
    display(fig)

    return best_lifetimes, best_amplitudes, R²
end


function MonoExpFitIRF_Global_GPTv2(file, columns::Vector{Int}, t0, uppertimebound, shared_params::Union{Nothing, Vector{Float64}}, initial_amplitudes::Union{Nothing, Vector{Float64}}, bound1, guess, num_optimizations, noise_level, size)

    # Resize function using in-place operations to avoid unnecessary memory allocations
    function Resize(file, column, size)
        x_out = collect(range(file[1, 1], file[end, 1], length=size))
        Interp = LinearInterpolation(file[:, 1], file[:, column])
        y_out = Interp(x_out)
        return DataFrame(x=x_out, y=y_out)
    end

    # Mono-exponential model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁)
        exp1 = @. exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])
        erf1 = @. erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2))
        return @. (p[4] + A₁ / p[3] * exp1) * (erf1 + 1) / 2
    end

    # Cut and resize data in a single step to minimize memory allocations
    t = file[:, 1]
    cut = t0 .< t .< uppertimebound
    t_cut = t[cut]
    KineticTraces = [Resize(hcat(t_cut, file[cut, col]), 2, size) for col in columns]

    # Automatically generate reasonable guesses if not provided
    if isnothing(shared_params)
        t0_guess = 0.01 # Guess the time shift 
        FWHM_guess = 1.5 # Guess FWHM 
        τ₁_guess = 0.2 * (maximum(t_cut) - minimum(t_cut)) # Medium lifetime guess
        A_inf_guess = 0.0 # Guess asymptotic amplitude as 0 (can be adjusted)
        shared_params = [t0_guess, FWHM_guess, τ₁_guess, A_inf_guess]
        println("Initialized shared_params with: $shared_params")
    end

    if isnothing(initial_amplitudes)
        initial_amplitudes = []
        for trace in KineticTraces
            amplitude_guess = 100 *(maximum(trace[:, 2]) + minimum(trace[:, 2])) # Difference between max and min for each trace
            push!(initial_amplitudes, amplitude_guess) # Single amplitude guess per trace
        end
        println("Initialized initial_amplitudes with: $initial_amplitudes")
    end

    # Global fitting function with vectorized operations
    function global_fit(params)
        shared_coeffs = params[1:4] # Shared coefficients (t0, FWHM, τ₁, A_inf)
        amplitudes = params[5:end] # Independent amplitudes for each trace
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            model_values = Y(shared_coeffs, x, amplitudes[i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess
        # Plot initial guess for each trace
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            display(Plots.plot!(x, KineticTrace, label="Data Trace $i"))
            display(Plots.plot!(x, Y(shared_params, x, initial_amplitudes[i]), label="Initial Guess"))
        end
    else
        # Optimization bounds
        lower_lifetimes = shared_params .* (1 .- bound1 / 100)
        upper_lifetimes = shared_params .* (1 .+ bound1 / 100)
        lower_amplitudes = fill(-Inf, length(initial_amplitudes))
        upper_amplitudes = fill(Inf, length(initial_amplitudes))

        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        # Use a regular array for results storage (replacing SharedArray)
        results = Array{Float64}(undef, num_optimizations, length(initial_params))
        Threads.@threads for i in 1:num_optimizations
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer), Optim.Options(time_limit=60, g_tol=1e-12))
            results[i, :] = Optim.minimizer(result)
        end

        # Average and standard deviation of results
        avg_values = round.(mean(results, dims=1), digits=3)
        std_values = round.(std(results, dims=1), digits=3)

        best_lifetimes = avg_values[1:4]
        best_amplitudes = avg_values[5:end]

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = round(1 - total_ss_residual / total_ss, digits=3)

        println("Best Fit Lifetimes: ", best_lifetimes)
        println("Best Fit Amplitudes: ", best_amplitudes)
        println("Stdev: ", std_values)
        println("Global R²: ", R²)

        # Plot data and fitted values
        Fit_Fig = Figure()
        width = 3
        ax1 = Axis(Fit_Fig[1, 1], title="IRF Fit", xlabel="Time")
        ax2 = Axis(Fit_Fig[2, 1], title="Residuals", xlabel="Time")

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[i])
            residuals = KineticTrace .- fit_values

            lines!(ax1, x, KineticTrace, linewidth=width, label="Data $i")
            lines!(ax1, x, fit_values, linewidth=width/2, linestyle=:dash, color=:black)
            lines!(ax2, x, residuals, linewidth=width)

            linkxaxes!(ax1, ax2)
        end
        axislegend(ax1, position=:rb)
        display(Fit_Fig)
    end
end

function BiExpFitIRF_Global_GPTv1(file, columns::Vector{Int}, t0, uppertimebound, shared_params, initial_amplitudes, bound1, bound2, guess, num_optimizations, noise_level, size)

    # Resize function using in-place operations to avoid unnecessary memory allocations
    function Resize(file, column, size)
        x_out = collect(range(file[1, 1], file[end, 1], length=size))
        Interp = LinearInterpolation(file[:, 1], file[:, column])
        y_out = Interp(x_out)
        return DataFrame(x=x_out, y=y_out) # directly return DataFrame without intermediate steps
    end

    # Bi-exponential model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁, A₂)
        exp1 = @. exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])
        exp2 = @. exp(0.5 * (p[2] / p[4])^2 - (x - p[1]) / p[4])
        erf1 = @. erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2))
        erf2 = @. erf(((x - p[1]) / p[2] - p[2] / p[4]) / sqrt(2))
        return @. (p[5] + A₁ / p[3] * exp1) * (erf1 + 1) / 2 + (A₂ / p[4] * exp2) * (erf2 + 1) / 2
    end

    # Cut and resize data in a single step to minimize memory allocations
    t = file[:, 1]
    cut = t0 .< t .< uppertimebound
    t_cut = t[cut]
    KineticTraces = [Resize(hcat(t_cut, file[cut, col]), 2, size) for col in columns]

    # Global fitting function with vectorized operations
    function global_fit(params)
        shared_coeffs = params[1:5]
        amplitudes = params[6:end]
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            model_values = Y(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
            residuals = KineticTrace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess
        # Plot initial guess for each trace
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            display(Plots.plot!(x, KineticTrace, label="Data Trace $i"))
            display(Plots.plot!(x, Y(shared_params, x, initial_amplitudes[2i-1], initial_amplitudes[2i]), label="Initial Guess"))
        end
    else
        # Optimization bounds
        lower_lifetimes = shared_params .* (1 .- bound1 / 100)
        upper_lifetimes = shared_params .* (1 .+ bound1 / 100)
        lower_amplitudes = fill(-Inf, length(initial_amplitudes))
        upper_amplitudes = fill(Inf, length(initial_amplitudes))
        
        lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
        upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

        # Parallel optimization to speed up multiple runs
        results = SharedArray{Float64}(num_optimizations, length(initial_params))
        Threads.@threads for i in 1:num_optimizations
            noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
            inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
            result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer), Optim.Options(time_limit=60, g_tol=1e-12))
            results[i, :] = Optim.minimizer(result)
        end

        # Average and standard deviation of results
        avg_values = round.(mean(results, dims=1), digits=3)
        std_values = round.(std(results, dims=1), digits=3)

        best_lifetimes = avg_values[1:5]
        best_amplitudes = avg_values[6:end]

        # Calculate global R² and residuals
        total_ss = 0.0
        total_ss_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
            residuals = KineticTrace .- fit_values
            total_ss_residual += sum(residuals.^2)
            total_ss += sum((KineticTrace .- mean(KineticTrace)).^2)
        end
        R² = round(1 - total_ss_residual / total_ss, digits=3)

        println("Best Fit Lifetimes: ", best_lifetimes)
        println("Best Fit Amplitudes: ", best_amplitudes)
        println("Stdev: ", std_values)
        println("Global R²: ", R²)

        # Plot data and fitted values
        Fit_Fig = Figure()
        width = 3
        ax1 = Axis(Fit_Fig[1, 1], title="IRF Fit", xlabel="Time")
        ax2 = Axis(Fit_Fig[2, 1], title="Residuals", xlabel="Time")

        for (i, trace) in enumerate(KineticTraces)
            x = trace[:, 1]
            KineticTrace = trace[:, 2]
            fit_values = Y(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
            residuals = KineticTrace .- fit_values

            lines!(ax1, x, KineticTrace, linewidth=width, label="Data $i")
            lines!(ax1, x, fit_values, linewidth=width/2, linestyle=:dash, color=:black)
            lines!(ax2, x, residuals, linewidth=width)

            linkxaxes!(ax1, ax2)
        end
        axislegend(ax1, position=:rb)
        display(Fit_Fig)
    end
end

function GlobalIRF_Fit(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing,initial_amplitudes=nothing, bound1=10, bound2=10, guess=false, num_optimizations=10, noise_level=0.1, size=1000)

    # Resize function to interpolate and resize the kinetic trace
    function Resize(x, y, size)
    """Resize x and y arrays to desired length of Data points"""
    x_out = collect(range(first(x), last(x), length=size))
    interp_func = LinearInterpolation(x, y)
    y_out = interp_func(x_out)
    return x_out, y_out
    end

    # Model functions for bi-exponential and mono-exponential models
    function Y_bi(p, x, A₁, A₂)
    @. ((p[5] + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
    (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
    (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) * 
    (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))
    end

    function Y_mono(p, x, A₁)
    @. ((p[4] + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
    (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2)
    end

    # Cut and resize data based on time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)
    t = file[cut, 1]

    # Preallocate kinetic traces
    KineticTraces = [DataFrame() for _ in columns]

    # Resize data for each column
    for (i, col) in enumerate(columns)
    trace_data = file[cut, col]
    t_resized, trace_resized = Resize(t, trace_data, size)
    KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)
    end

    # Initialize guess parameters if not provided
    if shared_params === nothing
    t0_guess = 0.01 # Guess the time shift 
    FWHM_guess = 1.5 # Guess FWHM 
    τ₁_guess = 0.2 * (maximum(t) - minimum(t)) # Medium lifetime guess
    A_inf_guess = 0.0 # Guess asymptotic amplitude as 0 (can be adjusted)
    if model_type == :bi
    τ₂_guess = 0.8 * (maximum(t) - minimum(t)) # Longer lifetime guess (for bi-exponential)
    shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess, A_inf_guess]
    else
    shared_params = [t0_guess, FWHM_guess, τ₁_guess, A_inf_guess]
    end
    println("Initialized shared_params with: $shared_params")
    end

    if initial_amplitudes === nothing
    # Default initial amplitudes: mean of each trace
    if model_type == :bi
    initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces]) # Two amplitudes per trace (for bi-exponential)
    else
    initial_amplitudes = [mean(trace.y) for trace in KineticTraces] # One amplitude per trace (for mono-exponential)
    end
    println("Auto-initialized initial_amplitudes: ", initial_amplitudes)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
    shared_coeffs = if model_type == :bi params[1:5] else params[1:4] end  # Shared lifetimes
    amplitudes = params[if model_type == :bi 6 else 5 end:end]   # Independent amplitudes
    total_residual = 0.0

    for (i, trace) in enumerate(KineticTraces)
    x = trace.x
    kinetic_trace = trace.y
    model_values = if model_type == :bi
    Y_bi(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
    else
    Y_mono(shared_coeffs, x, amplitudes[i])
    end
    residuals = kinetic_trace .- model_values
    total_residual += sum(residuals.^2)
    end
    return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess
    # Plot each trace and initial guess
    Plots.plot()
    for (i, trace) in enumerate(KineticTraces)
    x = trace.x
    kinetic_trace = trace.y
    Plots.plot!(x, kinetic_trace, label="Data Trace $i")
    if model_type == :bi
    Plots.plot!(x, Y_bi(shared_params, x, initial_amplitudes[2i-1], initial_amplitudes[2i]), label="Initial Guess $i")
    else
    Plots.plot!(x, Y_mono(shared_params, x, initial_amplitudes[i]), label="Initial Guess $i")
    end
    end
    return
    end

    # Set optimization bounds for shared lifetimes and amplitudes
    lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] - abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] - abs(shared_params[3]) * (bound2 / 100)), -Inf]

    upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] + abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] + abs(shared_params[3]) * (bound2 / 100)), Inf]

    if model_type == :bi
    push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
    push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
    end

    lower_amplitudes = fill(-Inf, length(columns) * (model_type == :bi ? 2 : 1))  # Allow amplitudes to vary
    upper_amplitudes = fill(Inf, length(columns) * (model_type == :bi ? 2 : 1))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    # Preallocate results
    results = zeros(num_optimizations, length(initial_params))

    # Optimize the model using multiple initial guesses with noise
    for i in 1:num_optimizations
    noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))

    # Perform global fit using the optimizer
    inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
    result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer))
    results[i, :] = Optim.minimizer(result)
    end

    # Calculate mean and standard deviation of the optimized parameters
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    # Extract the best-fit parameters
    best_lifetimes = round.(avg_values[1:(model_type == :bi ? 5 : 4)], digits=3)
    best_amplitudes = round.(avg_values[(model_type == :bi ? 6 : 5):end], digits=3)

    # Calculate global R² and residuals
    total_ss = 0.0
    total_ss_residual = 0.0

    for (i, trace) in enumerate(KineticTraces)
    x = trace.x
    kinetic_trace = trace.y
    fit_values = if model_type == :bi
    Y_bi(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
    else
    Y_mono(best_lifetimes, x, best_amplitudes[i])
    end
    residuals = kinetic_trace .- fit_values
    total_ss_residual += sum(residuals.^2)
    total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)
    end

    R² = round(1 - total_ss_residual / total_ss, digits=3)

    # Output the fit coefficients and R²
    println("Best Fit Lifetimes: ", best_lifetimes)
    println("Best Fit Amplitudes: ", best_amplitudes)
    println("Global R²: ", R²)

    # Plot the data and final fit for each trace
    fig = Figure(fontsize=20)
    ax1 = Axis(fig[1, 1], title="IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])", xlabel="Time", ylabel="Intensity")
    ax2 = Axis(fig[2, 1], title="Residuals", xlabel="Time", ylabel="Residuals")

    for (i, trace) in enumerate(KineticTraces)
    x = trace.x
    kinetic_trace = trace.y
    fit_values = if model_type == :bi
    Y_bi(best_lifetimes, x, best_amplitudes[2i-1], best_amplitudes[2i])
    else
    Y_mono(best_lifetimes, x, best_amplitudes[i])
    end
    residuals = kinetic_trace .- fit_values

    lines!(ax1, x, kinetic_trace, label="Data $i", linewidth=2)
    lines!(ax1, x, fit_values, label="Fit $i", linewidth=1, linestyle=:dash)
    lines!(ax2, x, residuals, linewidth=2)
    end

    axislegend(ax1)
    display(fig)

    return best_lifetimes, best_amplitudes, R²
end


function GlobalIRF_AllFit_v1(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing, initial_amplitudes=nothing, bound1=10, bound2=10, guess=false, num_optimizations=10, noise_level=0.1, size=1000)

    # Resize function to interpolate and resize the kinetic trace
    function Resize(x, y, size)
        """Resize x and y arrays to desired length of Data points"""
        x_out = collect(range(first(x), last(x), length=size))
        interp_func = LinearInterpolation(x, y)
        y_out = interp_func(x_out)
        return x_out, y_out
    end

    # Model functions for tri-exponential, bi-exponential, and mono-exponential models
    function Y_tri(p, x, A₁, A₂, A₃, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
        (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
        (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) * 
        (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2) +
        (A₃ / p[5] * exp(0.5*(p[2]/p[5])^2-(x-p[1])/p[5]) * 
        (erf(((x-p[1])/p[2] - p[2]/p[5])/sqrt(2))+1)/2))
    end

    function Y_bi(p, x, A₁, A₂, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
        (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
        (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) * 
        (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))
    end

    function Y_mono(p, x, A₁, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
        (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2)
    end

    # Cut and resize data based on time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)
    t = file[cut, 1]

    # Preallocate kinetic traces
    KineticTraces = [DataFrame() for _ in columns]

    # Resize data for each column
    for (i, col) in enumerate(columns)
        trace_data = file[cut, col]
        t_resized, trace_resized = Resize(t, trace_data, size)
        KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)
    end

    # Initialize guess parameters if not provided
    if shared_params === nothing
        t0_guess = 0.01 # Guess the time shift 
        FWHM_guess = 1.0 # Guess FWHM 
        τ₁_guess = 0.2 * (maximum(t) - minimum(t)) # Medium lifetime guess
        if model_type == :bi
            τ₂_guess = 0.8 * (maximum(t) - minimum(t)) # Longer lifetime guess (for bi-exponential)
            shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess]
        elseif model_type == :tri
            τ₂_guess = 0.6 * (maximum(t) - minimum(t)) # Medium-long lifetime guess
            τ₃_guess = 1.0 * (maximum(t) - minimum(t)) # Longest lifetime guess
            shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess, τ₃_guess]
        else
            shared_params = [t0_guess, FWHM_guess, τ₁_guess]
        end
        println("Initialized shared_params with: $shared_params")
    end

    if initial_amplitudes === nothing
        # Default initial amplitudes: mean of each trace and initial guesses for infinite offsets
        if model_type == :bi
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces]) # Two amplitudes and infinite offset per trace
        elseif model_type == :tri
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces]) # Three amplitudes and infinite offset per trace
        else
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces]) # One amplitude and infinite offset per trace
        end
        println("Auto-initialized initial_amplitudes: ", initial_amplitudes)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = if model_type == :bi params[1:4] elseif model_type == :tri params[1:5] else params[1:3] end  # Shared lifetimes
        amplitudes = params[if model_type == :bi 5 elseif model_type == :tri 6 else 4 end:end]   # Independent amplitudes
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            model_values = if model_type == :bi
                Y_bi(shared_coeffs, x, amplitudes[3i-2], amplitudes[3i-1], amplitudes[3i])
            elseif model_type == :tri
                Y_tri(shared_coeffs, x, amplitudes[4i-3], amplitudes[4i-2], amplitudes[4i-1], amplitudes[4i])
            else
                Y_mono(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
            end
            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes (including infinite offsets)
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess
        # Plot each trace and initial guess
        Plots.plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            Plots.plot!(x, kinetic_trace, label="Data Trace $i")
            if model_type == :bi
                Plots.plot!(x, Y_bi(shared_params, x, initial_amplitudes[3i-2], initial_amplitudes[3i-1], initial_amplitudes[3i]), label="Initial Guess $i")
            elseif model_type == :tri
                Plots.plot!(x, Y_tri(shared_params, x, initial_amplitudes[4i-3], initial_amplitudes[4i-2], initial_amplitudes[4i-1], initial_amplitudes[4i]), label="Initial Guess $i")
            else
                Plots.plot!(x, Y_mono(shared_params, x, initial_amplitudes[2i-1], initial_amplitudes[2i]), label="Initial Guess $i")
            end
        end
        return
    end

    # Set optimization bounds for shared lifetimes and amplitudes (amplitudes and infinite offset are unbounded)
    lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] - abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] - abs(shared_params[3]) * (bound2 / 100))]

    upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] + abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] + abs(shared_params[3]) * (bound2 / 100))]

    if model_type == :bi
        push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
    elseif model_type == :tri
        push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
        push!(lower_lifetimes, (shared_params[5] - abs(shared_params[5]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[5] + abs(shared_params[5]) * (bound2 / 100)))
    end

    lower_amplitudes = fill(-Inf, length(initial_amplitudes))  # Amplitudes and infinite offsets are unbounded
    upper_amplitudes = fill(Inf, length(initial_amplitudes))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    # Preallocate results
    results = zeros(num_optimizations, length(initial_params))

    # Optimize the model using multiple initial guesses with noise
    for i in 1:num_optimizations
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))

        # Perform global fit using the optimizer
        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer))
        results[i, :] = Optim.minimizer(result)
    end

    # Calculate mean and standard deviation of the optimized parameters
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    # Extract the best-fit parameters
    best_lifetimes = round.(avg_values[1:(model_type == :bi ? 4 : model_type == :tri ? 5 : 3)], digits=3)
    best_amplitudes = round.(avg_values[(model_type == :bi ? 5 : model_type == :tri ? 6 : 4):end], digits=3)

    # Calculate global R² and residuals
    total_ss = 0.0
    total_ss_residual = 0.0

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
        residuals = kinetic_trace .- fit_values
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)
    end

    R² = round(1 - total_ss_residual / total_ss, digits=3)

    # Output the fit coefficients and R²
    println("Best Fit Lifetimes: ", best_lifetimes)
    println("Best Fit Amplitudes: ", best_amplitudes)
    println("Global R²: ", R²)

    # Create a figure with padding and specific font size
        Fit_Fig = Figure(font = "", figure_padding = 25, fontsize = 20)
        width = 3  # Line width for plots

        # Define the color palette (using 4 colors from the HSV color scheme)
        color_palette = palette(ColorSchemes.hsv, 5)

        # Set up the first axis for the fit data
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

        # Plot the kinetic traces and fit values
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

            residuals = kinetic_trace .- fit_values

            # Plot the data trace
            lines!(Fit_Fig[1, 1], x, kinetic_trace, linewidth = width, label = "Trace $i", color = color_palette[i])
            
            # Plot the fit trace
            lines!(Fit_Fig[1, 1], x, fit_values, linewidth = width*0.5, linestyle = :dash, color = :black)
        end

        # Add a legend to the axis
        axislegend(ax1, position = :rb,nbanks = 3,  framevisible = false, fontsize = 15)

        # Set up the second axis for the residuals
        ax2 = CairoMakie.Axis(Fit_Fig[2, 1], title = "Residuals",
            xlabel = "", xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false, xminorticksvisible = true,
            ylabel = "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false, yminorticksvisible = true)

        # Plot the residuals for each trace
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

            # Plot residuals
            lines!(Fit_Fig[2, 1], x, residuals, linewidth = width, color = :black)
        end

        # Link the x-axes of the two subplots (fit and residuals)
        linkxaxes!(ax1, ax2)

        # Set aspect ratio and size of rows and columns
        colsize!(Fit_Fig.layout, 1, Aspect(1, 1.5))  # Aspect ratio of y vs. x
        rowsize!(Fit_Fig.layout, 1, 300)  # Size for the fit plot row
        rowsize!(Fit_Fig.layout, 2, 130)  # Size for the residuals row

        # Display and save the figure
        display(Fit_Fig)
        save("JULIA/Outputs/TriFit_Fig.png", Fit_Fig)

        

    return best_lifetimes, best_amplitudes, R²
end

function GlobalIRF_AllFit_v2(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing, initial_amplitudes=nothing, bound1=10, bound2=10, num_optimizations=10, noise_level=0.1, size=1000)

    # Resize function to interpolate and resize the kinetic trace
    function Resize(x, y, size)
        """Resize x and y arrays to desired length of Data points"""
        x_out = collect(range(first(x), last(x), length=size))
        interp_func = LinearInterpolation(x, y)
        y_out = interp_func(x_out)
        return x_out, y_out
    end

    # Model functions for tri-exponential, bi-exponential, and mono-exponential models
    function Y_tri(p, x, A₁, A₂, A₃, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
        (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
        (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) * 
        (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2) +
        (A₃ / p[5] * exp(0.5*(p[2]/p[5])^2-(x-p[1])/p[5]) * 
        (erf(((x-p[1])/p[2] - p[2]/p[5])/sqrt(2))+1)/2))
    end

    function Y_bi(p, x, A₁, A₂, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
        (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2 +
        (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) * 
        (erf(((x-p[1])/p[2] - p[2]/p[4])/sqrt(2))+1)/2))
    end

    function Y_mono(p, x, A₁, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) * 
        (erf(((x-p[1])/p[2] - p[2]/p[3])/sqrt(2))+1)/2)
    end

    # Cut and resize data based on time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)
    t = file[cut, 1]

    # Preallocate kinetic traces
    KineticTraces = [DataFrame() for _ in columns]

    # Resize data for each column
    for (i, col) in enumerate(columns)
        trace_data = file[cut, col]
        t_resized, trace_resized = Resize(t, trace_data, size)
        KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)
    end

    # Initialize guess parameters if not provided
    if shared_params === nothing
        t0_guess = 0.01 # Guess the time shift 
        FWHM_guess = 1.5 # Guess FWHM 
        τ₁_guess = 0.2 * (maximum(t) - minimum(t)) # Medium lifetime guess
        if model_type == :bi
            τ₂_guess = 0.8 * (maximum(t) - minimum(t)) # Longer lifetime guess (for bi-exponential)
            shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess]
        elseif model_type == :tri
            τ₂_guess = 0.6 * (maximum(t) - minimum(t)) # Medium-long lifetime guess
            τ₃_guess = 1.0 * (maximum(t) - minimum(t)) # Longest lifetime guess
            shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess, τ₃_guess]
        else
            shared_params = [t0_guess, FWHM_guess, τ₁_guess]
        end
    end

    if initial_amplitudes === nothing
        # Default initial amplitudes: mean of each trace and initial guesses for infinite offsets
        if model_type == :bi
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces]) # Two amplitudes and infinite offset per trace
        elseif model_type == :tri
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces]) # Three amplitudes and infinite offset per trace
        else
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces]) # One amplitude and infinite offset per trace
        end
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = if model_type == :bi params[1:4] elseif model_type == :tri params[1:5] else params[1:3] end  # Shared lifetimes
        amplitudes = params[if model_type == :bi 5 elseif model_type == :tri 6 else 4 end:end]   # Independent amplitudes
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            model_values = if model_type == :bi
                Y_bi(shared_coeffs, x, amplitudes[3i-2], amplitudes[3i-1], amplitudes[3i])
            elseif model_type == :tri
                Y_tri(shared_coeffs, x, amplitudes[4i-3], amplitudes[4i-2], amplitudes[4i-1], amplitudes[4i])
            else
                Y_mono(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
            end
            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes (including infinite offsets)
    initial_params = vcat(shared_params, initial_amplitudes)

    # Set optimization bounds for shared lifetimes and amplitudes (amplitudes and infinite offset are unbounded)
    lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] - abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] - abs(shared_params[3]) * (bound2 / 100))]

    upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] + abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] + abs(shared_params[3]) * (bound2 / 100))]

    if model_type == :bi
        push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
    elseif model_type == :tri
        push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
        push!(lower_lifetimes, (shared_params[5] - abs(shared_params[5]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[5] + abs(shared_params[5]) * (bound2 / 100)))
    end

    lower_amplitudes = fill(-Inf, length(initial_amplitudes))  # Amplitudes and infinite offsets are unbounded
    upper_amplitudes = fill(Inf, length(initial_amplitudes))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    # Preallocate results
    results = zeros(num_optimizations, length(initial_params))

    # Optimize the model using multiple initial guesses with noise
    for i in 1:num_optimizations
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))

        # Perform global fit using the optimizer
        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer))
        results[i, :] = Optim.minimizer(result)
    end

    # Calculate mean and standard deviation of the optimized parameters
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    # Extract the best-fit parameters
    best_lifetimes = round.(avg_values[1:(model_type == :bi ? 4 : model_type == :tri ? 5 : 3)], digits=3)
    best_amplitudes = round.(avg_values[(model_type == :bi ? 5 : model_type == :tri ? 6 : 4):end], digits=3)

    # Calculate global R² and residuals
    total_ss = 0.0
    total_ss_residual = 0.0

    # Data for DataFrame collection
    df_list = []
    
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
        residuals = kinetic_trace .- fit_values
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)

        # Add trace data to the DataFrame
        trace_df = DataFrame(
            Time = x,
            Raw_Data = kinetic_trace,
            Fit = fit_values,
            Residuals = residuals
        )

        push!(df_list, trace_df)
    end

    R² = round(1 - total_ss_residual / total_ss, digits=3)

    # Create a combined DataFrame of parameters and amplitudes
    param_labels = ["t₀", "FWHM", "τ₁"]
    if model_type == :bi
        push!(param_labels, "τ₂")
    elseif model_type == :tri
        push!(param_labels, "τ₂", "τ₃")
    end

    # Add amplitude labels
    amplitude_labels = []
    amp_per_trace = if model_type == :bi 3 else model_type == :tri ? 4 : 2 end
    for i in 1:length(columns)
        for j in 1:(amp_per_trace - 1)  # Amplitudes per trace except the infinite offset
            push!(amplitude_labels, "A" * string(j) * "_Trace_" * string(i))
        end
        push!(amplitude_labels, "A_inf_Trace_" * string(i))
    end

    param_labels = vcat(param_labels, amplitude_labels)

    # Combine lifetimes and amplitudes into one DataFrame
    param_amp_df = DataFrame(
        Parameter = param_labels,
        Value = vcat(best_lifetimes, best_amplitudes)
    )

    # Concatenate all trace DataFrames into one
    final_df = vcat(df_list...)

    ## PLOTTING USING MAKIE ##

    Fit_Fig = Figure(font="", figure_padding=25, fontsize=20)
    width = 3  # Line width for plots

    # Define the color palette (using HSV color scheme)
    color_palette = palette(ColorSchemes.hsv, 5)

    # Set up the first axis for the fit data
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

    # Plot the kinetic traces and fit values
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

        residuals = kinetic_trace .- fit_values

        # Plot the data trace
        lines!(Fit_Fig[1, 1], x, kinetic_trace, linewidth = width, label = "Trace $i", color = color_palette[i])
        
        # Plot the fit trace
        lines!(Fit_Fig[1, 1], x, fit_values, linewidth = width*0.5, linestyle = :dash, color = :black)
    end

    # Add a legend to the axis
    axislegend(ax1, position = :rb, nbanks = 3, framevisible = false, fontsize = 15)

    # Set up the second axis for the residuals
    ax2 = CairoMakie.Axis(Fit_Fig[2, 1], title = "Residuals",
        xlabel = "", xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false, xminorticksvisible = true,
        ylabel = "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false, yminorticksvisible = true)

    # Plot the residuals for each trace
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

        # Plot residuals
        lines!(Fit_Fig[2, 1], x, residuals, linewidth = width, color = :black)
    end

    # Link the x-axes of the two subplots (fit and residuals)
    linkxaxes!(ax1, ax2)

    # Set aspect ratio and size of rows and columns
    colsize!(Fit_Fig.layout, 1, Aspect(1, 1.5))  # Aspect ratio of y vs. x
    rowsize!(Fit_Fig.layout, 1, 300)  # Size for the fit plot row
    rowsize!(Fit_Fig.layout, 2, 130)  # Size for the residuals row

    # Display and save the figure
    display(Fit_Fig)
    save("JULIA/Outputs/TriFit_Fig.png", Fit_Fig)


    # Return the final DataFrame with data, fit, residuals, and combined parameters and amplitudes
    return final_df, param_amp_df, best_lifetimes, best_amplitudes, R²
end

function GlobalIRF_AllFit_v2_optimized(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing, initial_amplitudes=nothing, bound1=10, bound2=10, num_optimizations=10, noise_level=0.1, size=1000)

    # Precompute constants
    sqrt2_inv = 1 / sqrt(2)

    # Resize function to interpolate and resize the kinetic trace
    function Resize(x, y, size)
        x_out = collect(range(first(x), last(x), length=size))
        interp_func = LinearInterpolation(x, y)
        y_out = interp_func(x_out)
        return x_out, y_out
    end

    # Model functions for tri-exponential, bi-exponential, and mono-exponential models
    @inline function Y_tri(p, x, A₁, A₂, A₃, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) *
        (erf(((x-p[1])/p[2] - p[2]/p[3])*sqrt2_inv)+1)/2 +
        (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) *
        (erf(((x-p[1])/p[2] - p[2]/p[4])*sqrt2_inv)+1)/2) +
        (A₃ / p[5] * exp(0.5*(p[2]/p[5])^2-(x-p[1])/p[5]) *
        (erf(((x-p[1])/p[2] - p[2]/p[5])*sqrt2_inv)+1)/2))
    end

    @inline function Y_bi(p, x, A₁, A₂, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) *
        (erf(((x-p[1])/p[2] - p[2]/p[3])*sqrt2_inv)+1)/2 +
        (A₂ / p[4] * exp(0.5*(p[2]/p[4])^2-(x-p[1])/p[4]) *
        (erf(((x-p[1])/p[2] - p[2]/p[4])*sqrt2_inv)+1)/2))
    end

    @inline function Y_mono(p, x, A₁, A_inf)
        @. ((A_inf + A₁ / p[3] * exp(0.5*(p[2]/p[3])^2-(x-p[1])/p[3])) *
        (erf(((x-p[1])/p[2] - p[2]/p[3])*sqrt2_inv)+1)/2)
    end

    # Cut and resize data based on time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)
    t = file[cut, 1]

    # Preallocate kinetic traces
    KineticTraces = [DataFrame() for _ in columns]

    # Resize data for each column
    for (i, col) in enumerate(columns)
        trace_data = file[cut, col]
        t_resized, trace_resized = Resize(t, trace_data, size)
        KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)
    end

    # Initialize guess parameters if not provided
    if shared_params === nothing
        t0_guess = 0.01
        FWHM_guess = 1.5
        τ₁_guess = 0.2 * (maximum(t) - minimum(t))

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

    if initial_amplitudes === nothing
        if model_type == :bi
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces])
        elseif model_type == :tri
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces])
        else
            initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [0.0 for trace in KineticTraces])
        end
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = if model_type == :bi params[1:4] elseif model_type == :tri params[1:5] else params[1:3] end
        amplitudes = params[if model_type == :bi 5 elseif model_type == :tri 6 else 4 end:end]

        total_residual = 0.0
        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            model_values = if model_type == :bi
                Y_bi(shared_coeffs, x, amplitudes[3i-2], amplitudes[3i-1], amplitudes[3i])
            elseif model_type == :tri
                Y_tri(shared_coeffs, x, amplitudes[4i-3], amplitudes[4i-2], amplitudes[4i-1], amplitudes[4i])
            else
                Y_mono(shared_coeffs, x, amplitudes[2i-1], amplitudes[2i])
            end
            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    initial_params = vcat(shared_params, initial_amplitudes)

    # Set optimization bounds for shared lifetimes and amplitudes (amplitudes and infinite offset are unbounded)
    lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] - abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] - abs(shared_params[3]) * (bound2 / 100))]

    upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1 / 100)),
        (shared_params[2] + abs(shared_params[2]) * (bound1 / 100)),
        (shared_params[3] + abs(shared_params[3]) * (bound2 / 100))]

    if model_type == :bi
        push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
    elseif model_type == :tri
        push!(lower_lifetimes, (shared_params[4] - abs(shared_params[4]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[4] + abs(shared_params[4]) * (bound2 / 100)))
        push!(lower_lifetimes, (shared_params[5] - abs(shared_params[5]) * (bound2 / 100)))
        push!(upper_lifetimes, (shared_params[5] + abs(shared_params[5]) * (bound2 / 100)))
    end

    lower_amplitudes = fill(-Inf, length(initial_amplitudes))
    upper_amplitudes = fill(Inf, length(initial_amplitudes))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    results = zeros(num_optimizations, length(initial_params))

    for i in 1:num_optimizations
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer))
        results[i, :] = Optim.minimizer(result)
    end

    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    best_lifetimes = round.(avg_values[1:(model_type == :bi ? 4 : model_type == :tri ? 5 : 3)], digits=3)
    best_amplitudes = round.(avg_values[(model_type == :bi ? 5 : model_type == :tri ? 6 : 4):end], digits=3)

    total_ss = 0.0
    total_ss_residual = 0.0
    df_list = []

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
        residuals = kinetic_trace .- fit_values
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)

        trace_df = DataFrame(Time = x, Raw_Data = kinetic_trace, Fit = fit_values, Residuals = residuals)
        push!(df_list, trace_df)
    end

    R² = round(1 - total_ss_residual / total_ss, digits=3)

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

    final_df = vcat(df_list...)

    Fit_Fig = Figure(font="", figure_padding=25, fontsize=20)
    width = 3

    color_palette = palette(ColorSchemes.hsv, 5)

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

        residuals = kinetic_trace .- fit_values

        lines!(Fit_Fig[1, 1], x, kinetic_trace, linewidth = width, label = "Trace $i", color = color_palette[i])
        lines!(Fit_Fig[1, 1], x, fit_values, linewidth = width*0.5, linestyle = :dash, color = :black)
    end

    axislegend(ax1, position = :rb, nbanks = 3, framevisible = false, fontsize = 15)

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

    display(Fit_Fig)
    save("JULIA/Outputs/TriFit_Fig.png", Fit_Fig)

    return final_df, param_amp_df, best_lifetimes, best_amplitudes, R²
end

function GlobalIRF_AllFit_v3_optimized(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing, initial_amplitudes=nothing,bound1=10, bound2=10, num_optimizations=10, noise_level=0.1, size=1000)
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
        t0_guess = 0.01  # Initial guess for t₀
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

    # Calculate R² and Adjusted R²
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


function GlobalIRF_AllFit_v4_optimized(file, columns::Vector{Int}, t0, uppertimebound, model_type::Symbol = :bi, shared_params=nothing, initial_amplitudes=nothing, bound1=10, bound2=10, num_optimizations=10, noise_level=0.1, size=1000)
    sqrt2_inv = 1 / sqrt(2)  # Precompute for efficiency

    function Resize(x, y, size)
        x_out = collect(range(first(x), last(x), length=size))  # New x-axis with uniform size
        interp_func = LinearInterpolation(x, y)  # Linear interpolation
        y_out = interp_func(x_out)  # Interpolated y-values
        return x_out, y_out
    end

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
        t0_guess = 0.01  # Initial guess for t₀
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

            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Step 6: Optimization process with bounds
    initial_params = vcat(shared_params, initial_amplitudes)

    lower_lifetimes = [(shared_params[i] - abs(shared_params[i]) * (bound1 / 100)) for i in 1:length(shared_params)]
    upper_lifetimes = [(shared_params[i] + abs(shared_params[i]) * (bound1 / 100)) for i in 1:length(shared_params)]

    lower_amplitudes = fill(-Inf, length(initial_amplitudes))
    upper_amplitudes = fill(Inf, length(initial_amplitudes))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    results = zeros(num_optimizations, length(initial_params))

    optim_result = nothing  # Store the final optimization result for covariance analysis

    for i in 1:num_optimizations
        # Add small noise to initial parameters to avoid local minima
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))
        inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(inner_optimizer))
        results[i, :] = Optim.minimizer(result)
        if i == num_optimizations
            optim_result = result  # Save the final optimization result
        end
    end

    # Step 7: Calculate mean and std of the optimization results
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    best_lifetimes = round.(avg_values[1:(model_type == :bi ? 4 : model_type == :tri ? 5 : 3)], digits=3)
    best_amplitudes = round.(avg_values[(model_type == :bi ? 5 : model_type == :tri ? 6 : 4):end], digits=3)

    # Step 8: Compute R² and organize results
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

    # Step 10: Extract parameter uncertainties using the Hessian matrix
    # Define a function for Hessian computation using ForwardDiff
    hessian_fn = (params) -> ForwardDiff.hessian(global_fit, params)
    hessian_matrix = hessian_fn(Optim.minimizer(optim_result))

    # Regularize the Hessian matrix by adding a small value to the diagonal
    #hessian_matrix += 1e-9 * I

    covariance_matrix = inv(hessian_matrix)  # Covariance matrix is the inverse of the Hessian
    param_uncertainties = sqrt.(abs.(diag(covariance_matrix)))  # Standard errors of the parameters

    # Step 11: Prepare parameter labels and output DataFrames
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

    param_amp_df = DataFrame(
        Parameter = param_labels, 
        Value = vcat(best_lifetimes, best_amplitudes), 
        Stdev = param_uncertainties
    )

    # Step 12: Visualization using CairoMakie
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


    # Return the final DataFrame with traces, parameters, R², and Adjusted R²
    return final_df, param_amp_df, best_lifetimes, best_amplitudes, R², Adjusted_R²
end


##### COdes that work
function BiExpFitIRF_Global_GPTv3(file, columns::Vector{Int}, t0, uppertimebound, shared_params=nothing, initial_amplitudes=nothing, bound1=10, bound2=10, guess=false, num_optimizations=10, noise_level=0.1, size=1000)
    
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

    # Initialize guess parameters if not provided
    if shared_params === nothing
        t0_guess = 0.01 # Guess the time shift 
        FWHM_guess = 1.5 # Guess FWHM 
        τ₁_guess = 0.2 * (maximum(t) - minimum(t)) # Medium lifetime guess
        τ₂_guess = 0.8 * (maximum(t) - minimum(t)) # Longer lifetime guess
        A_inf_guess = 0.0 # Guess asymptotic amplitude as 0 (can be adjusted)
        shared_params = [t0_guess, FWHM_guess, τ₁_guess, τ₂_guess, A_inf_guess]
        println("Initialized shared_params with: $shared_params")
    end

    if initial_amplitudes === nothing
        # Default initial amplitudes: mean of each trace
        initial_amplitudes = vcat([mean(trace.y) for trace in KineticTraces], [mean(trace.y) for trace in KineticTraces]) # Two amplitudes per trace
        println("Auto-initialized initial_amplitudes: ", initial_amplitudes)
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

            lines!(ax2, x, residuals, linewidth = width)

            linkxaxes!(ax1, ax2)
        end
        axislegend(ax1, position = :rb, nbanks = 3, framevisible = false)

        colsize!(Fit_Fig.layout, 1, Aspect(1, 1.5)) # Set aspect ratio of the y vs. x axis
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout, 2, 130)
        resize_to_layout!(Fit_Fig)
        display(Fit_Fig)
    end
end

function MonoFitIRF_Global_GPTv3(file, columns::Vector{Int}, t0, uppertimebound, shared_params=nothing, initial_amplitudes=nothing, bound1=10, bound2=10, guess=false, num_optimizations=10, noise_level=0.1, size=1000)
    
    # Resize function to interpolate and resize the kinetic trace
    function Resize(x, y, size)
        """Resize x and y arrays to desired length of Data points"""
        x_out = collect(range(first(x), last(x), length=size))
        interp_func = LinearInterpolation(x, y)
        y_out = interp_func(x_out)
        return x_out, y_out
    end

    # Model function (lifetimes shared, amplitudes independent)
    function Y(p, x, A₁)
        @. ((p[4] + A₁ / p[3] * exp(0.5 * (p[2] / p[3])^2 - (x - p[1]) / p[3])) * 
            (erf(((x - p[1]) / p[2] - p[2] / p[3]) / sqrt(2)) + 1) / 2)
    end

    # Cut and resize data based on time bounds
    cut = (t0 .< file[:, 1]) .& (file[:, 1] .< uppertimebound)
    t = file[cut, 1]

    # Preallocate kinetic traces
    KineticTraces = [DataFrame() for _ in columns]

    # Resize data for each column
    for (i, col) in enumerate(columns)
        trace_data = file[cut, col]
        t_resized, trace_resized = Resize(t, trace_data, size)
        KineticTraces[i] = DataFrame(x=t_resized, y=trace_resized)
    end

    # Initialize guess parameters if not provided
    if shared_params === nothing
        t0_guess = 0.01 # Guess the time shift 
        FWHM_guess = 1.5 # Guess FWHM 
        τ₁_guess = 0.2 * (maximum(t) - minimum(t)) # Medium lifetime guess
        A_inf_guess = 0.0 # Guess asymptotic amplitude as 0 (can be adjusted)
        shared_params = [t0_guess, FWHM_guess, τ₁_guess, A_inf_guess]
        println("Initialized shared_params with: $shared_params")
    end

    if initial_amplitudes === nothing
        # Default initial amplitudes: mean of each trace
        initial_amplitudes = [mean(trace.y) for trace in KineticTraces]
        println("Auto-initialized initial_amplitudes: ", initial_amplitudes)
    end

    # Global fitting function with independent amplitudes for each trace
    function global_fit(params)
        shared_coeffs = params[1:4]  # Shared lifetimes
        amplitudes = params[5:end]   # Independent amplitudes
        total_residual = 0.0

        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            model_values = Y(shared_coeffs, x, amplitudes[i])
            residuals = kinetic_trace .- model_values
            total_residual += sum(residuals.^2)
        end
        return total_residual
    end

    # Initial guess parameters: shared lifetimes + independent amplitudes
    initial_params = vcat(shared_params, initial_amplitudes)

    if guess
        # Plot each trace and initial guess
        plot()
        for (i, trace) in enumerate(KineticTraces)
            x = trace.x
            kinetic_trace = trace.y
            plot!(x, kinetic_trace, label="Data Trace $i")
            plot!(x, Y(shared_params, x, initial_amplitudes[i]), label="Initial Guess $i")
        end
        return
    end

    # Set optimization bounds for shared lifetimes and amplitudes
    lower_lifetimes = [(shared_params[1] - abs(shared_params[1]) * (bound1 / 100)),
                       (shared_params[2] - abs(shared_params[2]) * (bound1 / 100)),
                       (shared_params[3] - abs(shared_params[3]) * (bound2 / 100)), -Inf]
    
    upper_lifetimes = [(shared_params[1] + abs(shared_params[1]) * (bound1 / 100)),
                       (shared_params[2] + abs(shared_params[2]) * (bound1 / 100)),
                       (shared_params[3] + abs(shared_params[3]) * (bound2 / 100)), Inf]

    lower_amplitudes = fill(-Inf, length(columns))  # Allow amplitudes to vary
    upper_amplitudes = fill(Inf, length(columns))

    lower_bounds = vcat(lower_lifetimes, lower_amplitudes)
    upper_bounds = vcat(upper_lifetimes, upper_amplitudes)

    # Preallocate results
    results = zeros(num_optimizations, length(initial_params))

    # Optimize the model using multiple initial guesses with noise
    for i in 1:num_optimizations
        noisy_coeffs = initial_params .* (1 .+ noise_level * (rand(length(initial_params)) .- 0.5))

        # Perform global fit using the optimizer
        result = optimize(global_fit, lower_bounds, upper_bounds, noisy_coeffs, Fminbox(LBFGS()))
        results[i, :] = Optim.minimizer(result)
    end

    # Calculate mean and standard deviation of the optimized parameters
    avg_values = mean(results, dims=1) |> vec
    std_values = std(results, dims=1) |> vec

    # Extract the best-fit parameters
    best_lifetimes = round.(avg_values[1:4], digits=3)
    best_amplitudes = round.(avg_values[5:end], digits=3)

    # Calculate global R² and residuals
    total_ss = 0.0
    total_ss_residual = 0.0

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x
        kinetic_trace = trace.y
        fit_values = Y(best_lifetimes, x, best_amplitudes[i])
        residuals = kinetic_trace .- fit_values
        total_ss_residual += sum(residuals.^2)
        total_ss += sum((kinetic_trace .- mean(kinetic_trace)).^2)
    end

    R² = round(1 - total_ss_residual / total_ss, digits=3)

    # Output the fit coefficients and R²
    println("Best Fit Lifetimes: ", best_lifetimes)
    println("Best Fit Amplitudes: ", best_amplitudes)
    println("Global R²: ", R²)

    # Plot the data and final fit for each trace
    fig = Figure(fontsize=20)
    ax1 = Axis(fig[1, 1], title="IRF: t₀ = $(best_lifetimes[1]); FWHM = $(best_lifetimes[2])", xlabel="Time", ylabel="Intensity")
    ax2 = Axis(fig[2, 1], title="Residuals", xlabel="Time", ylabel="Residuals")

    for (i, trace) in enumerate(KineticTraces)
        x = trace.x
        kinetic_trace = trace.y
        fit_values = Y(best_lifetimes, x, best_amplitudes[i])
        residuals = kinetic_trace .- fit_values

        lines!(ax1, x, kinetic_trace, label="Data $i", linewidth=2)
        lines!(ax1, x, fit_values, label="Fit $i", linewidth=1, linestyle=:dash)
        lines!(ax2, x, residuals, linewidth=2)
    end

    axislegend(ax1)
    display(fig)

    return best_lifetimes, best_amplitudes, R²
end

globalfit_test  = MonoFitIRF_Global_GPTv1(Zn_Chlorin_b2_Kinetics, [2, 3, 4], -10.0, 3500, [4.675, 1.998, 763.324, -0.0], [-76.729, -6.094, 7.29], 100, Inf, false,20,0.1, 200)
globalfit_test  = MonoExpFitIRF_Global_GPTv2(Zn_Chlorin_b2_Kinetics, [3], -10.0, 3500,[0.01, 1.5, 694.521, 0.0], [-70.916], 1000, false,3,0.01, 200)

globalfit_test  = BiExpFitIRF_Global_GPTv3(Zn_Chlorin_b2_Kinetics, [2,3,4], -10.0, 4500, nothing,nothing, 100, Inf, false,10,0.1, 200)
globalfit_test  = BiExpFitIRF_Global_GPTv3(Zn_Chlorin_b2_AIRFREE_Kinetics, [2, 3], -10.0, 3500, [5.512, 1.055, 563.081, 188.293, 0.0],[-34.3, -8.26,-5.0,-2.0], 1000, 100, false,3,0.5, 200)


globalfit_test  = BiExpFitIRF_Global_GPTv3(Zn_Chlorin_b2_AIRFREE_Kinetics, [2,5], -10.0, 4500, nothing,nothing, 100, 1000, false,10,0.01, 200)
globalfit_test  = MonoFitIRF_Global_GPTv3(Zn_Chlorin_b2_Kinetics, [2,3,4,5], -10.0, 3500, nothing, nothing, Inf, Inf, false,20,0.1, 200)


globalfit_test  = GlobalIRF_Fit(Zn_Chlorin_b2_Kinetics, [2,3,4], -10.0, 4000.0, :bi, nothing,nothing, 1000, 10000, false,10,0.1, 200)
globalfit_test  = GlobalIRF_AllFit_v1(Zn_Chlorin_b2_AIRFREE_Kinetics, [ 3, 4], -10.0, 4000.0, :bi, nothing,nothing, 100, 1000, false,1,0.1, 200)
globalfit_test  = GlobalIRF_AllFit_v2(Zn_Chlorin_b2_AIRFREE_Kinetics, [2, 3], -10.0, 4000.0, :bi, nothing,nothing, 100, 100,5,0.1, 200)
globalfit_test  = GlobalIRF_AllFit_v2_optimized(Zn_Chlorin_b2_AIRFREE_Kinetics, [2,3, 4], -100.0, 4000.0, :bi, nothing,nothing, 100, Inf,1,0.1, 200)
globalfit_test  = GlobalIRF(Zn_Chlorin_b2_Kinetics, [2, 3, 4], -100.0, 4000.0, :bi, nothing,nothing, 10, Inf,5,0.1, 200)


globalfit_test[2]

