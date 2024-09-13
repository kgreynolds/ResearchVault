using Interpolations,Plots,Statistics, DataFrames, CSV, Optim, ColorSchemes,MAT,LineSearches,SpecialFunctions, EasyFit,CairoMakie; CairoMakie.activate!()
# new comment
using MultivariateStats, DataFrames
plotly()
# Load your data (replace "data.csv" with your actual file)
data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_03/MesAcr_PhosphineOxide_TA_carpet.csv", header = false, DataFrame)

# Extract the absorbance data (excluding the first row and column)
absorbance_data = Matrix(data[2:end, 2:end])

# Perform PCA on the absorbance data
pca = fit(PCA, absorbance_data; maxoutdim=3)
PCAs = projection(pca)

x = data[2:end,1]
times = transpose(Matrix(data))[2:end,1]


Plots.plot(x, PCAs[:,1])
Plots.plot!(x, PCAs[:,2])
Plots.plot!(x, PCAs[:,3])

Yte = predict(pca, absorbance_data)
Yte[1,:]
Plots.plot(times, Yte[1,:])
Plots.plot!(times, Yte[2,:])
Plots.plot!(times, Yte[3,:])


Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Outside(),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)

function Fig(f = Figure())

    

    Data = absorbance_data # this needs to be a maatrix

    ### First Panel ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, 4),), ## or replace ncol with Integer 
                title = "carped data",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

            CairoMakie.heatmap!(Data)
            Colorbar(f[1,1,Right()])
            # axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #rt (right-top :rb (right bottom) etc)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)

            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        

    # ### Second Panel ##############################################################################

        ax2 =    CairoMakie.Axis(
                f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "PCA",
                xlabel = "Time (ns)", 
                ylabel = "ΔA (OD)", 
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

        lines!(x, PCAs[:,1],linewidth = width,label = "PCA 1")
        lines!(x, PCAs[:,2],linewidth = width, label = "PCA 2")
        # lines!(x, PCAs[:,3],linewidth = width, label = "PCA 3")

            axislegend("Time", position = :rt,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

    # ### Third Panel ##############################################################################

        ax3 =    CairoMakie.Axis(
                f[2,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "PCA",
                xlabel = "Time (ns)", 
                ylabel = "ΔA (OD)", 
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################


            lines!(times, Yte[1,:],linewidth = width,label = "PCA 1")
            lines!(times, Yte[2,:],linewidth = width, label = "PCA 2")
            # lines!(times, Yte[3,:],linewidth = width, label = "PCA 3")

            axislegend("Wavelength", position = :rb,nbanks = 3,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)


        Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
        Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
        colgap!(f.layout,1,20)
        rowsize!(f.layout,2,Aspect(1,0.75))
        colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
Fig()
save("JULIA/Outputs/Fig.pdf",Fig())
GLMakie.closeall()







using LinearAlgebra
using Statistics
using Plots

function pca_ta(data::Matrix{Float64})
    # Assume data is a matrix where each column is a spectrum at a different time point
    # and each row corresponds to a wavelength

    # Center the data
    data_centered = data .- mean(data, dims=2)

    # Perform SVD
    U, S, V = svd(data_centered)

    # Calculate explained variance ratio
    explained_variance_ratio = S.^2 / sum(S.^2)

    # Return principal components, singular values, and explained variance ratio
    return U, S, V, explained_variance_ratio
end

pca_ta(absorbance_data)


# Example usage
function plot_pca_results(wavelengths, times, U, S, V, explained_variance_ratio)
    # Plot first few principal components
    p1 = Plots.plot(wavelengths, U[:, 1:3], 
              title="Principal Components", 
              xlabel="Wavelength (nm)", 
              ylabel="Amplitude",
              label=["PC1" "PC2" "PC3"])

    # Plot singular values
    p2 = Plots.plot(1:length(S), S, 
              title="Singular Values", 
              xlabel="Component", 
              ylabel="Singular Value",
              yscale=:log10,
              seriestype=:scatter)

    # Plot explained variance ratio
    p3 = Plots.plot(1:length(explained_variance_ratio), cumsum(explained_variance_ratio), 
              title="Cumulative Explained Variance Ratio", 
              xlabel="Number of Components", 
              ylabel="Cumulative Explained Variance Ratio")

    # Plot time-dependent scores for first few PCs
    p4 = Plots.plot(times, V[:, 1:3], 
              title="Time-dependent Scores", 
              xlabel="Time", 
              ylabel="Score",
              label=["PC1" "PC2" "PC3"])

    # Combine plots
    Plots.plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))
end

# Main execution
function main()
    # Load your data here
    # data = load_ta_data()  # You need to implement this function
    # wavelengths = ...  # Vector of wavelengths
    # times = ...  # Vector of time points

    # Perform PCA
    U, S, V, explained_variance_ratio = pca_ta(absorbance_data)

    # Plot results
    plot_pca_results(data[2:end,1], times, U, S, V, explained_variance_ratio)
end

main()





using LinearAlgebra
using Statistics
using Plots
using Optim

function pca_ta(data::Matrix{Float64})
    # Assume data is a matrix where each column is a spectrum at a different time point
    # and each row corresponds to a wavelength

    # Center the data
    data_centered = data .- mean(data, dims=2)

    # Perform SVD
    U, S, V = svd(data_centered)

    # Calculate explained variance ratio
    explained_variance_ratio = S.^2 / sum(S.^2)

    # Return principal components, singular values, and explained variance ratio
    return U, S, V, explained_variance_ratio
end

function global_fit(data::Matrix{Float64}, times::Vector{Float64}, num_components::Int)
    # Define the model function (sum of exponentials)
    function model(t, p)
        return sum(p[i] * exp.(-t / p[i+num_components]) for i in 1:num_components)
    end

    # Define the objective function to minimize
    function objective(p)
        return sum((data .- model(times, p)).^2)
    end

    # Initial guess for parameters (amplitudes and time constants)
    initial_guess = vcat(ones(num_components), range(minimum(times), maximum(times), length=num_components))

    # Perform the optimization
    result = optimize(objective, initial_guess, LBFGS())

    # Extract the optimized parameters
    optimized_params = Optim.minimizer(result)

    return optimized_params
end

# Example usage
function plot_results(wavelengths, times, U, S, V, explained_variance_ratio, global_fit_params)
    # Plot first few principal components
    p1 = plot(wavelengths, U[:, 1:3], 
              title="Principal Components", 
              xlabel="Wavelength (nm)", 
              ylabel="Amplitude",
              label=["PC1" "PC2" "PC3"])

    # Plot singular values
    p2 = plot(1:length(S), S, 
              title="Singular Values", 
              xlabel="Component", 
              ylabel="Singular Value",
              yscale=:log10,
              seriestype=:scatter)

    # Plot explained variance ratio
    p3 = plot(1:length(explained_variance_ratio), cumsum(explained_variance_ratio), 
              title="Cumulative Explained Variance Ratio", 
              xlabel="Number of Components", 
              ylabel="Cumulative Explained Variance Ratio")

    # Plot time-dependent scores for first few PCs
    p4 = plot(times, V[:, 1:3], 
              title="Time-dependent Scores", 
              xlabel="Time", 
              ylabel="Score",
              label=["PC1" "PC2" "PC3"])

    # Plot global fit results
    num_components = length(global_fit_params) ÷ 2
    p5 = plot(times, [sum(global_fit_params[i] * exp.(-times / global_fit_params[i+num_components]) for i in 1:num_components)],
              title="Global Fit",
              xlabel="Time",
              ylabel="Amplitude",
              label="Fit")
    scatter!(p5, times, data[1, :], label="Data")

    # Combine plots
    plot(p1, p2, p3, p4, p5, layout=(3,2), size=(1000,800))
end

# Main execution
function main()
    # Load your data here
    # data = load_ta_data()  # You need to implement this function
    # wavelengths = ...  # Vector of wavelengths
    # times = ...  # Vector of time points

    # Perform PCA
    U, S, V, explained_variance_ratio = pca_ta(absorbance_data)

    # Perform global fit
    num_components = 3  # Choose the number of exponential components
    global_fit_params = global_fit(absorbance_data, times, num_components)

    # Plot results
    plot_results(data[2:end,1], times, U, S, V, explained_variance_ratio, global_fit_params)
end

main()

function model(t, p)
    return sum(p[i] * exp.(-t / p[i+num_components]) for i in 1:num_components)
end

num_components = 2
times

model(times,initial_guess)

# Define the objective function to minimize
function objective(p)
    return sum((data .- model(times, p)).^2)
end



# Initial guess for parameters (amplitudes and time constants)
initial_guess = vcat(ones(num_components), range(minimum(times), maximum(times), length=num_components))

result = optimize(objective, initial_guess, LBFGS())

# Extract the optimized parameters
optimized_params = Optim.minimizer(result)







###################################################################################################################
using LsqFit, LinearAlgebra, Plots

# Function to generate model matrix
function model_matrix(t, τ, n_components)
    A = zeros(length(t), n_components)
    for i in 1:n_components
        A[:, i] = exp.(-t ./ τ[i])
    end
    return A
end

# Global analysis function
function global_analysis(data, times, wavelengths, n_components)
    # Initial guess for time constants
    τ_init = 10.0 .^ range(0, stop=log10(maximum(times)), length=n_components)
    
    # Objective function for optimization
    function objective(params)
        τ = params[1:n_components]
        A = model_matrix(times, τ, n_components)
        C = data * A
        residuals = data - C * A'
        return residuals
    end

    # Perform the fit
    fit = fit = curve_fit(objective, τ_init, data, 3)
    τ_fit = fit.param
    
    # Calculate decay associated spectra (DAS)
    A_fit = model_matrix(times, τ_fit, n_components)
    DAS = data * A_fit

    return τ_fit, DAS
end

# Example usage
# Assuming 'data' is a matrix where rows are wavelengths and columns are time points
# 'times' is a vector of time points
# 'wavelengths' is a vector of wavelengths

n_components = 3  # Number of exponential components
τ_fit, DAS = global_analysis(absorbance_data, times, data[2:end,1], n_components)

# Plot results
plot(wavelengths, DAS, label=["Component $i" for i in 1:n_components], 
     xlabel="Wavelength (nm)", ylabel="Amplitude", title="Decay Associated Spectra")


     τ_init = 10.0 .^ range(0, stop=log10(maximum(times)), length=n_components)