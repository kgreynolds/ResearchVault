Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Outside(),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)


using Interpolations,Plots,Statistics, DataFrames, CSV, Optim, ColorSchemes,MAT,LineSearches,SpecialFunctions, EasyFit,CairoMakie; CairoMakie.activate!()

using MultivariateStats, DataFrames

using GLMakie
GLMakie.activate!() # hide

    # Load your data (replace "data.csv" with your actual file)
    data = Matrix(CSV.read("/Users/kris/Downloads/Test.csv", header = false, DataFrame))
    data[2:end,1]
    data[1,2:end]
    # Extract the absorbance data (excluding the first row and column)
    absorbance_data = data[2:end, 2:end]

    # Perform PCA on the absorbance data
    pca = fit(PCA, absorbance_data; maxoutdim=3)
    PCAs = projection(pca)

    x = data[2:end,1]
    times = transpose(Matrix(data))[2:end,1]

    Yte = predict(pca, absorbance_data)
    Yte[1,:]
    Plots.plot(times, Yte[1,:])
    Plots.plot!(times, Yte[2,:])
    Plots.plot!(times, Yte[3,:])

function PCA_Fig(f = Figure())

    
    Data = absorbance_data # this needs to be a maatrix

    ### First Panel ##############################################################################

        ax =    CairoMakie.Axis(
                f[1,1], palette = (color = palette(ColorSchemes.inferno, 4),), ## or replace ncol with Integer 
                title = "carped data",
                xlabel = "Wavelength (nm)", 
                ylabel = "Time (ns)", 
                # yscale = Makie.pseudolog10
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(nothing,0.04)

        ### Plot Data ##############################################################################

            CairoMakie.heatmap!(f[1,1],data[2:end,1],data[1,2:end] ,Data,colormap = Reverse(:rainbow))
            
            # Colorbar(f[1,1,Right()], hm, height = 350)
            
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
        lines!(x, PCAs[:,3],linewidth = width, label = "PCA 3")

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
            lines!(times, Yte[3,:],linewidth = width, label = "PCA 3")

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
PCA_Fig()
save("JULIA/Outputs/PCA_Fig.pdf",PCA_Fig())

########################################################################

using Optim, DataFrames

# Sample data (replace with your actual data)
data1 = DataFrame(times=1:10, value=exp.(-1*1:8))
data2 = DataFrame(times=1:10, value=exp.(-0.5*1:9))

value=exp.(-1*1:8)

# Define the model function
function monoexp_model(p, x)
    return p[1] .* exp.(-p[2] .* x)
end

# Create a function to calculate the objective function
function objective(p, x, y)
    for i in 1:length(x)
    return  (sum((monoexp_model(p, x[i]) - y[i]) .^ 2))
    end
end

# Fit the model to the data
x = [data1.times, data2.times]
y = [data1.value, data2.value]
p0 = [-1.0, 1.0] # Initial guess for parameters
result = optimize(p -> objective(p, x, y), p0)

# Extract fitted parameters
params = result.minimum
values = round.(Optim.minimizer(result),digits = 3)