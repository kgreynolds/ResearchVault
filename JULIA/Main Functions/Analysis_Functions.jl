using Interpolations,Plots,Statistics, DataFrames, CSV, Optim,LinearAlgebra, ColorSchemes,MAT,LineSearches,SpecialFunctions, EasyFit,CairoMakie; CairoMakie.activate!()
plotly()
set_theme!(theme_latexfonts())

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
function InterpOnRange(Data,start,stop, interval)

    x_out = collect(start:interval:stop)
    Interp  = LinearInterpolation(Data[:,1],Data[:,2])
    y_out = Interp(x_out)

    df = DataFrame(x = x_out,y = y_out)
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

