using Interpolations,Plots, DataFrames,Statistics, CSV, Optim, ColorSchemes,MAT,LineSearches,LaTeXStrings ,CairoMakie; CairoMakie.activate!()
plotly()


matfile = matopen("/Users/kris/Desktop/ResearchVault/JULIA/Inputs/Jane/24_0516_25 uM photobeta2-25 uM wt-alpha2_10 mM ruhex_410 nm.mat")
newmat = read(matfile,"data")
KineticTrace = hcat(newmat[:,1],newmat[:,2]*1000) # can use this to scale the y axis if too small
Smoothed_Kinetic = smoothMat(KineticTrace,2)

Plots.plot(KineticTrace[:,1],KineticTrace[:,2])
Plots.plot!(Smoothed_Kinetic[:,1],Smoothed_Kinetic[:,2])




BiFitIRF(Smoothed_Kinetic,2,-30.0,100000,[0.712, 7.231, 10.668, 23221.586, 9859.434, 242527.394, 311.445, 7348.662, 29745.925],1000, 0,0)
MonoFitIRF(Smoothed_Kinetic,2,-30.0,1000,[0.877, 7.491, 10.879, 4393.568, 300.116, 846.59, 30038.422],100000, 50,10)


save("JULIA/Outputs/Fig1.pdf",BiFitIRF(smoothedmat,2,0.0,100000,[0.219, 0.0, 10.451, 84.413, 13869.571, 376.454, 737.345, 35.513],false,1000))




function BiFitIRF(file,column,t0,uppertimebound,coeffs, size,bound1, bound2)

    function Fitfunc(x,coeffs)
    
        Gaussian(p) =  exp.(-((2 .* log(2) .* (x .- p[2])) ./ p[3]) .^2)
        Biexp(p) = (p[4] .* exp.((-(x .- p[2]) ./ p[5])) .+ (p[6] .* exp.((-(x .- p[2]) ./ p[7]))) .+ p[8])
    
        IRF = Gaussian(coeffs)
        Exp = Biexp(coeffs)
    
        function DirectLinearConvolution(f,g)
            N = length(f) # for This to work both f and g need to be the same size
            M = length(g)
        
            Conv = zeros(N)
            for n=1:N
                for m=1:M
                    if n-m+1 > 0
                        Conv[n] = Conv[n] + f[m] * g[n-m+1]
                    end
                    # n+1 <= m
                end
            end
            return Conv
        end

        Convolution = DirectLinearConvolution(IRF,Exp) .* (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp))) .-  (coeffs[9] .* exp.(-((2 .* log(2) .* (x .- coeffs[2])) ./ coeffs[3]) .^2))

    
        return Convolution
    end

    function Resize(file,column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1],file[end,1],size))
        Interp  = LinearInterpolation(file[:,1],file[:,column])
        y_out = Interp(x_out)
    
        ResizedDF = DataFrame(x = x_out,y = y_out)
        return ResizedDF
    end

    Data = Resize(file,column,size)

    t = Data[:,1]
    I = Data[:,2]

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)

    lower = [-Inf,(coeffs[2]-abs(coeffs[2])*(bound1/100)),(coeffs[3]-abs(coeffs[3]*(bound1/100))),-Inf, (coeffs[5]-abs(coeffs[5]*(bound2/100))),-Inf, (coeffs[7]-abs(coeffs[7]*(bound2/100))),-Inf, -Inf]
    upper = [Inf,(coeffs[2]+abs(coeffs[2])*(bound1/100)),(coeffs[3]+abs(coeffs[3]*(bound1/100))),Inf,(coeffs[5]+abs(coeffs[5]*(bound2/100))), Inf,(coeffs[7]+abs(coeffs[7]*(bound2/100))),Inf, Inf]

    inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking(order=3))
    res = optimize(peakfit0, lower, upper, coeffs, Fminbox(inner_optimizer),Optim.Options(time_limit = 200,g_tol = 1e-8))

    print("Fit Resuluts:", res)
    values = round.(Optim.minimizer(res),digits = 3) # Kinetics results
    residual = KineticTrace .- Fitfunc(time,values)
    print("Fit Coeff:", values)

    df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)

       
    
        Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
        width = 2

        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(values[2]); FWHM = $(values[3])" , subtitle = "τ₁ = $(values[5]); τ₂ = $(values[7])",
        palette = (color = palette(ColorSchemes.hsv, 4),),
        xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
        # xlims!(ax1,350,650)
        # ylims!(nothing,)

        lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data")
        lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = width, label = "Fit")
        lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = width, label = "Guess")
        axislegend(ax1, position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
        xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        # resize_to_layout!(Fit_Fig)

        display(Fit_Fig)
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
    return df
end

function MonoFitIRF(file,column,t0,uppertimebound,coeffs, size, bound1,bound2)

    function Fitfunction(x,coeffs)
    
        Gaussian(p) =  exp.(-((2 .* log(2) .* (x .- p[2])) ./ p[3]) .^2)
        mexp(p) = (p[4] .* exp.((-(x .- p[2]) ./ p[5])) .+ p[6])
    
        IRF = Gaussian(coeffs)
        Exp = mexp(coeffs)
    
        function DirectLinearConvolution(f,g)
            N = length(f) # for This to work both f and g need to be the same size
            M = length(g)
        
            Conv = zeros(N)
            for n=1:N
                for m=1:M
                    if n-m+1 > 0
                        Conv[n] = Conv[n] + f[m] * g[n-m+1]
                    end
                    # n+1 <= m
                end
            end
            return Conv
        end
    
        Convolution = (DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))) .-  (coeffs[7] .* exp.(-((2 .* log(2) .* (x .- coeffs[2])) ./ coeffs[3]) .^2))

        return Convolution
    end

    function Resize(file,column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1],file[end,1],size))
        Interp  = LinearInterpolation(file[:,1],file[:,column])
        y_out = Interp(x_out)
    
        ResizedDF = DataFrame(x = x_out,y = y_out)
        return ResizedDF
    end

    Data = Resize(file,column,size)
    
    t = Data[:,1]
    I = Data[:,2]

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunction(time,coeffs)).^2)

    lower = [-Inf,(coeffs[2]-abs(coeffs[2])*(bound1/100)),(coeffs[3]-abs(coeffs[3]*(bound1/100))),-Inf, (coeffs[5]-abs(coeffs[5]*(bound2/100))),-Inf,-Inf]
    upper = [Inf,(coeffs[2]+abs(coeffs[2])*(bound1/100)),(coeffs[3]+abs(coeffs[3]*(bound1/100))),Inf,(coeffs[5]+abs(coeffs[5]*(bound2/100))), Inf, Inf]

    inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
    res = optimize(peakfit0, lower, upper, coeffs, Fminbox(inner_optimizer),Optim.Options(time_limit = 200,g_tol = 1e-8))
    
    print("Fit Resuluts:", res)
    values = round.(Optim.minimizer(res),digits = 3) # Kinetics results
    residual = KineticTrace .- Fitfunction(time,values)
    print("Fit Coeff:", values)

    df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunction(time,values),Residual = residual)

    
    

        Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
        width = 2

        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title ="IRF: t₀ = $(values[2]); FWHM = $(values[3])" , subtitle = "A₁ = $(values[4]); τ₁ = $(values[5])",
        palette = (color = palette(ColorSchemes.hsv, 3),),
        xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
        # xlims!(ax1,350,650)
        # ylims!(nothing,)

        lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data")
        lines!(Fit_Fig[1,1],time, Fitfunction(time,values), linewidth = width, label = "Fit")
        # lines!(Fit_Fig[1,1],time,Fitfunction(time,coeffs), linewidth = width, label = "Guess")
        axislegend(ax1, position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
        xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        display(Fit_Fig)
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)


    return  df
    # return fitfunction(values)
end

function TriFitIRF(file,column,t0,uppertimebound,coeffs,size,bound1,bound2)

    function Fitfunc(x,coeffs)
    
        Gaussian(p) =  exp.(-((2 .* log(2) .* (x .- p[2])) ./ p[3]) .^2)
        Triexp(p) = (p[4] .* exp.((-(x .- p[2]) ./ p[5])) .+ (p[6] .* exp.((-(x .- p[2]) ./ p[7]))) .+ (p[8] .* exp.((-(x .- p[2]) ./ p[9]))) .+ p[10])
    
        IRF = Gaussian(coeffs)
        Exp = Triexp(coeffs)
    
        function DirectLinearConvolution(f,g)
            N = length(f) # for This to work both f and g need to be the same size
            M = length(g)
        
            Conv = zeros(N)
            for n=1:N
                for m=1:M
                    if n-m+1 > 0
                        Conv[n] = Conv[n] + f[m] * g[n-m+1]
                    end
                    # n+1 <= m
                end
            end
            return Conv
        end
    
        Convolution = (DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))) .-  (coeffs[11] .* exp.(-((2 .* log(2) .* (x .- coeffs[2])) ./ coeffs[3]) .^2))
    
        return Convolution
    end

    function Resize(file,column, size)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(file[1,1],file[end,1],size))
        Interp  = LinearInterpolation(file[:,1],file[:,column])
        y_out = Interp(x_out)
    
        ResizedDF = DataFrame(x = x_out,y = y_out)
        return ResizedDF
    end

    Data = Resize(file,column,size)

    t = Data[:,1]
    I = Data[:,2]

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)

    lower = [-Inf,(coeffs[2]-abs(coeffs[2])*(bound1/100)),(coeffs[3]-abs(coeffs[3]*(bound1/100))),-Inf, (coeffs[5]-abs(coeffs[5]*(bound2/100))),-Inf, (coeffs[7]-abs(coeffs[7]*(bound2/100))),-Inf, (coeffs[9]-abs(coeffs[9]*(bound2/100))), -Inf, -Inf]
    upper = [Inf,(coeffs[2]+abs(coeffs[2])*(bound1/100)),(coeffs[3]+abs(coeffs[3]*(bound1/100))),Inf,(coeffs[5]+abs(coeffs[5]*(bound2/100))), Inf,(coeffs[7]+abs(coeffs[7]*(bound2/100))),Inf, (coeffs[9]-abs(coeffs[9]*(bound2/100))), Inf,Inf]

    inner_optimizer = LBFGS(linesearch=LineSearches.BackTracking())
    res = optimize(peakfit0, lower, upper, coeffs, Fminbox(inner_optimizer),Optim.Options(time_limit = 200,g_tol = 1e-8))

    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual = KineticTrace .- Fitfunc(time,values)
    print("Fit Coeff:", values)

    df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)
  
    
        Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
        width = 2

        ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "IRF: t₀ = $(values[2]); FWHM = $(values[3])" , subtitle = "A₁ = $(values[4]); τ₁ = $(values[5]); A₂ = $(values[6]); τ₂ = $(values[7]); A₃ = $(values[8]); τ₃ = $(values[9])",
        palette = (color = palette(ColorSchemes.hsv, 3),),
        xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
        # xlims!(ax1,350,650)
        # ylims!(nothing,)

        lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
        lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = width, label = "Fit", color = :red)
        # lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = width, label = "Guess", color = :green)
        axislegend(ax1, position = :rt, framevisible = false)

        ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
        xlabel = "",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
        ylabel =  "", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)

        lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)

        colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
        rowsize!(Fit_Fig.layout, 1, 300)
        rowsize!(Fit_Fig.layout,2,100)
        resize_to_layout!(Fit_Fig)

        display(Fit_Fig)
        # save("JULIA/Outputs/MonoFit_Fig.png",Fit_Fig)
 
    return df

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


