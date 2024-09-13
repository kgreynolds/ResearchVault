## Funtions 
using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,MAT,LineSearches, Statistics,CairoMakie; CairoMakie.activate!()



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


    function LinearFit(x,y,residuals, ReturnFit)
    
        time = x
        KineticTrace = y
    
        ## FITTING FUNCTION
        fitfunction = mexp(p) = p[1] .* x .+ p[2]
        p0 = [5188.825, 0.00768]
    
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

    function MonoFit(file,column,t0,uppertimebound, residuals, ReturnFit, coeffs::Vector{Float64})

        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
        p0 = coeffs
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    

        # res = optimize(peakfit0, p0,LBFGS(linesearch=LineSearches.BackTracking()),Optim.Options(time_limit = 20))
        res = optimize(peakfit0, p0,LBFGS(linesearch=LineSearches.MoreThuente()),Optim.Options(time_limit = 20))

        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        residual =  KineticTrace .- fitfunction(values)
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values), residual)
        
    
        elseif residuals == true
            Fit_Fig = Figure(figure_padding = 25,fontsize =20, resolution = (600,600))
            width = 2
    
            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
            xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)
    
            lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
            lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
            lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
            axislegend(ax1, position = :rt, framevisible = false)
    
            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
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

    function BiFit(file,column,t0,uppertimebound, residuals, ReturnFit, coeffs::Vector{Float64})

        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
        p0 = coeffs
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        # res = optimize(peakfit0, p0,LBFGS(linesearch=LineSearches.BackTracking()),Optim.Options(time_limit = 20))
        res = optimize(peakfit0, p0,LBFGS(linesearch=LineSearches.MoreThuente()),Optim.Options(time_limit = 20))

        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        residual =  KineticTrace .- fitfunction(values)
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values), residual)
        
    
        elseif residuals == true
            Fit_Fig = Figure( figure_padding = 25,fontsize =20, resolution = (600,600))
            width = 2
    
            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
            xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)
    
            lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
            lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
            lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
            axislegend(ax1, position = :rt, framevisible = false)
    
            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    
            lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)
    
            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,100)
            resize_to_layout!(Fit_Fig)
    
            Fit_Fig
            # save("JULIA/Outputs/BiFit_Fig.png",Fit_Fig)
        else
    
            Plots.plot(time, KineticTrace, color = :blue, label = "Data")
            Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
            Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")
    
        end
    
        # return fitfunction(values)
    end
 
    function TriFit(file,column,t0,uppertimebound, residuals, ReturnFit)
    
        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
        p0 = [0.01, 2.19, 0.01, 8.7, 0.01, 19.0, 0.01]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        # res = optimize(peakfit0, p0,LBFGS(linesearch=LineSearches.BackTracking()),Optim.Options(time_limit = 20))
        res = optimize(peakfit0, p0,LBFGS(linesearch=LineSearches.MoreThuente()),Optim.Options(time_limit = 20))
        
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        residual = KineticTrace .- fitfunction(values)
        print("Fit Coeff:", values)
    
        df = DataFrame(Time = time,RawData = KineticTrace, Fit = Fitfunc(time,values),Residual = residual)

        if ReturnFit == true
            return df
        
    
        elseif residuals == true
            Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
            width = 2
    
            ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Fit",
            xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
            # xlims!(ax1,350,650)
            # ylims!(nothing,)
    
            lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "Data", color = :blue)
            lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "Fit", color = :red)
            lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
            axislegend(ax1, position = :rt, framevisible = false)
    
            ax2 = CairoMakie.Axis(Fit_Fig[2,1], title = "Residuals",
            xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
            ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    
            lines!(Fit_Fig[2,1],time, residual, linewidth = width, color = :black)
    
            colsize!(Fit_Fig.layout,1,Aspect(1, 1.5)) #Set aspect ration of the y vs. x axis)
            rowsize!(Fit_Fig.layout, 1, 300)
            rowsize!(Fit_Fig.layout,2,100)
            resize_to_layout!(Fit_Fig)
    
            Fit_Fig
            # save("JULIA/Outputs/TriFit_Fig.png",Fit_Fig)
        else
    
            Plots.plot(time, KineticTrace, color = :blue, label = "Data")
            Plots.plot!(time, fitfunction(p0), color = :green, label = "Guess")
            Plots.plot!(time, fitfunction(values), linewidth = 3, color = :red, label = "Fit")
    
        end
    
        # return fitfunction(values)
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

    # Resize Dataframe to new length
    function Resize(Data,column, length)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(Data[1,1],Data[end,1],length))
        Interp  = LinearInterpolation(Data[:,1],Data[:,column])
        y_out = Interp(x_out)

        df = DataFrame(x = x_out,y = y_out)
        return df
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

    function movingavg(x::Vector{Float64},y::Vector{Float64},n::Int64)
        outx::Vector{Float64} = []
        outy::Vector{Float64} = []
        
        for i in 1 : length(y) - n 
            localx = (mean(x[i:i+n]))*1000000 # take the average in the interval of n to the right of i 
            localy = (mean(y[i:i+n])) # take the average in the interval of n to the right of i
            push!(outx, localx)
            push!(outy,localy)
            
        end
        
        df = DataFrame(Time_us = outx, Intensity = outy) #become the headers for the CSV file
        

        return outx, outy, display(plot(outx, outy,xlabel = "Time", ylabel = "Abs")), CSV.write("/Users/kris/Documents/JULIA/Outputs/movingAvg_out.csv", df) ## This is the path and name of your file
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
        return sMat

    end










#Normalization to x=z
function Norm(x,y,z)
    x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
    y_interp = x_interp(x) # compute y values for the x vector
    y_norm = y_interp/x_interp(z) #normalize to desired value
    return y_norm
end

## Original Moving Average code from Ryan
function movingav(x::Vector{Float64},y::Vector{Float64},n::Int64)
    outx::Vector{Float64} = []
    outy::Vector{Float64} = []
    for i in 1 : length(y) - n 
        localx = sum(x[i:i+n]) / (n + 1)
        localy = sum(y[i:i+n]) / (n + 1) #std(x[i:i+n])
        push!(outx, localx)
        push!(outy,localy)
        
    end
    return outx, outy
end

function movingavTA(x::Vector{Float64},y::Vector{Float64},n::Int64)
    outx::Vector{Float64} = []
    outy::Vector{Float64} = []
    for i in 1 : length(y) - n 
        localx = (sum(x[i:i+n]) / (n + 1))*1000000 # take the average in the interval of n to the right of i and convert from seconds to μ seconds
        localy = (sum(y[i:i+n]) / (n + 1)) # take the average in the interval of n to the right of i
        push!(outx, localx)
        push!(outy,localy)
        
    end

    outy_t0 = sum(outy[1:20]/20) # take the average of the first 20 pre time zero points (make sure )

    deltaAbsorbance = real(log.((outy_t0)./outy))

    df = DataFrame(Time = outx, dAbs = deltaAbsorbance) #become the headers for the CSV file

    return outx, deltaAbsorbance, display(plot(outx, deltaAbsorbance, size=(1000, 800), xlabel = "Time(us)", ylabel = "dA(OD)")), CSV.write("/Users/kris/Documents/JULIA/Outputs/movingAvgTA.csv", df)
end

function movingav(x::Vector{Float64},y::Vector{Float64},n::Int64)
    outx::Vector{Float64} = []
    outy::Vector{Float64} = []
    for i in 1 : length(y) - n 
        localx = (mean(x[i:i+n])) # take the average in the interval of n to the right of i 
        localy = (mean(y[i:i+n])) # take the average in the interval of n to the right of i
        push!(outx, localx)
        push!(outy,localy)
        
    end
    
    df = DataFrame(Wavelength = outx, Abs = outy) #become the headers for the CSV file

    return outx, outy, display(plot!(outx, outy,xlims = [380,650], ylims = [0,0.0001], size=(1000, 800),xlabel = "Wavelength", ylabel = "A(OD)")), CSV.write("/Users/kris/Documents/JULIA/Outputs/movingAvgTA_out.csv", df)
end

function write_csv(path::String,headers::Vector{String},data::Vector{Vector{Float64}})
    input::String = ""
    for h in headers
        input *= "$h, "
    end
    input *= "\n"
    block = hcat(data...)
    for i in 1:size(block,1)
        line::String = ""
        for j in 1:size(block,2)
            line *= "$(block[i,j]),"
        end
        input *= "$line \n"
    end
    open(path,"w") do file #actually write the file
        write(file,input)
    end
end

function DirectLinearConvolution(f,g)
    N = length(f)
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

function Fitfunc(x,coeffs)
    
    Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
    mexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ p[6])

    IRF = Gaussian(coeffs)
    Exp = mexp(coeffs)

    function DirectLinearConvolution(f,g)
        N = length(f)
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

    Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))

    return Convolution
end

function heaveside(x::Vector{Float64},t0::Float64) # heaveside step function to esnure exponential is zero beofre the t0 
    
    # mexp(p) = (p[4] .* exp.((-(x .- t0) ./ p[5])) .+ p[6]) ## Chose your fitting function (parameters stored in the coff vector)

    biexp(p) = (p[4] .* exp.((-(x .- t0) ./ p[5]))) .+ (p[6] .* exp.((-(x .- t0) ./ p[7])) .+ p[8])

    cut = Time .>= t0 # generate bit vector for t0
    out = biexp(coeff) .* cut # sets values for the exponetnial to 0 before the cutoff

    return out
end

