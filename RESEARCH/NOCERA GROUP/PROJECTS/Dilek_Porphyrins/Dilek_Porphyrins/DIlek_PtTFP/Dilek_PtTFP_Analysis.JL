using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,Statistics,MAT ,CairoMakie; CairoMakie.activate!()
plotly()

## Photoluminesence Kinetics

PtTFP_665nm_PL_Airfree = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 665nm PL Kinetics 262us airfree FPT.csv", header = true, DataFrame)
PtTFP_665nm_PL_Airfree_fit = MonoFitIRF(Resize(PtTFP_665nm_PL_Airfree,2,1000),2,-500,260000,[0.5657113594191863, -11.324, 8.104, 219.958, 65792.762, 3.497],true,true,1) #τ₁ = 65.79 μs


PtTFP_740nm_PL_Airfree = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 740nm PL Kinetics 262us airfree FPT.csv", header = true, DataFrame)

PtTFP_665nm_PLKinetics_Air = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 665nm PL Kinetics 65us opened to air.csv", header = true, DataFrame)
PtTFP_665nm_PLKinetics_Air_fit = MonoFitIRF(Resize(PtTFP_665nm_PLKinetics_Air,2,2000),2,-50,10000,[0.498, 7.560, 11.604, 57.064, 1363.066, 0.203],true,true,1) #τ₁ = 1.36 μs


PtTFP_740nm_PLKinetics_Air = CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP 740nm PL Kinetics 65us opened to air.csv", header =true, DataFrame)

## TA Kinetics
PtTFP_N2_TA_spectra = Matrix(CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP N2 bubbled TA Spectrum 420-600nm, 5nm, 2nm step- photodegraded.csv", header = false, DataFrame))
PtTFP_N2_TA_spectra_T = smoothMat(transpose(PtTFP_N2_TA_spectra),3)[2:end, :]
PtTFP_N2_TA_Kinetics = smoothMat(PtTFP_N2_TA_spectra,3)[2:end, :]


PtTFP_air_TA_Spectra = Matrix(CSV.read("/Users/kris/Desktop/mAGNITUDE/PtTFP air TA Spectrum 420-600 nm, 5 nm resolution, 35 mins.csv", header = false,DataFrame))
PtTFP_air_TA_Spectra_T = smoothMat(transpose(PtTFP_air_TA_Spectra),3)[2:end, :]
PtTFP_air_TA_Kinetics = smoothMat(PtTFP_air_TA_Spectra,3)[2:end,:]

PtTFP_air_TA_Spectra[1802,1]

PtTFP_air_455nmTA_Kinetics = MonoFitIRF(Resize(PtTFP_air_TA_Kinetics,9,1000),2,-50,3400,[0.029, -11.325, 8.105, 67.444, 1185.809, -0.488],true,true,1) #τ₁ = 1.185 μs
PtTFP_air_515nmTA_Kinetics = MonoFitIRF(Resize(PtTFP_air_TA_Kinetics,21,1000),2,-50,3400,[-0.0062, -11.559, 7.781, 67.452, 1164.536, -0.179],true,true,1) #τ₁ = 1.164 μs
PtTFP_air_560nmTA_Kinetics = MonoFitIRF(Resize(PtTFP_air_TA_Kinetics,32,1000),2,-50,3400,[0.0043, -11.467, 5.567, 67.538, 1135.820, 0.201],true,true,1) #τ₁ = 1.135 μs


function PtTFP_PL_Kinetics_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pt-TFP PL 665 nm PL kinetics N2",
            palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (μs)", 
            ylabel = "Intensity (A.U.)"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 225:25:800

        # CairoMakie.xlims!(-50,3400)
        # CairoMakie.ylims!(nothing,)

        lines!(PtTFP_665nm_PL_Airfree_fit[:,1]/1000,PtTFP_665nm_PL_Airfree_fit[:,2],linewidth = width,label = "665 nm")
        lines!(PtTFP_665nm_PL_Airfree_fit[:,1]/1000,PtTFP_665nm_PL_Airfree_fit[:,3],linewidth = 2,label = "Fit, 65.8 us ",color = :black)


        # lines!(PtTFP_665nm_PLKinetics_Air_fit[:,1]/1000,PtTFP_665nm_PLKinetics_Air_fit[:,2],linewidth = width,label = "665 nm")
        # lines!(PtTFP_665nm_PLKinetics_Air_fit[:,1]/1000,PtTFP_665nm_PLKinetics_Air_fit[:,3],linewidth = 2,label = "Fit, 1.3 us ", color = :black)

        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_PL_Kinetics_plot()
save("JULIA/Outputs/PtTFP_PL_Kinetics_plot_N2.pdf",PtTFP_PL_Kinetics_plot())

function PtTFP_air_TA_Kinetics_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = GLMakie.Axis(f[1,1], title = "Pt-TFP TA Kinetics opened to air",
            palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (ns)", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 225:25:800

        GLMakie.xlims!(-50,3400)
        # CairoMakie.ylims!(nothing,)


        # lines!(PtTFP_air_TA_Spectra[2:end,1],MaxNorm(PtTFP_air_TA_Spectra[2:end,9]),linewidth = width,label = "455 nm ")
        # lines!(PtTFP_air_TA_Spectra[2:end,1],MaxNorm(PtTFP_air_TA_Spectra[2:end,21]),linewidth = width,label = "515 nm ")
        # lines!(PtTFP_air_TA_Spectra[2:end,1],MaxNorm(PtTFP_air_TA_Spectra[2:end,32]),linewidth = width,label = "560 nm ")

        lines!(PtTFP_air_TA_Spectra[2:end,1],PtTFP_air_TA_Spectra[2:end,9],linewidth = width,label = "455 nm ")
        lines!(PtTFP_air_455nmTA_Kinetics[:,1],PtTFP_air_455nmTA_Kinetics[:,3],linewidth = 2,label = "455 nm fit, 1.19 μs ",color = :black)


        lines!(PtTFP_air_TA_Spectra[2:end,1],PtTFP_air_TA_Spectra[2:end,21],linewidth = width,label = "515 nm ")
        lines!(PtTFP_air_515nmTA_Kinetics[:,1],PtTFP_air_515nmTA_Kinetics[:,3],linewidth = 2,label = "515 nm fit, 1.16 μs ", color = :black)


        lines!(PtTFP_air_TA_Spectra[2:end,1],PtTFP_air_TA_Spectra[2:end,32],linewidth = width,label = "560 nm ")
        lines!(PtTFP_air_560nmTA_Kinetics[:,1],PtTFP_air_560nmTA_Kinetics[:,3],linewidth = 2,label = "560 nm fit, 1.14 μs ", color = :black)


        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Wavelength", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_air_TA_Kinetics_plot()
save("JULIA/Outputs/PtTFP_air_TA_Kinetics_plot.eps",PtTFP_air_TA_Kinetics_plot())

function PtTFP_air_TA_Spectra_plot(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pt-TFP TA Spectra opened to air 3 pt smoothed",
            palette = (color = palette(ColorSchemes.cool, 6),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "ΔA"
        )
        # ax.yticks= 0:0.5:10
        # ax.xticks= 225:25:800

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)

        lines!(PtTFP_air_TA_Spectra_T[:,1],PtTFP_air_TA_Spectra_T[:,312],linewidth = width,label = "20 ns")
        lines!(PtTFP_air_TA_Spectra_T[:,1],PtTFP_air_TA_Spectra_T[:,427],linewidth = width,label = "250 ns")
        lines!(PtTFP_air_TA_Spectra_T[:,1],PtTFP_air_TA_Spectra_T[:,552],linewidth = width,label = "500 ns")
        lines!(PtTFP_air_TA_Spectra_T[:,1],PtTFP_air_TA_Spectra_T[:,802],linewidth = width,label = "1000 ns")
        lines!(PtTFP_air_TA_Spectra_T[:,1],PtTFP_air_TA_Spectra_T[:,1302],linewidth = width,label = "2000 ns")
        lines!(PtTFP_air_TA_Spectra_T[:,1],PtTFP_air_TA_Spectra_T[:,1802],linewidth = width,label = "3000 ns")



        hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
PtTFP_air_TA_Spectra_plot()
save("JULIA/Outputs/PtTFP_air_TA_Spectra_plot.eps",PtTFP_air_TA_Spectra_plot())










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

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
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

function MonoFitIRF(file,column,t0,uppertimebound,coeffs, residuals, ReturnFit, reduce)

    function Fitfunc(x,coeffs)
    
        Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
        mexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ p[6])
    
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
    
        Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))
    
        return Convolution
    end

    Data = file

    rebin = round(length(Data[:,1])/reduce)
    stepsize = (Data[end,1]-Data[1,1])/rebin
    InterpData = LinearInterpolation(Data[:,1],Data[:,column])
    
    t = collect(Data[1,1]:stepsize:Data[end,1])
    I = InterpData(t)

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)


    res = optimize(peakfit0, coeffs,LBFGS(),Optim.Options(time_limit = 200))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual = KineticTrace .- Fitfunc(time,values)
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
        lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = 2, label = "Fit", color = :red)
        lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = 2, label = "Guess", color = :green)
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
        Plots.plot!(time, Fitfunc(time,coeffs), color = :green, label = "Guess")
        Plots.plot!(time, Fitfunc(time,values), linewidth = 3, color = :red, label = "Fit")

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

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
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

function BiFitIRF(file,column,t0,uppertimebound,coeffs, residuals, ReturnFit, reduce)

    function Fitfunc(x,coeffs)
    
        Gaussian(p) = exp.(-4 .* log(2) .* (((x .- p[2]) ./ p[3]) .^2))
        Biexp(p) = (p[4] .* exp.((-(x) ./ p[5])) .+ (p[6] .* exp.((-(x) ./ p[7]))) .+ p[8])
    
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
    
        Convolution = DirectLinearConvolution(IRF,Exp) * (coeffs[1]/maximum(DirectLinearConvolution(IRF,Exp)))
    
        return Convolution
    end

    Data = file

    rebin = round(length(Data[:,1])/reduce)
    stepsize = (Data[end,1]-Data[1,1])/rebin
    InterpData = LinearInterpolation(Data[:,1],Data[:,2])
    
    t = collect(Data[1,1]:stepsize:Data[end,1])
    I = InterpData(t)

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]

    ## FITTING FUNCTION
    peakfit0(coeffs) = sum((KineticTrace .- Fitfunc(time,coeffs)).^2)


    res = optimize(peakfit0, coeffs,LBFGS(),Optim.Options(time_limit = 200))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual = KineticTrace .- Fitfunc(time,values)
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
        lines!(Fit_Fig[1,1],time, Fitfunc(time,values), linewidth = 2, label = "Fit", color = :red)
        lines!(Fit_Fig[1,1],time,Fitfunc(time,coeffs), linewidth = 2, label = "Guess", color = :green)
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
        Plots.plot!(time, Fitfunc(time,coeffs), color = :green, label = "Guess")
        Plots.plot!(time, Fitfunc(time,values), linewidth = 3, color = :red, label = "Fit")

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

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
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


# SOME FUNCTIONS

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
    y_norm = y/(maximum(abs.(y))) #normalize to desired value
    return y_norm
end

function Resize(Data,column, length)
    """Resize a column of a dataframe to desired length of Data points """
    x_out = collect(range(Data[1,1],Data[end,1],length))
    Interp  = LinearInterpolation(Data[:,1],Data[:,column])
    y_out = Interp(x_out)

    df = DataFrame(x = x_out,y = y_out)
    return df
end

function InterpOnRange(Data,start,stop, interval)

    x_out = collect(start:interval:stop)
    Interp  = LinearInterpolation(Data[:,1],Data[:,2])
    y_out = Interp(x_out)

    df = DataFrame(x = x_out,y = y_out)
    return df
end

function SpectralCorrection(Wavelength,Spectrum,Time,Kinetic,CorrWL,CorrTime)
    """Correct a TA spectrum to match with a given kinetic trace by scaling the spectrum by a factor"""
    SpectrumInterp = LinearInterpolation(Wavelength,Spectrum)
    Spectrum_Y = SpectrumInterp(Wavelength)

    KineticInterp = LinearInterpolation(Time,Kinetic)
    
    CorrFactor = KineticInterp[CorrTime]/SpectrumInterp[CorrWL]

    CorrSpectrum = Spectrum_Y .* CorrFactor

    print("factor is:", CorrFactor)
    return CorrSpectrum
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

function smoothMat(Mat,SmoothPts)
    """This function takes and input matrix Mar and generates a new matrix sMat with the moving average over
     some number of points defined by the SmoothPts argument"""
    sMat = zeros(size(Mat,1)-SmoothPts,size(Mat,2))

    sMat[1,:] = Mat[1,:]

    for i in 2:size(Mat,1)-SmoothPts
        for j in 1:size(Mat,2)
            sMat[i,j] = mean(skipmissing(Mat[i:i+SmoothPts,j])) # This works! calculates the mean of the i to i+2th row in the jth columns and replaces the E[i,j]th value with that mean
        end
    end
    return sMat

end 