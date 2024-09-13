using Interpolations,Plots, DataFrames,EasyFit, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA
#UvVis of Ir(III)K
IrK_Uv_Vis = CSV.read("RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_09_20_IrK_Exctinction_Coeff/IrK_ExcCoeffData.csv", header = true, DataFrame)
SV_IrK_PhSSPh = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_09_22_SternVolmer_IrK_PhSSPh/SternVolmer_IrK_200uM_PhSSPh.csv", header = true, DataFrame)
Uncorr_SV_IrK_PhSSPh = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_09_22_SternVolmer_IrK_PhSSPh/Uncorrected_SternVolmer_IrK_200uM_PhSSPh.csv", header = true, DataFrame)

SV_IrK_PhSSPh_norm = SV_IrK_PhSSPh[:,2:end] ./ maximum(abs.(SV_IrK_PhSSPh[:,2])) # Normalize the Array to inital intensity
Uncorr_SV_IrK_PhSSPh_norm = Uncorr_SV_IrK_PhSSPh[:,2:end] ./ maximum(abs.(Uncorr_SV_IrK_PhSSPh[:,2])) # Normalize the Array to inital intensity

Plots.plot(IrK_Uv_Vis[:,1], Norm(IrK_Uv_Vis[:,1],IrK_Uv_Vis[:,3],380))
Plots.plot!(SV_IrK_PhSSPh[:,1], Uncorr_SV_IrK_PhSSPh_norm[:,1])

# Stern Volmer analysis of 200uM Ir(III)K [(Ir[dF(CF3)ppy]2(dtbpy))PF6] and Diphenyl Disulfide in MeCN prepared in Glovebox

Plots.plot(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh[:,2])



SV_IrK_PhSSPh[60,1]
SV = SV_IrK_PhSSPh[84,2] ./ collect(SV_IrK_PhSSPh[84,2:end])
conc = [0.0,1.0,2.0,5.0,10.0] ./ 1000

fit = fitlinear(conc,SV)
round(fit.a ;digits = 3)
Plots.scatter(conc,SV)
Plots.plot!(fit.x,fit.y)
fit.R


### Make Figures 

function EmissionSpectra(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ir(III)K emission quenching by PhSSPh in MeCN",
            palette = (color = palette(ColorSchemes.cool, 6),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "Normalized Intensity"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        # lines!(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh_norm[:,1],linewidth = width,label = "0 mM")
        # lines!(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh_norm[:,2],linewidth = width,label = "1 mM")
        # lines!(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh_norm[:,3],linewidth = width,label = "2 mM")
        # lines!(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh_norm[:,4],linewidth = width,label = "5 mM")
        # lines!(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh_norm[:,5],linewidth = width,label = "10 mM")

        lines!(Uncorr_SV_IrK_PhSSPh[:,1],Uncorr_SV_IrK_PhSSPh_norm[:,1],linewidth = width,label = "0 mM")
        lines!(Uncorr_SV_IrK_PhSSPh[:,1],Uncorr_SV_IrK_PhSSPh_norm[:,2],linewidth = width,label = "1 mM")
        lines!(Uncorr_SV_IrK_PhSSPh[:,1],Uncorr_SV_IrK_PhSSPh_norm[:,3],linewidth = width,label = "2 mM")
        lines!(Uncorr_SV_IrK_PhSSPh[:,1],Uncorr_SV_IrK_PhSSPh_norm[:,4],linewidth = width,label = "5 mM")
        lines!(Uncorr_SV_IrK_PhSSPh[:,1],Uncorr_SV_IrK_PhSSPh_norm[:,5],linewidth = width,label = "10 mM")

        axislegend("[Ph-SS-Ph]", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, labelsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
EmissionSpectra()
save("JULIA/Outputs/EmissionSpectra.pdf",EmissionSpectra())

function E00(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ir(III)K electronic spectea in MeCN",
            palette = (color = palette(ColorSchemes.hsv, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "Normalized Intensity"
        )
        ax.xticks= 300:50:560
        CairoMakie.xlims!(300,560)
        # CairoMakie.ylims!(nothing,)
        
        lines!(IrK_Uv_Vis[:,1],IrK_Uv_Vis[:,3],linewidth = width,label = "Absorption")
        lines!(SV_IrK_PhSSPh[:,1],SV_IrK_PhSSPh_norm[:,1],linewidth = width,label = "Emission")
        lines!(Uncorr_SV_IrK_PhSSPh[:,1],Uncorr_SV_IrK_PhSSPh_norm[:,1],linewidth = width,label = "Emission_uncorr")


        vlines!(ax,455.5, ymax = 0.4, linewidth = 2,color = :black,linestyle = :dash)


        axislegend(position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        text!(455.5,1.1,text = "E₀₀ = 2.7 eV", color = :black, align = (:center, :center), labelsize = 25, font = :bold)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
E00()
save("JULIA/Outputs/E00.pdf",E00())

function SternVolmer(f = Figure())
    width = 3
    Io_I = Uncorr_SV_IrK_PhSSPh_norm[84,1] ./ collect(Uncorr_SV_IrK_PhSSPh_norm[84,1:end])
    Concentrations = [0.0,1.0,2.0,5.0,10.0] ./ 1000
    
    fit = fitlinear(Concentrations,Io_I, b = 1.0)
    a = round(fit.a; digits = 3)
    b = round(fit.b; digits = 3)
    R² = round(fit.R^2; digits = 3)

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Stern Volmer of Ir(III)K (200 μM) + PhSSPh in MeCN",
            palette = (color = palette(ColorSchemes.cool, 6),), xautolimitmargin = (0.05, 0.05),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "[Ph-SS-Ph] (mM)", 
            ylabel = "I₀/I"
        )
        # ax.xticks= 0:2:10
        ax.yticks= 0:1:10

        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        CairoMakie.scatter!(Concentrations,Io_I,linewidth = width,label = "Data", color = "red", markersize = 20, marker = :rect)
        CairoMakie.lines!(fit.x,fit.y,linewidth = width,label = "fit", color = "blue")



        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend(position = :lt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        text!(0.004,2,text = "y =  $a x + $b; R² = $R² ", color = :blue, textsize = 25)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
    # return fit
end
SternVolmer()
save("JULIA/Outputs/SternVolmer.eps",SternVolmer())



### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_450nm_Kinetics_180us_COnverted.csv",header = true,DataFrame)

    ## FITTING FUNCTIONS

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
        y_norm = y/maximum(abs.(y)) #normalize to desired value
        return y_norm
    end

    function MaxNormArray(y, column)
        """Normalize vector to its absolute largest value"""
        y_norm = y/maximum(abs.(y)) #normalize to desired value
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

