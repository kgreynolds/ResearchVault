using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,DelimitedFiles,Peaks ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

# Bulk Electrolysis Reduction of Ir(III)K in MeCN at -1.2V Absolute. 
IrK_Red = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_10_04_SpecEchem_IrK_MeCN/Reductive-Ir(III)_0.2mM_80mM_TBAPF6_AgAgCl -1.2V Reductive.csv", header = true, DataFrame)

# Bulk Electrolysis Re-Oxidation of Ir(III)K in MeCN at -1.2V Absolute. 
IrK_Red_Ox = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_10_04_SpecEchem_IrK_MeCN/Oxidative-Ir(III)_0.2mM_80mM_TBAPF6_AgAgCl +2V Oxidative.csv", header = true, DataFrame)


# CVs of Ferrocene and Ir(III)K in MeCN Electrolyte [TBA][PF6]
Fc_CV = readdlm("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_10_04_SpecEchem_IrK_MeCN/Fc_0p1M_TBAPF6_MeCN_AgAgCl.txt", skipstart = 34)
Ir_CV = readdlm("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_10_04_SpecEchem_IrK_MeCN/Ir(III)K_0p2uM_0p08M_TBAPF6_CV_AgAgCl_100mVps.txt", skipstart = 42)

Plots.plot(Ir_CV[:,1],Ir_CV[:,2])


Plots.plot()
for i in 1:50
    data = IrK_Red[i]
   display(Plots.plot!(data[:,1],data[:,2],xlims = (350,800)))
end

### Make Figures 

function Ferrocene_CV(f = Figure())
    width = 3

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "FC CV in MeCN",
            palette = (color = palette(ColorSchemes.hsv, 11),), xautolimitmargin = (0.05, 0.05),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Potneital vs. Ag/AgCl (V)", 
            ylabel = "Current (μA)"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(350,800)
        # CairoMakie.ylims!(-0.1,2.5)

        lines!(Fc_CV[615:end,1],Fc_CV[615:end,2]*(10^6),linewidth = width,label = "100 mV/s, 80 mM [TBA][PF6]")
        # lines!(IrK_Red,Data[:,2],linewidth = width,label = "10 ns")


        vlines!(ax,0.579, linewidth = 2,color = :black,linestyle = :dash)

        axislegend( position = :lb,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        text!(0.35,-10,text = "E(1/2) = 0.579 V", color = :black, textsize = 20, align = (:center, :center))

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Ferrocene_CV()
save("JULIA/Outputs/Ferrocene_CV.pdf",Ferrocene_CV())



function IrK_CV(f = Figure())
   """Reduction and oxidaion potentials are referenceed 
   to the E(1/2) of the CV of Fc/Fc+ couple in MeCN from above of 0.579 V """
   
    width = 3
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Ir(III)K CV in MeCN",
            palette = (color = palette(ColorSchemes.winter, 5),), xautolimitmargin = (0.05, 0.05),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Potneital vs. Fc/Fc+ (V)", 
            ylabel = "Current (μA)"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(350,800)
        # CairoMakie.ylims!(-0.1,2.5)

        lines!(Ir_CV[1700:end,1] .- 0.579,Ir_CV[1700:end,2]*(10^6),linewidth = width,label = "100 mV/s, 80 mM [TBA][PF6]")
        # lines!(IrK_Red,Data[:,2],linewidth = width,label = "10 ns")


        vlines!(ax,-1.731,ymin = 0.2,ymax = 1, linewidth = 2,color = :black,linestyle = :dash)
        vlines!(ax,1.323,ymin = 0,ymax = 0.8, linewidth = 2,color = :black,linestyle = :dash)

        axislegend( position = :lb,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        text!(-1,3,text = "E(red) = -1.731 V", color = :black, textsize = 20, align = (:center, :center))
        text!(0.8,-5,text = "E(ox) = 1.323 V", color = :black, textsize = 20, align = (:center, :center))

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
IrK_CV()
save("JULIA/Outputs/IrK_CV.pdf",IrK_CV())

function SpecEchem(f = Figure())
    width = 3

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Spec-Echem of Ir(III)K in MeCN",
            palette = (color = palette(ColorSchemes.winter, 10),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "A"
        )
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(350,800)
        CairoMakie.ylims!(-0.1,2.5)
        
        # # Plot the Reducing side
        for i in 2:ncol(IrK_Red)
            header = names(IrK_Red)[i]
            lines!(IrK_Red[:,1],IrK_Red[:,i],linewidth = width,label = "$header")
        end

        # Plot the Re-Oxidation
        # for i in 2:ncol(IrK_Red_Ox)
        #     header = names(IrK_Red_Ox)[i]
        #     lines!(IrK_Red_Ox[:,1],IrK_Red_Ox[:,i],linewidth = width,label = "$header")
        # end


        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(600,1.2,text = "-1.2 V vs. Ag/AgCl; 200 μM Ir(III)K, 80 mM TBAPF₆", color = :black, textsize = 20, align = (:center, :center))
        text!(600,1.2,text = "+2.0 V vs. Ag/AgCl; 200 μM Ir(III)K, 80 mM TBAPF₆", color = :black, textsize = 20, align = (:center, :center))

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
SpecEchem()
save("JULIA/Outputs/SpecEchem_red.png",SpecEchem())

### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA_GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_450nm_Kinetics_180us_COnverted.csv",header = true,DataFrame)

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
            push!(files, readdlm(joinpath(root, filename),skipstart = 13))
            end
        end
        end
    
        return files
    end
