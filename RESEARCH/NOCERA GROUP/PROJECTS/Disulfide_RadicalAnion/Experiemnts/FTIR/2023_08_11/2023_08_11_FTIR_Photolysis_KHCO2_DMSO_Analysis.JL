using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

DMSO_UVvis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_11_FTIR_Photolysis_KHCO2_DMSO/DMSO_390um_IRCell_UvVis.csv", header = true, DataFrame)
KHCO2_DMSO_UvVis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_11_FTIR_Photolysis_KHCO2_DMSO/10mM_KHCO2_DMSO_390um_IRCell_UvVis.csv", header = true, DataFrame)
KHCO2_DMSO_UvVis_postPhot = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_11_FTIR_Photolysis_KHCO2_DMSO/10mM_KHCO2_DMSO_390um_IRCell_UvVis_postPhotolysis.csv", header = true, DataFrame)

FTIR_KHCO2_Phtolysis = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_11_FTIR_Photolysis_KHCO2_DMSO/10mM_KHCO2_DMSO_390um_Photolysis.csv", header = true,DataFrame)

### Make Figures 

function Fig1(f = Figure())
    interval = 260:1:800
    Data1 = Interpolate(DMSO_UVvis[:,1],DMSO_UVvis[:,2],interval)
    Data2 = Interpolate(KHCO2_DMSO_UvVis[:,1],KHCO2_DMSO_UvVis[:,2], interval)
    Data3 = Interpolate(KHCO2_DMSO_UvVis_postPhot[:,1],KHCO2_DMSO_UvVis_postPhot[:,2],interval)

    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "UV-Vis of photolysis of KHCO₂ (10 mM) in DMSO",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "A", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0))
    # CairoMakie.xlims!(1250,1790)
    # CairoMakie.ylims!(nothing,)
    width = 4
    x = collect(interval)

    lines!(x,Data1[:,2],linewidth = width,label = "DMSO")
    lines!(x,Data2[:,2],linewidth = width,label = "10 mM KHCO2")
    lines!(x,Data3[:,2],linewidth = width,label = "10 mM KHCO2 Post Photolysis")

    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

Fig1()

save("JULIA/Outputs/Fig1.png",Fig1())

function Fig2(f = Figure())
    interval = 1500:8:2400
    Data1 = Interpolate2(FTIR_KHCO2_Phtolysis,2,interval)
    Data2 = Interpolate2(FTIR_KHCO2_Phtolysis,3,interval)
    Data3 = Interpolate2(FTIR_KHCO2_Phtolysis,4,interval)


    #Make Cairo Mackie figure
    font = :Arial
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "Steady State Photolysis FTIR of 10 mM KHCO₂ in DMSO, λₑ = 355nm", palette = (color = palette(ColorSchemes.winter, 3),),
    xlabel = "Wavenumber (cm⁻¹)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "T", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true, xautolimitmargin = (0.0, 0.0))
    # CairoMakie.xlims!(1250,1790)
    # CairoMakie.ylims!(nothing,)
    width = 2
    x = collect(interval)

    lines!(x,Norm(x,Data1[:,2],1800),linewidth = width,label = "0 sec")
    lines!(x,Norm(x,Data2[:,2],1800),linewidth = width,label = "60 sec")
    lines!(x,Norm(x,Data3[:,2],1800),linewidth = width,label = "90 sec")

    # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rb,nbanks = 1,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end

Fig2()
save("JULIA/Outputs/Fig2.png",Fig2())

### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = ""

    ## FITTING FUNCTIONS

    function MonoFit(file,column,t0,uppertimebound, residuals, ReturnFit)

        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
        p0 = [0.01, 11.0, 0.0]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values))
        
    
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
    
            lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)
    
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
    
    function BiFit(file,column,t0,uppertimebound, residuals, ReturnFit)

        Data = file
        t = Data[:,1]
        I = Data[:,column]
    
        # Create Bit-vectors for desrired time range and wavelength range for integration
        cut = t0 .< t .< uppertimebound
    
        time = t[cut]
        KineticTrace = I[cut]
    
        ## FITTING FUNCTION
        fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
        p0 = [0.01, 9.0, 0.02, 15.0, 0.0]
    
        peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,KineticTrace,fitfunction(values))
        
    
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
    
            lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)
    
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
    
        res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
        print("Fit Resuluts:", res)
        values = Optim.minimizer(res) # Kinetics results
        print("Fit Coeff:", values)
    
        if ReturnFit == true
            return hcat(time,fitfunction(values))
        
    
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
    
            lines!(Fit_Fig[2,1],time, KineticTrace .- fitfunction(values), linewidth = width, color = :black)
    
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
    #LinearInterpolationf
    function Interpolate(x,y,interval)

        range = collect(interval)
        interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = interp(range) # compute y values for the x vector
    
        Out = hcat(range,y_interp)

        return Out
    end


    function Interpolate2(file, column, interval)

        range = collect(interval)
        sorDF = sort!(file, rev =false)
        interp = LinearInterpolation(collect(sorDF[:,1]),collect(sorDF[:,column])) #interpolate function for x and y vectors
        y_interp = interp(range) # compute y values for the x vector
    
        Out = hcat(range,y_interp)

        return Out
    end

    #Normalization to x=z
    function Norm(x,y,z)
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/x_interp(z) #normalize to desired value
        return y_norm
    end
    # Normalize to Maximum
    function MaxNorm(x,y)
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/maximum(y_interp) #normalize to desired value
        return y_norm
    end

