using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

PhSSPh_20mM_475nm_Kinetics_NoFlow = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_noflow_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_025mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_0p25mlpmin_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_05mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_0p5mlpmin_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_1mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_1mlpmin_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_2mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_2mlpmin_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_3mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_3mlpmin_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_4mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_4mlpmin_18us_Converted.csv", header = true, DataFrame)
PhSSPh_20mM_475nm_Kinetics_5mlmin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_31_NS_TA_PhSSPh_DMSO_Flow_Kinetics/PhSSPh_20mM_DMSO_475nm_TAKinetics_5mlpmin_18us_Converted.csv", header = true, DataFrame)


### Make Figures 

function Fig1(f = Figure())
    Data1 = PhSSPh_20mM_475nm_Kinetics_NoFlow 
    Data2 = PhSSPh_20mM_475nm_Kinetics_025mlmin
    Data3 = PhSSPh_20mM_475nm_Kinetics_05mlmin
    Data4 = PhSSPh_20mM_475nm_Kinetics_1mlmin 
    Data5 = PhSSPh_20mM_475nm_Kinetics_2mlmin 
    Data6 = PhSSPh_20mM_475nm_Kinetics_3mlmin
    Data7 = PhSSPh_20mM_475nm_Kinetics_4mlmin
    Data8 = PhSSPh_20mM_475nm_Kinetics_5mlmin

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)

    #Make Cairo Mackie figure
    
    ax = CairoMakie.Axis(f[1,1], titlefont = :bold, title = "TA Kinetics of Ph-SS-Ph (20mM) in DMSO at various Flowrates", titlesize = 20,
    xlabel = "time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Norm ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(0,15)
    CairoMakie.ylims!(-0.05,nothing)
    width = 4
    x1 = Data1[:,1]
    x2 = Data2[:,1]
    x3 = Data3[:,1]
    x4 = Data4[:,1]
    x5 = Data5[:,1]
    x6 = Data6[:,1]
    x7 = Data7[:,1]
    x8 = Data8[:,1]

    lines!(x1,MaxNorm(x1,Data1[:,2]),linewidth = width,label = "No Flow")
    lines!(x2,MaxNorm(x2,Data2[:,2]),linewidth = width,label = "0.25")
    lines!(x3,MaxNorm(x3,Data3[:,2]),linewidth = width,label = "0.5")
    lines!(x4,MaxNorm(x4,Data4[:,2]),linewidth = width,label = "1")
    lines!(x5,MaxNorm(x5,Data5[:,2]),linewidth = width,label = "2")
    lines!(x6,MaxNorm(x6,Data6[:,2]),linewidth = width,label = "3")
    lines!(x7,MaxNorm(x7,Data7[:,2]),linewidth = width,label = "4")
    lines!(x8,MaxNorm(x8,Data8[:,2]),linewidth = width,label = "5")


    axislegend("475 nm", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    
    # Label(f[1,1,TopLeft()],"(a)")
    # colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    f
end


Fig1()

save("JULIA/Outputs/Fig1.png",Fig1())


### OPTIMIZER USING OPTIM 
    # Kinetics Fitting using Optim

    file = "/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_27_NS_TA_PhSSPh_KHCO2_DMSO_flow/PhSSPh_20mM_KHCO2_5mM_455nm_45us_TA_Kinetics_100avg_Converted.csv"
    Data = (CSV.read(file,header = true,DataFrame))

    Data = PhSSPh_20mM_KHCO2_5mM_455nm_Kinetics

    t = Data[:,1]
    I = Data[:,2]

    # Set start and end time for dataset
    uppertimebound = 41
    t0 = 0.5

    # Create Bit-vectors for desrired time range and wavelength range for integration
    cut = t0 .< t .< uppertimebound

    time = t[cut]
    KineticTrace = I[cut]


    ## FITTING FUNCTIONS

    # fitfunction = mexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3]
    # p0 = [0.01, 11.0, 0.0]

    fitfunction = bexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] 
    p0 = [0.01, 9.0, 0.02, 15.0, 0.0]

    # fitfunction = triexp(p) = p[1] .* exp.((-time) ./ p[2]) .+ p[3] .* exp.((-time) ./ p[4]) .+ p[5] .* exp.((-time) ./ p[6]).+ p[7] 
    # p0 = [0.05, 1.19, 0.04, 8.7, 0.03, 39.0, 0.011]

    Plots.plot(time, KineticTrace, xaxis = :log)
    Plots.plot!(time, fitfunction(p0))

    peakfit0(p) = sum((KineticTrace .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,BFGS(),Optim.Options(time_limit = 2000))
    values = Optim.minimizer(res) # Kinetics results
    mexp(values)
    Plots.plot!(time, fitfunction(values), linewidth = 3, color = :black)

    # Fit Plots  
    Fit_Fig = Figure(font = "", figure_padding = 25,fontsize =20, resolution = (600,600))
    width = 2

    ax1 = CairoMakie.Axis(Fit_Fig[1,1], title = "Ph-SS-Ph (20 mM) + KHCO₂ (5 mM) TA Kinetics", #xscale = log10,
    xlabel = "time",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "A.U.", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # xlims!(ax1,350,650)
    # ylims!(nothing,)

    lines!(Fit_Fig[1,1],time, KineticTrace,linewidth = width, label = "455 nm", color = :blue)
    lines!(Fit_Fig[1,1],time, fitfunction(values), linewidth = 2, label = "τ₁ = 2.9 μs, τ₂ = 9.7 μs", color = :red)
    # lines!(Fit_Fig[1,1],time,fitfunction(p0), linewidth = 2, label = "Guess", color = :green)
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

    save("JULIA/Outputs/Fit_Fig.png",Fit_Fig)

# SOME FUNCTIONS

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

    function GatherAllignedPoints(coordinate_sets::Vector{Vector{Vector{Float64}}},x_interval)
        interps::Vector{Any} = []
        out::Vector{Vector{Float64}} = []
    
        for v in coordinate_sets
            x,y = v
            push!(interps,LinearInterpolation(x,y))
        end
    
        for i in interps
            push!(out,i(collect(x_interval)))
        end
    
        return out
    end
    