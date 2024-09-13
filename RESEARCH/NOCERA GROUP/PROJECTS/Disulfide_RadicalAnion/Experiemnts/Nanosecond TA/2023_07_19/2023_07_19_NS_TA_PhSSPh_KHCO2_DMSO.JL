using Interpolations,Plots, DataFrames, CSV, Tables, Statistics, EasyFit,CairoMakie; CairoMakie.activate!()
### Moving Average Calculator, Used for nanosecond TA Kinetics
## Update path to file and appropriate Header info (NOTE: make sure the right rows and collums are dropped depending on input file structure)

plotly()

#Normalization to x=z
function Norm(x,y,z)
    x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
    y_interp = x_interp(x) # compute y values for the x vector
    y_norm = y_interp/x_interp(z) #normalize to desired value
    return y_norm
end

function MaxNorm(x,y)
    x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
    y_interp = x_interp(x) # compute y values for the x vector
    y_norm = y_interp/maximum(y_interp) #normalize to desired value
    return y_norm
end


PhssPh_450_180 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_450nm_Kinetics_180us_COnverted.csv", DataFrame)
Plots.plot(PhssPh_450_180[:,1],PhssPh_450_180[:,2], xlimits = (-10,180))



HCO2_PhssPh_1800 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_KHCO2_450nm_Kinetics_1800us_100avg_Converted.csv", DataFrame)
Plots.plot!(HCO2_PhssPh_1800[:,1],HCO2_PhssPh_1800[:,2], xlimits = (-10,180))

# Normalizesd TA Kinetics
PhssPh_450_180_norm = MaxNorm(PhssPh_450_180[:,1], PhssPh_450_180[:,2])
HCO2_PhssPh_norm = MaxNorm(HCO2_PhssPh_1800[:,1],HCO2_PhssPh_1800[:,2])

Plots.plot(PhssPh_450_180[:,1],PhssPh_450_180_norm, xlimits = (-10,180))
Plots.plot!(HCO2_PhssPh_1800[:,1], HCO2_PhssPh_norm, xlimits = (-10,180))

## TA Kinetics

    #Make Cairo Mackie figure
    font = :Arial
    f = Figure(font = font,figure_padding = 25,fontsize =20)
    ax1 = CairoMakie.Axis(f[1,1], titlefont = font, title = "450 nm TA Kinetics",
    xlabel = "Time (μs)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "Normalized ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(ax1,-10,180)
    CairoMakie.ylims!(-0.1, nothing)
    width = 3

    lines!(ax1, PhssPh_450_180[:,1],PhssPh_450_180_norm,linewidth = width,label = "Ph-SS-Ph",color = :blue)
    lines!(ax1, HCO2_PhssPh_1800[:,1],HCO2_PhssPh_norm,linewidth = width,label = "Ph-SS-Ph + KHCO₂",color = :red)


    hlines!(ax1,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend( position = :rt,nbanks = 1,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f

    save("JULIA/Outputs/2023_07_19_PhSSPh_TA_Kinetics.png",f, pix_per_unit = 8)


#TA spectra
    #Import Data to Plot
    PhSSPh_150ns_TA_Spectrum = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_DMSO_150nsdelay_100avg_TASpectrum.csv",DataFrame)
    PhSSPh_KHCO2_150ns_TA_Spectrum = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_19_NS_TA_PhSSPh_KHCO2_DMSO/Ph_SS_Ph_KHCO2_DMSO_150nsdelay_100avg_TASpectrum.csv",DataFrame)


    #Make Cairo Mackie figure
    f = Figure(figure_padding = 25,fontsize =20)
    ax = CairoMakie.Axis(f[1,1], title = " TA spectra at 150 ns",
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    CairoMakie.xlims!(ax,370,610)
    # ylims!(nothing,)
    width = 3
    x = PhSSPh_150ns_TA_Spectrum[:,1]
    x[750]
    baseline = 800

    lines!(PhSSPh_150ns_TA_Spectrum[:,1],PhSSPh_150ns_TA_Spectrum[:,2] .- PhSSPh_150ns_TA_Spectrum[750,2],linewidth = width,label = "Ph-SS-Ph",color = :red)
    lines!(PhSSPh_KHCO2_150ns_TA_Spectrum[:,1],PhSSPh_KHCO2_150ns_TA_Spectrum[:,2] .- PhSSPh_KHCO2_150ns_TA_Spectrum[750,2],linewidth = width,label = "Ph-SS-Ph + HCO₂",color = :blue)


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend(position = :lb,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
    save("JULIA/Outputs/2023_07_19_PhSSPh_TA_Spectra.png",f, pix_per_unit = 8)

