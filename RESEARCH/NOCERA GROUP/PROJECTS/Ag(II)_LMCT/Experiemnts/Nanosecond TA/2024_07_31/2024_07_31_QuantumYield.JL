## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

function Flux(vector,Wavelength)

    h = 6.62607015e-34
    c = 2.99792458e+8
    NA = 6.023e23 
    
    E_photon = h*c/(Wavelength*10^(-9))

    Avg_energy = mean(vector) # J per pulse
    Avg_Power = Avg_energy*50 # power in Watts at 50 Hz rep rate
    Avg_Flux = Avg_Power/E_photon/NA # mol photons /s

end

function QY(ΔA, ϵ, volume, flux, time)

    Δn = ΔA/(ϵ)*volume
    Quantum_yield = Δn / (flux * time) * 100

    return Quantum_yield
end

function QY_method_2(k, ϵ, volume, init_flux, init_A)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value
    I_abs = init_flux*(1-10^(-init_A))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end



#### FerriOxalate #########################################################################################################

    ## FerriOxalate UvVis
    FeOxalate_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/Ferrioxalate_Absorbance__0__10-31-59-317.csv",header = false, DataFrame))
    FeOxalate_UvVis_Kinetics_400nm = hcat(FeOxalate_UvVis[5:end,1], FeOxalate_UvVis[5:end,134])

    FeOxalate_UvVis[1,134]

    FeOxalate_Absorbances = transpose(FeOxalate_UvVis)

    Integer(round(size(FeOxalate_Absorbances[2:end,2:end],2)/11, digits = 0))

    Plots.plot(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,2:end])
    Plots.plot(FeOxalate_UvVis[15:end,1],FeOxalate_UvVis[15:end,134])

    ## Per-pulse energy after the Blank
    Post_stage_wBlank_FeOxalate = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/afterstage-blank.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_wBlank_FeOxalate[:,1],Post_stage_wBlank_FeOxalate[:,2])

    Post_stage_wBlank_FeOxalate_avg = mean(Post_stage_wBlank_FeOxalate[:,2])



    ## Per-pulse energy after the Sample for 20000 shots
    Post_stage_Sample_20000_FeOxalate = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-009/afterstage-sample-20000.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_20000_FeOxalate[:,1],Post_stage_Sample_20000_FeOxalate[:,2])

    ## calculated absorbance from the ΔP 
    Post_stage_Sample_20000_FeOxalate_abs = hcat(Post_stage_Sample_20000_FeOxalate[:,1], - log.(Post_stage_Sample_20000_FeOxalate[:,2] ./Post_stage_wBlank_FeOxalate_avg))
    mean(Post_stage_Sample_20000_FeOxalate_abs[1:10,2])
    Plots.plot(Post_stage_Sample_20000_FeOxalate_abs[:,1],Post_stage_Sample_20000_FeOxalate_abs[:,2])



    ## Fitting for rate constant
    FeOxalate_Fit = MonoFitIRF_New(FeOxalate_UvVis_Kinetics_400nm,2,0.0, 400, [0.08, 0.0001, 132.659, 210.211, 0.031], 10,Inf,false,10000)
    k = 1/210.211

    ## Calculation of the Flux and Quantum yields
    FLux_Blank = Flux(Post_stage_wBlank_FeOxalate_avg,436)

    QY_FeOxalate = QY_method_2(k, 50.26, 0.002,FLux_Blank,0.15078)


    function FeOx_Abs_decay(f = Figure())
        width = 3
        #Make Cairo Mackie figure
            ax = CairoMakie.Axis(f[1,1], title = "Ferrioxalate, 3mM QY",
                palette = (color = palette(ColorSchemes.inferno, 7),), xautolimitmargin = (0.0, 0.0),
                xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
                ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA"
            )
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(nothing,600)
            # CairoMakie.ylims!(-0.001,nothing)
            
            
            lines!(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,2],linewidth = width, label = "1 s")
            lines!(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,100],linewidth = width, label = "100 s")
            lines!(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,200],linewidth = width, label = "200 s")
            lines!(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,300],linewidth = width, label = "300 s")
            lines!(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,400],linewidth = width, label = "400 s")
            lines!(FeOxalate_Absorbances[2:end,1],FeOxalate_Absorbances[2:end,422],linewidth = width, label = "end")



            CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :blue)


            # vlines!(ax, 436, linewidth = 1,color = :blue,linestyle = :dash)
            # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

            axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

            # text!(x,y,text = "text", color = :red, textsize = 20)
            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            resize_to_layout!(f)
       f
            
    end
    FeOx_Abs_decay()
    save("JULIA/Outputs/FeOx_Abs_decay.pdf",Fe_Ox_AbFeOx_Abs_decays_decay())


#### Ag(II) low power #########################################################################################################

    ## Ag(II) low power
    Ag_LP_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-010/500 uj/Absorbance__0__11-28-10-069.csv",header = false, DataFrame))
    Ag_LP_UvVis_Kinetics_400nm = hcat(Ag_LP_UvVis[5:end,1], Ag_LP_UvVis[5:end,134])

    Ag_LP_UvVis[1,134]
    Ag_LP_UvVis[1,197]

    Ag_LP_Absorbances = transpose(Ag_LP_UvVis)
    Ag_LP_Absorbances[197,2]
    Plots.plot(Ag_LP_Absorbances[2:end,1],Ag_LP_Absorbances[2:end,2:end])
    Plots.plot(Ag_LP_UvVis[15:end,1],Ag_LP_UvVis[15:end,134])

    ## Per-pulse energy after the Sample for 20000 shots
    Post_stage_Sample_20000_Ag_LP = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-010/500 uj/afterstage-sample-500uj.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_20000_Ag_LP[:,1],Post_stage_Sample_20000_Ag_LP[:,2])

    ## Per-pulse energy after the Blank
    Post_Blank_Ag_LP = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-010/500 uj/afterstage-blank-500uj.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_Blank_Ag_LP[:,1],Post_Blank_Ag_LP[:,2])

    Post_Blank_Ag_LP_avg = mean(Post_Blank_Ag_LP[:,2])

    ## Fitting for rate constant
    Ag_LP_Fit = MonoFitIRF_New(Ag_LP_UvVis_Kinetics_400nm,2,0.0, 800, [0.08, 0.001, 32.343, 404.806, 0.002], 10,Inf,false,10000)
    k_lp = 1/404.806

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_LP = Flux(Post_Blank_Ag_LP_avg,450)

    QY_Ag_LP = QY_method_2(k_lp, 1297, 0.003,FLux_Blank_LP,0.0641)

#### Ag(II) High power #########################################################################################################

    ## Ag(II) high power
    Ag_HP_UvVis = Matrix(CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-010/2500 uj/Absorbance__0__11-47-46-064_TrueTime.csv",header = false, DataFrame))
    Ag_HP_UvVis_Kinetics_400nm = hcat(Ag_HP_UvVis[5:end,1], Ag_HP_UvVis[5:end,134])

    Ag_HP_UvVis[1,134]
    Ag_HP_UvVis[1,197]

    Ag_HP_Absorbances = transpose(Ag_HP_UvVis)
    Ag_HP_Absorbances[197,2]
    Plots.plot(Ag_HP_Absorbances[2:end,1],Ag_HP_Absorbances[2:end,2:end])
    Plots.plot(Ag_HP_UvVis[15:end,1],Ag_HP_UvVis[15:end,134])

    ## Per-pulse energy after the Sample for 20000 shots
    Post_stage_Sample_20000_Ag_HP = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-010/2500 uj/afterstage-sample-2500uj.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_stage_Sample_20000_Ag_HP[:,1],Post_stage_Sample_20000_Ag_HP[:,2])

    ## Per-pulse energy after the Blank
    Post_Blank_Ag_HP = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Experiemnts/Nanosecond TA/2024_07_31/BC-4-010/2500 uj/afterstage-blank-2500uj.csv", skipto = 6, header = true, DataFrame)
    Plots.plot(Post_Blank_Ag_HP[:,1],Post_Blank_Ag_HP[:,2])

    Post_Blank_Ag_HP_avg = mean(Post_Blank_Ag_HP[:,2])

    ## Fitting for rate constant
    Ag_HP_Fit = MonoFitIRF_New(Ag_HP_UvVis_Kinetics_400nm,2,0.0, 400, [0.08, 0.0, 14.786, 118.806, 0.002], 10,Inf,false,10000)
    k_Hp = 1/118.806

    ## Calculation of the Flux and Quantum yields
    FLux_Blank_HP = Flux(Post_Blank_Ag_HP_avg,450)

    QY_Ag_HP = QY_method_2(k_Hp, 1297, 0.003,FLux_Blank_HP,0.0987)



    # I_abs_0p5 = FLux_Blank_LP*(1-10^(-Ag_LP_Absorbances[197,2]))/0.03
    # I_abs_1 = FLux_Blank_3*(1-10^(-0.115))/0.03
    # I_abs_1p5 = FLux_Blank_AgTFA_1p5mJ_425nm*(1-10^(-AgTFA_1p5mJ_425nm_initial_Abs))/0.03
    # I_abs_2p5 = FLux_Blank_HP*(1-10^(-Ag_HP_Absorbances[197,2]))/0.03

    # fit = fitlinear([I_abs_1,I_abs_1p5,I_abs_2p5],[1/306.866, 1/155.194, 1/118.806])
    # Plots.plot(fit.x,fit.y)
    # Plots.scatter!([I_abs_0p5,I_abs_1,I_abs_1p5,I_abs_2p5],[1/404.806,1/306.866, 1/155.194, 1/118.806])