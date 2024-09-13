## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Pre_stage = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/450nm-prestage-1000shots-Agbpy2TFA-s1-1.csv", skipto = 6, header = true, DataFrame)
Post_stage = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/450nm-afterstage-1000shots-Agbpy2TFA-s1.csv", skipto = 6, header = true, DataFrame)
Post_stage_wBlank = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/450nm-afterstage-withblank-1000shots-Agbpy2TFA-s1.csv", skipto = 6, header = true, DataFrame)

Post_stage_Sample_1000 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/450nm-afterstage-withSample-1000shots-Agbpy2TFA-s1.csv", skipto = 6, header = true, DataFrame)
Post_stage_Sample_2000 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/450nm-afterstage-withSample-2000shots-Agbpy2TFA-s1.csv", skipto = 6, header = true, DataFrame)
Post_stage_Sample_3000 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/450nm-afterstage-withSample-3000shots-Agbpy2TFA-s1.csv", skipto = 6, header = true, DataFrame)

t1 = Post_stage_Sample_1000[end,1]
t2 = Post_stage_Sample_2000[end,1] - 10.9046
t3 = Post_stage_Sample_3000[end,1]

Plots.plot(Post_stage_Sample_1000[:,1],Post_stage_Sample_1000[:,2])
Plots.plot!(Post_stage_Sample_2000[2:end,1] .+ t1 .- 11,Post_stage_Sample_2000[2:end,2])
Plots.plot!(Post_stage_Sample_3000[:,1] .+ t3 .+ t1,Post_stage_Sample_3000[:,2])


Uv_Vis_0 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/4-007/Agbpy2TFA-0p6mM-before_Absorbance__0__00.txt", skipto = 15,header = false, DataFrame)
Uv_Vis_1000 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/4-007/Agbpy2TFA-0p6mM-after-1000shots-2p4mJ_Absorbance__0__01.txt", skipto = 15,header = false, DataFrame)
Uv_Vis_2000 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/4-007/Agbpy2TFA-0p6mM-after-2000shots-2p4mJ_Absorbance__0__02.txt", skipto = 15,header = false, DataFrame)
Uv_Vis_3000 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/Summary/2024_07_23_Qunatum Yield/BC-4-007/4-007/Agbpy2TFA-0p6mM-after-3000shots-2p4mJ_Absorbance__0__03.txt", skipto = 15,header = false, DataFrame)


#### 2024_07_30

Pre_stage_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/QY/BC-4-008/prestage-1000shots.csv", skipto = 6, header = true, DataFrame)
Post_stage_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/QY/BC-4-008/Afterstage-1000shots.csv", skipto = 6, header = true, DataFrame)
Post_stage_wBlank_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/QY/BC-4-008/Afterstage-blank-1000shots.csv", skipto = 6, header = true, DataFrame)
Post_stage_wBlank_avg = mean(Post_stage_wBlank_2024_07_30[:,2])

Plots.plot(Pre_stage_2024_07_30[1:1000,1],Pre_stage_2024_07_30[1:1000,2])
Plots.plot!(Post_stage_2024_07_30[:,1],Post_stage_2024_07_30[:,2])
Plots.plot!(Post_stage_wBlank_2024_07_30[:,1],Post_stage_wBlank_2024_07_30[:,2])


Post_stage_Sample_20000_2024_07_30 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Ag(II)_LMCT/QY/BC-4-008/Afterstage-sample-20000shots.csv", skipto = 6, header = true, DataFrame)

Plots.plot!(Post_stage_Sample_20000_2024_07_30[:,1],Post_stage_Sample_20000_2024_07_30[:,2])


Post_stage_Sample_20000_abs = hcat(Post_stage_Sample_20000[:,1], - log.(Post_stage_Sample_20000[:,2] ./Post_stage_wBlank_avg))
mean(Post_stage_Sample_20000_abs[1:10,2])
Plots.plot(Post_stage_Sample_20000_abs[:,1],Post_stage_Sample_20000_abs[:,2])
Plots.plot!(fit[:,1],fit[:,2])



fit = MonoFitIRF_New(Post_stage_Sample_20000_abs,2,0.0, 400, [0.008, 0.0001, 67.642, 281.581, 0.001], 10,50,false,10000)
k = 1/281.581
k = 1/246.888


FLux_Blank = Flux(Post_stage_wBlank_avg,450)
Flux_Sample = Flux(mean(Post_stage_Sample_20000_2024_07_30[1:10,2]),450)
Initial_Flux_abs = FLux_Blank - Flux_Sample

Proj_abs = FLux_Blank * (1-10^(-0.1))

QY_linear_Sample_Flux_20000 = QY_method_2(k, 1297, 0.003,FLux_Blank,mean(Post_stage_Sample_20000_abs[1:10,2]))
QY_linear_Sample_Flux_20000_2 = QY_method_2(k, 1297, 0.003,FLux_Blank,0.087)


# Absorbance values after each run of exposures
Uv_Vis_0[196,2]
Uv_Vis_1000[196,2]
Uv_Vis_2000[196,2]
Uv_Vis_3000[196,2]

# CONSTANTS
h = 6.62607015e-34
c = 2.99e+8
NA = 6.023e23
ϵ_450 =  1297


function Flux(vector,Wavelength)

    h = 6.62607015e-34
    c = 2.99e+8
    NA = 6.023e23 
    
    E_photon = h*c/(Wavelength*10^(-9))

    Avg_energy = mean(vector) # J per pulse
    Avg_Power = Avg_energy*50 # power in Watts at 50 Hz rep rate
    Avg_Flux = Avg_Power/E_photon/NA # mol photons /s

end

Flux_pre = Flux(Pre_stage[:,2],450)
Flux_postwBlank = Flux(Post_stage_wBlank[:,2],450)

# INITAL FLux transmitted through Sample
ΔFlux = Flux_postwBlank - Flux(Post_stage_Sample_1000[1:10,2], 450) 

ΔFlux_1 = Flux_postwBlank - Flux(Post_stage_Sample_1000[1:10,2], 450) 


# FLux transmitted through Sample after exposures
ΔFlux_1000 = Flux_postwBlank - Flux(Post_stage_Sample_1000[:,2], 450) 
ΔFlux_2000 = Flux_postwBlank - Flux(Post_stage_Sample_2000[:,2], 450)
ΔFlux_3000 = Flux_postwBlank - Flux(Post_stage_Sample_3000[:,2], 450)

# Projected intensity absorbed based on sample O.D.
I_abs1 = Flux_postwBlank*(1-10^(-Uv_Vis_0[196,2]))
I_abs2= Flux_pre*(1-10^(-Uv_Vis_0[196,2]))


function QY(ΔA, ϵ, volume, flux, time)

    Δn = ΔA/(ϵ)*volume
    Quantum_yield = Δn / (flux * time) * 100

    return Quantum_yield
end

ΔA_1 = Uv_Vis_0[196,2] - Uv_Vis_1000[196,2]
ΔA_2 = Uv_Vis_1000[196,2] - Uv_Vis_2000[196,2]
ΔA_3 = Uv_Vis_2000[196,2] - Uv_Vis_3000[196,2]

QY_Flux_pre = QY(ΔA_1,1297, 0.003,Flux_postwBlank, t1)
QY_Flux_postBlank = QY(ΔA_1,1297, 0.003,Flux_postwBlank, t1)

## These are the only feasible ones for this method because need to consider flux absorbed 
QY_InitSample_Flux = QY(ΔA_1,1297, 0.003,ΔFlux, t1)
QY_linear_Sample_Flux_1000 = QY(ΔA_1,1297, 0.003,ΔFlux_1000, t1)
QY_linear_Sample_Flux_1000 = QY(ΔA_2,1297, 0.003,ΔFlux_2000, t2)
QY_linear_Sample_Flux_1000 = QY(ΔA_3,1297, 0.003,ΔFlux_3000, t3)


function QY_method_2(k, ϵ, volume, init_flux, init_A)
    # here the flux is the incident flux, the absorbed Flux is calculated from the initial Absorbance value
    I_abs = init_flux*(1-10^(-init_A))/volume

    Quantum_yield = k * init_A / (I_abs * ϵ) * 100

    return Quantum_yield
end

fit_linear = fitlinear([0, t1,t1+t2,t1+t2+t3],[Uv_Vis_0[196,2],  Uv_Vis_1000[196,2],Uv_Vis_2000[196,2],Uv_Vis_3000[196,2]])
fit_exp = fitexp([0, t1,t1+t2,t1+t2+t3],[Uv_Vis_0[196,2],  Uv_Vis_1000[196,2],Uv_Vis_2000[196,2],Uv_Vis_3000[196,2]],n =1)

Plots.scatter([0, t1,t1+t2,t1+t2+t3],[Uv_Vis_0[196,2],  Uv_Vis_1000[196,2],Uv_Vis_2000[196,2],Uv_Vis_3000[196,2]])
Plots.plot!(fit_linear.x,fit_linear.y)
Plots.plot!(fit_exp.x,fit_exp.y)


## Rate constants determined by Linear fit of A vs time
QY_linear_Flux_pre =  QY_method_2(abs(fit_linear.a), 1297, 0.003,Flux_pre,Uv_Vis_0[196,2])
QY_linear_Flux_postBlank = QY_method_2(abs(fit_linear.a), 1297, 0.003,Flux_postwBlank,Uv_Vis_0[196,2])


# Inital Absorbance value taken to be the O.D. of the previous exposure for subsequent runs
QY_linear_Sample_Flux_2000 = QY_method_2(abs(fit_linear.a), 1297, 0.003,Flux_postwBlank,Uv_Vis_1000[196,2])
QY_linear_Sample_Flux_3000 = QY_method_2(abs(fit_linear.a), 1297, 0.003,Flux_postwBlank,Uv_Vis_2000[196,2]) 


## Rate constants determined by exponential fit of A vs time
QY_Exp_Flux_pre = QY_method_2(abs(1/fit_exp.b), 1297, 0.003,Flux_pre,Uv_Vis_0[196,2])
QY_Exp_Flux_postBlank = QY_method_2(abs(1/fit_exp.b), 1297, 0.003,Flux_postwBlank,Uv_Vis_0[196,2])

# Inital Absorbance value taken to be the O.D. of the previous exposure for subsequent runs
QY_Exp_Sample_Flux_2000 = QY_method_2(abs(1/fit_exp.b), 1297, 0.003,Flux_postwBlank,Uv_Vis_1000[196,2])
QY_Exp_Sample_Flux_3000 = QY_method_2(abs(1/fit_exp.b), 1297, 0.003,Flux_postwBlank,Uv_Vis_2000[196,2])


function Abs_decay(f = Figure())
    width = 3
    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Nitronbenzene + Phosphine TA Spectra in CME",
            palette = (color = palette(ColorSchemes.inferno, 7),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength (nm)", 
            ylabel = "ΔA"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(395,850)
        # CairoMakie.ylims!(-0.001,nothing)
        
        lines!(Uv_Vis_0[:,1],Uv_Vis_0[:,2],linewidth = width,label = "0 shots")
        lines!(Uv_Vis_1000[:,1],Uv_Vis_1000[:,2],linewidth = width,label = "1000 shots")
        lines!(Uv_Vis_2000[:,1],Uv_Vis_2000[:,2],linewidth = width,label = "2000 shots")
        lines!(Uv_Vis_3000[:,1],Uv_Vis_3000[:,2],linewidth = width,label = "3000 shots")



        vlines!(ax, 450, linewidth = 1,color = :black,linestyle = :dash)
       # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        axislegend("Exposures", position = :rt,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 
        # ax2 = CairoMakie.Axis(f[1,2], 
        #     palette = (color = palette(ColorSchemes.cool, 5),), xautolimitmargin = (0.0, 0.0),
        #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        #     xlabel = "Wavenumber", 
        #     ylabel = "ΔA")
        # linkxaxes!(ax, ax2)

        # lines!(ax2,Ph_SS_Ph_Kinetics_450nm[:,1],(Ph_SS_Ph_Kinetics_450nm[:,2] -Ph_SS_Ph_Kinetics_450nm[:,3]), color = :green )
        # Label(f[1,1,TopLeft()],"(a)", textsize = 20)

        # rowsize!(f.layout,1, 400)
        # colgap!(f.layout,1,10)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Abs_decay()
save("JULIA/Outputs/Abs_decay.pdf",Abs_decay())