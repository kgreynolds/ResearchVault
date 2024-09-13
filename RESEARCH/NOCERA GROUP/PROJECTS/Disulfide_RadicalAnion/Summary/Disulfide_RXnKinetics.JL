using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,MAT,LineSearches,DifferentialEquations,CairoMakie; CairoMakie.activate!()
plotly()

### INITIAL CONDITIONS
PC_ES₀ = 1 * 10^(-4)       #Photocatalyst excited state concentration after laser pulse
RSSR = 2 * 10^(-3)         # Disulfide starting conentration
HCO2 = 10 * 10^(-3)         # Formate
RSSR_EC₀ = 0.0              # Disulfide excited state
RS_Dot₀ = 0.0               # Thyil Radical 
CO2_RA = 0.0                # CO2 radical anion
RSSR_RA  = 0.0              # Disulfide radical anion
Thiolate = 0.0              # Thiolate
Sub = 0 *10^(-3)            # Substrate
CO2 = 0.0                   # CO2
Ox = 0.0                    # Oxalate

### RATE CONSTANTS
kEN = 2.7 * 10^(9)          # Rate of Energy transfer from PC to RSSR
kDec = 1 / (176*10^(-9))    # Rate of PC instrinsic decay
kDis1 = 1 * 10^(12)         # Rate of dissociation of disulfide excited state
kDis2 = 1 * 10^(6)          # Rate of dissociation of disulfide radical anion
knr = 0 * 10^(10)           # rate of non radiative relaxation of disulfide excited state
kGc = 1 * 10^(10)           # Rate og deminal radical recombination of thyil radical
kHAT = 1 * 10^(10)          # Rate of HAT from Formate
kET1 = 1*10^(10)            # Rate of ET 1
kET2 = 1*10^(10)            # Rate of ET 2
kET3 = 1*10^(10)            # Rate of ET 3
kET4 = 0*10^(10)            # Rate of ET 4
kET5 = 1*10^(10)            # Rate of ET 5

kDim = 1*10^(10)
 


# function Reacton!(du,u,p,t)
#     kEN,kDec,kDis1,kGc,kHAT,kET1,kET2,kET3,kDis2, kDim = p
#     du[1] = -kEN* u[1]*u[2] -kET2*u[6]*u[1] + kGc*((u[4])^2)                                                    # Dislfide
#     du[2] = -kEN* u[1]*u[2] - kDec * u[2]                                                                       # Photocatalyst excited state
#     du[3] = -kDis1 * u[3] + kEN *u[1]*u[2]                                                                      # Disulfide excited state
#     du[4] = -kGc*((u[4])^2)-kHAT*u[4]*u[5] - kET1*u[4]*u[6]  + 2*kDis1*u[3]+kDis2*u[7]                          # thyil radical 
#     du[5] = -kHAT*u[4]*u[5]                                                                                     # formate
#     du[6] = -kET1*u[6]*u[4] -kET2*u[6]*u[1] -kET3*u[6]*u[9] -kDim*((u[6])^2) + kHAT*u[4]*u[5]                   # CO2 radical anion
#     du[7] = -kDis2*u[7] + kET2*u[6]*u[1]                                                                        # disulfide radical anion
#     du[8] = kET1*u[4]*u[6]                                                                                      # Thiolate
#     du[9] = -kET3*u[9]*u[6]                                                                                     # Substrate
#     du[10] = kET1*u[6]*u[4] + kET2*u[6]*u[1] +kET3*u[6]*u[9]                                                    # CO2
#     du[11] = kDim*((u[6])^2)                                                                                    # Oxalate
# end

function Reacton!(du,u,p,t)
    kEN,kDec,kDis1,kGc,kHAT,kET1,kET2,kET3,kET4,kET5,kDis2, kDim = p
    du[1] = -kEN* u[1]*u[2] -kET2*u[6]*u[1] + kGc*((u[4])^2) +kET4*u[9]*u[7]+kET5*u[7]*u[4]                     # Dislfide
    du[2] = -kEN* u[1]*u[2] - kDec * u[2]                                                                       # Photocatalyst excited state
    du[3] = -kDis1 * u[3] + kEN *u[1]*u[2]                                                                      # Disulfide excited state
    du[4] = -kGc*((u[4])^2)-kHAT*u[4]*u[5] - kET1*u[4]*u[6] -kET5*u[7]*u[4]  + 2*kDis1*u[3]+kDis2*u[7]          # thyil radical 
    du[5] = -kHAT*u[4]*u[5]                                                                                     # formate
    du[6] = -kET1*u[6]*u[4] -kET2*u[6]*u[1] -kET3*u[6]*u[9] -kDim*((u[6])^2) + kHAT*u[4]*u[5]                   # CO2 radical anion
    du[7] = -kDis2*u[7]-kET4*u[9]*u[7] -kET5*u[7]*u[4] + kET2*u[6]*u[1]                                         # disulfide radical anion
    du[8] = kET1*u[4]*u[6]+kET5*u[7]*u[4]                                                                       # Thiolate
    du[9] = -kET3*u[9]*u[6] -kET4*u[9]*u[7]                                                                     # Substrate
    du[10] = kET1*u[6]*u[4] + kET2*u[6]*u[1] +kET3*u[6]*u[9]                                                    # CO2
    du[11] = kDim*((u[6])^2)                                                                                    # Oxalate
end

    u0 = [RSSR,PC_ES₀,RSSR_EC₀, RS_Dot₀,HCO2,CO2_RA,RSSR_RA,Thiolate, Sub, CO2, Ox]

    p = (kEN,kDec,kDis1,kGc,kHAT,kET1,kET2,kET3,kET4,kET5,kDis2, kDim) # we could also make this an array, or any other sequence type!
    tspan = (0.0,100*10^(-6))
    prob = ODEProblem(Reacton!,u0,tspan,p)
    sol = solve(prob)

    matrix = sol.u
    newsol = stack(sol,dims=1).* 1000 #converts vector of vectors in sol to a matrix


# function Reacton!(du,u,p,t)
#     kEN,kDec,kDis1,kGc = p
#     du[1] = -kEN* u[1]*u[2]+ kGc*((u[4])^2)                                                                     # Dislfide
#     du[2] = -kEN* u[1]*u[2] - kDec * u[2]                                                                       # Photocatalyst excited state
#     du[3] = -kDis1 * u[3] + kEN *u[1]*u[2]                                                                      # Disulfide excited state
#     du[4] = -kGc*((u[4])^2)+ 2*kDis1*u[3]                                                                       # thyil radical 

# end

#     u0 = [RSSR,PC_ES₀,RSSR_EC₀, RS_Dot₀]

#     p = (kEN,kDec,kDis1,kGc) # we could also make this an array, or any other sequence type!
#     tspan = (0.0,5*10^(-6))
#     prob = ODEProblem(Reacton!,u0,tspan,p)
#     sol = solve(prob)

#     matrix = sol.u
#     newsol = stack(sol,dims=1).* 1000 #converts vector of vectors in sol to a matrix

### Make Figures 
function Fig1(f = Figure())
    width = 3
    x = sol.t *10^6



    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], 
            palette = (color = palette(ColorSchemes.hsv, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (us)", 
            ylabel = "[C] (mM)"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(0,2)
        # CairoMakie.ylims!(-0.000005,0.0001)
        

        lines!(ax,x,newsol[:,1],linewidth = width,label = "RSSR", color = :black)
        lines!(ax,x,newsol[:,5],linewidth = width,label = "HCO2", color = :lightblue)
        lines!(ax,x,newsol[:,10],linewidth = width,label = "CO2", color = :green)
        # lines!(ax,x,newsol[:,9],linewidth = width,label = "Substrate", color = :orange)





        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        # axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        f[1, 2] = Legend(f, ax, framevisible = false)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 

        ax2 = CairoMakie.Axis(f[2,1], 
            palette = (color = palette(ColorSchemes.hsv, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (us)", 
            ylabel = "[C] (mM)"
        )
        
        
        lines!(ax2,x,newsol[:,2],linewidth = width,label = "PC ES", color = :violet)
        # lines!(ax2,x,newsol[:,3],linewidth = width,label = "RSSR EC")
        # lines!(ax2,x,newsol[:,6],linewidth = width,label = "CO2 RA")
        # lines!(ax2,x,newsol[:,4],linewidth = width,label = "RS Dot")
        lines!(ax2,x,newsol[:,7],linewidth = width,label = "RSSR RA", color = :gold)


        # axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        f[2, 2] = Legend(f, ax2, framevisible = false)

        ax3 = CairoMakie.Axis(f[3,1], 
        palette = (color = palette(ColorSchemes.hsv,5),), xautolimitmargin = (0.0, 0.0),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
        xlabel = "Time (us)", 
        ylabel = "[C] (mM)"
    )
        lines!(ax3,x,newsol[:,4],linewidth = width,label = "RS Dot")
        lines!(ax3,x,newsol[:,6],linewidth = width,label = "CO2 RA")
        lines!(ax3,x,newsol[:,8],linewidth = width,label = "Thiolate")
        lines!(ax3,x,newsol[:,11],linewidth = width,label = "Oxalate")

                #  CairoMakie.ylims!(-0,0.003)

        # axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        f[3, 2] = Legend(f, ax3, framevisible = false)


    # Label(f[1,1,TopLeft()],"(a)", textsize = 20)
    linkxaxes!(ax, ax2)
    # rowsize!(f.layout,1, 100)
    # colgap!(f.layout,1,10)
    # colsize!(f.layout,1,Aspect(1, 4)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Fig1()
save("JULIA/Outputs/RSSR_RxN_2xHCO2_300us_v2.pdf",Fig1())

function Fig2(f = Figure())
    width = 3
    x = sol.t *10^6



    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], 
            palette = (color = palette(ColorSchemes.hsv, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (us)", 
            ylabel = "[C] (mM)"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(0,2)
        # CairoMakie.ylims!(-0.000005,0.0001)
        

        lines!(ax,x,newsol[:,1],linewidth = width,label = "RSSR", color = :black)
        # lines!(ax,x,newsol[:,2],linewidth = width,label = "PC", color = :lightblue)
        # lines!(ax,x,newsol[:,3],linewidth = width,label = "RS EC", color = :green)
        # lines!(ax,x,newsol[:,4],linewidth = width,label = "RS Dot", color = :orange)





        # hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

        # axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        f[1, 2] = Legend(f, ax, framevisible = false)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    # Second Panel of figure 

        ax2 = CairoMakie.Axis(f[2,1], 
            palette = (color = palette(ColorSchemes.hsv, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
            xlabel = "Time (us)", 
            ylabel = "[C] (mM)"
        )
        
        
        lines!(ax2,x,newsol[:,2],linewidth = width,label = "PC")
        # lines!(ax2,x,newsol[:,3],linewidth = width,label = "RS EC")
        lines!(ax2,x,newsol[:,4],linewidth = width,label = "RS Dot")


        # axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
        f[2, 2] = Legend(f, ax2, framevisible = false)

    #     ax3 = CairoMakie.Axis(f[3,1], 
    #     palette = (color = palette(ColorSchemes.hsv,5),), xautolimitmargin = (0.0, 0.0),
    #     xgridvisible = false, xminorticksvisible = true, xtickalign = 0, xticksize = 10, xlabelsize = 20,
    #     ygridvisible = false, yminorticksvisible = true, ytickalign = 0, yticksize = 10, ylabelsize = 20,
    #     xlabel = "Time (us)", 
    #     ylabel = "[C] (mM)"
    # )
    #     lines!(ax3,x,newsol[:,4],linewidth = width,label = "RS Dot")
    #     lines!(ax3,x,newsol[:,6],linewidth = width,label = "CO2 RA")
    #     lines!(ax3,x,newsol[:,8],linewidth = width,label = "Thiolate")
    #     lines!(ax3,x,newsol[:,11],linewidth = width,label = "Oxalate")

    #             #  CairoMakie.ylims!(-0,0.003)

    #     # axislegend( position = :rt,nbanks = 2,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)
    #     f[3, 2] = Legend(f, ax3, framevisible = false)


    # Label(f[1,1,TopLeft()],"(a)", textsize = 20)
    # linkxaxes!(ax, ax2)
    # rowsize!(f.layout,1, 100)
    # colgap!(f.layout,1,10)
    # colsize!(f.layout,1,Aspect(1, 4)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Fig2()
save("JULIA/Outputs/RSSROnly_EnT_5us.pdf",Fig2())

# function Reacton!(du,u,p,t)
#     kEN,kDec,kDis1,kGc = p
#     # du[1] = -kEN* u[1] - kDec * u[1] # Photocatalyst excited state
#     du[1] = -kEN* u[1] - kDec * u[1] # Photocatalyst excited state
#     du[2] = -kDis1 * u[2] + kEN * u[1] #Disulfide excited state
#     du[3] = -kGc * u[3] + 2*kDis1 * u[2] #Disulfide radical 
# end


# u0 = [PC_ES₀,RSSR₀, RS_Dot₀]

# p = (kEN,kDec,kDis1,kGc) # we could also make this an array, or any other sequence type!
# tspan = (0.0,5*10^(-8))
# prob = ODEProblem(Reacton!,u0,tspan,p)
# sol = solve(prob)

# Plots.plot(sol)