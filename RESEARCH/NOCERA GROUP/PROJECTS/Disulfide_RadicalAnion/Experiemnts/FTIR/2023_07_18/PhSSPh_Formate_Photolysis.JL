using CSV, DataFrames, ColorSchemes,CairoMakie; CairoMakie.activate!()


KHCO2 = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_18_FTIR_PhssPh_KHCO2_DMSO_Photolysis/KHCO2_20mM_DMSO_390um_abs.0.dpt",header = false,DataFrame)
KHCO2_PhssPh = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_18_FTIR_PhssPh_KHCO2_DMSO_Photolysis/KHCO2_10mM_PhssPh_2p5mM_DMSO_390um_Abs.0.dpt",header = false,DataFrame) 
KHCO2_PhssPh_370nm = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_18_FTIR_PhssPh_KHCO2_DMSO_Photolysis/KHCO2_10mM_PhssPh_2p5mM_DMSO_390um_post370_1min_Abs.0.dpt",header = false,DataFrame)

PhssPh = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_18_FTIR_PhssPh_KHCO2_DMSO_Photolysis/PhssPh_5mM_DMSO_390um.1.dpt",header = false,DataFrame)
PhssPh_370nm = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_18_FTIR_PhssPh_KHCO2_DMSO_Photolysis/PhssPh_5mM_DMSO_390um_post370.0.dpt",header = false,DataFrame)



#Make Cairo Mackie figure

f = Figure(font = font,figure_padding = 25, fontsize = 20)
width = 3
ax = CairoMakie.Axis(f[1,1],
xlabel = "Wavenumber (cm⁻¹)", xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
ylabel =  "A", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
# CairoMakie.xlims!(1500,2400)
# hidespines!(ax,:r)

x = KHCO2[:,1]
# lines!(x,(KHCO2[:,2]),linewidth = width,label = "KHCO2",color = :blue)
lines!(ax,x, KHCO2_PhssPh[:,2],linewidth = width,label = "KHCO2+Ph-SS-Ph",color = :blue)
lines!(ax,x, KHCO2_PhssPh_370nm[:,2],linewidth = width,label = "KHCO2+Ph-SS-Ph post photolysis",color = :red)

lines!(ax,x, PhssPh[:,2],linewidth = width,label = "Ph-SS-Ph",color = :green)
lines!(ax,x, PhssPh_370nm[:,2],linewidth = width,label = "Ph-SS-Ph post photolysis",color = :orange)

hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

axislegend(ax,position = :rt,nbanks = 1,framevisible = false)


# ax2 = CairoMakie.Axis(f[2,1], titlefont = font,
# xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
# ylabel =  "A", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
# # CairoMakie.xlims!(1500,2400)
# # hidespines!(ax,:r)

# x = KHCO2[:,1]
# # lines!(x,(KHCO2[:,2]),linewidth = width,label = "KHCO2",color = :blue)
# lines!(ax2,x, PhssPh[:,2],linewidth = width,label = "Ph-SS-Ph",color = :blue)
# lines!(ax2,x, PhssPh_370nm[:,2],linewidth = width,label = "Ph-SS-Ph post photolysis",color = :red)

# hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

# axislegend(ax2,position = :rt,nbanks = 1,framevisible = false)


# colsize!(f.layout,1,Aspect(1, 1.25))
#Set aspect ration of the y vs. x axis)
resize_to_layout!(f)
f
save("JULIA/Outputs/Figure.png",f, pix_per_unit = 32)






f = Figure(font = font,figure_padding = 25, fontsize = 20)
width = 3
ax = CairoMakie.Axis(f[1,1],
xlabel = "Wavenumber (cm⁻¹)", xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
ylabel =  "A", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
# CairoMakie.xlims!(1500,2400)
# hidespines!(ax,:r)

x = KHCO2[:,1]
# lines!(x,(KHCO2[:,2]),linewidth = width,label = "KHCO2",color = :blue)
lines!(ax,x, KHCO2_PhssPh_370nm[:,2]-KHCO2_PhssPh[:,2],linewidth = width,label = "KHCO2+Ph-SS-Ph",color = :blue)
# lines!(ax,x, KHCO2_PhssPh_370nm[:,2],linewidth = width,label = "KHCO2+Ph-SS-Ph post photolysis",color = :red)

lines!(ax,x, PhssPh_370nm[:,2]-PhssPh[:,2],linewidth = width,label = "Ph-SS-Ph",color = :green)
# lines!(ax,x, PhssPh_370nm[:,2],linewidth = width,label = "Ph-SS-Ph post photolysis",color = :orange)

hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

axislegend(ax,position = :rt,nbanks = 1,framevisible = false)


# ax2 = CairoMakie.Axis(f[2,1], titlefont = font,
# xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
# ylabel =  "A", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
# # CairoMakie.xlims!(1500,2400)
# # hidespines!(ax,:r)

# x = KHCO2[:,1]
# # lines!(x,(KHCO2[:,2]),linewidth = width,label = "KHCO2",color = :blue)
# lines!(ax2,x, PhssPh[:,2],linewidth = width,label = "Ph-SS-Ph",color = :blue)
# lines!(ax2,x, PhssPh_370nm[:,2],linewidth = width,label = "Ph-SS-Ph post photolysis",color = :red)

# hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

# axislegend(ax2,position = :rt,nbanks = 1,framevisible = false)


# colsize!(f.layout,1,Aspect(1, 1.25))
#Set aspect ration of the y vs. x axis)
resize_to_layout!(f)
f
save("JULIA/Outputs/Figure.png",f, pix_per_unit = 32)