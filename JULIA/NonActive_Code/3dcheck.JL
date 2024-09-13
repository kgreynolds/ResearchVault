# quick 3-D check

using GLMakie, CSV, DataFrames

GLMakie.activate!()
gl = GLMakie

x = collect(0:0.01:100) 

y = x .* cos.(x).^2

z =  y.*cos.(x) .+ 0.3y

gl.scatter(x,y,z, color = :green, markersize = 0.1)

Ph_SS_Ph_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_21_NSTA_FTIR_PhSSPh_TBAHCO2_MeCN/NS_TA/20mM_PhSSPh_MeCN_ArPurge_420nm_TAKinetics_500avg.csv", header = true, DataFrame)

gl.lines(Ph_SS_Ph_Kinetics[:,1],Ph_SS_Ph_Kinetics[:,2])