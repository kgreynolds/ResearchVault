using CSV, DataFrames, ColorSchemes,CairoMakie; CairoMakie.activate!()
using Images


    #Import Data to Plot
    Data1 = CSV.read("RESEARCH/NOCERA GROUP.nosync/PROJECTS/Ag(II)Project/Experiments/2023_07_13_TRIR_Ag(Bpy)(TFA)2_MeCN/Ag(Bpy)(TFA)2_TRIR_Spectra_1100_1300.csv",DataFrame)
    Data2 = CSV.read("RESEARCH/NOCERA GROUP.nosync/PROJECTS/Ag(II)Project/Experiments/2023_07_13_TRIR_Ag(Bpy)(TFA)2_MeCN/Ag(Bpy)(TFA)2_TRIR_spectra_1600_2200Abs.csv",DataFrame)
    Data3 = CSV.read("RESEARCH/NOCERA GROUP.nosync/PROJECTS/Ag(II)Project/Experiments/2023_07_13_TRIR_Ag(Bpy)(TFA)2_MeCN/Ag(Bpy)(TFA)2_TRIR_Spectra_2320_2450.csv",DataFrame)
    a = collect(names(Data3))
    ncol(Data1)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data1)-1)
    #Make Cairo Mackie figure
    font = :Arial
    f = Figure(font = font,figure_padding = 25,fontsize =20)
    ax = CairoMakie.Axis(f[1,1], titlefont = font,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    # CairoMakie.xlims!(ax,350,650)
    # ylims!(nothing,)


    width = 3
    size = 20
    x1 = Data1[:,1]
    # bv1 = 1125 .< x1 .< 1275

    x2 = Data2[:,1]
    x3 = Data3[:,1]



#First Figure
function Fig1(f = Figure(fontsize = size,figure_padding = 25))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber (cm⁻¹)",xlabelsize = size, xtickalign = 0, xticksize = 5, xminorticksvisible = true,xgridvisible = false,xreversed = true,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false, xautolimitmargin = (0.0, 0.0))
    # xlims!(ax,1125,1275)
    
    lines!(f[1,1],x1[bv1],Data1[:,2][bv1],linewidth = width,label = "0 ns", color = cm[n])
    # lines!(f[1,1],x1,Data1[:,3],linewidth = width,label = "10 ns", color = cm[n])
    # lines!(f[1,1],x1,Data1[:,4],linewidth = width,label = "20 ns", color = cm[3*n])
    # lines!(f[1,1],x1,Data1[:,5],linewidth = width,label = "30 ns", color = cm[4*n])
    # lines!(f[1,1],x1,Data1[:,6],linewidth = width,label = "40 ns", color = cm[5*n])
    # lines!(f[1,1],x1,Data1[:,7],linewidth = width,label = "50 ns", color = cm[6*n])
    # lines!(f[1,1],x1,Data1[:,8],linewidth = width,label = "60 ns", color = cm[7*n])
    # lines!(f[1,1],x1,Data1[:,9],linewidth = width,label = "80 ns", color = cm[8*n])
    # lines!(f[1,1],x1,Data1[:,10],linewidth = width,label = "100 ns", color = cm[9*n])
    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :rb, nbanks = 2,framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end

Fig1()

save("JULIA/Outputs/Fig1.eps",Fig1())

#Second Figure
function Fig2(f = Figure(fontsize = size,figure_padding = 25))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber (cm⁻¹)",xlabelsize = size, xtickalign = 0, xticksize = 5,xminorticksvisible = true, xgridvisible = false,xreversed = true,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false)
    # xlims!(ax,1250,1790)
    lines!(f[1,1],x2,Data2[:,2],linewidth = width,label = "0 ns", color = cm[n])
    lines!(f[1,1],x2,Data2[:,3],linewidth = width,label = "10 ns", color = cm[2*n])
    lines!(f[1,1],x2,Data2[:,4],linewidth = width,label = "20 ns", color = cm[3*n])
    lines!(f[1,1],x2,Data2[:,5],linewidth = width,label = "30 ns", color = cm[4*n])
    lines!(f[1,1],x2,Data2[:,6],linewidth = width,label = "40 ns", color = cm[5*n])
    lines!(f[1,1],x2,Data2[:,7],linewidth = width,label = "50 ns", color = cm[6*n])
    lines!(f[1,1],x2,Data2[:,8],linewidth = width,label = "60 ns", color = cm[7*n])
    lines!(f[1,1],x2,Data2[:,9],linewidth = width,label = "80 ns", color = cm[8*n])
    lines!(f[1,1],x2,Data2[:,10],linewidth = width,label = "100 ns", color = cm[9*n])
    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :lb, nbanks = 2,framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end
Fig2()

save("JULIA/Outputs/Fig2.eps",Fig2())

#Third Figure
function Fig3(f = Figure(fontsize = size,figure_padding = 25))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber (cm⁻¹)",xlabelsize = size, xtickalign = 0, xticksize = 5,xminorticksvisible = true, xgridvisible = false,xreversed = true,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false)
    # CairoMakie.xlims!(2325,2375)
    lines!(f[1,1],x3,Data3[:,2],linewidth = width,label = "0 ns", color = cm[n])
    lines!(f[1,1],x3,Data3[:,3],linewidth = width,label = "10 ns", color = cm[2*n])
    lines!(f[1,1],x3,Data3[:,4],linewidth = width,label = "20 ns", color = cm[3*n])
    lines!(f[1,1],x3,Data3[:,5],linewidth = width,label = "30 ns", color = cm[4*n])
    lines!(f[1,1],x3,Data3[:,6],linewidth = width,label = "40 ns", color = cm[5*n])
    lines!(f[1,1],x3,Data3[:,7],linewidth = width,label = "50 ns", color = cm[6*n])
    lines!(f[1,1],x3,Data3[:,8],linewidth = width,label = "200 ns", color = cm[7*n])

    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :lt, nbanks = 2,framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end
Fig3()

save("JULIA/Outputs/Fig3.eps",Fig3())


# Set up the Maine figure Panel with 4 quadrants
p = Figure(font = :Helvetica, fontsize = size, resolution = res,figure_padding = 25)

# Define Quadrants
p1 = p[1,1]
p2 = p[2,1]
p3 = p[3,1]


Fig1(p1)
Fig2(p2)
Fig3(p3)


#Set spacing/size between columns and rows
colsize!(p.layout,1,Aspect(1, 1.25))

rowsize!(p.layout,1,500)
rowsize!(p.layout,2,500)
rowsize!(p.layout,3,500)



## Figure Panel labels
Label(p[1, 1,TopLeft()], "(a)")
Label(p[2, 1,TopLeft()], "(b)")
Label(p[3, 1,TopLeft()], "(c)")


#reszie the figure content to the layout
resize_to_layout!(p)
p

save("JULIA/Outputs/TRIR_Figure_Ag(bpy)TFA2.png",p)
