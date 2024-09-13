using CSV, DataFrames, ColorSchemes, Interpolations, CairoMakie; CairoMakie.activate!()
using LaTeXStrings

size_in_inches = (3.25, 3.25)
dpi = 300
res = size_in_pixels = size_in_inches .* dpi

Data = CSV.read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/Summaries/Interpolated_Data/IrB_C12Bicarb_TRIR.csv",DataFrame)
a = collect(names(Data))
cm = ColorSchemes.winter #color Scheme
n = 1/ncol(Data) #index for colorscheme
size = 20 #fontsize
width = 4 #Linewidth for traces

x = Data[:,1]

#First Figure
function Fig1(f = Figure(fontsize = size,resolution = res))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber",xlabelsize = size, xtickalign = 0, xticksize = 5, xminorticksvisible = true,xgridvisible = false,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false)
    xlims!(ax,1250,1790)
    lines!(f[1,1],x,Data[:,3],linewidth = width,label = "10 ns", color = cm[n])
    lines!(f[1,1],x,Data[:,3],linewidth = width,label = "25 ns", color = cm[2*n])
    lines!(f[1,1],x,Data[:,4],linewidth = width,label = "50 ns", color = cm[3*n])
    lines!(f[1,1],x,Data[:,5],linewidth = width,label = "75 ns", color = cm[4*n])
    lines!(f[1,1],x,Data[:,6],linewidth = width,label = "100 ns", color = cm[5*n])
    lines!(f[1,1],x,Data[:,7],linewidth = width,label = "150 ns", color = cm[6*n])
    lines!(f[1,1],x,Data[:,8],linewidth = width,label = "200 ns", color = cm[7*n])
    lines!(f[1,1],x,Data[:,9],linewidth = width,label = "400 ns", color = cm[8*n])
    lines!(f[1,1],x,Data[:,10],linewidth = width,label = "1000 ns", color = cm[9*n])
    lines!(f[1,1],x,Data[:,11],linewidth = width,label = "2000 ns", color = cm[10*n])
    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :rt, nbanks = 2,framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end
Fig1()

#Second Figure
function Fig2(f = Figure(fontsize = size, resolution = res))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber",xlabelsize = size, xtickalign = 0, xticksize = 5,xminorticksvisible = true, xgridvisible = false,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false)
    xlims!(ax,1250,1790)
    lines!(f[1,1],x,Data[:,3],linewidth = width,label = "10 ns", color = cm[n])
    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :rt, nbanks = 2, framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end
Fig2()

#Third Figure
function Fig3(f = Figure(fontsize = size, resolution = res))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber",xlabelsize = size, xtickalign = 0, xticksize = 5,xminorticksvisible = true, xgridvisible = false,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false)
    xlims!(ax,1250,1790)
    lines!(f[1,1],x,Data[:,3],linewidth = width,label = "10 ns", color = cm[n])
    lines!(f[1,1],x,Data[:,4],linewidth = width,label = "50 ns", color = cm[3*n])
    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :rt, nbanks = 2, framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end
Fig3()

#Fourth Figure
function Fig4(f = Figure( fontsize = size, resolution = res))
    ax = CairoMakie.Axis(f[1,1],xlabel = "Wavenumber",xlabelsize = size, xtickalign = 0, xticksize = 5,xminorticksvisible = true, xgridvisible = false,
    ylabel =  "ΔA", ylabelsize = size, ytickalign = 0, yticksize = 5,yminorticksvisible = true, ygridvisible = false)
    xlims!(ax,1250,1790)
    lines!(f[1,1],x,Data[:,3],linewidth = width,label = "10 ns", color = cm[n])
    lines!(f[1,1],x,Data[:,4],linewidth = width,label = "50 ns", color = cm[3*n])
    lines!(f[1,1],x,Data[:,9],linewidth = width,label = "400 ns", color = cm[8*n])
    lines!(f[1,1],x,Data[:,11],linewidth = width,label = "2000 ns", color = cm[10*n])
    hlines!(ax,0,linewidth = 1,linestyle = :dash,color = :black)
    axislegend("Delay Time", position = :rt, nbanks = 2, framevisible = false)
    colsize!(f.layout,1,Aspect(1, 1.25))
    f
end
Fig4()


# Set up the Maine figure Panel with 4 quadrants
p = Figure(font = :Helvetica, fontsize = size, resolution = res)

# Define Quadrants
p1 = p[1,1]
p2 = p[1,2]
p3 = p[2,1]
p4 = p[2,2]

Fig1(p1)
Fig2(p2)
Fig3(p3)
Fig4(p4)

#Set spacing/size between columns and rows
colsize!(p.layout,1,Aspect(1, 1.25))
colsize!(p.layout,2,Aspect(1,1.25))
rowsize!(p.layout,1,500)
rowsize!(p.layout,2,500)
rowgap!(p.layout,1,10)

## Figure Panel labels
Label(p[1, 1,TopLeft()], "(a)")
Label(p[1, 2,TopLeft()], "(b)")
Label(p[2, 1,TopLeft()], "(c)")
Label(p[2, 2,TopLeft()], "(d)")

#reszie the figure content to the layout
resize_to_layout!(p)
p

save("/Users/kris/Documents/JULIA/Bicarb_Paper_Figures/Bicarb_TRIR_Panel.pdf",p)

