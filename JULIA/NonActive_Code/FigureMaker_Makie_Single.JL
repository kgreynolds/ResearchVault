using CSV, DataFrames, ColorSchemes,CairoMakie; CairoMakie.activate!()
using Images


    #Import Data to Plot
    Data = CSV.read("RESEARCH/NOCERA GROUP.nosync/LASERS/Femtosecond/DATA/2023_05_31/4/Zn-Tpp_ACN_400pump-4ns-test-averaged.csv",DataFrame)
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    font = :Arial
    f = Figure(font = font,figure_padding = 25,fontsize =20)
    ax = CairoMakie.Axis(f[1,1], titlefont = font,
    xlabel = "Wavelength (nm)",xlabelsize = 20, xtickalign = 0, xticksize = 10, xgridvisible = false,xminorticksvisible = true,
    ylabel =  "ΔA", ylabelsize = 20, ytickalign = 0, yticksize = 10, ygridvisible = false,yminorticksvisible = true)
    xlims!(ax,350,650)
    # ylims!(nothing,)
    width = 3
    x = Data[:,1]

    lines!(x,Data[:,2],linewidth = width,label = "-500 fs",color = cm[n])
    lines!(x,Data[:,3],linewidth = width,label = "500 fs",color = cm[2*n])
    lines!(x,Data[:,4],linewidth = width,label = "1 ps",color = cm[3*n])
    lines!(x,Data[:,5],linewidth = width,label = "2 ps",color = cm[4*n])
    lines!(x,Data[:,6],linewidth = width,label = "4 ps",color = cm[5*n])
    lines!(x,Data[:,7],linewidth = width,label = "6 ps",color = cm[6*n])
    lines!(x,Data[:,8],linewidth = width,label = "10 ps",color = cm[7*n])
    lines!(x,Data[:,9],linewidth = width,label = "100 ps",color= cm[8*n])
    lines!(x,Data[:,10],linewidth = width,label = "500 ps",color = cm[9*n])
    lines!(x,Data[:,11],linewidth = width,label = "1000 ps",color = cm[10*n])
    lines!(x,Data[:,12],linewidth = width,label = "2000 ps",color = cm[11*n])
    lines!(x,Data[:,13],linewidth = width,label = "3000 ps",color = cm[12*n])
    lines!(x,Data[:,14],linewidth = width,label = "4000 ps",color = cm[13*n])

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 2,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f

    save("JULIA/Outputs/Figure.png",f, pix_per_unit = 8) #save file as eps, png, pdf etc.
