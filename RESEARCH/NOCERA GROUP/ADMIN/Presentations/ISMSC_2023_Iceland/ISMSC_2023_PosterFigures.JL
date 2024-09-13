using CSV, DataFrames, ColorSchemes,CairoMakie; CairoMakie.activate!()
using Images


## Pd(II) TFP TA Spectra in THF
    #Import Data to Plot
    Data = CSV.read("JULIA/Inputs/Pd-THF 400nm Pump TA spectra.csv",DataFrame)
    a = collect(names(Data))
    ncol(Data)

    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    # cm = ColorSchemes.linear_kbc_5_95_c73_n256
    cm = ColorSchemes.linear_kry_0_97_c73_n256
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    font = :HelveticaNue
    f = Figure(font = font,figure_padding = 25,fontsize =25)
    ax = CairoMakie.Axis(f[1,1],titlefont = font,
    xlabel = "Wavelength (nm)",xlabelsize = 25, xtickalign = 1, xticksize = 10, xgridvisible = false,xminorticksvisible = true, xminortickalign = 1, xminorticksize = 10, xminorticks = IntervalsBetween(5),
    ylabel =  "ΔA", ylabelsize = 25, ytickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, yminortickalign = 1, yminorticksize = 10 , yminorticks = IntervalsBetween(5) )
    xlims!(ax,420,760)
    # ylims!(nothing,)
    width = 3
    x = Data[:,1]

    # lines!(x,Data[:,2],linewidth = width,label = "-500 fs",color = cm[n])
    lines!(x,Data[:,3],linewidth = width,label = "0.1 ps",color = cm[3*n])
    lines!(x,Data[:,4],linewidth = width,label = "0.2 ps",color = cm[4*n])
    lines!(x,Data[:,5],linewidth = width,label = "0.3 ps",color = cm[5*n])
    lines!(x,Data[:,6],linewidth = width,label = "0.4 ps",color = cm[6*n])
    lines!(x,Data[:,7],linewidth = width,label = "0.5 ps",color = cm[7*n])
  
    lines!(x,Data[:,9],linewidth = width,label = "1 ps",color= cm[8*n])
    lines!(x,Data[:,10],linewidth = width,label = "10 ps",color = cm[9*n])
    lines!(x,Data[:,11],linewidth = width,label = "50 ps",color = cm[10*n])
    lines!(x,Data[:,12],linewidth = width,label = "100 ps",color = cm[11*n])
    lines!(x,Data[:,13],linewidth = width,label = "200 ps",color = cm[12*n])
    lines!(x,Data[:,14],linewidth = width,label = "1000 ps",color = cm[13*n])
    lines!(x,Data[:,15],linewidth = width,label = "5000 ps",color = cm[14*n])

    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    f
    save("JULIA/Outputs/Figure.png",f, pix_per_unit = 8) #save file as eps, png, pdf etc.


## Mn(III) TFP TA Spectra in THF
    #Import Data to Plot
    Data = CSV.read("JULIA/Inputs/Mn-P in THF TA Spectra.csv",DataFrame)
    a = collect(names(Data))
    ncol(Data)
    #Set color Scheme (n increments over the colorscheme depends on dataframe size)
    # cm = ColorSchemes.linear_kbc_5_95_c73_n256
    cm = ColorSchemes.winter
    n = 1/(ncol(Data)-1)
    #Make Cairo Mackie figure
    font = :HelveticaNue
    f = Figure(font = font,figure_padding = 25, fontsize = 25)
    ax = CairoMakie.Axis(f[1,1], titlefont = font,
    xlabel = "Wavelength (nm)",xlabelsize = 25, xtickalign = 1, xticksize = 10, xgridvisible = false,xminorticksvisible = true, xminortickalign = 1, xminorticksize = 10, xminorticks = IntervalsBetween(5),
    ylabel =  "ΔA", ylabelsize = 25, ytickalign = 1, yticksize = 10, ygridvisible = false,yminorticksvisible = true, yminortickalign = 1, yminorticksize = 10 , yminorticks = IntervalsBetween(5) )
    xlims!(ax,420,760)
    # ylims!(nothing,)
    width = 3
    x = Data[:,1]

    # lines!(x,Data[:,2],linewidth = width,label = "-500 fs",color = cm[n])
    lines!(x,Data[:,3],linewidth = width,label = "0.1 ps",color = cm[3*n])
    lines!(x,Data[:,4],linewidth = width,label = "0.2 ps",color = cm[4*n])
    lines!(x,Data[:,5],linewidth = width,label = "0.3 ps",color = cm[5*n])
    lines!(x,Data[:,6],linewidth = width,label = "0.4 ps",color = cm[6*n])
    lines!(x,Data[:,7],linewidth = width,label = "0.5 ps",color = cm[7*n])
  
    lines!(x,Data[:,9],linewidth = width,label = "1 ps",color= cm[8*n])
    lines!(x,Data[:,10],linewidth = width,label = "10 ps",color = cm[9*n])
    lines!(x,Data[:,11],linewidth = width,label = "50 ps",color = cm[10*n])
    lines!(x,Data[:,12],linewidth = width,label = "100 ps",color = cm[11*n])
    lines!(x,Data[:,13],linewidth = width,label = "200 ps",color = cm[12*n])


    hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)

    axislegend("Delay Time", position = :rb,nbanks = 3,framevisible = false) #:rt (right-top :rb (right bottom) etc)
    # Label(f[1,1,TopLeft()],"(a)")
    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis)
    resize_to_layout!(f)
    
    f
    save("JULIA/Outputs/Figure.png",f, pix_per_unit = 8) 

