using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,MAT ,CairoMakie; CairoMakie.activate!()
plotly()

n1 = 1.52 # refractive index of Glass
θ1 = 0:0.001:(pi/2)
collect(θ1)
n2 = 1.00029 #refractive index of Air

    ref = []
for i in 1:length(θ1)   
    r = (n1 *cos(θ1[i]) - √(Complex((n2)^2 - (n1 * sin(θ1[i]))^ 2)))/(n1 * cos(θ1[i]) + √(Complex((n2)^2 -(n1 * sin(θ1[i]))^2)))
    push!(ref, r)
end
return R = abs.(ref)


function TE_Ref(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Glass-to-air transition",
            palette = (color = palette(ColorSchemes.hsv, 4),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Incident angle (θᵢ)", 
        )
        ax.xticks= 0:10:90
        CairoMakie.xlims!(nothing,60)
        # CairoMakie.ylims!(nothing,)
        
        CairoMakie.lines!(θ1* (180/pi),1 .- (R .^ 2),linewidth = width, label = "T")

        CairoMakie.lines!(θ1 * (180/pi),R .^ 2,linewidth = width, label = "R")



        axislegend( position = :rt,nbanks = 1,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
TE_Ref()
save("JULIA/Outputs/TE_Ref.pdf",TE_Ref())