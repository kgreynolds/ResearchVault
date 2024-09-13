using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes,MAT ,CairoMakie; CairoMakie.activate!()
plotly()

## IMPORT DATA

Pyridostigmine_UvVis_Calc = CSV.read("/Users/kris/Desktop/ResearchVault/COURSES/Chem165/Pyridostigmine/Pyridostigmine_uvvis.csv", header = true, DataFrame)


MehtylPyridinium_uvvis= CSV.read("/Users/kris/Desktop/ResearchVault/COURSES/Chem165/MehtylPyridinium_decomp/MehtylPyridinium_uvvis.csv", header = true, DataFrame)

### Make Figures 

function Pyridostigmine_UvVis(f = Figure())
    width = 3

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], title = "Pyridostigmine and Decomp product Calculated Absorption Spectra",
            palette = (color = palette(ColorSchemes.hsv, 3),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "Wavelength", 
            ylabel = "ϵ"
        )
        # ax.xticks= 0:100:400
        # CairoMakie.xlims!(200,350)
        # CairoMakie.ylims!(nothing,)
        
        CairoMakie.lines!(Pyridostigmine_UvVis_Calc[:,1],Pyridostigmine_UvVis_Calc[:,2],linewidth = width, label = "Pyridostigmine")
        CairoMakie.lines!(MehtylPyridinium_uvvis[:,1],MehtylPyridinium_uvvis[:,2],linewidth = width, label = "Methyl Pyridinium")



        axislegend( position = :rt,nbanks = 1,framevisible = false, labelsize = 15) #:rt (right-top :rb (right bottom) etc)

        # text!(x,y,text = "text", color = :red, textsize = 20)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Pyridostigmine_UvVis()
save("JULIA/Outputs/Pyridostigmine_UvVis.pdf",Pyridostigmine_UvVis())

# SOME FUNCTIONS

    #Normalization to x=z
    function Norm(x,y,z)
        """Normalize Data to a value Z """
        x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
        y_interp = x_interp(x) # compute y values for the x vector
        y_norm = y_interp/x_interp(z) #normalize to desired value
        return y_norm
    end

    # Normalize to Maximum
    function MaxNorm(y)
        """Normalize vector to its absolute largest value"""
        y_norm = y/maximum(abs.(y)) #normalize to desired value
        return y_norm
    end

    function Resize(Data,column, length)
        """Resize a column of a dataframe to desired length of Data points """
        x_out = collect(range(Data[1,1],Data[end,1],length))
        Interp  = LinearInterpolation(Data[:,1],Data[:,column])
        y_out = Interp(x_out)
    
        df = DataFrame(x = x_out,y = y_out)
        return df
    end

    function smoothMat(Mat,SmoothPts)
        """This function takes and input matrix Mar and generates a new matrix sMat with the moving average over
         some number of points defined by the SmoothPts argument"""
        sMat = zeros(size(Mat,1)-SmoothPts,size(Mat,2))
    
        for i in 1:size(Mat,1)-SmoothPts
            for j in 1:size(Mat,2)
                sMat[i,j] = mean(skipmissing(Mat[i:i+SmoothPts,j])) # This works! calculates the mean of the i to i+2th row in the jth columns and replaces the E[i,j]th value with that mean
            end
        end
        return sMat
    
    end
