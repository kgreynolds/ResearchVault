## IMPORT ANALYSIS FUNCTIONS
include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

Kris_Figure_Theme = merge(Theme(
    Axis = (
        xautolimitmargin = (0.01, 0.01),  yautolimitmargin = (0.05, 0.05), alignmode = Outside(),
        xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20, xlabelpadding = 10,
        ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20, ylabelpadding = 10,
    )
    ), theme_latexfonts())
set_theme!(Kris_Figure_Theme)


## IMPORT DATA
file1 = CSV.read("path", header = true, DataFrame)

# Data Workup
####################################################################################################
Data = MonoFitIRF_New(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017, 0.0], 10,Inf,false,10000)
Data = BiFitIRF_new(Data,2,0.0, 800, [0.008, 0.001, 50.592, 307.017,10.0, 100, 0.0], 10,Inf,false,10000)


### Make Figures 
####################################################################################################
Anthracene_TA_Spectra = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_18/Anthracene_MeCN_TA_Spectra_Surface_processed.csv", header =true, DataFrame)
Anthracene_TA_Kinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_18/Anthracene_MeCN_TA_Kinetics_Surface_processed.csv", header =true, DataFrame)

Anthracene_TA_Kinetics_367nm = BiFitIRF_new(Anthracene_TA_Kinetics,4,-10.0, 1000, [4.0808, 0.9464, 0.5103, 3.4901, 0.9754, 161.3797, -0.0001], 1,100,false,200)
Anthracene_TA_Kinetics_399nm = BiFitIRF_new(Anthracene_TA_Kinetics,5,-10.0, 1000, [4.122, 0.956, -0.1378, 3.7869, 7.2715, 146.7088, 0.0001], 1,1000,false,200)
Anthracene_TA_Kinetics_420nm = BiFitIRF_new(Anthracene_TA_Kinetics,2,-10.0, 1000, [4.122, 0.956, -0.6092, 3.947, 22.6853, 147.9505, 0.0001], 1,1000,false,200)
Anthracene_TA_Kinetics_601nm = MonoFitIRF_New(Anthracene_TA_Kinetics,3,0.0, 1000, [4.122, 0.956, 0.119, 3.692, -0.0], 10,10,false,200)

Singlet  = mean([3.4901,3.7869,3.947])
Singlet_stdev  = std([3.4901,3.7869,3.947,3.692])

Triplet  = mean([161.3797,146.7088,147.9505])
Triplet_stdev  = std([161.3797,146.7088,147.9505])


function Anthracene_TA_Fig(f = Figure())

    TA_Spectra = Anthracene_TA_Spectra
    TA_Kinetics = [
        (Anthracene_TA_Kinetics_399nm[1], " 399 nm"),
        (Anthracene_TA_Kinetics_420nm[1], "420 nm"),
        (Anthracene_TA_Kinetics_601nm[1], "601 nm")
    ]
    
    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, 20+1),), ## or replace ncol with Integer 
                title = "Anthracene Ns TA Spectra",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(360,700)
            # CairoMakie.ylims!(nothing,0.04)

            for i in [2,6,7,8,9,10,11,12,13,15,17,19,20,21,22,23,24,25,26,27] ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
            
            end

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false,labelsize = 15)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    #### Second Panel TA KINETICS ##############################################################################

        ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "Anthracene Ns TA Kinetics",
                xlabel = "Time (ns)", 
                ylabel = "ΔA (OD)", 
                # xscale = Makie.Symlog10(1.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                xminorticks = IntervalsBetween(2)
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(-10, 2000)
            # CairoMakie.ylims!(-0.001,0.02)

        for (data, label) in TA_Kinetics
            lines!(ax2, data[:, 1], data[:, 2],linewidth = 3,label = label) # Kinetic Trace
            lines!(ax2, data[:, 1], data[:, 3], linewidth = 2, linestyle = :dash, color = :black) # Fit 
        end

            axislegend("Wavelength", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                ### Add Inset axis into plot for kinetics
                # ax2_1 = Axis(f[1,2], width=Relative(0.5), height=Relative(0.5),halign=0.8, valign=0.5, backgroundcolor=:white, 
                #xlabel = "Time (ns)", 
                #ylabel = "ΔA (OD)", 
                # )

                #CairoMakie.xlims!(-10, 20)
                # CairoMakie.ylims!(-0.003,0.02)

                # for (data, label) in TA_Kinetics
                #     Makie.lines!(ax2_1, data[:, 1], data[:, 2],linewidth = 3,label = label)
                #     lines!(ax2_1, data[:, 1], data[:, 3],linewidth = 2, linestyle = :dash, color = :black)
                # end

            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
Anthracene_TA_Fig()
save("JULIA/Outputs/Anthracene_TA_Fig.png",Anthracene_TA_Fig())

[1,2,4]

function Anthracene_TA_Fig_2(f = Figure())

    TA_Spectra = Anthracene_TA_Spectra

    
    ### First Panel TA Spectra  ##############################################################################

        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.inferno, 9+1),), ## or replace ncol with Integer 
                title = "Anthracene early time TA Spectra",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                xminorticks = IntervalsBetween(2)
            )
        
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(360,700)
            # CairoMakie.ylims!(nothing,0.04)

            for i in [2,6,7,8,9,10,11,12,13]## or replace ncol(Data) with a vector of custom indecies
                lines!(ax, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
            
            end

            axislegend("Delay Time", position = :rt,nbanks = 3,framevisible = false,labelsize = 15)
            # Legend(f[1,2], ax, "Delay Time (ns)", framevisible = false, nbanks = 1 ) #:rt (right-top :rb (right bottom) etc)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            # text!(x,y,text = "text", color = :red, textsize = 20)
            # CairoMakie.vspan!(433, 439; ymin = 0.0, ymax = 1.0, color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    #### Second Panel TA KINETICS ##############################################################################

        ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.inferno, 12+1),),
                title = "Anthracene late time TA Spectra",
                xlabel = "Wavelength (nm)", 
                ylabel = "ΔA (OD)", 
                # xscale = Makie.Symlog10(1.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                xminorticks = IntervalsBetween(2)
            )

        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            CairoMakie.xlims!(360, 700)
            # CairoMakie.ylims!(-0.001,0.02)

            for i in [13,15,17,19,20,21,22,23,24,25,26,27] ## or replace ncol(Data) with a vector of custom indecies
                lines!(ax2, TA_Spectra[:,1],TA_Spectra[:,i],linewidth = width,label = " $(names(TA_Spectra)[i])")
            
            end

            axislegend("Delay Time", position = :rt, nbanks = 3, framevisible = false, labelsize = 15)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

                ### Add Inset axis into plot for kinetics
                # ax2_1 = Axis(f[1,2], width=Relative(0.5), height=Relative(0.5),halign=0.8, valign=0.5, backgroundcolor=:white, 
                #xlabel = "Time (ns)", 
                #ylabel = "ΔA (OD)", 
                # )

                #CairoMakie.xlims!(-10, 20)
                # CairoMakie.ylims!(-0.003,0.02)

                # for (data, label) in TA_Kinetics
                #     Makie.lines!(ax2_1, data[:, 1], data[:, 2],linewidth = 3,label = label)
                #     lines!(ax2_1, data[:, 1], data[:, 3],linewidth = 2, linestyle = :dash, color = :black)
                # end

            Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
            Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
            colgap!(f.layout,1,20)
            colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size

    resize_to_layout!(f)
    f
end
Anthracene_TA_Fig_2()
save("JULIA/Outputs/Anthracene_TA_Fig_2.png",Anthracene_TA_Fig_2())


####### PCA AcrPO ##########################################################################################
TA_data = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP/PROJECTS/Gabbai_Collaboration/Experiemnts/Nanosecond TA/2024_09_18/Anthracene_TA_carpet.csv", header = false, DataFrame) # TA carpet data with wavelength in rows and time points as columns

# Extract the absorbance data (excluding the first row and column)
absorbance_data = Matrix(TA_data[2:end, 2:end])

Wavelengths = TA_data[2:end,1]
Delay_times = transpose(Matrix(TA_data))[2:end,1]

function pca_ta(data::Matrix{Float64})
    # Assume data is a matrix where each column is a spectrum at a different time point
    # and each row corresponds to a wavelength

    # Center the data
    data_centered = data .- mean(data, dims=2)

    # Perform SVD
    U, S, V = svd(data_centered)
    ## U = principle components
    ## S = Singular Values
    ## V = time dependnet scores

    # Calculate explained variance ratio
    explained_variance_ratio = S.^2 / sum(S.^2)

    # Return principal components, singular values, and explained variance ratio
    return U, S, V, explained_variance_ratio
end

U, S, V, explained_variance_ratio = pca_ta(absorbance_data) .* (-1)
Anthracene_PCA_Spectra = DataFrame(Wavelength_nm = Wavelengths, PC_1 = U[:,1], PC_2 = U[:,2], PC_3 = U[:,3])
Anthracene_PCA_Kinetics = DataFrame(time_µs = Delay_times,Kinetic_PC1 = V[:,1],Kinetic_PC2 = V[:,2],Kinetic_PC3 = V[:,2])

PCA_1_Kinetics = BiFitIRF_new(PCA_Kinetics,2,-10.0, 200,[0.0018, 0.0009, -0.0209, 0.3853, -0.0012, 0.0016, 0.0274], 10,10,false,400)
PCA_2_Kinetics = BiFitIRF_new(PCA_Kinetics,3,-10.0, 200,[0.0086, 0.001, -0.0666, 0.3682, 0.0085, 0.1588, 0.0536], 10,10,false,400)


function PCA_Fig(f = Figure())

    U, S, V, explained_variance_ratio = pca_ta(absorbance_data) .* (-1)
    Wavelengths = TA_data[2:end,1]
    Delay_times = transpose(Matrix(TA_data))[2:end,1]
    components = collect(1:2)
    ### 1st Panel Principal Components  ##############################################################################
        ax =    CairoMakie.Axis(f[1,1], palette = (color = palette(ColorSchemes.hsv, 4),), ## or replace ncol with Integer
                title = "Principal Components",
                xlabel = "Wavelength (nm)", 
                ylabel = "Amplitude", 
                xminorticks = IntervalsBetween(2)
            )
        width = 3 

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(375,910)
            # CairoMakie.ylims!(-0.02,nothing)

        ### Plot Data ##############################################################################

        for i in components #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
            lines!(ax, Wavelengths,U[:,i],linewidth = width,label = "PCA $i")
        end

            axislegend(position = :rt,nbanks = 3, framevisible = false, labelsize = 15)
            hlines!(ax,0, linewidth = 1,color = :black,linestyle = :dash)
            vlines!(ax,436, linewidth = 1,color = :lightblue)

            colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
            rowsize!(f.layout,1, 400)
        
    # ### 2nd Panel SVD Kinetics ##############################################################################
        ax2 =    CairoMakie.Axis(f[1,2], palette = (color = palette(ColorSchemes.hsv, 4),),
                title = "Time-dependent Scores",
                xlabel = "Time (ns)", 
                ylabel = "Score", 
                # xscale = Makie.Symlog10(10.0), 
                # xticks = [-5, 0, 5, 10,50, 100, 1000], 
                # xminorticks = IntervalsBetween(2)
        )
        width = 3

        ### Custom Axis limits ##############################################################################
            # ax.xticks= 1:2:17
            # CairoMakie.xlims!(-10, 1500)
            # CairoMakie.ylims!(-0.01,0.03)

        ### Plot Kinetic Data ##############################################################################

        for i in components #2:ncol(TA_Spectra) ## or replace ncol(Data) with a vector of custom indecies
            lines!(ax2, Delay_times,V[:,i],linewidth = width,label = "PCA $i")
        end
            axislegend( position = :rt, nbanks = 3, framevisible = false,labelsize = 15)
            hlines!(ax2,0, linewidth = 1,color = :black,linestyle = :dash)

    # ### 3rd Panel: Explained Varience ##############################################################################

        ax3 =    CairoMakie.Axis(f[2,1], palette = (color = palette(ColorSchemes.hsv, 4),),
            title = "Explained varience",
            xlabel = "Principal Component", 
            ylabel = "Cumulative Explained Varience Ratio", 
            # xscale = Makie.Symlog10(10.0), 
            # xticks = [-5, 0, 5, 10,50, 100, 1000], 
            # xminorticks = IntervalsBetween(2)
        )
        width = 3

        ### Custom Axis limits ##############################################################################
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(0, 10)
        # CairoMakie.ylims!(-0.01,0.03)

        ### Plot Kinetic Data ##############################################################################

        CairoMakie.scatter!(ax3,1:length(explained_variance_ratio), cumsum(explained_variance_ratio),markersize = 10, color = :black)
    

    # ### 4th Panel: Singular Values ##############################################################################

        ax3 =    CairoMakie.Axis(f[2,2], palette = (color = palette(ColorSchemes.hsv, 4),),
        title = "Singular Values",
        xlabel = "Principal Component", 
        ylabel = "Singular Value", 
        # xscale = Makie.Symlog10(10.0), 
        # xticks = [-5, 0, 5, 10,50, 100, 1000], 
        # xminorticks = IntervalsBetween(2)
            )
        width = 3

        ### Custom Axis limits ##############################################################################
        # ax.xticks= 1:2:17
        CairoMakie.xlims!(0, 10)
        # CairoMakie.ylims!(-0.01,0.03)

        ### Plot Kinetic Data ##############################################################################

        CairoMakie.scatter!(ax3,1:length(S), S,markersize = 10, color = :black)




    Label(f[1,2,TopLeft()],"(b)", fontsize = 20)
    Label(f[1,1,TopLeft()],"(a)", fontsize = 20)
    Label(f[2,1,TopLeft()],"(c)", fontsize = 20)
    Label(f[2,2,TopLeft()],"(d)", fontsize = 20)
    colgap!(f.layout,1,20)
    colsize!(f.layout,2,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    rowsize!(f.layout,2, 400)

    resize_to_layout!(f)
    f
end


PCA_Fig()
save("JULIA/Outputs/ PCA_AcrPO_TA_Fig.png",PCA_Fig())





