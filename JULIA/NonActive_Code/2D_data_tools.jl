# written for plotting and subtracting 2D intensity data
# Author: Ryan James Evenson, G5 Nocera Research Group, evenson@g.harvard.edu

using CSV, DataFrames, ColorSchemes, Interpolations, CairoMakie; CairoMakie.activate!()

mutable struct Intensity2D
    matrix::Matrix{Float64}
    wavenumbers::Vector{Float64}
    xlabel::String
    ylabel::String
    title::String
end

mutable struct UserDefinedKinetics
    τ::Float64
    A0::Float64
end

currenttype = Intensity2D

function read2D(file::String, x::String, y::String, title::String) #read common formats
    m = Matrix(CSV.read(file,DataFrame))
    wn = m[:,1]
    m = m[:,begin+1:end]
    return Intensity2D(m,wn,x,y,title)
end

function write2D() #write to commmon format

end

function GatherAllignedPoints(coordinate_sets::Vector{Vector{Vector{Float64}}},x_interval::Vector{Float64})
    interps::Vector{Any} = []
    out::Vector{Vector{Float64}} = []

    for v in coordinate_sets
        x,y = v
        print("\nx = $x\ny = $y")
        push!(interps,LinearInterpolation(reverse(x),reverse(y)))
    end

    for i in interps
        push!(out,i(x_interval))
    end

    return out
end


function GenerateKineticSingal(data2D::Intensity2D,signal1D::Vector{Vector{Float64}},signal_kinetics::UserDefinedKinetics) #signal1D as [x,y] for interpolation reasons
    #create a 2D data set of a 1D signal with specific kinetics, such that the created data set matches the dimensions of
    # the provided 2D reference data. This allows you to in-place add or subtract the input and output 2D data sets

    kinetics_matrix = deepcopy(data2D.matrix)
    sample_matrix = deepcopy(data2D.matrix)

    kinetics_matrix .= 1.0 #clean up matrix

    gapinput::Vector{Vector{Vector{Float64}}} = [[data2D.wavenumbers,Vector(data2D.matrix[:,1])],signal1D]

    signal1D_properspacing = GatherAllignedPoints(gapinput,data2D.wavenumbers)[2]

    A0 = signal_kinetics.A0
    τ = signal_kinetics.τ

    time = range(0,stop = 5.5*length(data2D.matrix[1,:]), length = length(data2D.matrix[1,:]))

    kinetics_matrix = [abs*A0*exp(-t/τ) for abs in signal1D_properspacing, t in time]

    sample_matrix = hcat(repeat([signal1D_properspacing],length(time))...)

    return sample_matrix .* kinetics_matrix

end




function nice2Dplot(data::currenttype,xpx::Int64,ypx::Int64)

    #data.matrix = data.matrix ./ maximum(data.matrix)

    figure = ( 
        resolution=(xpx, ypx),
        #font="C:\\Users\\R. Evenson\\Downloads\\nimbus-roman-no9-l.regular.otf",
        fontsize = 20)

    axis = (
        xlabel = data.xlabel,
        ylabel = data.ylabel,
        title = data.title,
        yminorticksvisible = true,
        xminorticksvisible = true,
        xreversed = true,
        xminorticks = IntervalsBetween(4),
        yminorticks = IntervalsBetween(4)
        )

        cm = ColorSchemes.roma

        lim = (minimum(data.matrix),maximum(data.matrix))

        x = data.wavenumbers
        y = range(0,stop = 5.5*length(data.matrix[:,1]), length = length(data.matrix[:,1]))

    fig, ax, lin = CairoMakie.heatmap(x,y,data.matrix; colormap = cm, axis=axis, figure=figure)
    CairoMakie.Colorbar(fig[1,2], limits = lim, label = "Itensity (A.U.)",colormap = cm)
    
    fig

end

csv1 = "C:\\Users\\R. Evenson\\Documents\\GitHub\\Evenson\\screwing around\\bullshitdata.csv"
csv1 = "C:\\Users\\R. Evenson\\Downloads\\IrB-2mM-Bicarb-10mM-dACN-8res-2mm-1VAC-5VDC-A-Flow-Absorbance.0.DPT"

data = read2D(csv1,"Wavenumber (cm⁻¹)","time (ns)","IrB-2mM-Bicarb-10mM-dACN-8res-2mm-1VAC-5VDC-A-Flow-Absorbance")
data.matrix[1,:]
data.matrix[:,1]
data.matrix
nice2Dplot(data,900,900) # the two numbers are the aspects of the resulting plot in pixels, x pixles, y pixels



kin = UserDefinedKinetics(50.0,0.02)
out = GenerateKineticSingal(data,[data.wavenumbers,data.matrix[:,1] ],kin)

subtracted_data = deepcopy(data)
subtracted_data.matrix = subtracted_data.matrix .- 500out
subtracted_data.title = "adjustments made"

nice2Dplot(data, 900, 1000)
nice2Dplot(subtracted_data, 900, 1000)