## Data Files
using CSV, Plots, DataFrames,Interpolations, StatsBase


#interpolate overlap function 
#x interval example input : 1.0:0.01:10.0
#coordinate_sets example input : [[x1,y1],[x2,y2],[x3,y3]...]
function GatherAllignedPoints(coordinate_sets::Vector{Vector{Vector{Float64}}},x_interval)
    interps::Vector{Any} = []
    out::Vector{Vector{Float64}} = []

    for v in coordinate_sets
        x,y = v
        push!(interps,LinearInterpolation(x,y))
    end

    for i in interps
        push!(out,i(collect(x_interval)))
    end

    return out
end


#brute force coefficient search function
# OD is the vector being fitted, and components is a vector of component vectors
# the function returns a vector of weights for each component in the order they were input

function bfcs(OD::Vector{Float64},components::Vector{Vector{Float64}})

    N = length(components) #setting up for N D combination space searching
    grid::Int64 = 1 # how fine the search mesh is, 2grid+1 is the side length of each dimension
    V = repeat([2grid+1],N)
    dims = Tuple(x for x in V) #just some hacky shit that works for generalizing this code to N components being fit
    # for some reason I had trouble writing it in a saner looking way. 

    coeffs::Array{Float64,N} = Array{Float64,N}(undef,dims) # this is the combination space area/volume/hypervolume/whatever we are searching over
    #all of these undefined values are rapidly overwritten. No undefined behavior should result. 

    cycle::Int64 = 1 # a starting point, don't change this. This increments as cycles complete. 
    converged::Bool = false # keeps track of whether or not convergance has been reached. Checked at the end of every cycle.


    #*******************
    # change me if this isn't working
    weight_range::Float64 = 2 #⨦ This SHOULD be changed. This is the upper bound on the magnetiude of the coefficients you are searching for. Make it larger than you expect by a bit. The loops will narrow it down. 
    thresh::Float64 = -1e-10 # convergance is determined by the slope of the residual. When that negative slope flattens out, convergance is reached. 
    resthresh::Float64 = 0.00001 # a goal residual once you know your data. If you make thresh too small to hit, then this can be your exit condition. 

    #less recommended for changing, but may help 
    warm_up_cycles::Int64 = 50 # number of cycles that must elapse before convergance is considered. 
    limit::Int64 = 8000 # just a high number, change it how ever long you want this loop to run as an upper bound. 

    #*******************

    
    reses::Vector{Float64} = Vector{Float64}(undef,15) # this vector containts the numeric residuals. A derivative is taken, and averaged. This average value is used to judge convergance.
    reses .= -500.0 # just arbitrary and large in magnitude.
    
    weights::Vector{Float64} = Vector{Float64}(undef,N) # for use in weighting the components at a particular point in the combination space. 
    minima_weights::Vector{Float64} = Vector{Float64}(undef,N) # for exporting the best fit weights. 
    minima_res::Float64 = 1e30 #arbitrary and large, so that it is replaced by smaller values. If your system every has a residual this big, I would worry. If your scale is natually larger, perhaps you should incerase this number is there are issues. 
    
    w::Vector{Float64} = range(-weight_range, stop = weight_range, length = 2grid+1) # this is the initial weight range before we enter the main loop. 
    ws::Vector{Vector{Float64}} = repeat([w],N) #each component's coefficient dimension is independantly searched.
    # if you didn't independantly search them, you would not be able to simultaneously refine coefficients. 

    #in case the components are different lengths, Interpolations.jl will still let us work with the subset they overlap with.....


    while !converged && cycle <= limit # main loop, could easily be a big for-loop too.

        for ci in CartesianIndices(coeffs) # this loop calculates the residual at each possible coefficient combination within the search space
            # and notes the lowest residuals. 
    
            i = Tuple(ci) #just parses the coordinates into something iterable 

            for (q,d) in enumerate(i) # for each corodinate, assemble the weights their location in space represents.
                weights[q] = ws[q][d] # with those in a vector, they can easily be multipled element-wise
            end

            simulation = sum(components .* weights)  # the math is prettier now. This weights the components and assembles them

            res = sum((OD .- simulation).^2) # a residual is calculated for the simulated data a this point in coefficient space

            coeffs[ci] = copy(res) # in retrospect, I don't know why these were being saved. I had some idea about plotting it or taking fancy gradients or something, but you can literally delete this line and it changes nothing. 

            if res < minima_res # really all that matters is  the displacement of the old best fit by the new best fit. 

                #print("\n$cycle  $(mean(diff(reses)))")
                reses[begin:end-1] = reses[begin+1:end]
                reses[end] = deepcopy(res)
                minima_res = deepcopy(res)
                minima_weights = deepcopy(weights)

            end 
            

        end # if the residual improvement has stagnated, then convergance has been achieved. 

        
        if  minima_res <= resthresh # this is where we check for convergance. 
            converged = true
            print("\nGood enough.\n")
            printstyled("\n\n*************\n  CONVERGED on specific residual   \n\nalternative condition\n~$(((2grid+1)^N)*(cycle+1)) linear combinations tried\nresidual = $minima_res \ncomponent weights = $minima_weights\n*************\n\n", color = :green)
            break
        end
        
        
        if cycle > warm_up_cycles && mean(diff(reses)) >= thresh # this is where we check for convergance. 
            converged = true
            print("\nGood enough.\n")
            printstyled("\n\n*************\n  CONVERGED   \n\n~$(((2grid+1)^N)*(cycle+1)) linear combinations tried\nresidual = $minima_res \ncomponent weights = $minima_weights\n*************\n\n", color = :green)
            break
        end
        

        cycle += 1 #increment the cycle limiter. 

        #TODO: need a new way of doing weights for each dimension seperately, so they can zoom independantly
        #done, this was a development note. 
        
        for (i,d_weight) in enumerate(ws) #this refines the search range to be finer and finer as more cycles accumulate. 

            if minima_weights[i] == 0.0
                minima_weights[i] += rand()
            end

            ws[i] = range(minima_weights[i] - 1.5minima_weights[i]/0.1cycle,
             stop = minima_weights[i] + 1.5minima_weights[i]/0.1cycle, length = 2grid+1) 
            
        end

        if cycle % 100 == 0 
            print("\nworking on it ... residual = $minima_res\n")
        end

        if cycle == limit
            printstyled("\n\n*************\nNOT CONVERGED\n\nresidual = $minima_res \ncomponent weights = $minima_weights\n*************\n\n", color = :yellow)
        end

    end

    return minima_weights # really the only thing we are interested in finding. Other things could exported, feel free to insert them. 
    
end

function write_csv(path::String,headers::Vector{String},data::Vector{Vector{Float64}})
    input::String = ""
    for h in headers
        input *= "$h, "
    end
    input *= "\n"
    block = hcat(data...)
    for i in 1:size(block,1)
        line::String = ""
        for j in 1:size(block,2)
            line *= "$(block[i,j]),"
        end
        input *= "$line \n"
    end
    open(path,"w") do file #actually write the file
        write(file,input)
    end
end

# demo code
# if you haven't used @time before, it's a macro that times how long the code after it takes. 
# I used it to optimize some of the settings within bfcs(), like the grid size.

# Adds y to x 

r = 1250:1:2500
Wvn = collect(r)

##FILES
IrB_TRIR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_TRIR.csv"))
IrB_TRIR_times = GatherAllignedPoints([[IrB_TRIR["Wavenumber"],IrB_TRIR["100ns"]],[IrB_TRIR["Wavenumber"],IrB_TRIR["200ns"]],[IrB_TRIR["Wavenumber"],IrB_TRIR["300ns"]],[IrB_TRIR["Wavenumber"],IrB_TRIR["400ns"]],[IrB_TRIR["Wavenumber"],IrB_TRIR["500ns"]],[IrB_TRIR["Wavenumber"],IrB_TRIR["750ns"]],[IrB_TRIR["Wavenumber"],IrB_TRIR["950ns"]]],r)
IrB_TRIR_100ns = IrB_TRIR_times[1]
IrB_TRIR_200ns = IrB_TRIR_times[2]
IrB_TRIR_300ns = IrB_TRIR_times[3]
IrB_TRIR_400ns = IrB_TRIR_times[4]
IrB_TRIR_500ns = IrB_TRIR_times[5]
IrB_TRIR_750ns = IrB_TRIR_times[6]
IrB_TRIR_950ns = IrB_TRIR_times[7]
plot(r, IrB_TRIR_times[1:end])

IrB_C12Bicarb_TRIR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/C12_bicarb_TRIR Spectra.csv"))
IrB_C12Bicarb_TRIR_times = GatherAllignedPoints([[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t10ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t25ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t50ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t75ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t100ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t150ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t200ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t400ns"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t1us"]],[IrB_C12Bicarb_TRIR["Wavenumber"],IrB_C12Bicarb_TRIR["t2us"]]],r)
C12_10ns = IrB_C12Bicarb_TRIR_times[1]
C12_25ns = IrB_C12Bicarb_TRIR_times[2]
C12_50ns = IrB_C12Bicarb_TRIR_times[3]
C12_75ns = IrB_C12Bicarb_TRIR_times[4]
C12_100ns = IrB_C12Bicarb_TRIR_times[5]
C12_150ns = IrB_C12Bicarb_TRIR_times[6]
C12_200ns = IrB_C12Bicarb_TRIR_times[7]
C12_400ns = IrB_C12Bicarb_TRIR_times[8]
C12_1us = IrB_C12Bicarb_TRIR_times[9]
C12_2us = IrB_C12Bicarb_TRIR_times[10]
plot(r, IrB_C12Bicarb_TRIR_times[1:end])

IrB_C13Bicarb_TRIR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/C13_Bicarb_TRIR_Spectra.csv"))
IrB_C13Bicarb_TRIR_times = GatherAllignedPoints([[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t10ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t25ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t50ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t75ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t100ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t150ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t200ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t400ns"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t1us"]],[IrB_C13Bicarb_TRIR["Wavenumber"],IrB_C13Bicarb_TRIR["t2us"]]],r)
C13_10ns = IrB_C13Bicarb_TRIR_times[1]
C13_25ns = IrB_C13Bicarb_TRIR_times[2]
C13_50ns = IrB_C13Bicarb_TRIR_times[3]
C13_75ns = IrB_C13Bicarb_TRIR_times[4]
C13_100ns = IrB_C13Bicarb_TRIR_times[5]
C13_150ns = IrB_C13Bicarb_TRIR_times[6]
C13_200ns = IrB_C13Bicarb_TRIR_times[7]
C13_400ns = IrB_C13Bicarb_TRIR_times[8]
C13_1us = IrB_C13Bicarb_TRIR_times[9]
C13_2us = IrB_C13Bicarb_TRIR_times[10]
plot(r, IrB_C13Bicarb_TRIR_times[1:end])

IRSEC_data = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IRSEC.csv"))
IRSEC_Wavenumber = IRSEC_data["Wavenumber"]
IRSEC_30s = IRSEC_data["-1.3V 30s"]
IRSEC_90s = IRSEC_data["-1.3V 90s"]
IRSEC_0s = IRSEC_data["0V 0sec"]

IRSEC = GatherAllignedPoints([[IRSEC_Wavenumber,IRSEC_0s],[IRSEC_Wavenumber,IRSEC_30s],[IRSEC_Wavenumber,IRSEC_90s]],r)
plot(r,IRSEC[1:2],linewidth = 3)

dIr = IRSEC[2]-IRSEC[1]

IrB_bicarb_FTIR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_Bicarb_FTIR.csv"))
    IrB_Bicarb_FTIR_Ampl = GatherAllignedPoints([[IrB_bicarb_FTIR["Wavenumber"],IrB_bicarb_FTIR["Ampl"]]],r)
    plot(r, IrB_Bicarb_FTIR_Ampl[1])

IrB_bicarb_postFTIR1 = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_bicarb_FTIR_post1.csv"))
    IrB_bicarb_postFTIR1_Ampl = GatherAllignedPoints([[IrB_bicarb_postFTIR1["Wavenumber"],IrB_bicarb_postFTIR1["Ampl"]]],r)
    plot!(r,IrB_bicarb_postFTIR1_Ampl[1])

PPN_Bicarb_FITR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/PPN_Bicarb_FTIR.csv"))
    PPN_Bicarb = GatherAllignedPoints([[PPN_Bicarb_FITR["Wavenumber"],PPN_Bicarb_FITR["A"]]],r)    
    Bicarb_Ampl = PPN_Bicarb[1]
    plot!(r, Bicarb_Ampl)

Bicarb_coeff = bfcs(IrB_Bicarb_FTIR_Ampl[1],[Bicarb_Ampl])
    PPN_Bicarb_Fit = Bicarb_Ampl*Bicarb_coeff[1]
    plot!(r,PPN_Bicarb_Fit)


    Fitcoeff = Array([])

PPNHCO_FTIRs = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/PPNHCO3_FTIRs.csv"))
PPNHCO3 = GatherAllignedPoints([[PPNHCO_FTIRs["Wavenumber"],PPNHCO_FTIRs["50mM"]],[PPNHCO_FTIRs["Wavenumber"],PPNHCO_FTIRs["10mM"]],[PPNHCO_FTIRs["Wavenumber"],PPNHCO_FTIRs["5mM"]]],r)
plot(r, PPNHCO3[1:end])
write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wavenumber","5mM","10mM","50mM"],[Wvn,PPNHCO3[3],PPNHCO3[2],PPNHCO3[1]])

# for i in 1:1:10
#     coeff = bfcs(IrB_C13Bicarb_TRIR_times[i],[dIr])
#     push!(Fitcoeff,coeff)
# end

# C13_diff = []
# for i in 1:1:10
#     diff = IrB_C13Bicarb_TRIR_times[i].-dIr.*Fitcoeff[i]
#     push!(C13_diff,diff)
# end

for i in 1:1:10
    coeff = bfcs(IrB_C12Bicarb_TRIR_times[i],[IrB_C13Bicarb_TRIR_times[i]])
    push!(Fitcoeff,coeff)
end

C12_13_diff = []
for i in 1:1:10
    diff = IrB_C12Bicarb_TRIR_times[i]-IrB_C13Bicarb_TRIR_times[i].*Fitcoeff[i]
    push!(C12_13_diff,diff)
end
plot(r,C12_13_diff[3:5],size=(1600, 1000), linewidth =3)




plot(r,C13_diff[1:end],size=(1600, 1000), linewidth =3, xlims = (1250,1790),ticks = 50)


write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn","IrB_TRIR_100ns","IrB_TRIR_200ns","IrB_TRIR_300ns","IrB_TRIR_400ns","IrB_TRIR_500ns","IrB_TRIR_750ns","IrB_TRIR_950ns"],[Wvn,IrB_TRIR_100ns,IrB_TRIR_200ns,IrB_TRIR_300ns,IrB_TRIR_400ns,IrB_TRIR_500ns,IrB_TRIR_750ns,IrB_TRIR_950ns])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn", "C12_10ns","C12_25ns","C12_50ns","C12_75ns","C12_100ns","C12_150ns","C12_200ns","C12_400ns","C12_1000ns","C12_2000ns"],[Wvn, C12_10ns,C12_25ns,C12_50ns,C12_75ns,C12_100ns,C12_150ns,C12_200ns,C12_400ns,C12_1us,C12_2us])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn", "C13_10ns","C13_25ns","C13_50ns","C13_75ns","C13_100ns","C13_150ns","C13_200ns","C13_400ns","C13_1000ns","C13_2000ns"],[Wvn, C13_10ns,C13_25ns,C13_50ns,C13_75ns,C13_100ns,C13_150ns,C13_200ns,C13_400ns,C13_1us,C13_2us])

    write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn","IrB_Bicarb_FTIR_Ampl","IrB_bicarb_postFTIR_Amp","PPN_Bicarb_Fit"],[Wvn,IrB_Bicarb_FTIR_Ampl[1],IrB_bicarb_postFTIR1_Ampl[1],PPN_Bicarb_Fit])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn","IRSEC_0s","IRSEC_30s","IRSEC_90s","IRSEC_diff"],[Wvn,IRSEC[1],IRSEC[2],IRSEC[3],dIr])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn", "C12_10ns_diff","C12_25ns_diff","C12_50ns_diff","C12_75n_diff","C12_100ns_diff","C12_150ns_diff","C12_200ns_diff","C12_400n_diff","C12_1000ns_diff","C12_2000ns_diff"],[Wvn,C12_diff[1],C12_diff[2],C12_diff[3],C12_diff[4],C12_diff[5],C12_diff[6],C12_diff[7],C12_diff[8],C12_diff[9],C12_diff[10]])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wvn", "C13_10ns_diff","C13_25ns_diff","C13_50ns_diff","C13_75n_diff","C13_100ns_diff","C13_150ns_diff","C13_200ns_diff","C13_400n_diff","C13_1000ns_diff","C13_2000ns_diff"],[Wvn,C13_diff[1],C13_diff[2],C13_diff[3],C13_diff[4],C13_diff[5],C13_diff[6],C13_diff[7],C13_diff[8],C13_diff[9],C13_diff[10]])

k = 7

bv = IrB_C12Bicarb_TRIR_times[k] .> 0

IrB_C12Bicarb_TRIR_times[k][bv]
wvn1 = Wvn[bv]
plot(wvn1,IrB_C12Bicarb_TRIR_times[k][bv])

bv1 = IrB_C13Bicarb_TRIR_times[k] .> 0

IrB_C13Bicarb_TRIR_times[k][bv1]
wvn2 = Wvn[bv1]
plot!(wvn2,IrB_C13Bicarb_TRIR_times[k][bv1])