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
    weight_range::Float64 = 5.0 #⨦ This SHOULD be changed. This is the upper bound on the magnetiude of the coefficients you are searching for. Make it larger than you expect by a bit. The loops will narrow it down. 
    thresh::Float64 = -1e-9 # convergance is determined by the slope of the residual. When that negative slope flattens out, convergance is reached. 
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

r = 1250:1:1790
Wvn = collect(r)



IrB_bicarb_FTIR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_Bicarb_FTIR.csv"))
    IrB_bicarb_FTIR_Wavenumber = IrB_bicarb_FTIR["Wavenumber"]
    IrB_Bicarb_FTIR_Ampl = IrB_bicarb_FTIR["Ampl"]
    plot(IrB_bicarb_FTIR_Wavenumber, IrB_Bicarb_FTIR_Ampl)

IrB_bicarb_postFTIR1 = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_bicarb_FTIR_post1.csv"))
    IrB_bicarb_postFTIR1_Wavenumber = IrB_bicarb_postFTIR1["Wavenumber"]
    IrB_bicarb_postFTIR1_Ampl = IrB_bicarb_postFTIR1["Ampl"]
    plot!(IrB_bicarb_postFTIR1_Wavenumber,IrB_bicarb_postFTIR1_Ampl)


IrB_bicarb_postFTIR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_bicarb_postkessel_FTIR.csv"))
    IrB_bicarb_postFTIR_Wavenumber = IrB_bicarb_postFTIR["Wavenumber"]
    IrB_bicarb_postFTIR_Ampl = IrB_bicarb_postFTIR["Ampl"]
    plot!(IrB_bicarb_postFTIR_Wavenumber,IrB_bicarb_postFTIR_Ampl)

Bicarb_FTIR = GatherAllignedPoints([[IrB_bicarb_FTIR_Wavenumber,IrB_Bicarb_FTIR_Ampl],[IrB_bicarb_postFTIR_Wavenumber, IrB_bicarb_postFTIR_Ampl],[IrB_bicarb_postFTIR1_Wavenumber,IrB_bicarb_postFTIR1_Ampl]],r)
    IrB_Bicarb_pre = Bicarb_FTIR[1]
    IrB_Bicarb_post = Bicarb_FTIR[2]
    IrB_Bicarb_post1 = Bicarb_FTIR[3]
    IrB_bicarb_diff = (IrB_Bicarb_post-IrB_Bicarb_pre)
    plot(r,[IrB_Bicarb_pre,IrB_Bicarb_post1])
    plot!(r,IrB_bicarb_diff,linewidth = 3)


PPN_Bicarb_FITR = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/PPN_Bicarb_FTIR.csv"))
    PPN_Bicarb = GatherAllignedPoints([[PPN_Bicarb_FITR["Wavenumber"],PPN_Bicarb_FITR["A"]]],r)    
    Bicarb_Ampl = PPN_Bicarb[1]
    plot!(r, Bicarb_Ampl)

Bicarb_coeff = bfcs(IrB_Bicarb_pre,[Bicarb_Ampl])
    PPN_Bicarb_Fit = Bicarb_Ampl*Bicarb_coeff[1]
    plot(r,PPN_Bicarb_Fit)

plot(r,IrB_Bicarb_pre-PPN_Bicarb_Fit)

IrB_bicarb_spectra = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IrB_bicarb_spectra.csv"))
    IrB_bicarb_Wavenumber = IrB_bicarb_spectra["Wavenumber"]
    _50ns = IrB_bicarb_spectra["51.1ns"]
    _70ns = IrB_bicarb_spectra["70.3ns"]
    _120ns = IrB_bicarb_spectra["102ns"]
    _150ns = IrB_bicarb_spectra["148ns"]
    _290ns = IrB_bicarb_spectra["289ns"]
    _620ns = IrB_bicarb_spectra["619ns"]
    _2000ns = IrB_bicarb_spectra["2000ns"]

TRIR = GatherAllignedPoints([[IrB_bicarb_Wavenumber,_50ns],[IrB_bicarb_Wavenumber,_70ns],[IrB_bicarb_Wavenumber,_120ns],[IrB_bicarb_Wavenumber,_150ns],[IrB_bicarb_Wavenumber,_290ns],[IrB_bicarb_Wavenumber,_620ns],[IrB_bicarb_Wavenumber,_2000ns]], r)
plot!(r,TRIR[1:end])



IRSEC_data = CSV.File(read("/Users/kris/Desktop/RESEARCH/HARVARD/NOCERA GROUP.nosync/PROJECTS/Bicarbonate Oxidaion/SCRATCH/IRSEC.csv"))
IRSEC_Wavenumber = IRSEC_data["Wavenumber"]
IRSEC_30s = IRSEC_data["-1.3V 30s"]
IRSEC_90s = IRSEC_data["-1.3V 90s"]
IRSEC_0s = IRSEC_data["0V 0sec"]

IRSEC = GatherAllignedPoints([[IRSEC_Wavenumber,IRSEC_0s],[IRSEC_Wavenumber,IRSEC_30s],[IRSEC_Wavenumber,IRSEC_90s]],r)
plot!(r,IRSEC[1:1],linewidth = 3)

dIr = IRSEC[2]-IRSEC[1]
plot(r,dIr)



fit = dIr

_50nscoeff = bfcs(TRIR[1],[fit])
_50nsfit = _50nscoeff[1]*fit
d50 = TRIR[1]-_50nsfit
plot(r,[_50nsfit,TRIR[1],d50])

_70nscoeff = bfcs(TRIR[2],[fit])
_70nsfit = _70nscoeff[1]*fit
d70 = TRIR[2]-_70nsfit
plot(r,[_70nsfit,TRIR[2],d70])

_120nscoeff = bfcs(TRIR[3],[fit])
_120nsfit = _120nscoeff[1]*fit
d120 = TRIR[3]-_120nsfit
plot(r,[_120nsfit,TRIR[3],d120])

_150nscoeff = bfcs(TRIR[4],[fit])
_150nsfit = _150nscoeff[1]*fit
d150 = TRIR[4]-_150nsfit
plot(r,[_150nsfit,TRIR[4],d150])

_290nscoeff = bfcs(TRIR[5],[fit])
_290nsfit = _290nscoeff[1]*fit
d290 = TRIR[5]-_290nsfit
plot(r,[_290nsfit,TRIR[5],d290])

_620nscoeff = bfcs(TRIR[6],[fit])
_620nsfit = _620nscoeff[1]*fit
d620 = TRIR[6]-_620nsfit
plot(r,[_620nsfit,TRIR[6],d620])


_2000nscoeff = bfcs(TRIR[7],[fit])
_2000nsfit = _2000nscoeff[1]*fit
d2000 = TRIR[7]-_2000nsfit
plot(r,[_2000nsfit,TRIR[7],d2000])

plot(r,[d50,d70,d120,d150,d290,d620,d2000],size=(1600, 1000), linewidth =3, title = "best fit differences")
plot(r,[d120,d150,d290],size=(1600, 1000), linewidth =3, title = "best fit differences")


function Norm(x,y,z)
    x_interp = LinearInterpolation(x,y) #interpolate function for x and y vectors
    y_interp = x_interp(x) # compute y values for the x vector
    y_norm = y_interp/x_interp(z) #normalize to desired value
    return y_norm
end

z = 1492

norm_dIr = Norm(r, fit,z)
norm_70ns = Norm(r,TRIR[2],z)
norm_120ns = Norm(r,TRIR[3],z)
norm_150ns = Norm(r,TRIR[4],z)
norm_290ns = Norm(r,TRIR[5],z)
norm_620ns = Norm(r,TRIR[6],z)
norm_2000ns = Norm(r,TRIR[7],z)

plot(r, [norm_dIr,norm_70ns,norm_150ns,norm_290ns,norm_620ns,norm_2000ns], size=(1600, 1000), linewidth = 3)
plot(r, [ norm_150ns-norm_dIr,norm_290ns-norm_dIr,norm_620ns-norm_dIr,norm_2000ns-norm_dIr], size=(1600, 1000), linewidth = 3)

plot(r, [norm_120ns-norm_70ns, norm_150ns-norm_70ns,norm_dIr], linewidth = 3)


plot(r, [IrB_Bicarb_pre*1.25,IRSEC[1],IRSEC[2],IrB_Bicarb_post*1.25])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wavenumber","d50","d70","d120","d150","d290","d620","d2000"],[Wvn,d50,d70,d120,d150,d290,d620,d2000])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wavenumber","PPN_Bicarb_Fit","IrB_Bicarb_pre","IrB_bicarb_FTIR_diff","IRSEC_0s","IRSEC_90s","IRSEC_diff"],[Wvn,PPN_Bicarb_Fit,IrB_Bicarb_pre,IrB_bicarb_diff,IRSEC[1],IRSEC[3],dIr])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wavenumber","50ns","70ns","120ns","150ns","290ns","620ns","2000ns"],[Wvn,TRIR[1],TRIR[2],TRIR[3],TRIR[4],TRIR[5],TRIR[6],TRIR[7]])

write_csv("/Users/kris/Documents/JULIA/Outputs/newCSV.csv",["Wavenumber","IrB_Bicarb_pre","IrB_Bicarb_post","PPN_Bicarb_Fit"],[Wvn,IrB_Bicarb_pre,IrB_Bicarb_post1,PPN_Bicarb_Fit])