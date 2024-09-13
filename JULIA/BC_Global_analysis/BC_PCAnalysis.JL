# Brute Force Coefficient Search . J L 
# Author: Ryan James Evenson; G5 Nocera Research Group, Harvard CCB evenson@g.harvard.edu

using CSV, Plots, DataFrames,Interpolations, StatsBase #,StaticArrays Haven't yet done static arrays yet. 

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
    grid::Int64 = 10 # how fine the search mesh is, 2grid+1 is the side length of each dimension
    V = repeat([2grid+1],N)
    dims = Tuple(x for x in V) #just some hacky shit that works for generalizing this code to N components being fit
    # for some reason I had trouble writing it in a saner looking way. 

    coeffs::Array{Float64,N} = Array{Float64,N}(undef,dims) # this is the combination space area/volume/hypervolume/whatever we are searching over
    #all of these undefined values are rapidly overwritten. No undefined behavior should result. 

    cycle::Int64 = 1 # a starting point, don't change this. This increments as cycles complete. 
    converged::Bool = false # keeps track of whether or not convergance has been reached. Checked at the end of every cycle.


    #*******************
    # change me if this isn't working
    weight_range::Float64 = 10.0 #⨦ This SHOULD be changed. This is the upper bound on the magnetiude of the coefficients you are searching for. Make it larger than you expect by a bit. The loops will narrow it down. 
    thresh::Float64 = -1e-6 # convergance is determined by the slope of the residual. When that negative slope flattens out, convergance is reached. 
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

Agbpy2_Solution = CSV.File(read("/Users/kris/Documents/JULIA/BC_Global_analysis/Agbpy2_solution.csv"))
W_Agbpy2_Solution = Agbpy2_Solution["w"]
A_Agbpy2_Solution = Agbpy2_Solution["A"]

Agbpy2TFA_dr = CSV.File(read("/Users/kris/Documents/JULIA/BC_Global_analysis/Agbpy2TFA_diffusereflectance.csv"))  
W_Agbpy2TFA_dr = Agbpy2TFA_dr["w"]
F_Agby2TFA_dr = Agbpy2TFA_dr["F"]

Agbpy2TFA_solution = CSV.File(read("/Users/kris/Documents/JULIA/BC_Global_analysis/Agbpy2TFA_solution.csv"))
W_Agbpy2TFA_solution = Agbpy2TFA_solution["w"]
A_Agbpy2TFA_solution = Agbpy2TFA_solution["Abs"]

AgbpyTFA2_dr = CSV.File(read("/Users/kris/Documents/JULIA/BC_Global_analysis/AgbpyTFA2_diffusereflectance.csv"))
W_AgbpyTFA2_dr = AgbpyTFA2_dr["w"]
F_AgbpyTFA2_dr = AgbpyTFA2_dr["F"]

AgbpyTFA2_Solution = CSV.File(read("/Users/kris/Documents/JULIA/BC_Global_analysis/AgbpyTFA2_solution.csv"))
w_AgbpyTFA2_Solution = AgbpyTFA2_Solution["Wavelength"]
A_AgbpyTFA2_Solution = AgbpyTFA2_Solution["Abs"]

Range = 600:1:800
gp = GatherAllignedPoints([[W_Agbpy2_Solution,A_Agbpy2_Solution], [W_Agbpy2TFA_dr, F_Agby2TFA_dr],[W_Agbpy2TFA_solution,A_Agbpy2TFA_solution],[W_AgbpyTFA2_dr,F_AgbpyTFA2_dr],[w_AgbpyTFA2_Solution,A_AgbpyTFA2_Solution]],Range)

AgBPY2_S = gp[1]
AgBpy2TFA_dr = gp[2]
AgBpy2TFA_S = gp[3]
AgBpyTFA2_dr = gp[4]
AgBpyTFA2_S = gp[5]



Coeff = bfcs(AgBpy2TFA_S,[AgBpy2TFA_dr,AgBPY2_S, AgBpyTFA2_dr])

fit = Coeff[1]*AgBpy2TFA_dr+Coeff[2]*AgBPY2_S +Coeff[3]*AgBpyTFA2_dr



plot(Range, AgBpy2TFA_S)
plot!(Range, fit)