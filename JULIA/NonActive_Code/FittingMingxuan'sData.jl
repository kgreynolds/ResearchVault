# further fitting of mingxuan's data on mac

#for Mingxuan
using Plots, Statistics, Optim, Peaks, CSV, CairoMakie, LaTeXStrings, DataFrames, ColorSchemes, GeometryBasics

mutable struct PVoight
    μ::Float64 #the center of the peak
    Γ::Float64 # the peak width parameter       
    η::Float64 # the ratio of gaussian to lorentzian contribution to the peak 
    scale::Float64 # peak intensity
end

#functions for fitting convolutions

function G(x,f)
    return exp(-x^2 / 2f^2) / (f*√2π)
end

function L(x,f)
    return f / (π*(x^2 + f^2))
end

function pv(x,fg,fl)
    return η(fg,fl)*L(x,fwhm(fg,fl)) + (1 - η(fg,fl))*G(x,fwhm(fg,fl))
end

function pv(x,η)
    return η*L(x,1) + (1 - η)*G(x,1)
end

function npv(x,η)
    result = pv.(x,η)
    return result / maximum(result)
end

function pv(x,η,Γ)
    return η*L(x,Γ) + (1 - η)*G(x,Γ)
end

function npv(x,η,Γ)
    result = pv.(x,η,Γ)
    return result / maximum(result)
end

function peak(x::Vector,p::PVoight)
    return npv(x .- p.μ, p.η, p.Γ) * p.scale
end

function peak(x::Float64,p::PVoight)
    return npv(x - p.μ, p.η, p.Γ) * p.scale
end

function con(t::Vector{Float64},τ::Float64,irf::PVoight)
    out::Vector{Float64} = []
    for (i,v) in enumerate( peak(t,irf) )
        push!(out, v*sum(  lt(t,t[i],τ) ))
    end
    return out / length(t)
end

function lt(x::Float64,t::Int64)
    if x <= t
        return 0.0
    else
        return exp(-(x)/t) 
    end
end

function lt(x::Vector{Float64},t0::Float64,t::Float64)
    out::Vector{Float64} = []
    for v in x
        if v <= t0
            push!(out, 0.0)
        else
            push!(out, exp(-(v - t0)/t) )
        end
    end
    return out
end

function con(t::Vector{Float64},τ::Float64,irf::PVoight)
    out::Vector{Float64} = []
    #peak(t,irf) 
    #lt(t,delay,τ) 

    #make the set of lifetimes, then weight them, and sum them
    lifetimes::Vector{Vector{Float64}} = []#[lt(t,delay,τ) for delay in t,time in t]
    #=
    for delay in t
        push!(lifetimes, lt(t,delay,τ))
    end
    =#
    N = length(t)


    for (i,weight) in enumerate(peak(t,irf))
        delay = t[i]
        push!(lifetimes, weight * lt(t,delay,τ))
    end
    
    #=
    for (i,weight) in enumerate(peak(t,irf))
        delay = t[i]
        push!(out, sum(weight * lt(t,delay,τ)))
    end
    =#

    m = hcat(lifetimes...)

    for i in collect(axes(m)[1])
        push!(out,sum(m[i,:]))
    end

    return out ./ maximum(out)#(sum(m)/N)

end

function rconv(t::Vector{Float64},τ::Float64,irf::Vector{Float64},scale::Float64) # this convolution method is written to take any curve as the "IRF"
    out::Vector{Float64} = []
    lifetimes::Vector{Vector{Float64}} = []
 
    for (i,weight) in enumerate(irf)
        delay = t[i]
        push!(lifetimes, weight * lt(t,delay,τ))
    end

    m = hcat(lifetimes...)

    for i in collect(axes(m)[1])
        push!(out,sum(m[i,:]))
    end

    return out * (scale/maximum(out))

end




mxd = "/Users/evenson/Downloads/QD560 lifetime.csv"
msd_ns = "/Users/evenson/Downloads/20230125 nanosecond laser.csv"
mxd2 = "/Users/evenson/Downloads/QD540SRB_Streak camera.csv"

#data = readRJE.readGrII(mxd2)

m = Matrix(CSV.read(mxd2,DataFrame))
names(CSV.read(mxd2,DataFrame))
t = m[:,1]; i = m[:,2]

fx = t#py100n.time
fy = i

Plots.plot(fx,fy, label = "data")#,yscale = :log)
i1 = 1:300
i2 = 220:479
#Plots.scatter(fx,fy, label = "data", yscale = :log)
Plots.scatter!(fx[i1],fy[i1], label = "data")
#Plots.scatter!(fx[i2],fy[i2], label = "data")

nx = fx[i1]#vcat(fx[i1],fx[i2])
ny = fy[i1]#vcat(fy[i1],fy[i2])
using Peaks
pks, vals = findmaxima(ny)

pks, proms = peakproms(pks, ny; strict = true, minprom  = 200, maxprom = nothing)

pks, widths, leftedge, rightedge = peakwidths(pks, ny, proms)


Plots.scatter(nx,ny, label = "data")
Plots.scatter!(nx[pks],ny[pks], label = "peak", marker = :x, markersize = 10)

#modelpeak(p) = con(fx,p[5],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).*p[6] .+ p[4]
modelpeak(p) = con(nx,p[5],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).*p[4] .+ p[6] .+ con(nx,p[7],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).*p[8] #biexp model
modelpeak(p) = peak(nx,PVoight(p[10],abs(p[5]),abs(sin(p[9])),1)).*p[4] .+ p[6] .+ con(nx,p[7],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).*p[8] # monoexp and scattering
peakfit0(p) = sum((ny .- modelpeak(p)).^2)

#modelpeak(p) = con(nx,p[5],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).*abs(p[6]) .+ con(nx,p[4],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).* abs(p[7]) .+ 25


p0 = rand(8) 
p0[1] = nx[pks][1] - 0.2 
p0[2] *= 0.5#rand()
p0[4] = 500.
p0[5] = 20 # t1
p0[6] = 20 # constant
p0[7] = 0.2 # t2
p0[8] = 8000 # scale

p0 = rand(10) 
p0[1] = nx[pks][1] - 0 
p0[2] = 0.1
p0[4] = 990. #scale scat
p0[5] = 0.29 # width scat
p0[6] = 0 # constant
p0[7] = 10. # t2
p0[8] = 300 # scale decay
p0[9] = 0.0#rand() # LG scat
p0[10] = nx[pks][1]- 0.015# location of scattering peak



#p0[6] = 150.#700. * rand()
#p0[7] = 0#150.#700. * rand()
#p0[8] = 20.
p0

Plots.plot!(nx,modelpeak(p0), label = "initial guess")
Plots.plot!(nx,ny .- modelpeak(p0), label = "initial guess diff")

res3 = optimize(peakfit0,p0,LBFGS(),Optim.Options(time_limit = 120))
gpside = Optim.minimizer(res)
gp2 = deepcopy(Optim.minimizer(res3))
initres = deepcopy(Optim.minimizer(res3))

gp2[7] = 0

irf20ns = deepcopy(PVoight(gp2[1],abs(sin(gp2[2])),abs(sin(gp2[3])),1) )
irf20ns

Plots.scatter(nx,ny, label = false,marker = :+)
#Plots.plot!(x,modelpeak(p0), label = "initial guess")
Plots.plot!(nx,modelpeak(gp2), 
label = "fitted peak\nμ = $(round(gp2[1], digits = 2)) ns\nτ1 = $(round(gp2[7],digits = 4)) ns\nτ2 = $(round(gp2[5],digits = 4)) ns\nη = $(round(abs(sin(gp2[3])),digits = 2)), Γ = $(round(abs(gp2[2]),digits = 2))")
Plots.plot!(nx,ny .- modelpeak(gp2) .- 1000, label = "diff", color = :black)
Plots.plot!(nx,peak(nx,PVoight(gp2[1],abs(gp2[2]),abs(sin(gp2[3])),1) ).* gp2[4] .- 2500, label = "fitted IRF")
Plots.plot!(nx,lt(nx,gp2[1],gp2[7]).* gp2[8] .- 000, label = "decay 1")
Plots.plot!(nx,lt(nx,gp2[1],gp2[5]).* gp2[4] .- 000, label = "decay 2")
lt

Plots.xlims!(0,25)
Plots.ylims!(-1000,2800)

Plots.scatter(nx,ny, label = false,marker = :+)
#Plots.plot!(x,modelpeak(p0), label = "initial guess")
Plots.plot!(nx,modelpeak(gp2), 
label = "fitted peak\nμ = $(round(gp2[1], digits = 2)) ns\nη = $(round(abs(sin(gp2[3])),digits = 2)), Γ = $(round(abs(gp2[2]),digits = 2))")
Plots.plot!(nx,ny .- modelpeak(gp2) .- 1000, label = "diff", color = :black)
Plots.plot!(nx,peak(nx,PVoight(gp2[1],abs(gp2[2]),abs(sin(gp2[3])),1) ).* gp2[8] .- 2500, label = "fitted IRF")
Plots.plot!(nx,lt(nx,gp2[1],gp2[7]).* gp2[8] .- 000, label = "decay, τ1 = $(round(gp2[7],digits = 4)) ns")
Plots.plot!(nx,peak(nx,PVoight(gp2[10],abs(gp2[5]),abs(sin(gp2[9])),1) ).* gp2[4] .- 12000, label = "scattering")
Plots.xlims!(0,30)

[gp2[1],abs(gp2[2]),abs(sin(gp2[3]))]
[gp2[10],abs(gp2[5]),abs(sin(gp2[9]))]
gp2[1] - gp2[10]

plotpeak(x,p) = peak(x,PVoight(p[10],abs(p[5]),abs(sin(p[9])),1)).*p[4] .+ p[6] .+ con(x,p[7],PVoight(p[1],abs(p[2]),abs(sin(p[3])),1)).*p[8]

function mxp(xpx,ypx,data,gp,head)
    ix, iy = data
    
    fnt = "Serif Roman"
    figure = ( 
        resolution=(xpx, ypx),
        font = fnt,#font="C:\\Users\\R. Evenson\\Downloads\\nimbus-roman-no9-l.regular.otf",
        fontsize = 20)

    axis = (
        xlabel = L"Time (ns)",
        ylabel = L"\int{I(nm,ns)} \thinspace \partial nm", #L"∫ PL  Intensity _{a}^{b}   (counts)",
        title = "Quantum Dot Emission Modeling",#data.title,
        titlefont = fnt,
        yminorticksvisible = true,
        xminorticksvisible = true,
        #yscale = log10,
        #xscale = log10,
        xminorticks = IntervalsBetween(10),
        yminorticks = IntervalsBetween(4)
        )

        cm = :viridis#reverse(ColorSchemes.roma)
        c = ColorSchemes.Reds_3[3]
        #lim = (minimum(data.matrix),maximum(data.matrix))

        
        #x = data.time#range(0,stop = (5.5)*length(data.matrix[:,1]), length = length(data.matrix[:,1]))
        #y = sum(data.matrix,dims = 2)[:,1]
        #y /= maximum(y) 
        #x .-= 2
        #return x , y
    #fig, ax, lin = CairoMakie.heatmap(x,y,data.matrix; colormap = cm, axis=axis, figure=figure)
    lw = 1
    c = ColorSchemes.Reds[6]
    
    fig, ax, lin = CairoMakie.scatter(ix,iy; marker =  '▢',markersize = 20, color = c, axis=axis, figure=figure)
    lin2 = CairoMakie.lines!(ix,plotpeak(ix,gp); color = :gray,linewidth = lw*3, linestyle = :solid)
    lin3 = CairoMakie.lines!(ix,iy .- plotpeak(ix,gp) .- 500; color = c,linewidth = lw)
    iv = 1:50
    lin4 = CairoMakie.lines!(ix[iv],peak(ix[iv],PVoight(gp[1],abs(gp[2]),abs(sin(gp[7])),1) ).* gp[8] .- 1200; color = :green,linewidth = lw)
    lin5 = CairoMakie.lines!(ix,lt(ix,gp[1],gp[7]).* gp[8] .- 000; color = :black,linewidth = lw*2)
    lin6 = CairoMakie.lines!(ix[iv],peak(ix[iv],PVoight(gp[10],abs(gp[5]),abs(sin(gp[9])),1) ).* gp[4] .- 2400; color = :green,linewidth = lw, linestyle = :dash)

    m1 = "$head, streak camera data"
    m2 = "Fitted model (scattering + kinetics)"
    m3 = "Residual"
    m4 = "Fitted IRF"
    m5 = "Kinetics τ = $(round(gp[7],digits = 2)) ns"
    m6 = "Scattering contribution"
    
    ms = [m1,m2,m3,m4,m5,m6]
    CairoMakie.axislegend(ax, [lin,lin2,lin3,lin4,lin5,lin6],ms,position = :rb)

    fig

end

mxp(850,600)

# next fit using the coef of the first fit 

t = m[:,1]; i = m[:,3]

fx = t#py100n.time
fy = i

Plots.plot(fx,fy, label = "data")#,yscale = :log)
i1 = 1:300
i2 = 1:300
#Plots.scatter(fx,fy, label = "data", yscale = :log)
Plots.scatter!(fx[i1],fy[i1], label = "data")
Plots.plot!(nx,modelpeak(gp2), label = "initial guess")

nx = fx[i1]#vcat(fx[i1],fx[i2])
ny = fy[i1]#vcat(fy[i1],fy[i2])
using Peaks
pks, vals = findmaxima(ny)

pks, proms = peakproms(pks, ny; strict = true, minprom  = 200, maxprom = nothing)

pks, widths, leftedge, rightedge = peakwidths(pks, ny, proms)


Plots.scatter(nx,ny, label = "data")
Plots.scatter!(nx[pks],ny[pks], label = "peak", marker = :x, markersize = 10)

res4 = optimize(peakfit0,gp2,LBFGS(),Optim.Options(time_limit = 120))
gp3 = deepcopy(Optim.minimizer(res4))

heads = names(CSV.read(mxd2,DataFrame))
[1:300]

mxp(850,600,[m[:,1][1:300],m[:,2][1:300]],gp2,heads[2])
mxp(850,600,[m[:,1][1:300],m[:,3][1:300]],gp3,heads[3])
mxp(850,600,[m[:,1][1:300],m[:,4][1:300]],gp4,heads[4])
mxp(850,600,[m[:,1][1:300],m[:,5][1:300]],gp5,heads[5])
mxp(850,600,[m[:,1][1:300],m[:,6][1:300]],gp6,heads[6])
mxp(850,600,[m[:,1][1:300],m[:,7][1:300]],gp7,heads[7])

t = m[:,1]; i = m[:,4]

fx = t#py100n.time
fy = i

Plots.plot(fx,fy, label = "data")#,yscale = :log)
i1 = 1:300
i2 = 1:300
#Plots.scatter(fx,fy, label = "data", yscale = :log)
Plots.scatter!(fx[i1],fy[i1], label = "data")
Plots.plot!(nx,modelpeak(gp2), label = "initial guess")

nx = fx[i1]#vcat(fx[i1],fx[i2])
ny = fy[i1]#vcat(fy[i1],fy[i2])
using Peaks
pks, vals = findmaxima(ny)

pks, proms = peakproms(pks, ny; strict = true, minprom  = 200, maxprom = nothing)

pks, widths, leftedge, rightedge = peakwidths(pks, ny, proms)


Plots.scatter(nx,ny, label = "data")
Plots.scatter!(nx[pks],ny[pks], label = "peak", marker = :x, markersize = 10)

res5 = optimize(peakfit0,gp3,LBFGS(),Optim.Options(time_limit = 120))
gp4 = deepcopy(Optim.minimizer(res5))

t = m[:,1]; i = m[:,5]

fx = t#py100n.time
fy = i

Plots.plot(fx,fy, label = "data")#,yscale = :log)
i1 = 1:300
i2 = 1:300
#Plots.scatter(fx,fy, label = "data", yscale = :log)
Plots.scatter!(fx[i1],fy[i1], label = "data")
Plots.plot!(nx,modelpeak(gp3), label = "initial guess")

nx = fx[i1]#vcat(fx[i1],fx[i2])
ny = fy[i1]#vcat(fy[i1],fy[i2])
using Peaks
pks, vals = findmaxima(ny)

pks, proms = peakproms(pks, ny; strict = true, minprom  = 200, maxprom = nothing)

pks, widths, leftedge, rightedge = peakwidths(pks, ny, proms)


Plots.scatter(nx,ny, label = "data")
Plots.scatter!(nx[pks],ny[pks], label = "peak", marker = :x, markersize = 10)

res6 = optimize(peakfit0,gp6,LBFGS(),Optim.Options(time_limit = 120)) # inital fitting is susepct
gp5 = deepcopy(Optim.minimizer(res6))

t = m[:,1]; i = m[:,6]

fx = t#py100n.time
fy = i

Plots.plot(fx,fy, label = "data")#,yscale = :log)
i1 = 1:300
i2 = 1:300
#Plots.scatter(fx,fy, label = "data", yscale = :log)
Plots.scatter!(fx[i1],fy[i1], label = "data")
Plots.plot!(nx,modelpeak(gp2), label = "initial guess")

nx = fx[i1]#vcat(fx[i1],fx[i2])
ny = fy[i1]#vcat(fy[i1],fy[i2])
using Peaks
pks, vals = findmaxima(ny)

pks, proms = peakproms(pks, ny; strict = true, minprom  = 200, maxprom = nothing)

pks, widths, leftedge, rightedge = peakwidths(pks, ny, proms)


Plots.scatter(nx,ny, label = "data")
Plots.scatter!(nx[pks],ny[pks], label = "peak", marker = :x, markersize = 10)

res7 = optimize(peakfit0,gp4,LBFGS(),Optim.Options(time_limit = 120))
gp6 = deepcopy(Optim.minimizer(res7))

t = m[:,1]; i = m[:,7]

fx = t#py100n.time
fy = i

Plots.plot(fx,fy, label = "data")#,yscale = :log)
i1 = 1:300
i2 = 1:300
#Plots.scatter(fx,fy, label = "data", yscale = :log)
Plots.scatter!(fx[i1],fy[i1], label = "data")
Plots.plot!(nx,modelpeak(gp6), label = "initial guess")

nx = fx[i1]#vcat(fx[i1],fx[i2])
ny = fy[i1]#vcat(fy[i1],fy[i2])
using Peaks
pks, vals = findmaxima(ny)

pks, proms = peakproms(pks, ny; strict = true, minprom  = 200, maxprom = nothing)

pks, widths, leftedge, rightedge = peakwidths(pks, ny, proms)


Plots.scatter(nx,ny, label = "data")
Plots.scatter!(nx[pks],ny[pks], label = "peak", marker = :x, markersize = 10)

res8 = optimize(peakfit0,gp6,LBFGS(),Optim.Options(time_limit = 120))
gp7 = deepcopy(Optim.minimizer(res8))
