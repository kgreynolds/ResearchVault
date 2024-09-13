using Plots, DataFrames, CSV,Roots, ColorSchemes,Polynomials,SpecialFunctions,LaTeXStrings,CairoMakie; CairoMakie.activate!()
plotly()


# PLOT_PROB=false plots the wavefunction, psi; PLOT_PROB=true plots |psi|^2
PLOT_PROB = true

# Maximum vibrational quantum number to calculate wavefunction for
VMAX = 15
α = 1
# Some appearance settings
# Pad the q-axis on each side of the maximum turning points by this fraction
QPAD_FRAC = 0.1
# Scale the wavefunctions by this much so they don't overlap
SCALING = 1
# Colours of the positive and negative parts of the wavefunction
COLOUR1 = RGBA(0.6196, 0.0039, 0.2588, 1.0)
COLOUR2 = RGBA(0.3686, 0.3098, 0.6353, 1.0)

# Normalization constant and energy for vibrational state v
N(v) = 1 ./ sqrt(2^v * factorial(v)) * sqrt(sqrt(α/pi))
get_E(v) = v + 0.5


function make_Hr()
    """Return a list of np.poly1d objects representing Hermite polynomials."""

    # Define the Hermite polynomials up to order VMAX by recursion:
    # H_[v] = 2qH_[v-1] - 2(v-1)H_[v-2]
    Hr = Vector{Any}(undef, VMAX + 1)
    Hr[1] = Polynomial([1.0])
    Hr[2] = Polynomial([0, 2])
    for v in 3:VMAX+1
        Hr[v] = Hr[2]*Hr[v-1] - 2*(v-2)*Hr[v-2]
    end
    return Hr
end

Hr = make_Hr()

function get_psi(v, x)
    """Return the harmonic oscillator wavefunction for level v on grid q."""
    return N(v) * Hr[v+1].(sqrt(α)*x) .* exp.(-α .* x.^2/2)
end

function get_potential(x)
    """Return potential energy on scaled oscillator displacement grid x."""
    return α*x.^2 / 2
end

function MorsePotential(D,β,x)
    """Construct a morse potential over x"""
    return D * (1 .- exp.(-β * x)).^2
end

function Anharmonicity(χ,v)
    """Return the anharmonicity factor(χ) correcion for the vibrational energies"""
    return Polynomial([1,-χ,χ,-χ])(v+0.5)
end



x = collect(range(-25,25,length = 500))
V = get_potential(x)

# V = MorsePotential(10,α/pi,x)

MorsePotential(1,1,x)
morse = Plots.plot!(x,MorsePotential(10,α/pi,x))

p= Plots.plot()

p = Plots.plot!(x,V,color=:black,label="",linewidth=1.5,ylims = (-1,get_E(VMAX)+0.5))

for v in 0:VMAX
    psi_v = get_psi(v, x)
    E_v = get_E(v)
    n = v+1
    if PLOT_PROB == true
     display(Plots.plot!(p,x,psi_v.^2 .*SCALING.*1.5 .+ E_v,color=COLOUR1,linewidth = 2,label=""))
     display(Plots.plot!(p,x,E_v .+ zeros(length(x)),fillrange=E_v,color=COLOUR2,alpha=0.5,label=""))
    else
        display(Plots.plot!(x,psi_v.*SCALING .+ E_v,color=COLOUR1,linewidth = 2,label=""))
        display(Plots.plot!(x,E_v .+ zeros(length(x)),fillrange=E_v,color=COLOUR2,alpha=0.2,label=""))
        # display(Plots.plot!(x,(Anharmonicity(0.2,v) + E_v -1) .+ zeros(length(x)),color = :black,linewidth = 2))
    end
    
end


ymax = get_E(VMAX)+0.5

if PLOT_PROB == true
    ylabel = "|ψ(q)|²"
else
    ylabel = "ψ(q)"
end


savefig(p,"sho-psi$(PLOT_PROB+1)-$(VMAX).pdf")
display(p)


