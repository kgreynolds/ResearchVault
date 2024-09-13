using Interpolations,Plots, DataFrames, CSV, Optim, ColorSchemes, Statistics,LineSearches ,CairoMakie; CairoMakie.activate!()
plotly()

mₑ = 9.1093837*10^(-31)
h = 6.62607015*10^(-34)
ħ = h/(2*π)
e = 1.60217663*10^(-19)
c = 299792458
λ = 6000*10^(-9)
ν = c/(λ)
ϵ₀ = 8.85418782*10^(-12)
N = 10^6

Δν = 20
ϵ = 1000


int = Δν*ϵ
# oscillator strength from linewidth
f = 4.32*10^(-9)*int

eps = 0.2/(4.32*10^(-9)*5000)


# Prefactor for oscillator strength QM
cont = (8*π^2*ν*mₑ)/(3*h*e^2)

# Transition dipole moment Units: Cm
μ = √((f/cont))

# Convert to Debye
D = μ/(3.33564*10^(-30)) # in DEbye

# Vaccum Cavity mode electric field
E = √(h*ν/(2*ϵ₀*(λ^3)))

# coupling strength converted to eV ( μ*E*√(N) is in units of J)
g = μ*E*√(N)*6.24*10^(18)

splitting = round(2*g;digits = 4)

# Calculate extinction coeff for IR absorption intesnities in units of KM/mol
Intensity = 23
linewidth = 2
epsilonIR = 27.648*(Intensity/linewidth)

# number of molecules in some concentration for a mode-volume
N = 0.01*(6.23e23)*(λ^3/0.001)

function Polariton_Coupling(N,λnm,f)
    #N = number of molecules coupled to
    #λnm is the wavelength of the transition
    #f is the oscillator strength of the transition

    # mₑ = 9.1093837*10^(-31)
    # h = 6.62607015*10^(-34)
    # ħ = h/(2*π)
    # e = 1.60217663*10^(-19)
    # c = 299792458
    # λ = λnm*10^(-9)
    # ν = c/(λ)
    # ϵ₀ = 8.85418782*10^(-12)

    mₑ = 9.1093837e-31
    h = 6.62607015e-34
    ħ = h/(2*π)
    e = 1.60217663e-19
    c = 299792458
    λ = λnm*10^(-9)
    ν = c/(λ)
    ϵ₀ = 8.85418782e-12
    
    # Prefactor for oscillator strength QM
    cont = (8*π^2*ν*mₑ)/(3*h*e^2)

    # Transition dipole moment Units: Cm
    μ = √((f/cont))

    # Convert to Debye
    D = μ/(3.33564e-30) # in DEbye

    # Vaccum Cavity mode electric field
    E = √(h*ν/(2*ϵ₀*(λ^3)))

    # coupling strength converted to meV ( μ*E*√(N) is in units of J)
    g = μ*E*√(N)

    # Rabi Splitting in meV
    Ω = 2*g*6.24e18 *1000

    return μ ,round(D;digits = 3),round(E;digits =3),g,round(Ω;digits = 3)
end

Polariton_Coupling(10^7,300,1)

Polariton_Coupling(10^7,2941,0.1)

function CavityDispersion(n,λ,m)
    h = 6.62607015*10^(-34)
    ħ = h/(2*π)
    c = 299792458
    θ = (-π*0.25):0.001:(π*0.25)

    kₛ = n*m*π/λ
    kₚ = (2*π/λ) .* sin(θ)

    E = (ħ*c/n) .* √(kₛ^2 .* kₚ^2)
end

n = 1
λ = 10000*10^(-9)
m = 1
h = 6.62607015*10^(-34)
ħ = h/(2*π)
c = 299792458
θ = (-π*0.25):0.0001:(π*0.25)
Eₘ = 0.06/  6.242e18
Γ = 0

kₛ = n*m*π/λ
kₚ = []
for i in θ
    out = (2*π/λ) .* sin(i)
    push!(kₚ,out)
end

kₚ

E = (ħ*c/n) * .√(kₛ .^2 .+ kₚ .^ 2) 



gₙ = Polariton_Coupling(10^8,10000,0.5)[4]


Eup = (E .+ Eₘ .- im .* ħ .* Γ) ./ 2 .+ (.√((gₙ .^2) .+ 0.25 .* ((Eₘ .- E) .+ im .* ħ .* Γ) .^2 ))
Elp = (E .+ Eₘ .- im .* ħ .* Γ) ./ 2 .- (.√((gₙ .^2) .+ 0.25 .* ((Eₘ .- E) .+ im .* ħ .* Γ) .^2 ))

Plots.plot(θ,E* 6.242e18)
Plots.plot!(θ,real(Eup* 6.242e18))
Plots.plot!(θ,real(Elp* 6.242e18))

function Fig1(f = Figure())
    width = 4

    #Make Cairo Mackie figure
        ax = CairoMakie.Axis(f[1,1], 
            palette = (color = palette(ColorSchemes.hsv, 5),), xautolimitmargin = (0.0, 0.0),
            xgridvisible = false, xminorticksvisible = true, xtickalign = 1, xminortickalign = 1, xticksize = 10, xlabelsize = 20,
            ygridvisible = false, yminorticksvisible = true, ytickalign = 1, yminortickalign = 1, yticksize = 10, ylabelsize = 20,
            xlabel = "θ", 
            ylabel = "E (meV)"
        )
        # ax.xticks= 1:2:17
        # CairoMakie.xlims!(1250,1790)
        # CairoMakie.ylims!(nothing,)
        
        lines!(θ,E* 6.242e18,linewidth = width,label = "Cavity")
        lines!(θ,real(Eup* 6.242e18),linewidth = width,label = "UP")
        lines!(θ,real(Elp* 6.242e18),linewidth = width,label = "LP")

        hlines!(ax,0.06, linewidth = 1,color = :black,linestyle = :dash)

        axislegend( position = :ct,nbanks = 1,framevisible = false, labelsize = 20) #:rt (right-top :rb (right bottom) etc)

    colsize!(f.layout,1,Aspect(1, 1.25)) #Set aspect ration of the y vs. x axis) can also be set to Relative(2/3) for leative column size
    resize_to_layout!(f)
    f
end
Fig1()
save("JULIA/Outputs/Dispersion.pdf",Fig1())