include("/Users/kris/Desktop/ResearchVault/JULIA/Main Functions/Analysis_Functions.jl")

#### See Robben, K. C.; Cheatum, C. M. Increasing Pump-Probe Signal toward Asymptotic Limits. J. Phys. Chem. B 2023, 127 (21), 4694–4707.
#### for detials 
  


h = 6.62607015e-34 #Plancks Constant
c = 2.99792458e+8 # Speed of light
NA = 6.0221408e23 # Avogadros number
Wavelength = 450

E_photon = h*c/(Wavelength*10^(-9)) # Photon Energy in J
Spotsize = 0.5 # Diameter of pump beam in cm
ϵ = 1305 # molar absorptivity in M⁻¹ cm⁻¹
OD_init = 0.1
E_pump = 0.002 # energy of pump pulse in J


Vobs = (π*(Spotsize/2)^2*1) /1000 ## observation volume in L
n = OD_init/ϵ * Vobs * NA # mols of initial absorber
N_photon = E_pump/E_photon
# N_photon = E_pump/(E_photon*NA)*50

r = n/N_photon # Ratio of absorbers to photons 

F = (E_pump/(π*(Spotsize/2)^2)) # Avg Fluence in J/ cm²
σ = log(10)*1000*ϵ/NA #Absorption crosssection of absorber
Fsat = E_photon/σ  # Saturation Fluence of the absorber

Saturation = F/Fsat

## Plot of OD vs Saturation ratio 
E_pumps = collect(0:0.0001:2)
F_1 = @. (E_pumps/(π*(Spotsize/2)^2)) # Avg Fluence in J/ cm²
σ_1 = log(10)*1000*ϵ/NA #Absorption crosssection of absorber

Fsat_1 = E_photon/σ

Saturation_1 = F_1/Fsat_1

OD_sat = exp.(-(Saturation_1))
OD_sat[1:5]
Plots.plot(log.(Saturation_1), OD_sat)
Plots.scatter!(log.(Saturation_1[1:5]), OD_sat[1:5])