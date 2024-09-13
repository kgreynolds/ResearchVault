## Femtosecond laser Probe Light Calibration
## Use any CSV file containing the calibration data points as Pixel vs. Wavelenght 

using Plots, DataFrames, CSV, CurveFit
plotly()
function Makefile(x,y)
    df = DataFrame(Pixel = x, Wavelength = y)

    return(CSV.write("JULIA/TA_Codes/FS_WhitelightCalibration/Femtosecond_Calibration_Output.csv",df))
end

## Path to the Calibration points to be fitted
Input = CSV.File(read("JULIA/TA_Codes/FS_WhitelightCalibration/Femtosecond_Calibration_Input_2023_07_09.csv"))

Pixel = Input["Pixel"]
Wavelength = Input["Wavelength"]

## Fit the calibration points to a nth degree polynomial
fit = curve_fit(Polynomial, Pixel, Wavelength, 3)
ran = collect(1:1:600)

Plots.plot(ran, fit.(ran), label = "Fit", size=(900, 600), xlabel = "pixel", ylabel = "Wavelength") 
Plots.scatter!(Pixel, Wavelength, label = "Raw", markersize = 5)


Makefile(ran, fit.(ran))