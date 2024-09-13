using Interpolations,Plots, DataFrames, CSV, Tables, Statistics, Optim
### Moving Average Calculator, Used for nanosecond TA Kinetics
## Update path to file and appropriate Header info (NOTE: make sure the right rows and collums are dropped depending on input file structure)

plotly()

function movingavTA(x::Vector{Float64},y::Vector{Float64},bin::Int64) ## moving average and calculate absorbance for kinetic traces obtain from NS oscilloscope
    outx::Vector{Float64} = []
    outy::Vector{Float64} = []

    n = floor(Int, bin*10^(-9)/(x[3]-x[2]))

    for i in 1 : length(y) - n 
        localx = (mean(x[i:i+n]))*1000000 # take the average in the interval of n to the right of i and convert from seconds to μ seconds
        localy = (mean(y[i:i+n])) # take the average in the interval of n to the right of i
        push!(outx, localx)
        push!(outy,localy)
        
    end

    outy_t0 = mean(outy[1:500])
    Absorbance = real(log10.(complex(((outy_t0)./outy))))

    df = DataFrame(Time_μs = outx, Abs = Absorbance)

    return outx, Absorbance, display(Plots.plot(outx, Absorbance, size=(1000, 800),xlabel = "Time(μs)", ylabel = "ΔA")), CSV.write("JULIA/Outputs/movingAvgTA_out.csv", df)
end



File_to_Convert = CSV.File(read("/Users/kris/Downloads/001_Trp-nh2_TA_Kinetics_450nm_500-sweeps.trc"), skipto = 6,drop = [3], header = ["Time", "Ampl"])
x = File_to_Convert["Time"]
y = File_to_Convert["Ampl"]

x[500] # Check this for time zero in function


movingavTA(x,y,0)
movingavTA(x,y,10)




### Bullshit Emission subtraction
IrK_PhssPh_450TAKin = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/RAW_Data/NS/2023_08_24/IrK_PhSSPh_450nm_TAKinetics_500avg_18us.trc", skipto = 6,drop = [3], header = ["Time", "Ampl"], DataFrame)
Plots.plot(IrK_PhssPh_450TAKin[:,1],IrK_PhssPh_450TAKin[:,2])

IrK_PhSSPh_480nm_EmissionKinetics = CSV.read("/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_08_24_NS_TA_PhSSPh_IrK_MeCN/Kinetics/IrK_PhSSPh_480nm_EmissionKinetics_25avg_900ns.csv", header = true, DataFrame)
IrK_PhSSPh_480nm = DataFrame(time = (IrK_PhSSPh_480nm_EmissionKinetics[:,1] ./ 1000000), I = IrK_PhSSPh_480nm_EmissionKinetics[:,2])

Plots.plot!(IrK_PhSSPh_480nm[:,1],(IrK_PhSSPh_480nm[:,2] .* 0.265).- 0.0768)

x = collect(range(IrK_PhSSPh_480nm[1,1],IrK_PhSSPh_480nm[end,1],2000))

# x = collect(range(IrK_PhssPh_450TAKin[1,1],IrK_PhssPh_450TAKin[end,1],1000))

IrK_PhSSPh_480nm_Interp = LinearInterpolation(IrK_PhSSPh_480nm[:,1],IrK_PhSSPh_480nm[:,2])
IrK_PhSSPh_480nm_Interp_Y = IrK_PhSSPh_480nm_Interp(x)

IrK_PhssPh_450TAKin_Inerp = LinearInterpolation(IrK_PhssPh_450TAKin[:,1],IrK_PhssPh_450TAKin[:,2])
IrK_PhssPh_450TAKin_Inerp_Y = IrK_PhssPh_450TAKin_Inerp(x)

Plots.plot(x,IrK_PhSSPh_480nm_Interp_Y)
Plots.plot!(x,IrK_PhssPh_450TAKin_Inerp_Y)
Plots.plot!(x,fitfunction(values))
Plots.plot(x,IrK_PhssPh_450TAKin_Inerp_Y .- fitfunction(values))
    ## FITTING FUNCTION
    fitfunction(p) = p[1] .* IrK_PhSSPh_480nm_Interp_Y .+ p[2]
    p0 = [0.265, -0.0768]

    peakfit0(p) = sum((IrK_PhssPh_450TAKin_Inerp_Y .- fitfunction(p)).^2)

    res = optimize(peakfit0, p0,LBFGS(),Optim.Options(time_limit = 2000))
    print("Fit Resuluts:", res)
    values = Optim.minimizer(res) # Kinetics results
    residual =  KineticTrace .- fitfunction(values)
    print("Fit Coeff:", values)