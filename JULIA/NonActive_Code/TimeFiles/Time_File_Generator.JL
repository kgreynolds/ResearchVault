using CSV, DataFrames

## Time file generator for Femtosecond TA LabView Program

# Maximum Time delay in picoseconds
Δt = 6000

# Number of time points in the linear interval
N = 20

# Total number of time points in the file
M = round(N*(1+log10(Δt)))

# Ranges for linear and exponential parts
i1 = 0:1:(N-1)
i2 = N:1:M

linear_part = round.((collect(-1 .+ (2*i1/(N)))),digits = 1)
exponential_part = round.(10 .^(-1 .+ i2/N),digits = 1)

timefile = vcat(linear_part,exponential_part) .*1000
timefile_df = DataFrame(timefile=timefile)


CSV.write("2023_05_12_TA_timefile_6ns.txt",timefile_df, header = false)


