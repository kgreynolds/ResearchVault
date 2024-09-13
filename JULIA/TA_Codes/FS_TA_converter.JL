using CSV
using DataFrames

# Read the CSV file into a DataFrame
file_to_convert = "/Users/kris/Desktop/ResearchVault/RESEARCH/NOCERA GROUP.nosync/PROJECTS/DisulfideRadical Anion/Experiemnts/2023_07_26_FS_TA_PhSSPh_KHCO2_5mM_DMSO/3. PhSSPh_20mM_formate_5mM_DMSO_360pump_3ns_flow_averaged"
df = DataFrame(CSV.File(file_to_convert, header=false))

    # Transpose the DataFrame
    df = permutedims(df)
    N = nrow(df)

    # Add a new column to the left with integers from 1 to N
    insertcols!(df, 1, :new_column => 1:1:N)

    # Sort the new column in descending order
    sort!(df, :new_column, rev=true)

        # Read the calibration CSV file into a DataFrame
        calibration_df = DataFrame(CSV.File("JULIA/TA_Codes/FS_WhitelightCalibration/Femtosecond_Calibration_Output_2023_07_09.csv", header=true))
        sort!(calibration_df, :Pixel,rev = true)

    # Replace the first column of df with the second column of calibration_df
    df[!, 1] = calibration_df[!, 2]
df


## Select Time file for conversion
timefile = DataFrame(CSV.File("JULIA/TA_Codes/TimeFiles/2023_05_12_TA_timefile_3ns.txt",header = false))

    # Convert from fs to ps
    timefile = timefile ./1000

    # Add a zero to the top of the first column of timefile
    insert!(timefile, 1, Dict(:Column1 => 0))

    #Transpose timefile
    timefile = permutedims(timefile)

    # Insert the transposed timefile rows into the df
    for i in ncol(timefile)
        insert!(df, 1, [ timefile[1,i] for i in 1:ncol(timefile)])
    end
df

# Write the resulting DataFrame to a new CSV file
CSV.write("JULIA/Outputs/FS_TA_Converter_output.csv", df, header = false)
