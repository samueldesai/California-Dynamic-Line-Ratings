import Pkg

#Pkg.add("DataFrames")

#Pkg.add("Statistics")

#Pkg.add("CSV")

using DataFrames
using Statistics
using CSV

path = "C:\\Users\\gogre\\Downloads\\Load_Agg_Post_Assignment_v3_latest.csv"

# Read CSV data (no headers)
raw_data = CSV.File(path, header=false) |> DataFrame
rows, cols = size(raw_data)

# Initialize matrices for real and imaginary parts
real_matrix = Array{Float64}(undef, rows, cols)
imag_matrix = Array{Float64}(undef, rows, cols)

for i in 1:rows
    for j in 1:cols
        val = raw_data[i, j]  # Fixed: Use DataFrame indexing properly

        # Convert to string in case it's not already
        val_str = string(val)
        
        plus_idx = findfirst(==('+'), val_str)
        minus_idx = findfirst(==('-'), val_str)
        
        # Decide split index
        split_idx = nothing
        if plus_idx !== nothing && plus_idx > 1  # Ensure + is not the first character
            split_idx = plus_idx
        elseif minus_idx !== nothing && minus_idx > 1  # Ensure - is not the first character
            split_idx = minus_idx
        else
            error("Invalid format: $val_str")
        end

        # Extract real and imaginary parts as substrings
        real_str = val_str[1:split_idx-1]
        imag_str = val_str[split_idx:end-1]  # remove the trailing 'i'

        # Parse to float
        real_val = parse(Float64, real_str)
        imag_val = parse(Float64, imag_str)

        # Assign to matrices
        real_matrix[i, j] = real_val
        imag_matrix[i, j] = imag_val
    end
end

function timeavg_day(data)
    n_rows, n_cols = size(data)
    
    if n_cols % 24 != 0
        error("Each row must have a number of columns divisible by 24. Found $n_cols columns in $(n_rows)×$(n_cols) array.")
    end

    chunk_count = n_cols ÷ 24
    output = zeros(Float64, n_rows, 24)

    # Process each row
    for i in 1:n_rows
        row = data[i, :]
        for j in 1:24
            indices = j:24:(j + 24*(chunk_count-1))
            output[i, j] = mean(row[indices])
        end
    end
    return output
end

function timeavg_monthly(data)
    n_rows, n_cols = size(data)
    
    if n_cols != 8760
        error("Expected 8760 hourly values per row. Got $n_cols.")
    end

    # Days per month (non-leap year)
    days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    hours_per_month = days_in_month .* 24  # array of hours in each month

    output = zeros(Float64, n_rows, 12 * 24)  # 12 months × 24 hourly values

    for i in 1:n_rows
        row = data[i, :]
        idx = 1  # Index into original row
        for m in 0:11  # m = 0 to 11
            month_hours = hours_per_month[m+1]
            month_data = row[idx : idx + month_hours - 1]
            for h in 0:23
                # Collect all h-th hours across all days in this month
                hour_values = month_data[h+1:24:end]
                output[i, m*24 + h + 1] = mean(hour_values)
            end
            idx += month_hours
        end
    end

    return output
end

# Create header names: "january, 0", ..., "december, 23"
months = ["january", "february", "march", "april", "may", "june",
          "july", "august", "september", "october", "november", "december"]



col_names = String[]

for month in months
    for hour in 0:23
        push!(col_names, "$(month), $(hour)")
    end
end

# Process and save results
real_output = timeavg_monthly(real_matrix)
reactive_output = timeavg_monthly(imag_matrix)

real_df = DataFrame(real_output, Symbol.(col_names))
imag_df = DataFrame(reactive_output, Symbol.(col_names))

CSV.write("C:\\Users\\gogre\\Desktop\\real_output_288.csv", real_df)
CSV.write("C:\\Users\\gogre\\Desktop\\reactive_output_288.csv", imag_df)
