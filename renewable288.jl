using CSV
using DataFrames
using Statistics
using Plots

# Load the data
df = CSV.read("C:\\Users\\gogre\\Downloads\\Expanded_Profiles.csv", DataFrame, header = 3)

# Rename for easier access
rename!(df, Dict(
    :Month => :month,
    :Hour => :hour_index,
    Symbol("Solar CF") => :solar_cf,
    Symbol("Wind CF") => :wind_cf
))

# Convert hour index to hour of day (1 to 24)
df.hour_of_day = mod.(df.hour_index .- 1, 24) .+ 1

# Group by month and hour of day, then compute means
grouped = combine(groupby(df, [:month, :hour_of_day]), 
    :solar_cf => mean => :avg_solar_cf,
    :wind_cf => mean => :avg_wind_cf
)

# Sort by month and hour of day
sort!(grouped, [:month, :hour_of_day])

# Optional: Save to CSV
CSV.write("C:\\Users\\gogre\\Desktop\\Expanded_Profiles.csv", grouped)

# Show result
println(grouped)

# Reshape into 24x12 matrices (24 hours, 12 months)
solar_matrix = reshape(grouped.avg_solar_cf, 24, 12)
wind_matrix = reshape(grouped.avg_wind_cf, 24, 12)

# Set up ticks
xticks = 1:1:12
yticks = 1:1:24

# Solar CF Heatmap
solar_plot = heatmap(xticks, yticks, solar_matrix;
    xlabel="Month", ylabel="Hour of Day", title="Average Solar Capacity Factor",
    colorbar_title="Solar CF", yflip=false, xticks=xticks, yticks=yticks)

# Save solar plot
savefig(solar_plot, "C:\\Users\\gogre\\Desktop\\avg_solar_cf.png")

# Wind CF Heatmap
wind_plot = heatmap(xticks, yticks, wind_matrix;
    xlabel="Month", ylabel="Hour of Day", title="Average Wind Capacity Factor",
    colorbar_title="Wind CF", yflip=false, xticks=xticks, yticks=yticks)

# Save wind plot
savefig(wind_plot, "C:\\Users\\gogre\\Desktop\\avg_wind_cf.png")