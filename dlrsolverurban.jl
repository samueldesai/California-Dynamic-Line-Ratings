########################################################################
# 0.  Packages
########################################################################
using PowerModels              # network parsing + matrices
using JuMP, Gurobi, SparseArrays
using LinearAlgebra: Diagonal
using CSV, DataFrames
using Plots
using Dates
using Statistics

########################################################################
# 1.  Load case & build “basic” matrices
########################################################################
casefile  = "CaliforniaTestSystem.m"                   
raw_data  = parse_file(casefile)
data      = make_basic_network(raw_data)    # guarantees matrix form :contentReference[oaicite:0]{index=0}

# Matrices & parameters from the documentation ------------------------

A_t = calc_basic_incidence_matrix(data)                 # 𝑨ᵀ (branch × bus)
A = A_t'                                                # 𝑨 (bus × branch)
bz  = calc_basic_branch_series_impedance(data)          # complex z_ℓ in pu
bvec = -imag.(inv.(bz))                                 # susceptance b_ℓ in pu
Yd  = Diagonal(bvec)                                    # diag(b)
nb, nl = size(A_t,2), size(A_t,1)
baseMVA = data["baseMVA"]                            # base MVA

# Getting Information for model --------------------------------------

# Find the column indices where the first row of A is 1
start = now()
time_steps = 24
months_name = ["Jan", "Apr", "Jul", "Oct"]

folder = "DLRoutput/" #TODO CHANGE
mkpath(folder)

# finding slack bus index
slack_bus_idx = 1951 # manual assignment for debugging
for bus in values(data["bus"])
    bus_type = bus["bus_type"]
    if bus_type == 3 # slack bus
        global slack_bus_idx = bus["bus_i"]
        break
    end
end

Pd_jan_path  = "finaldata/january_load.csv"
Pd_jan_df = CSV.read(Pd_jan_path, DataFrame) 
Pd_apr_path  = "finaldata/april_load.csv"
Pd_apr_df = CSV.read(Pd_apr_path, DataFrame) 
Pd_july_path  = "finaldata/july_load.csv"
Pd_july_df = CSV.read(Pd_july_path, DataFrame) 
Pd_oct_path  = "finaldata/october_load.csv"
Pd_oct_df = CSV.read(Pd_oct_path, DataFrame) 

Pd = [Pd_jan_df, Pd_apr_df, Pd_july_df, Pd_oct_df]

Pg_jan_path  = "finaldata/january_gen.csv"
Pg_jan_df = CSV.read(Pg_jan_path, DataFrame) 
Pg_apr_path  = "finaldata/april_gen.csv"
Pg_apr_df = CSV.read(Pg_apr_path, DataFrame) 
Pg_july_path  = "finaldata/july_gen.csv"
Pg_july_df = CSV.read(Pg_july_path, DataFrame) 
Pg_oct_path  = "finaldata/october_gen.csv"
Pg_oct_df = CSV.read(Pg_oct_path, DataFrame) 

Pg = [Pg_jan_df, Pg_apr_df, Pg_july_df, Pg_oct_df]

# Branch MW limits, careful that line ratings are in the same order as A

slr_line_limits_df = CSV.read("SLR_Ratings_6_2.csv", DataFrame) #TODO change 
CSV.write(folder * "slr_line_limits.csv", slr_line_limits_df)
slr_line_limit = slr_line_limits_df[!, :Rating]

dlr_limits_imported_df = CSV.read("DLR_Ratings_6_2.csv", DataFrame)
dlr_line_limit = Array{Float64}(undef, nl, 24, 4)
for hour in 0:23
    for month in 1:4
        column = months_name[month] * " " * string(hour)
        dlr_line_limit[:, hour+1, month] = dlr_limits_imported_df[:, column]
    end 
end

# Read solar_available CSV and create a dictionary
solar_df = CSV.read("finaldata/solar_available.csv", DataFrame)
solar_dict = Dict(row.Bus_number => row.Solar_Potential_MW for row in eachrow(solar_df))
# Create a vector of solar availability per bus (default 0 if not present)
solar_avail_vec = [get(solar_dict, bus, 0.0) for bus in 1:nb] ./ baseMVA
CSV.write(folder * "solar_pot.csv", DataFrame(SolarAvailability = solar_avail_vec))

# Read windpotential CSV and create a dictionary
wind_df = CSV.read("finaldata/windpotential.csv", DataFrame)
wind_dict = Dict(row.Bus_number => row.wind_potential for row in eachrow(wind_df))
# Create a vector of wind availability per bus (default 0 if not present)
wind_avail_vec = [get(wind_dict, bus, 0.0) for bus in 1:nb] ./ baseMVA
CSV.write(folder * "wind_pot.csv", DataFrame(WindAvailability = wind_avail_vec))

jan_cf_df = CSV.read("finaldata/renewable_cf/jan_renewables.csv", DataFrame)
apr_cf_df = CSV.read("finaldata/renewable_cf/apr_renewables.csv", DataFrame)
july_cf_df = CSV.read("finaldata/renewable_cf/july_renewables.csv", DataFrame)
oct_cf_df = CSV.read("finaldata/renewable_cf/oct_renewables.csv", DataFrame)

renewable_cfs = [jan_cf_df, apr_cf_df, july_cf_df, oct_cf_df]

load_inc_df = CSV.read("finaldata/exportedbuses.csv", DataFrame, select=["bus_i"])

# Get allowed buses from load_inc_df

allowed_buses = Set(load_inc_df.bus_i)


println("done reading data, now building model ... ")

# ...existing package and data loading code...

months = 1:4
time_steps = 24

# --- Storage for combined results ---
all_solar_flows = Array{Float64}(undef, nl, time_steps * 4)
all_solar_gen = Array{Float64}(undef, nb, time_steps * 4)
all_solar_extra_demand = Array{Float64}(undef, nb, time_steps * 4)
all_solar_load = Array{Float64}(undef, nb, time_steps * 4)
all_solar_thermal_gen = Array{Float64}(undef, nb, time_steps * 4)
all_solar_capacity = zeros(nb, 4)

all_wind_flows = Array{Float64}(undef, nl, time_steps * 4)
all_wind_gen = Array{Float64}(undef, nb, time_steps * 4)
all_wind_extra_demand = Array{Float64}(undef, nb, time_steps * 4)
all_wind_load = Array{Float64}(undef, nb, time_steps * 4)
all_wind_thermal_gen = Array{Float64}(undef, nb, time_steps * 4)
all_wind_capacity = zeros(nb, 4)

for (mi, month) in enumerate(months)
    # --- SOLAR ---
    m_solar = Model(Gurobi.Optimizer)
    set_optimizer_attribute(m_solar, "LogToConsole", 0)
    set_optimizer_attribute(m_solar, "log_file", folder * "solar_solve_$mi.log")
    println("building solar model for month $month")

    @variable(m_solar, SolarCapacity[1:nb] >= 0)
    @variable(m_solar, θ[1:nb, 1:time_steps])
    @variable(m_solar, F[1:nl, 1:time_steps])
    @variable(m_solar, SolarGen[1:nb, 1:time_steps])
    @variable(m_solar, ExtraDemandSolar[1:nb, 1:time_steps] >= 0)

    Pg_matrix = Matrix(Pg[month])
    Pd_matrix = Matrix(Pd[month])
    renewable_cf_df = renewable_cfs[month]

    for i in 1:time_steps
        @constraint(m_solar, θ[slack_bus_idx, i] == 0)
        @constraint(m_solar, F[:, i] .== Diagonal(bvec) * A_t * θ[:, i])
        @constraint(m_solar, A * F[:, i] .== Pg_matrix[:, i] + SolarGen[:, i] - Pd_matrix[:, i] - ExtraDemandSolar[:, i])
    end

    for line in 1:nl
        for i in 1:time_steps
            if (slr_line_limit[line] > dlr_line_limit[line, i, month])
                @constraint(m_solar, F[line, i] <= slr_line_limit[line])
                @constraint(m_solar, -slr_line_limit[line] <= F[line, i])
            else
                @constraint(m_solar, F[line, i] <= dlr_line_limit[line, i, month])
                @constraint(m_solar, -dlr_line_limit[line, i, month] <= F[line, i])                
            end
        end 
    end
    
    @constraint(m_solar, [bus=1:nb, i=1:time_steps], SolarGen[bus, i] <= SolarCapacity[bus] * renewable_cf_df[i, :avg_solar_cf])
    @constraint(m_solar, [bus=1:nb, i=1:time_steps], SolarGen[bus, i] >= 0)
    @constraint(m_solar, [bus=1:nb], SolarCapacity[bus] <= solar_avail_vec[bus])
    for bus in 1:nb
        if solar_avail_vec[bus] > 0
            for t in 1:time_steps
                @constraint(m_solar, ExtraDemandSolar[bus, t] == 0)
            end
        end
        if bus ∉ allowed_buses
            for t in 1:time_steps
                @constraint(m_solar, ExtraDemandSolar[bus, t] == 0)
            end
        end
    end

    @objective(m_solar, Max, sum(SolarGen) - sum(SolarCapacity) * 0.1)
    println("optimizing solar model for month $mi")
    optimize!(m_solar)

    idx_start = (mi-1)*time_steps + 1
    idx_end = mi*time_steps

    all_solar_flows[:, idx_start:idx_end] = value.(F)
    all_solar_gen[:, idx_start:idx_end] = value.(SolarGen)
    all_solar_extra_demand[:, idx_start:idx_end] = value.(ExtraDemandSolar)
    all_solar_load[:, idx_start:idx_end] = Pd_matrix
    all_solar_thermal_gen[:, idx_start:idx_end] = Pg_matrix
    all_solar_capacity[:, mi] = value.(SolarCapacity)

    # --- WIND ---
    m_wind = Model(Gurobi.Optimizer)
    set_optimizer_attribute(m_wind, "LogToConsole", 0)
    set_optimizer_attribute(m_wind, "log_file", folder * "wind_solve_$mi.log")
    println("building wind model for month $month")

    @variable(m_wind, WindCapacity[1:nb] >= 0)
    @variable(m_wind, θ[1:nb, 1:time_steps])
    @variable(m_wind, F[1:nl, 1:time_steps])
    @variable(m_wind, WindGen[1:nb, 1:time_steps])
    @variable(m_wind, ExtraDemandWind[1:nb, 1:time_steps] >= 0)

    for i in 1:time_steps
        @constraint(m_wind, θ[slack_bus_idx, i] == 0)
        @constraint(m_wind, F[:, i] .== Diagonal(bvec) * A_t * θ[:, i])
        @constraint(m_wind, A * F[:, i] .== Pg_matrix[:, i] + WindGen[:, i] - Pd_matrix[:, i] - ExtraDemandWind[:, i])
    end

    for line in 1:nl
        for i in 1:time_steps
            if (slr_line_limit[line] > dlr_line_limit[line, i, month])
                @constraint(m_wind, F[line, i] <= slr_line_limit[line])
                @constraint(m_wind, -slr_line_limit[line] <= F[line, i])
            else
                @constraint(m_wind, F[line, i] <= dlr_line_limit[line, i, month])
                @constraint(m_wind, -dlr_line_limit[line, i, month] <= F[line, i])                
            end
        end 
    end

    @constraint(m_wind, [bus=1:nb, i=1:time_steps], WindGen[bus, i] <= WindCapacity[bus] * renewable_cf_df[i, :avg_wind_cf])
    @constraint(m_wind, [bus=1:nb, i=1:time_steps], WindGen[bus, i] >= 0)
    @constraint(m_wind, [bus=1:nb], WindCapacity[bus] <= wind_avail_vec[bus])
    for bus in 1:nb
        if wind_avail_vec[bus] > 0
            for t in 1:time_steps
                @constraint(m_wind, ExtraDemandWind[bus, t] == 0)
            end
        end
        if bus ∉ allowed_buses
            for t in 1:time_steps
                @constraint(m_wind, ExtraDemandWind[bus, t] == 0)
            end
        end
    end

    @objective(m_wind, Max, sum(WindGen) - sum(WindCapacity) * 0.1)
    println("optimizing wind model for month $mi")
    optimize!(m_wind)
    
    all_wind_flows[:, idx_start:idx_end] = value.(F)
    all_wind_gen[:, idx_start:idx_end] = value.(WindGen)
    all_wind_extra_demand[:, idx_start:idx_end] = value.(ExtraDemandWind)
    all_wind_load[:, idx_start:idx_end] = Pd_matrix
    all_wind_thermal_gen[:, idx_start:idx_end] = Pg_matrix
    all_wind_capacity[:, mi] = value.(WindCapacity)
end

# --- Save combined results ---
CSV.write(folder * "solar_combined_line_flows.csv", DataFrame(all_solar_flows, :auto))
CSV.write(folder * "solar_combined_solar_gen.csv", DataFrame(all_solar_gen, :auto))
CSV.write(folder * "solar_combined_extra_demand.csv", DataFrame(all_solar_extra_demand, :auto))
CSV.write(folder * "solar_combined_load.csv", DataFrame(all_solar_load, :auto))
CSV.write(folder * "solar_combined_thermal_gen.csv", DataFrame(all_solar_thermal_gen, :auto))
CSV.write(folder * "solar_capacity_by_month.csv", DataFrame(all_solar_capacity, :auto))

CSV.write(folder * "wind_combined_line_flows.csv", DataFrame(all_wind_flows, :auto))
CSV.write(folder * "wind_combined_wind_gen.csv", DataFrame(all_wind_gen, :auto))
CSV.write(folder * "wind_combined_extra_demand.csv", DataFrame(all_wind_extra_demand, :auto))
CSV.write(folder * "wind_combined_load.csv", DataFrame(all_wind_load, :auto))
CSV.write(folder * "wind_combined_thermal_gen.csv", DataFrame(all_wind_thermal_gen, :auto))
CSV.write(folder * "wind_capacity_by_month.csv", DataFrame(all_wind_capacity, :auto))

# --- Combined plots (example for solar) ---
time = 1:(time_steps*4)
total_solar = sum(all_solar_gen, dims=1)[:]
total_thermal = sum(all_solar_thermal_gen, dims=1)[:]
total_load = sum(all_solar_load, dims=1)[:]
total_new_load = sum(all_solar_extra_demand, dims=1)[:] + total_load

plot(time, total_thermal, label="Thermal Gen", lw=2)
plot!(time, total_solar, label="Solar Gen", lw=2)
plot!(time, total_load, label="Original Load", lw=2)
plot!(time, total_new_load, label="Total Load (with Extra Demand)", lw=2)
xlabel!("Hour")
ylabel!("MW (per unit $baseMVA)")
title!("Solar: Total Generation and Load Over 4 Months (96 hours)")
savefig(folder * "solar_combined_generation_load_timeseries.png")

# --- Combined wind plots ---
time = 1:(time_steps*4)
total_wind = sum(all_wind_gen, dims=1)[:]
total_thermal = sum(all_wind_thermal_gen, dims=1)[:]
total_load = sum(all_wind_load, dims=1)[:]
total_new_load = sum(all_wind_extra_demand, dims=1)[:] + total_load

plot(time, total_thermal, label="Thermal Gen", lw=2)
plot!(time, total_wind, label="Wind Gen", lw=2)
plot!(time, total_load, label="Original Load", lw=2)
plot!(time, total_new_load, label="Total Load (with Extra Demand)", lw=2)
xlabel!("Hour")
ylabel!("MW (per unit $baseMVA)")
title!("Wind: Total Generation and Load Over 4 Months (96 hours)")
savefig(folder * "wind_combined_generation_load_timeseries.png")

bar_vals = [sum(total_thermal), sum(total_wind), sum(total_load), sum(total_new_load)]
bar_labels = ["Thermal Gen", "Wind Gen", "Original Load", "Total Load"]
bar_plot = bar(bar_labels, bar_vals, legend=false, ylabel="MW·h", title="Wind: Total Generation over 4 Months")
savefig(bar_plot, folder * "wind_combined_total_generation_bar.png")