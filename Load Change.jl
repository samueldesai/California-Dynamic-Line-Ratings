using PowerModels, JuMP, KNITRO, SparseArrays
using LinearAlgebra: Diagonal
using CSV, DataFrames
using Dates
using Plots

# Load case file and build network matrices
casefile = "CaliforniaTestSystem.m"
raw_data = parse_file(casefile)
data = make_basic_network(raw_data)

A_t = calc_basic_incidence_matrix(data)
A = A_t'
bz = calc_basic_branch_series_impedance(data)
bvec = -imag.(inv.(bz))
nb, nl = size(A_t, 2), size(A_t, 1)
baseMVA = data["baseMVA"]
Yd = spdiagm(bvec)

# Output folder
folder = "finaldata"
mkpath(folder)

# Load Pd and generator data
Pd_df = CSV.read(folder * "\\real_output_288.csv", DataFrame)
Pd_df = Matrix{Float64}(Pd_df) ./ baseMVA

line_limits_df = CSV.read("test.csv", DataFrame)
line_limit = line_limits_df[!, :Rating]

# Slack bus index
slack_bus_idx = 1951
for bus in values(data["bus"])
    if bus["bus_type"] == 3
        slack_bus_idx = bus["bus_i"]
        break
    end
end

# Generator data
gen_df = DataFrame(Pg_min=zeros(nb), Pg_max=zeros(nb), c2=zeros(nb), c1=zeros(nb), c0=zeros(nb))
for gen in values(data["gen"])
    bus = gen["gen_bus"]
    gen_df[bus, :Pg_min] += gen["pmin"]
    gen_df[bus, :Pg_max] += gen["pmax"]
    gen_df[bus, :c2] += gen["cost"][3]
    gen_df[bus, :c1] += gen["cost"][2]
    gen_df[bus, :c0] += gen["cost"][1]
end

# Run each 24-hour block
month_names = ["january", "april", "july", "october"]
time_steps = 24

for block in 0:3 
    month = month_names[block + 1]
    println(">>> Solving for $month")

    col_start = block * time_steps + 1
    col_end = col_start + time_steps - 1
    Pd_block = Pd_df[:, col_start:col_end]

    m = Model(KNITRO.Optimizer)
    set_optimizer_attribute(m, "honorbnds", 1)
    set_optimizer_attribute(m, "algorithm", 2)
    
    # Check for warm-start files first
    gen_file = folder * "\\$(month)_gen.csv"
    load_file = folder * "\\$(month)_load.csv"
    flow_file = folder * "\\$(month)_flow.csv"
    adj_file  = folder * "\\$(month)_adj.csv"
    
    has_warmstart = isfile(gen_file) && isfile(load_file) && isfile(flow_file) && isfile(adj_file)
    
    if has_warmstart
        println("[$month] Using warm-start from previous solution...")
        # Use algorithm 4 (interior-point) for warm-start
        set_optimizer_attribute(m, "algorithm", 1)
    else
        println("[$month] Cold start - no previous solution found")
        # Use algorithm 1 (interior-point with crossover) for cold start
        set_optimizer_attribute(m, "algorithm", 1)
    end

    @variable(m, θ[1:nb, 1:time_steps])
    @variable(m, Pg[1:nb, 1:time_steps])
    @variable(m, F[1:nl, 1:time_steps])
    @variable(m, Pd_adjustment[1:nb, 1:time_steps] >= 0)
    @variable(m, Pd[1:nb, 1:time_steps])

    @constraint(m, [b=1:nb, t=1:time_steps], Pd[b, t] == Pd_block[b, t] * Pd_adjustment[b, t])
    @constraint(m, [b=1:nb, t=1:time_steps], Pd_adjustment[b, t] <= 1)

    for t in 1:time_steps
        @constraint(m, θ[slack_bus_idx, t] == 0)
        @constraint(m, A * F[:, t] .== Pg[:, t] - Pd[:, t])
        @constraint(m, F[:, t] .== Yd * A_t * θ[:, t])
    end

    @constraint(m, [b=1:nb, t=1:time_steps], gen_df.Pg_min[b] <= Pg[b, t] <= gen_df.Pg_max[b])
    @constraint(m, [l=1:nl, t=1:time_steps], -line_limit[l] <= F[l, t] <= line_limit[l])

    @expression(m, power_cost[b=1:nb, t=1:time_steps],
        gen_df.c2[b] * Pg[b, t]^2 + gen_df.c1[b] * Pg[b, t] + gen_df.c0[b])
    @expression(m, penalty[b=1:nb, t=1:time_steps], 1e20 * (1 - Pd_adjustment[b, t]))

    @objective(m, Min, sum(power_cost) + sum(penalty))

    # Apply warm-start values AFTER model is built
    if has_warmstart
        try
            Pg_init = Matrix(CSV.read(gen_file, DataFrame))
            Pd_init = Matrix(CSV.read(load_file, DataFrame))
            F_init  = Matrix(CSV.read(flow_file, DataFrame))
            adj_init = Matrix(CSV.read(adj_file, DataFrame))

            # Set starting values for all variables
            for b in 1:nb, t in 1:time_steps
                if size(Pg_init, 1) >= b && size(Pg_init, 2) >= t
                    set_start_value(Pg[b, t], Pg_init[b, t])
                end
                if size(Pd_init, 1) >= b && size(Pd_init, 2) >= t
                    set_start_value(Pd[b, t], Pd_init[b, t])
                end
                if size(adj_init, 1) >= b && size(adj_init, 2) >= t
                    set_start_value(Pd_adjustment[b, t], adj_init[b, t])
                end
            end
            
            for l in 1:nl, t in 1:time_steps
                if size(F_init, 1) >= l && size(F_init, 2) >= t
                    set_start_value(F[l, t], F_init[l, t])
                end
            end
            
            # Calculate θ values from the flow values using F = Yd * A_t * θ
            # θ = (A_t' * Yd^-1) * F (approximately, for starting values)
            for t in 1:time_steps
                try
                    # Simple approach: set θ to zero except at non-slack buses
                    for b in 1:nb
                        if b != slack_bus_idx
                            set_start_value(θ[b, t], 0.0)
                        end
                    end
                catch e
                    println("Warning: Could not set θ starting values: $e")
                end
            end
            
        catch e
            println("Warning: Error loading warm-start files: $e")
            println("Proceeding with cold start...")
        end
    end

    optimize!(m)
    
    # Check solution status
    if termination_status(m) == MOI.OPTIMAL
        println("[$month] Optimal solution found. Obj = ", objective_value(m))
    elseif termination_status(m) == MOI.LOCALLY_SOLVED
        println("[$month] Locally optimal solution found. Obj = ", objective_value(m))
    else
        println("[$month] Solution status: ", termination_status(m))
        println("[$month] Obj = ", objective_value(m))
    end

    # Save results
    mkpath(folder * "\\6_2_load")
    CSV.write(folder * "\\6_2_load\\$(month)_gen.csv", DataFrame(value.(Pg), :auto))
    CSV.write(folder * "\\6_2_load\\$(month)_load.csv", DataFrame(value.(Pd), :auto))
    CSV.write(folder * "\\6_2_load\\$(month)_flow.csv", DataFrame(value.(F), :auto))
    CSV.write(folder * "\\6_2_load\\$(month)_adj.csv", DataFrame(value.(Pd_adjustment), :auto))

    # Plot generation and load over time
    gen = value.(Pg)
    load = value.(Pd)
    orig_load = Pd_block

    total_gen = sum(gen, dims=1)'
    total_load = sum(orig_load, dims=1)'
    total_new_load = sum(load, dims=1)'

    time = 1:time_steps
    plot(time, total_gen[1:time_steps], label="Total Generation", lw=2)
    plot!(time, total_load[1:time_steps], label="Total Load", lw=2)
    plot!(time, total_new_load[1:time_steps], label="Total New Load", lw=2)
    xlabel!("Hour")
    ylabel!("MW (per unit $baseMVA)")
    title!("Total Generation and Load Over 24 Hours ($month)")
    savefig(folder * "\\finaldata\\gen_load_timeseries_$(month).png")

    println("[$month] Done.\n")
end