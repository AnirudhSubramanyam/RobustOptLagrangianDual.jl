using JuMP
using TimerOutputs
using Printf, DelimitedFiles, Distances

const SOLVER = "Gurobi"
const THREADLIM = 12
if SOLVER == "Gurobi" && !@isdefined(GUROBI_ENV)
    using Gurobi
    const GUROBI_ENV = Gurobi.Env()
end
if SOLVER == "Mosek"
    using Mosek, MosekTools
end
if SOLVER == "CPLEX"
    using CPLEX
end

include("problem.jl")
include("facility.jl")
include("network.jl")
include("rostering.jl")
include("algorithms.jl")

function test_facility_location(instance::String, budget_range::Vector{Int})
    time_limit = 3600.0
    info_d = Tuple{Int, Float64, Float64, Float64}[]
    info_p = Tuple{Int, Float64, Float64, Float64}[]
    
    reset_timer!()
    for (i, budget) in enumerate(budget_range)
        problem = FacilityLocation(instance, budget)
        it_d, lb_d, ub_d, ttime_d = run_ccg(problem, LinearizedKKT, time_limit)
        it_p, lb_p, ub_p, ttime_p = run_ccg(problem, PenaltyDual, time_limit)
        # it_d, lb_d, ub_d, ttime_d = run_benders(problem, LinearizedDual, time_limit)
        # it_p, lb_p, ub_p, ttime_p = run_benders(problem, PenaltyDual, time_limit)
        push!(info_d, (it_d, lb_d, ub_d, ttime_d))
        push!(info_p, (it_p, lb_p, ub_p, ttime_p))
    end

    for (i, budget) in enumerate(budget_range)
        @printf("budget = %d. OPT-default = [%.2f, %.2f] (time = %.2f, iter = %d)\n", budget, info_d[i][2], info_d[i][3], info_d[i][4], info_d[i][1])
    end
    for (i, budget) in enumerate(budget_range)
        @printf("budget = %d. OPT-penalty = [%.2f, %.2f] (time = %.2f, iter = %d)\n", budget, info_p[i][2], info_p[i][3], info_p[i][4], info_p[i][1])
    end
    print_timer()

    return nothing
end

function test_network_design(instance::String, budget_range::Vector{Int})
    time_limit = 600.0
    info_d = Tuple{Int, Float64, Float64, Float64}[]
    info_p = Tuple{Int, Float64, Float64, Float64}[]
    det = solve_deterministic_problem(NetworkDesign(instance, 1.0, 0))
    
    reset_timer!()
    for (i, budget) in enumerate(budget_range)
        problem = NetworkDesign(instance, 1.0, budget)
        it_d, lb_d, ub_d, ttime_d = run_ccg(problem, IndicatorDual, time_limit)
        it_p, lb_p, ub_p, ttime_p = run_ccg(problem, PenaltyDual, time_limit)
        push!(info_d, (it_d, lb_d, ub_d, ttime_d))
        push!(info_p, (it_p, lb_p, ub_p, ttime_p))
    end

    @printf("OPT-deterministic = %.2f\n", det)
    for (i, budget) in enumerate(budget_range)
        @printf("budget = %d. OPT-default = [%.2f, %.2f] (time = %.2f, iter = %d)\n", budget, info_d[i][2], info_d[i][3], info_d[i][4], info_d[i][1])
    end
    for (i, budget) in enumerate(budget_range)
        @printf("budget = %d. OPT-penalty = [%.2f, %.2f] (time = %.2f, iter = %d)\n", budget, info_p[i][2], info_p[i][3], info_p[i][4], info_p[i][1])
    end
    print_timer()

    return nothing
end

function test_rostering(scale::Int, budget_range::Vector{Int})
    scale_time = scale_demand = scale
    time_limit = 100.0
    reset_timer!()
    for (i, budget) in enumerate(budget_range) #3:3:18
        problem = Rostering(budget, scale_time, scale_demand)
        println("starting LinearizedKKT-$budget")
        run_ccg(problem, LinearizedKKT, time_limit)
        println("starting PenaltyDual-$budget")
        run_ccg(problem, PenaltyDual, time_limit)
    end
end

# test_facility_location("data/CFLP/Cap_F20_C25.txt", collect(1:2))
# test_network_design("data/SNDLIB/dfn-bwin.txt", collect(1:5))
# test_rostering(1, collect(12:12))

