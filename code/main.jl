using JuMP
using TimerOutputs
using Printf, DelimitedFiles, Distances

const solver = "Gurobi"
const THREADLIM = 12
if solver == "Gurobi" && !@isdefined(GUROBI_ENV)
    using Gurobi
    const GUROBI_ENV = Gurobi.Env()
end
if solver == "Mosek"
    using Mosek, MosekTools
end

include("problem.jl")
include("facility.jl")
include("network.jl")

function test_facility_location(instance::String)
    budget = 2
    time_limit = 120.0
    timer = TimerOutput()
    reset_timer!(timer)

    problem = FacilityLocation(instance, budget)
    # run_default(problem, KKT, CCG, time_limit, timer)
    # run_default(problem, IndicatorDual, CCG, time_limit, timer)
    # run_default(problem, IndicatorDual, BD, time_limit, timer)
    run_adaptive_penalty_outer(problem, CCG, time_limit, timer)
    # run_adaptive_penalty_outer(problem, BD, time_limit, timer)
    @show(timer)

    return nothing
end

function test_network_design(instance::String)
    budget = 1
    time_limit = 120.0
    timer = TimerOutput()
    reset_timer!(timer)

    # det = solve_deterministic_problem(NetworkDesign("data/SNDLIB/dfn-gwin.txt", 1.0, 0)); @show(det);

    problem = NetworkDesign(instance, 1.0, budget)
    run_default(problem, IndicatorDual, CCG, time_limit, timer)
    run_adaptive_penalty_outer(problem, CCG, time_limit, timer)
    @show(timer)

    return nothing
end

# test_facility_location("data/CFLP/Cap_F20_C25.txt")
# test_network_design("data/SNDLIB/dfn-bwin.txt")
# test_network_design("data/SNDLIB/dfn-gwin.txt")
