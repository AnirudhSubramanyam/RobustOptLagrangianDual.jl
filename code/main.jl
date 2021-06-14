include("problem.jl")
include("facility.jl")

function test_facility_location(instance::String)
    budget = 10
    time_limit = 7200.0
    timer = TimerOutput()
    reset_timer!(timer)

    problem = FacilityLocation(instance, budget)
    run_default(problem, KKT, CCG, time_limit, timer)
    run_default(problem, LinearizedDual, CCG, time_limit, timer)
    run_default(problem, IndicatorDual, CCG, time_limit, timer)
    run_default(problem, IndicatorDual, BD, time_limit, timer)
    #=
    run_default(problem, IndicatorDual, BD, time_limit, timer)
    run_default(problem, Penalty, CCG, time_limit, timer)
    run_default(problem, Penalty, BD, time_limit, timer)
    run_adaptive_penalty(problem, CCG, time_limit, timer)
    run_adaptive_penalty(problem, BD, time_limit, timer)
    =#

    @show(timer)

    return 0
end

test_facility_location("data/CFLP/Cap_F10_C10.txt")
