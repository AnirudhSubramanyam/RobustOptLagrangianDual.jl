"""
    run_ccg(problem::AbstractProblem, subproblemtype::SubproblemType, time_limit::Float64, opt_tol::Float64 = 1e-4, feas_tol::Float64 = 1e-5)

Solve `problem` using the column-and-constraint generation algorithm using the
`subproblemtype` subproblem with specified `time_limit` (in seconds) and
optimality and feasibility tolerances of `opt_tol` and `feas_tol`, respectively.

Returns `(num_iter, lb, ub, total_time, iter_inner)`
"""
function run_ccg(
    problem::AbstractProblem,
    subproblemtype::SubproblemType,
    time_limit::Float64,
    opt_tol::Float64=1e-4,
    feas_tol::Float64=1e-5,
)
    if mixed_integer_recourse(problem)
        if !complete_recourse(problem)
            @error("to do - CCG algorithm for mixed-integer instances without complete recourse")
            return nothing
        end
        if subproblemtype == LagrangianDual && indicator_uncertainty(problem)
            @error("to do - CCG algorithm for mixed-integer instances with indicator uncertainties")
            return nothing
        end
        return run_ccg_mixed_integer_recourse(problem, subproblemtype, time_limit, opt_tol, feas_tol)
    end
    return run_iterative_continuous_recourse(problem, CCG, subproblemtype, time_limit, opt_tol, feas_tol)
end

"""
    run_benders(problem::AbstractProblem, subproblemtype::SubproblemType, time_limit::Float64, opt_tol::Float64 = 1e-4, feas_tol::Float64 = 1e-5)

Solve `problem` using the benders decomposition algorithm using the
`subproblemtype` subproblem with specified `time_limit` (in seconds) and
optimality and feasibility tolerances of `opt_tol` and `feas_tol`, respectively.

Returns `(num_iter, lb, ub, total_time, nothing)`
"""
function run_benders(
    problem::AbstractProblem,
    subproblemtype::SubproblemType,
    time_limit::Float64,
    opt_tol::Float64=1e-4,
    feas_tol::Float64=1e-5,
)
    if mixed_integer_recourse(problem)
        @error("Benders algorithm invalid for problems with mixed-integer recourse")
        return nothing
    end
    if subproblemtype ∈ [LinearizedKKT, IndicatorKKT]
        @error("Benders-$(subproblemtype) configuration is invalid")
        return nothing
    end
    return run_iterative_continuous_recourse(problem, Benders, subproblemtype, time_limit, opt_tol, feas_tol)
end

function run_iterative_continuous_recourse(
    problem::AbstractProblem,
    masterproblemtype::MasterType,
    subproblemtype::SubproblemType,
    time_limit::Float64,
    opt_tol::Float64,
    feas_tol::Float64,
)
    @assert !mixed_integer_recourse(problem)
    @assert masterproblemtype != Benders || subproblemtype ∉ [LinearizedKKT, IndicatorKKT]

    LB = -Inf
    UB = +Inf
    iter = 0
    start_t = time()

    SP = JuMP.Model()
    MP = init_master(problem)
    scenario_list = Dict()

    λ = nothing
    if subproblemtype == LagrangianDual
        λ = 1.0
    end

    @timeit "$(subproblemtype)-$(masterproblemtype)" begin
        while time() - start_t <= time_limit
            iter += 1

            lb = solve_MP(problem, MP, time_limit - (time() - start_t))
            !isnan(lb) && (LB = lb) # normal termination (optimal or infeasible)
            !isfinite(lb) && break  # infeasible or non-normal termination

            # Feasibility subproblem
            found_infeasible_scenario = false
            if !complete_recourse(problem)
                SP = build_feasibility_sp(problem, MP, subproblemtype)
                ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                isnan(ub) && break  # non-normal termination
                if termination_status(SP) != MOI.OBJECTIVE_LIMIT && objective_value(SP) > feas_tol
                    found_infeasible_scenario = true
                end
            end

            # Optimality subproblem
            if !found_infeasible_scenario
                if subproblemtype == LagrangianDual
                    normal_termination = true
                    while true
                        SP = build_sp(problem, MP, subproblemtype, λ)
                        ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                        if !isfinite(ub) # non-normal termination
                            normal_termination = false
                            break
                        end
                        @timeit "solve_second_stage_problem_lagrangian" begin
                            step = solve_second_stage_problem_lagrangian(problem, MP, SP, λ)
                        end
                        if step <= feas_tol
                            UB = min(UB, ub)
                            break
                        end
                        λ *= 2.0
                    end
                    normal_termination || break
                else
                    SP = build_sp(problem, MP, subproblemtype)
                    ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                    !isfinite(ub) && break  # non-normal termination
                    UB = min(UB, ub)
                end
            end

            # Check optimality of the Lagrangian multiplier
            if subproblemtype == LagrangianDual && !found_infeasible_scenario && gap(UB, LB) <= opt_tol
                @timeit "check_lagrangian_parameter" begin
                    try
                        SP = build_checking_sp(problem, MP)
                        zb = solve_SP(problem, SP, time_limit - (time() - start_t))
                        λbar = init_lagrangian_coefficient(problem, SP)
                        if gap(zb, UB) > opt_tol && gap(λbar, λ) > opt_tol
                            UB = zb
                            λ = λbar
                        end
                    catch
                        @warn("Unable to build and/or solve indicator constrained model. You may want to try a different solver.")
                    end
                end
            end

            # Print progress
            print_progress(iter, LB, UB, time() - start_t, λ)

            # Terminate or update master
            if found_infeasible_scenario || gap(UB, LB) > opt_tol
                debug_repeated_scenarios(scenario_list, masterproblemtype, problem, SP, LB, UB, found_infeasible_scenario)
                update_master_continuous(problem, MP, SP, masterproblemtype, subproblemtype)
            else
                break
            end
        end
    end
    return iter, LB, UB, time() - start_t, nothing
end

function run_ccg_mixed_integer_recourse(
    problem::AbstractProblem,
    subproblemtype::SubproblemType,
    time_limit::Float64,
    opt_tol::Float64,
    feas_tol::Float64,
)
    @assert complete_recourse(problem)
    @assert subproblemtype != LagrangianDual || !indicator_uncertainty(problem)

    masterproblemtype = CCG
    algname = "$(subproblemtype)-$(masterproblemtype)"
    LB = -Inf
    UB = +Inf
    iter = 0
    start_t = time()
    iter_inner = Vector{Int}()

    MP_outer = init_master(problem)
    scenario_list = Dict()

    λ = nothing
    if subproblemtype == LagrangianDual
        λ = 1.0
    end

    @timeit algname begin
        while time() - start_t <= time_limit
            iter += 1
            push!(iter_inner, 0)

            lb = solve_MP(problem, MP_outer, time_limit - (time() - start_t))
            !isnan(lb) && (LB = lb) # normal termination (optimal or infeasible)
            !isfinite(lb) && break  # infeasible or non-normal termination

            @timeit "$algname-InnerLevel" begin
                LB_inner = -Inf
                UB_inner = +Inf
                MP_inner = init_master_inner_level(problem)
                MP_inner_opt = nothing
                discrete_decision_list = Dict()

                if subproblemtype == LagrangianDual
                    @timeit "compute_lagrangian_parameter" begin
                        λ = max(λ, compute_lagrangian_coefficient(problem, MP_outer))
                    end
                end

                normal_termination = true
                while time() - start_t <= time_limit
                    iter_inner[end] += 1

                    ub = solve_MP(problem, MP_inner, time_limit - (time() - start_t))
                    if isnan(ub) # non-normal termination
                        normal_termination = false
                        break
                    end
                    UB_inner = ub
                    UB = min(UB, UB_inner)

                    SP = build_second_stage_problem(problem, MP_outer, MP_inner)
                    lb = solve_SP(problem, SP, time_limit - (time() - start_t))
                    if !isfinite(lb) # non-normal termination
                        normal_termination = false
                        break
                    end
                    if lb > LB_inner
                        LB_inner = lb
                        MP_inner_opt = init_master_inner_level(problem, MP_inner)
                        # if gap(LB_inner, LB) > 10.0*opt_tol
                        #     break
                        # end
                    end
                    if subproblemtype == LagrangianDual
                        record_discrete_second_stage_decision(problem, SP, discrete_decision_list)
                        while true
                            @timeit "solve_second_stage_problem_lagrangian" begin
                                step = solve_second_stage_problem_lagrangian(problem, MP_outer, MP_inner, λ)
                            end
                            if step <= feas_tol
                                break
                            end
                            λ *= 2.0
                        end

                        # Check optimality of the Lagrangian multiplier
                        if gap(UB_inner, LB_inner) <= opt_tol && gap(UB, LB) <= opt_tol
                            @timeit "check_lagrangian_parameter" begin
                                try
                                    MP_inner_check = build_master_inner_level_check(problem, MP_outer, discrete_decision_list)
                                    zb = solve_MP(problem, MP_inner_check, time_limit)
                                    λbar = init_lagrangian_coefficient(problem, MP_inner_check)
                                    if gap(zb, UB) > opt_tol && gap(λbar, λ) > opt_tol
                                        UB = zb
                                        λ = λbar
                                        MP_inner = init_master_inner_level(problem, MP_outer, discrete_decision_list, subproblemtype, λ)
                                        continue
                                    end
                                catch
                                    @warn("Unable to build and/or solve indicator constrained model. You may want to try a different solver.")
                                end
                            end
                        end
                    end

                    # Print progress
                    print_progress(iter_inner[end], LB_inner, UB_inner, time() - start_t, λ, true)

                    # Terminate or update master
                    if gap(UB_inner, LB_inner) > opt_tol
                        update_master_inner_level(problem, MP_outer, MP_inner, SP, subproblemtype, λ)
                    else
                        break
                    end
                end
                normal_termination || break
            end

            # Print progress
            print_progress(iter, LB, UB, time() - start_t, λ)

            # Terminate or update master
            if gap(UB, LB) > opt_tol
                debug_repeated_scenarios(scenario_list, masterproblemtype, problem, MP_inner_opt, LB, UB, false)
                update_master_mixed_integer(problem, MP_outer, MP_inner_opt, masterproblemtype)
            else
                break
            end
        end
    end
    return iter, LB, UB, time() - start_t, iter_inner
end
