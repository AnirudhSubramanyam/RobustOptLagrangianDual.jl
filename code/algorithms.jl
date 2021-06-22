include("utilities.jl")

function run_default(
    problem::Problem,
    subproblemtype::SubproblemType,
    masterproblemtype::MasterType,
    time_limit::Float64,
    opt_tol::Float64 = 1e-4,
    feas_tol::Float64 = 1e-5,
)
    algname = "$(subproblemtype)-$(masterproblemtype)"    
    if subproblemtype == LinearizedKKT && masterproblemtype == BD
        @error("$algname is not guaranteed to converge.")
        return 0, -Inf, +Inf, 0.0
    end
    if subproblemtype == Penalty
        @error("$subproblemtype-subproblem not supported in 'default' mode.")
        return 0, -Inf, +Inf, 0.0
    end

    LB = -Inf
    UB = +Inf
    iter = 0
    start_t = time()

    SP = JuMP.Model()
    MP = init_master(problem)
    scenario_list = Dict()

    @timeit TIMER algname begin
        while time() - start_t <= time_limit
            iter += 1

            lb = solve_MP(problem, MP, time_limit - (time() - start_t))
            !isnan(lb) && (LB = lb) # normal termination (optimal or infeasible)
            !isfinite(lb) && break  # infeasible or non-normal termination

            # Feasibility subproblem
            found_infeasible_scenario = false
            if !problem.CompleteRecourse
                SP = build_feasibility_sp(problem, MP, subproblemtype)
                if SOLVER == "Gurobi"
                    JuMP.set_optimizer_attribute(SP, "SolutionLimit", 1)
                    JuMP.set_optimizer_attribute(SP, "Cutoff", feas_tol)
                end
                ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                isnan(ub) && break  # non-normal termination
                if termination_status(SP) != MOI.OBJECTIVE_LIMIT && objective_value(SP) > feas_tol
                    found_infeasible_scenario = true
                end
            end

            # Optimality subproblem
            if !found_infeasible_scenario
                SP = build_sp(problem, MP, subproblemtype)
                ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                !isfinite(ub) && break  # non-normal termination
                UB = min(UB, ub)
            end

            # Print progress
            print_progress(iter, LB, UB, time() - start_t)

            # Terminate or update master
            if gap(UB, LB) > opt_tol
                debug_repeated_scenarios(scenario_list, masterproblemtype, problem, SP, LB, UB, found_infeasible_scenario)
                update_master(problem, MP, SP, masterproblemtype, subproblemtype)
            else
                break
            end
        end
    end
    return iter, LB, UB, time() - start_t
end

function run_lagrangian(
    problem::Problem,
    masterproblemtype::MasterType,
    time_limit::Float64,
    opt_tol::Float64 = 1e-4,
    feas_tol::Float64 = 1e-5,
)
    subproblemtype = Penalty
    algname = "Lagrangian-$(masterproblemtype)"

    LB = -Inf
    UB = +Inf
    iter = 0
    start_t = time()

    SP = JuMP.Model()
    MP = init_master(problem)
    scenario_list = Dict()

    rho = 1.0
    UB_inner_loop = UB

    @timeit TIMER algname begin
        while time() - start_t <= time_limit
            iter += 1

            lb = solve_MP(problem, MP, time_limit - (time() - start_t))
            !isnan(lb) && (LB = lb) # normal termination (optimal or infeasible)
            !isfinite(lb) && break  # infeasible or non-normal termination
            

            # Feasibility subproblem
            found_infeasible_scenario = false
            if !problem.CompleteRecourse
                SP = build_feasibility_sp(problem, MP, subproblemtype, rho)
                if SOLVER == "Gurobi"
                    JuMP.set_optimizer_attribute(SP, "SolutionLimit", 1)
                    JuMP.set_optimizer_attribute(SP, "Cutoff", feas_tol)
                end
                ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                isnan(ub) && break  # non-normal termination
                if termination_status(SP) != MOI.OBJECTIVE_LIMIT && objective_value(SP) > feas_tol
                    found_infeasible_scenario = true
                end
            end


            # Optimality subproblem
            if !found_infeasible_scenario
                SP = build_sp(problem, MP, subproblemtype, rho)
                ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                !isfinite(ub) && break  # non-normal termination
                UB_inner_loop = min(UB_inner_loop, ub)

                # Heuristic for actual upper bound -- this is optional for convergence
                step = solve_second_stage_problem_penalty(problem, MP, SP, rho)
                (step <= feas_tol) && (UB = min(UB, ub))

                # Exact method for actual upper bound - necessary for convergence
                if step > feas_tol && gap(UB_inner_loop, LB) <= opt_tol
                    flag = false
                    while true
                        rho *= 2.0
                        SP = build_sp(problem, MP, subproblemtype, rho)
                        ub = solve_SP(problem, SP, time_limit - (time() - start_t))
                        if !isfinite(ub)# non-normal termination
                            flag = true
                            break
                        end
                        step = solve_second_stage_problem_penalty(problem, MP, SP, rho)
                        (step <= feas_tol) && (UB = min(UB, ub))
                        (step <= feas_tol) && break
                    end
                    flag && break # non-normal termination
                    if gap(UB, LB) > opt_tol
                        UB_inner_loop = UB
                    end
                end
            end

            # Print progress
            print_progress(iter, LB, UB, time() - start_t, rho)

            # Terminate or update master
            if found_infeasible_scenario || gap(UB, LB) > opt_tol
                debug_repeated_scenarios(scenario_list, masterproblemtype, problem, SP, LB, UB, found_infeasible_scenario)
                update_master(problem, MP, SP, masterproblemtype, subproblemtype)
            else
                break
            end
        end 
    end
    return iter, LB, UB, time() - start_t
end
