
gap(ub, lb) = isinf(ub) ? Inf : ((ub - lb)/ub)

function initializeJuMPModel()
    if solver == "Mosek"
        return Model(optimizer_with_attributes(
            Mosek.Optimizer,
            "MSK_IPAR_LOG" => 0,
            "MSK_IPAR_NUM_THREADS" => THREADLIM))
    end
    if solver == "Gurobi"
        return Model(optimizer_with_attributes(
            with_optimizer(Gurobi.Optimizer, GUROBI_ENV),
            "OutputFlag" => 0,
            "MIPGap" => 0,
            "Threads" => THREADLIM))
    end
end

function solve_MP(problem::Problem, MP::JuMP.Model)
    optimize!(MP)
    status = termination_status(MP)
    if problem.CompleteRecourse
        if status != MOI.OPTIMAL
            @error("could not solve master problem to optimality. status = $(status)")
        end
    else
        if status ∉ [MOI.OPTIMAL, MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            @error("could not solve master problem to optimality or infeasibility. status = $(status)")
        end
    end
end

function solve_SP(problem::Problem, SP::JuMP.Model)
    optimize!(SP)
    status = termination_status(SP)
    if problem.CompleteRecourse
        if status != MOI.OPTIMAL
            @error("could not solve subproblem to optimality. status = $(status)")
        end
    else
        if status ∉ [MOI.OPTIMAL, MOI.SOLUTION_LIMIT, MOI.OBJECTIVE_LIMIT]
            @error("could not solve subproblem to optimality, solution_limit, or objective_limit. status = $(status)")
        end
    end
end

function run_default(
    problem::Problem,
    subproblemtype::SubproblemType,
    masterproblemtype::MasterType,
    time_limit::Float64,
    timer::TimerOutput = TimerOutputs.get_defaulttimer(),
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
    @timeit timer algname begin
        LB = -Inf
        UB = +Inf
        epsilon = 1e-4 # optimality tolerance
        start_t = time()

        @timeit timer "Build Master" begin
            MP = init_master(problem)
        end

        @timeit timer "Main loop" begin
            iter = 0
            while time() - start_t <= time_limit
                iter += 1
                @timeit timer "Optimize Master" begin
                    solve_MP(problem, MP)
                end
                if problem.CompleteRecourse || termination_status(MP) == MOI.OPTIMAL
                    LB = problem.ObjScale * objective_value(MP)
                else
                    LB = +Inf
                    @warn("infeasible problem instance. Quitting main loop.")
                    break
                end
                @timeit timer "Build Subproblem" begin
                    SP = build_sp(problem, MP, subproblemtype)
                end
                @timeit timer "Optimize Subproblem" begin
                    solve_SP(problem, SP)
                end
                if problem.CompleteRecourse
                    UB = min(UB, problem.ObjScale*objective_value(SP))
                elseif termination_status(SP) == MOI.OBJECTIVE_LIMIT || objective_value(SP) <= epsilon
                    @timeit timer "Build Subproblem" begin
                        SP = build_sp(problem, MP, subproblemtype, 1.0, false)
                    end
                    @timeit timer "Optimize Subproblem" begin
                        solve_SP(problem, SP)
                    end
                    UB = min(UB, problem.ObjScale*objective_value(SP))
                end

                optgap = gap(UB, LB)
                @printf("iter %3d: LB = %6.2e UB = %6.2e gap = %8.2f%% time=%8.1fsec\n", iter, LB, UB, optgap*100.0, time() - start_t)
                if optgap <= epsilon
                    break
                end
                @timeit timer "Build Master" begin
                    update_master(problem, MP, SP, masterproblemtype, subproblemtype)
                end
            end
        end
    end
    return iter, LB, UB, time() - start_t
end

function run_adaptive_penalty_outer(
    problem::Problem,
    masterproblemtype::MasterType,
    time_limit::Float64,
    timer::TimerOutput = TimerOutputs.get_defaulttimer(),
)
    function compute_upper_bound(MP::JuMP.Model, UB_initial::Float64)
        UB_true = UB_initial
        rho_UB = 1000.0
        while true
            @timeit timer "Build Subproblem" begin
                SP_UB = build_sp(problem, MP, Penalty, rho_UB)
            end
            @timeit timer "Optimize Subproblem" begin
                solve_SP(problem, SP_UB)
            end
            if !problem.CompleteRecourse
                if termination_status(SP_UB) != MOI.OBJECTIVE_LIMIT && objective_value(SP_UB) > epsilon
                    break
                end
            end
            @timeit timer "Build Subproblem" begin
                SP_UB = build_sp(problem, MP, Penalty, rho_UB, false)
            end
            @timeit timer "Optimize Subproblem" begin
                solve_SP(problem, SP_UB)
            end
            @timeit timer "Optimize Second-Stage Problem" begin
                step = solve_second_stage_problem_penalty(problem, MP, SP_UB, rho_UB)
            end
            if step <= epsilon
                UB_true = min(UB_true, problem.ObjScale*objective_value(SP_UB))
                break
            end
            rho_UB *= 2.0
        end
        return UB_true
    end
    subproblemtype = Penalty
    algname = "AdaptivePenaltyOuter-$(masterproblemtype)"
    @timeit timer algname begin
        LB = -Inf
        UB = +Inf
        epsilon = 1e-4
        rho = 1.0
        start_t = time()

        @timeit timer "Build Master Problem" begin
            MP = init_master(problem)
        end

        @timeit timer "Main loop" begin
            iter = 0
            UB_inner_loop = +Inf
            while time() - start_t <= time_limit
                iter += 1
                @timeit timer "Optimize Master" begin
                    solve_MP(problem, MP)
                end
                if problem.CompleteRecourse || termination_status(MP) == MOI.OPTIMAL
                    LB = problem.ObjScale*objective_value(MP)
                else
                    LB = +Inf
                    @warn("infeasible problem instance. Quitting main loop.")
                    break
                end
                @timeit timer "Build Subproblem" begin
                    SP = build_sp(problem, MP, subproblemtype, rho)
                end
                @timeit timer "Optimize Subproblem" begin
                    if !problem.CompleteRecourse && solver == "Gurobi"
                        JuMP.set_optimizer_attribute(SP, "SolutionLimit", 1)
                        JuMP.set_optimizer_attribute(SP, "Cutoff", epsilon)
                    end
                    solve_SP(problem, SP)
                end
                if problem.CompleteRecourse
                    UB_inner_loop = min(UB_inner_loop, problem.ObjScale*objective_value(SP))
                elseif termination_status(SP) == MOI.OBJECTIVE_LIMIT || objective_value(SP) <= epsilon
                    @timeit timer "Build Subproblem" begin
                        SP = build_sp(problem, MP, subproblemtype, rho, false)
                    end
                    @timeit timer "Optimize Subproblem" begin
                        solve_SP(problem, SP)
                    end
                    UB_inner_loop = min(UB_inner_loop, problem.ObjScale*objective_value(SP))
                end

                if gap(UB_inner_loop, LB) > epsilon
                    @timeit timer "Build Master" begin
                        update_master(problem, MP, SP, masterproblemtype, subproblemtype)
                    end

                    optgap = gap(UB, LB)
                    @printf("iter %3d: LB = %8.2e UB = %8.2e gap = %8.2f%% time=%8.1fsec rho = %8.2f\n", iter, LB, UB, optgap*100.0, time() - start_t, rho)
                    continue
                end
                
                @timeit timer "Compute upper bound" begin
                    UB = compute_upper_bound(MP, UB)
                end
                @timeit timer "Optimize Second-Stage Problem" begin
                    step = solve_second_stage_problem_penalty(problem, MP, SP, rho)
                end

                optgap = gap(UB, LB)
                @printf("iter %3d: LB = %8.2e UB = %8.2e gap = %8.2f%% time=%8.1fsec rho = %8.2f\n", iter, LB, UB, optgap*100.0, time() - start_t, rho)

                if optgap <= epsilon
                    break
                end
                rho *= 2.0
                UB_inner_loop = +Inf
            end

            
        end
    end
    return iter, LB, UB, time() - start_t
end
