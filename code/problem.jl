using JuMP
using Mosek, MosekTools
using TimerOutputs
using Printf

abstract type Problem end

function initializeJuMPModel()
    return Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_IPAR_LOG" => 0))
end

function solve_MP(MP::JuMP.Model)
    optimize!(MP)
    status = termination_status(MP)
    if status != MOI.OPTIMAL
        @error("could not solve master problem to optimality. status = $(status)")
    end
end

function solve_SP(SP::JuMP.Model)
    optimize!(SP)
    status = termination_status(SP)
    if status != MOI.OPTIMAL
        @error("could not solve subproblem to optimality. status = $(status)")
    end
end


gap(ub, lb) = (ub - lb)/ub


@enum(SubproblemType,
    KKT,
    Penalty,
    LinearizedDual,
    IndicatorDual,
)
@enum(MasterType,
    CCG,
    BD,
)

function run_default(
    problem::Problem,
    subproblemtype::SubproblemType,
    masterproblemtype::MasterType,
    time_limit::Float64,
    timer::TimerOutput = TimerOutputs.get_defaulttimer(),
)
    algname = "$(subproblemtype)-$(masterproblemtype)"    
    if subproblemtype == KKT && masterproblemtype == BD
        @error("$algname is not guaranteed to converge.")
        return 0, -Inf, +Inf, 0.0
    end
    @timeit timer algname begin
        LB = -Inf
        UB = +Inf
        epsilon = 0.01/100.0 # optimality tolerance
        start_t = time()

        @timeit timer "Compute penalty parameter" begin
            rho = 100.0
            if subproblemtype == Penalty
                rho = compute_penalty_dual(problem)
            end
        end
        @timeit timer "Build Master" begin
            MP = init_master(problem)
        end

        @timeit timer "Main loop" begin
            iter = 0
            while time() - start_t <= time_limit
                iter += 1
                @timeit timer "Optimize Master" begin
                    solve_MP(MP)
                    LB = objective_value(MP)
                end
                @timeit timer "Build Subproblem" begin
                    SP = build_sp(problem, MP, subproblemtype, rho)
                end
                @timeit timer "Optimize Subproblem" begin
                    solve_SP(SP)
                    # @timeit timer "Optimize SecondStage" begin
                    #     step = ccg_sp_penalty_solve_second_stage(problem, MP, SP, rho)
                    #     @show(step)
                    # end
                    UB = min(UB, objective_value(SP))
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


function run_adaptive_penalty(
    problem::Problem,
    masterproblemtype::MasterType,
    time_limit::Float64,
    timer::TimerOutput = TimerOutputs.get_defaulttimer(),
)
    subproblemtype = Penalty
    algname = "AdaptivePenalty-$(masterproblemtype)"
    @timeit timer algname begin
        LB = -Inf
        UB = +Inf
        epsilon = 0.01/100.0
        rho = 1.0
        start_t = time()

        @timeit timer "Build Master Problem" begin
            MP = init_master(problem)
        end

        @timeit timer "Main loop" begin
            iter = 0
            while time() - start_t <= time_limit
                iter += 1
                @timeit timer "Optimize Master" begin
                    solve_MP(MP)
                    LB = objective_value(MP)
                end
                @timeit timer "Build Subproblem" begin
                    SP = build_sp(problem, MP, subproblemtype, rho)
                end
                while true
                    @timeit timer "Optimize Subproblem" begin
                        solve_SP(SP)
                        # if gap(objective_value(SP), LB) > epsilon
                        #     break
                        # end
                    end
                    @timeit timer "Optimize Second-Stage Problem" begin
                        step = solve_second_stage_problem_penalty(problem, MP, SP, rho)
                    end
                    if step <= epsilon
                        UB = min(UB, objective_value(SP))
                        break
                    end
                    rho *= 2.0
                    @timeit timer "Build Subproblem" begin
                        SP = build_sp(problem, MP, subproblemtype, rho)
                    end
                end

                optgap = gap(UB, LB)
                @printf("iter %3d: LB = %6.2e UB = %6.2e gap = %8.2f%% time=%8.1fsec rho = %8.2f\n", iter, LB, UB, optgap*100.0, time() - start_t, rho)
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

