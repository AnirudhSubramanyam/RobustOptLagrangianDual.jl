
function gap(ub, lb)
    return isinf(ub) ? Inf : ((ub - lb) / ub)
end

function print_progress(iter, lb, ub, elapsed_time, λ = nothing, inner = false)
    str = inner ? "\tinner-iter" : "iter"
    if λ === nothing
        @printf("%s %3d: LB = %8.2e UB = %8.2e gap = %8.2f%% time=%8.1fsec\n", str, iter, lb, ub, gap(ub, lb) * 100.0, elapsed_time)
    else
        @printf("%s %3d: LB = %8.2e UB = %8.2e gap = %8.2f%% time=%8.1fsec λ = %8.2f\n", str, iter, lb, ub, gap(ub, lb) * 100.0, elapsed_time, λ)
    end
end

function set_optimizer_time_limit(m::JuMP.Model, time_limit::Float64)
    if SOLVER == "Mosek"
        JuMP.set_optimizer_attribute(m, "MSK_DPAR_OPTIMIZER_MAX_TIME", max(time_limit, 0.01))
    end
    if SOLVER == "Gurobi"
        JuMP.set_optimizer_attribute(m, "TimeLimit", max(time_limit, 0.01))
    end
    if SOLVER == "CPLEX"
        JuMP.set_optimizer_attribute(m, "CPXPARAM_TimeLimit", max(time_limit, 0.01))
    end
    if SOLVER == "Cbc"
        JuMP.set_optimizer_attribute(m, "seconds", max(time_limit, 0.01))
    end
end

function initializeJuMPModel()
    if SOLVER == "Mosek"
        return Model(optimizer_with_attributes(
            Mosek.Optimizer,
            "MSK_IPAR_LOG" => 0,
            "MSK_DPAR_MIO_TOL_REL_GAP" => 0,
            "MSK_IPAR_NUM_THREADS" => THREADLIM
        ))
    end
    if SOLVER == "Gurobi"
        return Model(optimizer_with_attributes(
            with_optimizer(Gurobi.Optimizer, GUROBI_ENV),
            "OutputFlag" => 0,
            "MIPGap" => 0,
            "Threads" => THREADLIM
        ))
    end
    if SOLVER == "CPLEX"
        return Model(optimizer_with_attributes(
            CPLEX.Optimizer,
            "CPXPARAM_ScreenOutput" => 0,
            "CPXPARAM_MIP_Tolerances_MIPGap" => 0,
            "CPXPARAM_Threads" => THREADLIM
        ))
    end
    if SOLVER == "Cbc"
        return Model(optimizer_with_attributes(
            Cbc.Optimizer,
            "logLevel" => 0,
            "ratioGap" => 0,
            "threads" => THREADLIM
        ))
    end
end

function solve_MP(problem::AbstractProblem, MP::JuMP.Model, time_limit::Float64)
    set_optimizer_time_limit(MP, time_limit)
    @timeit "Optimize Master" optimize!(MP)
    status = termination_status(MP)
    if complete_recourse(problem)
        if status != MOI.OPTIMAL
            @warn("could not solve master problem to optimality. status = $(status)")
            return NaN
        end
    else
        if status ∉ [MOI.OPTIMAL, MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            @warn("could not solve master problem to optimality or infeasibility. status = $(status)")
            return NaN
        end
    end
    if status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        @warn("potentially infeasible problem instance")
        return +Inf
    end
    return objective_scale(problem) * objective_value(MP)
end

function solve_SP(problem::AbstractProblem, SP::JuMP.Model, time_limit::Float64)
    set_optimizer_time_limit(SP, time_limit)
    @timeit "Optimize Subproblem" optimize!(SP)
    status = termination_status(SP)
    if complete_recourse(problem)
        if status != MOI.OPTIMAL
            @warn("could not solve subproblem to optimality. status = $(status)")
            return NaN
        end
    else
        if status ∉ [MOI.OPTIMAL, MOI.SOLUTION_LIMIT, MOI.OBJECTIVE_LIMIT]
            @warn("could not solve subproblem to optimality, solution_limit, or objective_limit. status = $(status)")
            return NaN
        end
    end
    if status == MOI.OBJECTIVE_LIMIT
        return -Inf
    end
    return objective_scale(problem) * objective_value(SP)
end

function debug_repeated_scenarios(scenario_list::Dict, masterproblemtype::MasterType, problem::AbstractProblem, SP::JuMP.Model, LB::Float64, UB::Float64, found_infeasible_scenario::Bool)
    in_list = record_scenario(problem, SP, scenario_list)
    if in_list && masterproblemtype == CCG
        found_infeasible_scenario && @warn("repeated scenario - claimed infeasible", objective_value(SP))
        found_infeasible_scenario || @warn("repeated scenario - claimed better", gap(UB, LB))
        @assert false
    end
end
