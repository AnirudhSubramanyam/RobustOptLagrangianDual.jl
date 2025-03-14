module RobustOptLagrangianDual

using JuMP
using TimerOutputs
using Printf, DelimitedFiles, Distances

THREADLIM = 8
"""
    set_num_threads(threadlim::Int)

Set the number of threads to use for parallel branch & bound.
"""
function set_num_threads(threadlim::Int)
    global THREADLIM = threadlim
end

SOLVER = "SCIP"
"""
    set_solver_SCIP()

Set the SCIP solver.
"""
function set_solver_SCIP()
    @eval using SCIP
    global SOLVER = "SCIP"
end

GUROBI_ENV = nothing
"""
    set_solver_Gurobi()

Set the Gurobi solver.
Gurobi must be separately installed and added via `pkg> add Gurobi`.
The package defaults to SCIP if this function throws an error.
"""
function set_solver_Gurobi()
    try
        @eval using Gurobi
        global SOLVER = "Gurobi"
        if isnothing(GUROBI_ENV)
            global GUROBI_ENV = Gurobi.Env()
        end
    catch
        @warn("Unable to use Gurobi. Is it properly installed? Have you called `pkg> add Gurobi`?\nDefaulting to SCIP")
        set_solver_SCIP()
    end
end

"""
    set_solver_CPLEX()

Set the CPLEX solver.
CPLEX must be separately installed and added via `pkg> add CPLEX`.
The package defaults to SCIP if this function throws an error.
"""
function set_solver_CPLEX()
    try
        @eval using CPLEX
        global SOLVER = "CPLEX"
    catch
        @warn("Unable to use CPLEX. Is it properly installed? Have you called `pkg> add CPLEX`?\nDefaulting to SCIP")
        set_solver_SCIP()
    end
end

"""
    set_solver_Mosek()

Set the Mosek solver.
Mosek must be separately installed and added via `pkg> add Mosek, MosekTools`.
The package defaults to SCIP if this function throws an error.
"""
function set_solver_Mosek()
    try
        @eval using Mosek, MosekTools
        global SOLVER = "Mosek"
    catch
        @warn("Unable to use Mosek. Is it properly installed? Have you called `pkg> add Mosek, MosekTools`?\nDefaulting to SCIP")
        set_solver_SCIP()
    end
end

include("problem.jl")
include("facility.jl")
include("network.jl")
include("rostering.jl")
include("utilities.jl")
include("algorithms.jl")

export run_ccg, run_benders, solve_deterministic_problem
export FacilityLocation, NetworkDesign, Rostering
export SubproblemType, LinearizedKKT, IndicatorKKT, LinearizedDual, IndicatorDual, LagrangianDual
export MasterType, CCG, Benders
export set_num_threads, set_solver_SCIP, set_solver_CPLEX, set_solver_Gurobi, set_solver_Mosek

end # module


