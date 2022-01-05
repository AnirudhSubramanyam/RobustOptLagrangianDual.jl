module RobustOptLagrangianDual

using JuMP
using TimerOutputs
using Printf, DelimitedFiles, Distances

const SOLVER = "CPLEX"
const THREADLIM = 8
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
include("utilities.jl")
include("algorithms.jl")

export run_ccg, run_benders, solve_deterministic_problem
export FacilityLocation, NetworkDesign, Rostering
export SubproblemType, LinearizedKKT, IndicatorKKT, LinearizedDual, IndicatorDual, LagrangianDual
export MasterType, CCG, Benders

end # module


