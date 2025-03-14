
"""
    SubproblemType

Enum representing the subproblem type. Possible values are `LinearizedKKT`,
`IndicatorKKT`, `LinearizedDual`, `IndicatorDual` and `LagrangianDual`.
"""
@enum(
    SubproblemType,
    LinearizedKKT,
    IndicatorKKT,
    LinearizedDual,
    IndicatorDual,
    LagrangianDual,
)
"""
    MasterType

Enum representing the master problem type. Possible values are `CCG` and `Benders`.
"""
@enum(
    MasterType,
    CCG,
    Benders,
)

"""
    AbstractProblem

Abstract supertype for the definition of specific problems.
"""
abstract type AbstractProblem end

"""
    mixed_integer_recourse(problem::AbstractProblem)

Return `true` if `problem` has mixed-integer recourse variables
and `false` otherwise.
"""
function mixed_integer_recourse end

"""
    complete_recourse(problem::AbstractProblem)

Return `true` if `problem` satisfies complete recourse
and `false` otherwise.
"""
function complete_recourse end

"""
    objective_scale(problem::AbstractProblem)

Return the objective scaling factor of `problem`.
"""
function objective_scale end

"""
    indicator_uncertainty(problem::AbstractProblem)

Return `true` if `problem` has indicator uncertainties only
and `false` otherwise.
"""
function indicator_uncertainty end

"""
    solve_deterministic_problem(problem::AbstractProblem)

Solve deterministic problem for nominal parameter values.
    
Return the optimal objective value.
"""
function solve_deterministic_problem end

"""
    solve_second_stage_problem_lagrangian(problem::AbstractProblem, MP::JuMP.Model, SP::JuMP.Model, λ::Float64)

Solve the second-stage Lagrangian problem with coefficient `λ`
based on the optimal values of the master model `MP` and subproblem `SP`.

Return the optimal value of the Lagrangian penalty term.
"""
function solve_second_stage_problem_lagrangian end

"""
    function compute_lagrangian_parameter(problem::AbstractProblem, MP::JuMP.Model)

Return an optimal Lagrangian coefficient based on the optimal
value of the (outer-level) master model `MP`.

Must be implemented if `indicator_uncertainty(problem) == false`.
"""
function compute_lagrangian_coefficient end

"""
    build_second_stage_problem(problem::AbstractProblem, MP_outer::JuMP.Model, MP_inner::JuMP.Model)

Return the JuMP model of the deterministic second-stage problem
based on the optimal values of the outer-level `MP_outer` and
inner-level `MP_inner` master models.

Must be implemented if `mixed_integer_recourse(problem) == true`.
"""
function build_second_stage_problem end

"""
    record_scenario(problem::AbstractProblem, SP::JuMP.Model, scenario_list::Dict)

This is only for debugging purposes. Record the optimal value of the
subproblem model `SP` in the `scenario_list` dictionary.

Return `true` if the scenario already existed in `scenario_list`
and `false` otherwise.
"""
function record_scenario end

"""
    record_discrete_second_stage_decision(problem::AbstractProblem, SP::JuMP.Model, discrete_decision_list::Dict)

Record the optimal discrete second-stage decisions from the solved
subproblem model `SP` in the `discrete_decision_list` dictionary.

Return `true` if the decision already existed in `discrete_decision_list`
and `false` otherwise.

Must be implemented if `mixed_integer_recourse(problem) == true`.
"""
function record_discrete_second_stage_decision end

"""
    init_master(problem::AbstractProblem)

Initialize the (outer-level) master model with variables and constraints.
"""
function init_master end

"""
    function init_lagrangian_coefficient(problem::AbstractProblem, P::JuMP.Model)

Return an optimal Lagrangian coefficient computing using the solved JuMP model `P`
built using either `build_checking_sp` if `mixed_integer_recourse(problem) == false`
or `build_master_inner_level_check` if `mixed_integer_recourse(problem) == true`.
"""
function init_lagrangian_coefficient end

"""
    update_master_continuous(problem::AbstractProblem, MP::JuMP.Model, SP::JuMP.Model, master::MasterType, subproblem::SubproblemType)

Update the master model `MP` based on the optimal value of the
continuous subproblem model `SP`, the `master` type and `subproblem` type.

Must be implemented if `mixed_integer_recourse(problem) == false`.
"""
function update_master_continous end

"""
    update_master_mixed_integer(problem::AbstractProblem, MP_outer::JuMP.Model, MP_inner::JuMP.Model, master::MasterType)

Update the outer-level master model `MP_outer` based on the optimal value
of the mixed-integer inner-level master model `MP_inner` and
the `master` type of the inner-level master model.

Must be implemented if `mixed_integer_recourse(problem) = true`.
"""
function update_master_mixed_integer end

"""
    build_sp(problem::AbstractProblem, MP::JuMP.Model, subproblem::SubproblemType, λ::Float64 = 1.0)

Return the JuMP model of the optimality type `subproblem` based on the
optimal value of the master model `MP` and the Lagrangian coefficient `λ`.

Must be implemented if `mixed_integer_recourse(problem) == false`.
"""
function build_sp end

"""
    build_feasibility_sp(problem::AbstractProblem, MP::JuMP.Model, subproblem::SubproblemType)

Return the JuMP model of the feasibility type `subproblem` based on the
optimal value of the master model `MP`.

Must be implemented if `mixed_integer_recourse(problem) == false && complete_recourse(problem) == false`.
"""
function build_feasibility_sp end

"""
    build_checking_sp(problem::AbstractProblem, MP::JuMP.Model)

Return the JuMP model of the optimality type subproblem, which also computes an
optimal Lagrangian coefficient based on the optimal values of the solved
master problem `MP`.

Must be implemented if `mixed_integer_recourse(problem) == false`.
"""
function build_checking_sp end

"""
    init_master_inner_level(problem::AbstractProblem)
    init_master_inner_level(problem::AbstractProblem, MP_inner::JuMP.Model)
    init_master_inner_level(problem::AbstractProblem, MP_outer::JuMP.Model, discrete_decision_list::Dict, master_inner::SubproblemType, λ = nothing)

Initialize the inner-level master model with variables and constraints.

Must be implemented if `mixed_integer_recourse(problem) == true`.

In the second version, the variable values must be fixed to
optimal values of the solved inner-level master model `MP_inner`.

In the third version, the model must be built using the
first-stage decisions in the solved outer-level master model `MP_outer`,
the discrete second-stage decisions stored in the `discrete_decision_list`
dictionary, and the `master_inner` type of the inner-level master model.
The Lagrangian coefficient `λ` must be `Float64` if `master_inner == LagrangianDual`.
"""
function init_master_inner_level end

"""
    update_master_inner_level(problem::AbstractProblem, MP_inner::JuMP.Model, MP_outer::JuMP.Model, SP::JuMP.Model, master_inner::SubproblemType, λ = nothing)

Update the inner-level master model `MP_inner` based on the optimal values
of the outer-level master model `MP_outer`, the subproblem model `SP`, and
the `master_inner` type of the inner-level master model.

The Lagrangian coefficient `λ` must be `Float64` if `master_inner == LagrangianDual`.

Must be implemented if `mixed_integer_recourse(problem) = true`.
"""
function update_master_inner_level end

"""
    build_master_inner_level_check(problem::AbstractProblem, MP_outer::JuMP.Model, discrete_decision_list::Dict)

Return the JuMP model of the inner-level master problem,
which also computes an optimal Lagrangian coefficient,
based on the optimal values of the solved outer-level master problem `MP_outer`.

Must be implemented if `mixed_integer_recourse(problem) = true`.
"""
function build_master_inner_level_check end

