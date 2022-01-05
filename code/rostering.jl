using Random
"""
    Rostering <: AbstractProblem

Default implementation of `Rostering` problem.
"""
struct Rostering <: AbstractProblem
    I::Int
    J::Int
    T::Int
    N::Int
    FixedCostRegular::Array{Float64,2}
    FixedCostPartTime::Array{Float64,2}
    HourlyCostPartTime::Array{Float64,2}
    MinShiftsRegular::Vector{Int}
    MaxShiftsRegular::Vector{Int}
    MinShiftsPartTime::Vector{Int}
    MaxShiftsPartTime::Vector{Int}
    PenaltyCost::Vector{Float64}
    DemandAvg::Vector{Float64}
    DemandDev::Vector{Float64}
    budget::Int

    function Rostering(budget::Int, scale_time::Int, scale_demand::Int, seed::Int = 1)
        Random.seed!(seed)
        I = 12*scale_demand
        J = 3*scale_demand
        T = 21*scale_time
        N = 8
        FixedCostRegular = zeros(I, T)
        FixedCostPartTime = zeros(J, T)
        HourlyCostPartTime = zeros(J, T)
        MinShiftsRegular = zeros(Int, I)
        MaxShiftsRegular = zeros(Int, I)
        MinShiftsPartTime = zeros(Int, J)
        MaxShiftsPartTime = zeros(Int, J)
        PenaltyCost = zeros(T)
        DemandAvg = zeros(T)
        DemandDev = zeros(T)
        for i=1:I, t=1:T
            FixedCostRegular[i,t] = rand(5:15)
        end
        for j=1:J, t=1:T
            FixedCostPartTime[j,t] = rand(20:30)
            HourlyCostPartTime[j,t] = rand(4:8)
        end
        for i=1:I
            MinShiftsRegular[i] = rand(4:8)*scale_time
            MaxShiftsRegular[i] = rand(8:14)*scale_time
        end
        for j=1:J
            MinShiftsPartTime[j] = rand(2:4)*scale_time
            MaxShiftsPartTime[j] = rand(4:6)*scale_time
        end
        for t=1:T
            PenaltyCost[t] = rand(40:50)
            DemandAvg[t] = rand(30:80)*scale_demand
            DemandDev[t] = 0.05*DemandAvg[t]
        end

        new(I,
            J,
            T,
            N,
            FixedCostRegular,
            FixedCostPartTime,
            HourlyCostPartTime,
            MinShiftsRegular,
            MaxShiftsRegular,
            MinShiftsPartTime,
            MaxShiftsPartTime,
            PenaltyCost,
            DemandAvg,
            DemandDev,
            budget)
    end
end

mixed_integer_recourse(R::Rostering) = true

complete_recourse(R::Rostering) = true

objective_scale(R::Rostering) = 1.0

indicator_uncertainty(R::Rostering) = false

function solve_deterministic_problem(R::Rostering, xval = nothing)
    m = initializeJuMPModel()

    @variable(m, x[1:R.I,1:R.T], Bin)
    @variable(m, y[1:R.J,1:R.T], Bin)
    @variable(m, z[1:R.J,1:R.T] >= 0)
    @variable(m, w[1:R.T] >= 0)

    @objective(m, Min,
        +sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
        +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
        +sum(R.HourlyCostPartTime[j,t]*z[j,t] for j in 1:R.J, t in 1:R.T)
        +sum(R.PenaltyCost[t]*w[t] for t in 1:R.T)
    )

    @constraint(m, [i in 1:R.I, t in 1:(R.T-2)],
        x[i,t] + x[i,t+1] + x[i,t+2] <= 2
    )
    @constraint(m, [i in 1:R.I],
        sum(x[i,t] for t in 1:R.T) <= R.MaxShiftsRegular[i]
    )
    @constraint(m, [i in 1:R.I],
        sum(x[i,t] for t in 1:R.T) >= R.MinShiftsRegular[i]
    )
    @constraint(m, [j in 1:R.J, t in 1:(R.T-1)],
        y[j,t] + y[j,t+1] <= 1
    )
    @constraint(m, [j in 1:R.J],
        sum(y[j,t] for t in 1:R.T) <= R.MaxShiftsPartTime[j]
    )
    @constraint(m, [j in 1:R.J],
        sum(y[j,t] for t in 1:R.T) >= R.MinShiftsPartTime[j]
    )
    @constraint(m, [j in 1:R.J, t in 1:R.T],
        z[j,t] <= R.N*y[j,t]
    )
    @constraint(m, [t in 1:R.T],
        sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] >= R.DemandAvg[t]
    )
    if xval !== nothing
        for i in 1:R.I, t in 1:R.T
            fix(x[i,t], xval[i,t]; force = true)
        end
    end

    optimize!(m)

    return objective_value(m)
end

function init_master(R::Rostering)
    m = initializeJuMPModel()

    @variable(m, x[1:R.I,1:R.T], Bin)
    @constraint(m, [i in 1:R.I, t in 1:(R.T-2)],
        x[i,t] + x[i,t+1] + x[i,t+2] <= 2
    )
    @constraint(m, [i in 1:R.I],
        sum(x[i,t] for t in 1:R.T) <= R.MaxShiftsRegular[i]
    )
    @constraint(m, [i in 1:R.I],
        sum(x[i,t] for t in 1:R.T) >= R.MinShiftsRegular[i]
    )
    @variable(m, s >= 0)
    @objective(m, Min, s)

    return m
end

function update_master_mixed_integer(R::Rostering, MP_outer::JuMP.Model, MP_inner::JuMP.Model, master::MasterType)
    s = MP_outer[:s]
    x = MP_outer[:x]

    if master == CCG
        g = JuMP.value.(MP_inner[:g])
        g = Float64.(Int.(round.(g))) # should be 0-1

        y = @variable(MP_outer, [1:R.J,1:R.T], Bin)
        z = @variable(MP_outer, [1:R.J,1:R.T], lower_bound = 0)
        w = @variable(MP_outer, [1:R.T], lower_bound = 0)

        @constraint(MP_outer, [j in 1:R.J, t in 1:(R.T-1)],
            y[j,t] + y[j,t+1] <= 1
        )
        @constraint(MP_outer, [j in 1:R.J],
            sum(y[j,t] for t in 1:R.T) <= R.MaxShiftsPartTime[j]
        )
        @constraint(MP_outer, [j in 1:R.J],
            sum(y[j,t] for t in 1:R.T) >= R.MinShiftsPartTime[j]
        )
        @constraint(MP_outer, [j in 1:R.J, t in 1:R.T],
            z[j,t] <= R.N*y[j,t]
        )
        @constraint(MP_outer, [t in 1:R.T],
            sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] >= R.DemandAvg[t] + (R.DemandDev[t]*g[t])
        )
        @constraint(MP_outer,
            s >= sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
                +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
                +sum(R.HourlyCostPartTime[j,t]*z[j,t] for j in 1:R.J, t in 1:R.T)
                +sum(R.PenaltyCost[t]*w[t] for t in 1:R.T)
        )
    end

    if master == Benders
        error("Master problem of type $master not supported by Rostering instances")
    end
end

function init_master_inner_level(R::Rostering)
    m = initializeJuMPModel()

    @variable(m, g[1:R.T], Bin)
    @variable(m, s,
        upper_bound = sum(R.FixedCostRegular) +
                      sum(R.FixedCostPartTime) +
                      (R.N*sum(R.HourlyCostPartTime)) +
                      sum((R.DemandAvg[t] + R.DemandDev[t])*R.PenaltyCost[t] for t in 1:R.T)
    )
    @constraint(m,
        sum(g[t] for t in 1:R.T) <= R.budget
    )
    @objective(m, Max, s)

    return m
end

function init_master_inner_level(R::Rostering, MP_inner::JuMP.Model)
    gval = JuMP.value.(MP_inner[:g])
    gval = Float64.(Int.(round.(gval))) # should be 0-1

    m = initializeJuMPModel()
    @variable(m, g[1:R.T])
    @objective(m, Max, 0)
    fix.(g, gval)
    optimize!(m)

    return m
end

function update_master_inner_level(R::Rostering, MP_outer::JuMP.Model, MP_inner::JuMP.Model, SP_inner::JuMP.Model, master_inner::SubproblemType, λ = nothing)
    s = MP_inner[:s]
    g = MP_inner[:g]
    x = JuMP.value.(MP_outer[:x])
    x = Float64.(Int.(round.(x))) # should be 0-1
    y = JuMP.value.(SP_inner[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1

    # dual variables (common to all)
    α = @variable(MP_inner, [1:R.J,1:R.T], lower_bound = 0)
    β = @variable(MP_inner, [1:R.T], lower_bound = 0)

    @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
        -α[j,t] + β[t] <= R.HourlyCostPartTime[j,t]
    )
    @constraint(MP_inner, [t in 1:R.T],
        β[t] <= R.PenaltyCost[t]
    )

    if master_inner ∈ [LinearizedKKT, IndicatorKKT]
        z = @variable(MP_inner, [1:R.J,1:R.T], lower_bound = 0)
        w = @variable(MP_inner, [1:R.T], lower_bound = 0)
        Iz = @variable(MP_inner, [1:R.J,1:R.T], Bin)
        Iα = @variable(MP_inner, [1:R.J,1:R.T], Bin)
        Iw = @variable(MP_inner, [1:R.T], Bin)
        Iβ = @variable(MP_inner, [1:R.T], Bin)

        @constraint(MP_inner,
            s <= sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
                +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
                +sum(R.HourlyCostPartTime[j,t]*z[j,t] for j in 1:R.J, t in 1:R.T)
                +sum(R.PenaltyCost[t]*w[t] for t in 1:R.T)
        )

        # primal feasibility
        @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
            z[j,t] <= R.N*y[j,t]
        )
        @constraint(MP_inner, [t in 1:R.T],
            sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] >= R.DemandAvg[t] + (R.DemandDev[t]*g[t])
        )

        # complementary slackness
        if master_inner == LinearizedKKT
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                α[j,t] <= R.PenaltyCost[t]*(1 - Iα[j,t])
            )
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                (R.N*y[j,t]) - z[j,t] <= R.N*Iα[j,t]
            )
            @constraint(MP_inner, [t in 1:R.T],
                β[t] <= R.PenaltyCost[t]*(1 - Iβ[t])
            )
            @constraint(MP_inner, [t in 1:R.T],
                sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] - (R.DemandAvg[t] + (R.DemandDev[t]*g[t])) <= R.I*R.N*Iβ[t]
            )
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                z[j,t] <= R.N*(1 - Iz[j,t])
            )
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                α[j,t] - β[t] + R.HourlyCostPartTime[j,t] <= (R.HourlyCostPartTime[j,t] + R.PenaltyCost[t])*Iz[j,t]
            )
            @constraint(MP_inner, [t in 1:R.T],
                w[t] <= (R.DemandAvg[t] + R.DemandDev[t])*(1 - Iw[t])
            )
            @constraint(MP_inner, [t in 1:R.T],
                R.PenaltyCost[t] - β[t] <= R.PenaltyCost[t]*Iw[t]
            )
        end

        if master_inner == IndicatorKKT
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                Iα[j,t] => {α[j,t] <= 0}
            )
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                !Iα[j,t] => {(R.N*y[j,t]) - z[j,t] <= 0}
            )
            @constraint(MP_inner, [t in 1:R.T],
                Iβ[t] => {β[t] <= 0}
            )
            @constraint(MP_inner, [t in 1:R.T],
                !Iβ[t] => {sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] - (R.DemandAvg[t] + (R.DemandDev[t]*g[t])) <= 0}
            )
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                Iz[j,t] => {z[j,t] <= 0}
            )
            @constraint(MP_inner, [j in 1:R.J, t in 1:R.T],
                !Iz[j,t] => {α[j,t] - β[t] + R.HourlyCostPartTime[j,t] <= 0}
            )
            @constraint(MP_inner, [t in 1:R.T],
                Iw[t] => {w[t] <= 0}
            )
            @constraint(MP_inner, [t in 1:R.T],
                !Iw[t] => {R.PenaltyCost[t] - β[t] <= 0}
            )
        end
    end

    if master_inner ∈ [LinearizedDual, IndicatorDual]
        βg = @variable(MP_inner, [1:R.T], lower_bound = 0)
        @constraint(MP_inner,
            s <= sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
                +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
                -sum(R.N*y[j,t]*α[j,t] for j in 1:R.J, t in 1:R.T)
                +sum((R.DemandAvg[t]*β[t]) + (R.DemandDev[t]*βg[t]) for t in 1:R.T)
                -sum(R.N*x[i,t]*β[t] for i in 1:R.I, t in 1:R.T)
        )
        if master_inner == LinearizedDual
            @constraint(MP_inner, [t in 1:R.T],
                βg[t] >= β[t] - (R.PenaltyCost[t]*(1 - g[t]))
            )
            @constraint(MP_inner, [t in 1:R.T],
                βg[t] <= R.PenaltyCost[t]*g[t]
            )
            @constraint(MP_inner, [t in 1:R.T],
                βg[t] <= β[t]
            )
        end
        if master_inner == IndicatorDual
            @constraint(MP_inner, [t in 1:R.T],
                g[t] => {βg[t] == β[t]}
            )
            @constraint(MP_inner, [t in 1:R.T],
                !g[t] => {βg[t] == 0}
            )
        end
    end

    if master_inner ∈ [LagrangianDual]
        #= default (unscaled)
        γ = @variable(MP_inner, [1:R.T], lower_bound = 0)
        @constraint(MP_inner,
            s <= sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
                +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
                -sum(R.N*y[j,t]*α[j,t] for j in 1:R.J, t in 1:R.T)
                +sum(R.DemandAvg[t]*β[t] for t in 1:R.T)
                -sum(R.N*x[i,t]*β[t] for i in 1:R.I, t in 1:R.T)
                -sum(γ[t] for t in 1:R.T)
                +sum(λ*g[t] for t in 1:R.T)
        )
        @constraint(MP_inner, [t in 1:R.T],
            -γ[t] - (R.DemandDev[t]*β[t]) <= λ*(1 - 2g[t])
        )
        =#
        ###= scale gamma
        γ = @variable(MP_inner, [1:R.T])
        @constraint(MP_inner,
            s <= sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
                +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
                -sum(R.N*y[j,t]*α[j,t] for j in 1:R.J, t in 1:R.T)
                +sum(R.DemandAvg[t]*β[t] for t in 1:R.T)
                -sum(R.N*x[i,t]*β[t] for i in 1:R.I, t in 1:R.T)
                +sum(λ*γ[t] for t in 1:R.T)
        )
        @constraint(MP_inner, [t in 1:R.T],
            γ[t] <= g[t]
        )
        @constraint(MP_inner, [t in 1:R.T],
            γ[t] <= (R.DemandDev[t]*β[t]/λ) + (1 - g[t])
        )
        ##=#
    end
end

function compute_lagrangian_coefficient(R::Rostering, MP_outer::JuMP.Model)
    x = JuMP.value.(MP_outer[:x])
    x = Float64.(Int.(round.(x))) # should be 0-1
    firstStageCost = sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)

    temparr = sort(R.PenaltyCost.*R.DemandDev, rev = true)
    
    L = U = 0
    L = solve_deterministic_problem(R, x)
    # L += firstStageCost
    # L += sum(R.MinShiftsPartTime[j]*minimum(R.FixedCostPartTime[j,:]) for j in 1:R.J)
    # L += sum(min(minimum(R.HourlyCostPartTime[:,t]), R.PenaltyCost[t])*max(R.DemandAvg[t] - (R.N*sum(x[i,t] for i in 1:R.I)), 0.0) for t in 1:R.T)
    U += firstStageCost
    U += sum(R.MaxShiftsPartTime[j]*maximum(R.FixedCostPartTime[j,:] .+ (R.N*R.HourlyCostPartTime[j,:])) for j in 1:R.J)
    U += sum(R.PenaltyCost[t]*R.DemandAvg[t] for t in 1:R.T)
    U += sum(temparr[t] for t in 1:min(R.budget, R.T))
    @assert L <= U

    return U - L
end

function build_second_stage_problem(R::Rostering, MP_outer::JuMP.Model, MP_inner::JuMP.Model)
    x = JuMP.value.(MP_outer[:x])
    x = Float64.(Int.(round.(x))) # should be 0-1
    g = JuMP.value.(MP_inner[:g])
    g = Float64.(Int.(round.(g))) # should be 0-1

    m = initializeJuMPModel()

    @variable(m, y[1:R.J,1:R.T], Bin)
    @variable(m, z[1:R.J,1:R.T] >= 0)
    @variable(m, w[1:R.T] >= 0)

    @objective(m, Min,
        +sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
        +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
        +sum(R.HourlyCostPartTime[j,t]*z[j,t] for j in 1:R.J, t in 1:R.T)
        +sum(R.PenaltyCost[t]*w[t] for t in 1:R.T)
    )

    @constraint(m, [j in 1:R.J, t in 1:(R.T-1)],
        y[j,t] + y[j,t+1] <= 1
    )
    @constraint(m, [j in 1:R.J],
        sum(y[j,t] for t in 1:R.T) <= R.MaxShiftsPartTime[j]
    )
    @constraint(m, [j in 1:R.J],
        sum(y[j,t] for t in 1:R.T) >= R.MinShiftsPartTime[j]
    )
    @constraint(m, [j in 1:R.J, t in 1:R.T],
        z[j,t] <= R.N*y[j,t]
    )
    @constraint(m, [t in 1:R.T],
        sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] >= R.DemandAvg[t] + (R.DemandDev[t]*g[t])
    )

    return m
end

function solve_second_stage_problem_lagrangian(R::Rostering, MP_outer::JuMP.Model, MP_inner::JuMP.Model, λ::Float64)
    x = JuMP.value.(MP_outer[:x])
    x = Float64.(Int.(round.(x))) # should be 0-1
    g = JuMP.value.(MP_inner[:g])
    g = Float64.(Int.(round.(g))) # should be 0-1

    m = initializeJuMPModel()

    @variable(m, y[1:R.J,1:R.T], Bin)
    @variable(m, z[1:R.J,1:R.T] >= 0)
    @variable(m, w[1:R.T] >= 0)
    @variable(m, 0 <= u[1:R.T] <= 1)

    @objective(m, Min,
        +sum(R.FixedCostRegular[i,t]*x[i,t] for i in 1:R.I, t in 1:R.T)
        +sum(R.FixedCostPartTime[j,t]*y[j,t] for j in 1:R.J, t in 1:R.T)
        +sum(R.HourlyCostPartTime[j,t]*z[j,t] for j in 1:R.J, t in 1:R.T)
        +sum(R.PenaltyCost[t]*w[t] for t in 1:R.T)
        +sum(λ*(((1 - 2g[t])*u[t]) + g[t]) for t in 1:R.T)
    )

    @constraint(m, [j in 1:R.J, t in 1:(R.T-1)],
        y[j,t] + y[j,t+1] <= 1
    )
    @constraint(m, [j in 1:R.J],
        sum(y[j,t] for t in 1:R.T) <= R.MaxShiftsPartTime[j]
    )
    @constraint(m, [j in 1:R.J],
        sum(y[j,t] for t in 1:R.T) >= R.MinShiftsPartTime[j]
    )
    @constraint(m, [j in 1:R.J, t in 1:R.T],
        z[j,t] <= R.N*y[j,t]
    )
    @constraint(m, [t in 1:R.T],
        sum(R.N*x[i,t] for i in 1:R.I) + sum(z[j,t] for j in 1:R.J) + w[t] >= R.DemandAvg[t] + (R.DemandDev[t]*u[t])
    )

    optimize!(m)

    step = sum(((1 - 2g[t])*value(u[t])) + g[t] for t in 1:R.T)

    return step
end

function record_scenario(R::Rostering, MP_inner::JuMP.Model, scenario_list::Dict)
    g = JuMP.value.(MP_inner[:g])
    g = Int.(round.(g)) # should be 0-1
    demandDeviations = Vector{Int}()
    for t in 1:R.T
        if g[t] == 1
            push!(demandDeviations, t)
        end
    end
    in_list = true
    if !haskey(scenario_list, demandDeviations)
        in_list = false
        scenario_list[demandDeviations] = true
    end

    return in_list
end
