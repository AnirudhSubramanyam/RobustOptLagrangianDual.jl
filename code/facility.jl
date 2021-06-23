struct FacilityLocation <: Problem
    NumFac::Int
    NumCust::Int
    Facilities::Vector{Int}
    Customers::Vector{Int}
    Long::Vector{Float64}
    Lat::Vector{Float64}
    Demand::Vector{Float64}
    FixedCost::Vector{Float64}
    Capacity::Vector{Float64}
    Distance::Dict{Tuple{Int, Int}, Float64}
    PenaltyCost::Vector{Float64}
    budget::Int
    Mα::Vector{Float64}
    Mβ::Vector{Float64}
    Mx::Dict{Tuple{Int, Int}, Float64}
    Mu::Vector{Float64}
    ObjScale::Float64
    CompleteRecourse::Bool

    function FacilityLocation(filename::String, budget::Int, omega_paper::Bool = true, pvalue = 1.0)
        (data, header) = readdlm(filename; header=true)
        NumFac = parse(Int, header[2])
        NumCust = parse(Int, header[4])
        Facilities = collect(1:NumFac)
        Customers = collect(1:NumCust)
        Long = data[:,2]
        Lat = data[:,3]
        Demand = data[1:NumCust,4]/1e5
        FixedCost = data[1:NumFac,5]
        Capacity = data[1:NumFac,6]
        Distance = Dict{Tuple{Int, Int}, Float64}()
        if omega_paper
            for i in Customers, j in Facilities
                xi = (Long[i], Lat[i])
                xj = (Long[j], Lat[j])
                Distance[(i,j)] = haversine(xi, xj, 3958.755866)
            end
            penalty = maximum(values(Distance))
        else
            Demand .= floor.(Demand)
            for i in Customers, j in Facilities
                xi = (Long[i], Lat[i])
                xj = (Long[j], Lat[j])
                Distance[(i,j)] = floor(euclidean(xi, xj)*20)
            end
            ds = sort(collect(Distance), by=x->x[2])
            penalty = ds[Int(ceil(pvalue*NumCust*NumFac))][2]
        end
        PenaltyCost = penalty*ones(NumCust)
        Mα = deepcopy(PenaltyCost)
        Mβ = [max(Capacity[j], maximum([Demand[i]*max(Distance[(i,j)], PenaltyCost[i]) for i in Customers])) for j in Facilities]
        Mx = deepcopy(Distance)
        for i in Customers, j in Facilities
            Mx[(i,j)] += max(Capacity[j], Demand[i]*Distance[(i,j)], Demand[i]*PenaltyCost[i])
        end
        Mu = [max(PenaltyCost[i], Demand[i]) for i in Customers]
        new(NumFac,
            NumCust,
            Facilities,
            Customers,
            Long,
            Lat,
            Demand,
            FixedCost,
            Capacity,
            Distance,
            PenaltyCost,
            budget,
            Mα,
            Mβ,
            Mx,
            Mu,
            1.0,
            true)
    end
end

function solve_deterministic_problem(FL::FacilityLocation)
    m = initializeJuMPModel()

    @variable(m, y[FL.Facilities], Bin)
    @variable(m, x[FL.Customers, FL.Facilities] >= 0)
    @variable(m, u[FL.Customers] >= 0)

    @constraint(m, [i in FL.Customers],
        sum(x[i,j] for j in FL.Facilities) + u[i] >= FL.Demand[i]
    )
    @constraint(m, [j in FL.Facilities],
        sum(x[i,j] for i in FL.Customers) <= FL.Capacity[j]*y[j]
    )
    @objective(m, Min,
        +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Distance[(i,j)]*x[i,j] for i in FL.Customers, j in FL.Facilities)
        +sum(FL.PenaltyCost[i]*u[i] for i in FL.Customers)
    )

    optimize!(m)

    return objective_value(m)
end

function init_master(FL::FacilityLocation)
    m = initializeJuMPModel()

    @variable(m, y[FL.Facilities], Bin)
    @variable(m, s >= 0)
    @objective(m, Min, s)

    return m
end

function update_master(FL::FacilityLocation, MP::JuMP.Model, SP::JuMP.Model, master::MasterType, subproblem::SubproblemType)
    s = MP[:s]
    y = MP[:y]

    if master == CCG
        z = JuMP.value.(SP[:z])
        z = Float64.(Int.(round.(z))) # should be 0-1
        x = @variable(MP, [FL.Customers, FL.Facilities], lower_bound = 0)
        u = @variable(MP, [FL.Customers], lower_bound = 0)

        @constraint(MP, [i in FL.Customers],
            sum(x[i,j] for j in FL.Facilities) + u[i] >= FL.Demand[i]
        )
        @constraint(MP, [j in FL.Facilities],
            sum(x[i,j] for i in FL.Customers) <= FL.Capacity[j]*y[j]*(1 - z[j])
        )
        @constraint(MP,
            s >= sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
                +sum(FL.Distance[(i,j)]*x[i,j] for i in FL.Customers, j in FL.Facilities)
                +sum(FL.PenaltyCost[i]*u[i] for i in FL.Customers)
        )
    end

    if master == Benders
        if subproblem == PenaltyDual
            @warn("optimality cut not finalized for $(master)-$(subproblem) for FacilityLocation instances")
        end

        if subproblem == PenaltyDual
            α = JuMP.value.(SP[:α])
            β = JuMP.value.(SP[:β])
            @constraint(MP,
                s >= sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
                    +sum(FL.Demand[i]*α[i] for i in FL.Customers)
                    -sum(FL.Capacity[j]*y[j]*β[j] for j in FL.Facilities)
            )
        end

        if subproblem ∈ [PenaltyDual, LinearizedDual, IndicatorDual]
            α = JuMP.value.(SP[:α])
            β = JuMP.value.(SP[:β])
            z = JuMP.value.(SP[:z])
            z = Float64.(Int.(round.(z))) # should be 0-1
            @constraint(MP,
                s >= sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
                    +sum(FL.Demand[i]*α[i] for i in FL.Customers)
                    -sum(FL.Capacity[j]*y[j]*β[j]*(1-z[j]) for j in FL.Facilities)
            )
        end
    end
end

function build_sp_linearized_kkt(FL::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1

    SP = initializeJuMPModel()
    @variable(SP, x[FL.Customers, FL.Facilities] >= 0)
    @variable(SP, u[FL.Customers] >= 0)
    @variable(SP, α[FL.Customers] >= 0)
    @variable(SP, β[FL.Facilities] >= 0)
    @variable(SP, wx[FL.Customers, FL.Facilities], Bin)
    @variable(SP, wu[FL.Customers], Bin)
    @variable(SP, wα[FL.Customers], Bin)
    @variable(SP, wβ[FL.Facilities], Bin)
    @variable(SP, z[FL.Facilities], Bin)

    # objective
    @objective(SP, Max,
        +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Distance[(i,j)]*x[i,j] for i in FL.Customers, j in FL.Facilities)
        +sum(FL.PenaltyCost[i]*u[i] for i in FL.Customers)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FL.Facilities) <= FL.budget)

    # primal feasibility
    @constraint(SP, [i in FL.Customers],
        sum(x[i,j] for j in FL.Facilities) + u[i] >= FL.Demand[i]
    )
    @constraint(SP, [j in FL.Facilities],
        sum(x[i,j] for i in FL.Customers) <= FL.Capacity[j]*y[j]*(1 - z[j])
    )
    # dual feasibility
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        α[i] - β[j] <= FL.Distance[(i,j)]
    )
    @constraint(SP, [i in FL.Customers],
        α[i] <= FL.PenaltyCost[i]
    )

    # complementarity - primal feas
    @constraint(SP, [i in FL.Customers],
        sum(x[i,j] for j in FL.Facilities) + u[i] <= FL.Demand[i] + (FL.Mα[i]*(1 - wα[i]))
    )
    @constraint(SP, [i in FL.Customers],
        α[i] <= FL.Mα[i]*wα[i]
    )
    @constraint(SP, [j in FL.Facilities],
        sum(x[i,j] for i in FL.Customers) >= ((FL.Capacity[j]*y[j]*(1 - z[j])) + (FL.Mβ[j]*(wβ[j] - 1)))
    )
    @constraint(SP, [j in FL.Facilities],
        β[j] <= (FL.Mβ[j]*wβ[j])
    )
    # complementarity - dual feas
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        α[i] - β[j] >= (FL.Distance[(i,j)] + (FL.Mx[(i,j)]*(wx[i,j] - 1)))
    )
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        x[i,j] <= (FL.Mx[(i,j)]*wx[i,j])
    )
    @constraint(SP, [i in FL.Customers],
        α[i] >= FL.PenaltyCost[i] + (FL.Mu[i]*(wu[i] - 1))
    )
    @constraint(SP, [i in FL.Customers],
        u[i] <= FL.Mu[i]*wu[i]
    )

    return SP
end

function build_sp_indicator_kkt(FL::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1

    SP = initializeJuMPModel()
    @variable(SP, x[FL.Customers, FL.Facilities] >= 0)
    @variable(SP, u[FL.Customers] >= 0)
    @variable(SP, α[FL.Customers] >= 0)
    @variable(SP, β[FL.Facilities] >= 0)
    @variable(SP, wx[FL.Customers, FL.Facilities], Bin)
    @variable(SP, wu[FL.Customers], Bin)
    @variable(SP, wα[FL.Customers], Bin)
    @variable(SP, wβ[FL.Facilities], Bin)
    @variable(SP, z[FL.Facilities], Bin)

    # objective
    @objective(SP, Max,
        +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Distance[(i,j)]*x[i,j] for i in FL.Customers, j in FL.Facilities)
        +sum(FL.PenaltyCost[i]*u[i] for i in FL.Customers)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FL.Facilities) <= FL.budget)

    # primal feasibility
    @constraint(SP, [i in FL.Customers],
        sum(x[i,j] for j in FL.Facilities) + u[i] >= FL.Demand[i]
    )
    @constraint(SP, [j in FL.Facilities],
        sum(x[i,j] for i in FL.Customers) <= FL.Capacity[j]*y[j]*(1 - z[j])
    )
    # dual feasibility
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        α[i] - β[j] <= FL.Distance[(i,j)]
    )
    @constraint(SP, [i in FL.Customers],
        α[i] <= FL.PenaltyCost[i]
    )

    # complementarity - primal feas
    @constraint(SP, [i in FL.Customers],
        wα[i] => {sum(x[i,j] for j in FL.Facilities) + u[i] <= FL.Demand[i]}
    )
    @constraint(SP, [i in FL.Customers],
        !wα[i] => {α[i] <= 0}
    )
    @constraint(SP, [j in FL.Facilities],
        wβ[j] => {sum(x[i,j] for i in FL.Customers) >= (FL.Capacity[j]*y[j]*(1 - z[j]))}
    )
    @constraint(SP, [j in FL.Facilities],
        !wβ[j] => {β[j] <= 0}
    )
    # complementarity - dual feas
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        wx[i,j] => {α[i] - β[j] >= FL.Distance[(i,j)]}
    )
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        !wx[i,j] => {x[i,j] <= 0}
    )
    @constraint(SP, [i in FL.Customers],
        wu[i] => {α[i] >= FL.PenaltyCost[i]}
    )
    @constraint(SP, [i in FL.Customers],
        !wu[i] => {u[i] <= 0}
    )

    return SP
end

function build_sp_linearized_dual(FL::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    SP = initializeJuMPModel()
    @variable(SP, α[FL.Customers] >= 0)
    @variable(SP, β[FL.Facilities] >= 0)
    @variable(SP, z[FL.Facilities], Bin)
    @variable(SP, βz[FL.Facilities] >= 0)

    # objective
    @objective(SP, Max,
        +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Demand[i]*α[i] for i in FL.Customers)
        -sum(FL.Capacity[j]*y[j]*(β[j] - βz[j]) for j in FL.Facilities)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FL.Facilities) <= FL.budget)

    # dual feasibility
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        α[i] - β[j] <= FL.Distance[(i,j)]
    )
    @constraint(SP, [i in FL.Customers],
        α[i] <= FL.PenaltyCost[i]
    )

    # bilinearities
    @constraint(SP, [j in FL.Facilities],
        βz[j] <= FL.Mβ[j]*z[j]
    )
    @constraint(SP, [j in FL.Facilities],
        βz[j] >= β[j] - (FL.Mβ[j]*(1 - z[j]))
    )
    @constraint(SP, [j in FL.Facilities],
        βz[j] <= β[j]
    )

    return SP
end

function build_sp_indicator_dual(FL::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    SP = initializeJuMPModel()
    @variable(SP, α[FL.Customers] >= 0)
    @variable(SP, β[FL.Facilities] >= 0)
    @variable(SP, z[FL.Facilities], Bin)
    @variable(SP, βz[FL.Facilities])

    # objective
    @objective(SP, Max,
        +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Demand[i]*α[i] for i in FL.Customers)
        -sum(FL.Capacity[j]*y[j]*(β[j] - βz[j]) for j in FL.Facilities)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FL.Facilities) <= FL.budget)

    # dual feasibility
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        α[i] - β[j] <= FL.Distance[(i,j)]
    )
    @constraint(SP, [i in FL.Customers],
        α[i] <= FL.PenaltyCost[i]
    )

    # indicators
    @constraint(SP, [j in FL.Facilities],
        z[j] => {βz[j] == β[j]}
    )
    @constraint(SP, [j in FL.Facilities],
        !z[j] => {βz[j] == 0}
    )

    return SP
end

function build_sp_fixed_penalty(FL::FacilityLocation, MP::JuMP.Model, rho::Float64)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    SP = initializeJuMPModel()
    @variable(SP, α[FL.Customers] >= 0)
    @variable(SP, β[FL.Facilities] >= 0)
    @variable(SP, z[FL.Facilities], Bin)

    # objective
    @objective(SP, Max,
        +sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Demand[i]*α[i] for i in FL.Customers)
        -sum(FL.Capacity[j]*y[j]*β[j] for j in FL.Facilities)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FL.Facilities) <= FL.budget)

    # dual feasibility
    @constraint(SP, [i in FL.Customers, j in FL.Facilities],
        α[i] - β[j] <= FL.Distance[(i,j)] + (rho*z[j])
    )
    @constraint(SP, [i in FL.Customers],
        α[i] <= FL.PenaltyCost[i]
    )

    return SP
end

function build_sp(FL::FacilityLocation, MP::JuMP.Model, subproblem::SubproblemType, rho::Float64 = 1.0)
    if subproblem == LinearizedKKT
        return build_sp_linearized_kkt(FL, MP)
    end

    if subproblem == IndicatorKKT
        return build_sp_indicator_kkt(FL, MP)
    end

    if subproblem == LinearizedDual
        return build_sp_linearized_dual(FL, MP)
    end

    if subproblem == IndicatorDual
        return build_sp_indicator_dual(FL, MP)
    end

    if subproblem == PenaltyDual
        return build_sp_fixed_penalty(FL, MP, rho)
    end
end

function solve_second_stage_problem_lagrangian(FL::FacilityLocation, MP::JuMP.Model, SP::JuMP.Model, rho::Float64)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    z = JuMP.value.(SP[:z])
    z = Float64.(Int.(round.(z))) # should be 0-1
    
    m = initializeJuMPModel()
    @variable(m, x[FL.Customers, FL.Facilities] >= 0)
    @variable(m, u[FL.Customers] >= 0)
    @constraint(m, [i in FL.Customers],
        sum(x[i,j] for j in FL.Facilities) + u[i] >= FL.Demand[i]
    )
    @constraint(m, cap[j in FL.Facilities],
        sum(x[i,j] for i in FL.Customers) <= FL.Capacity[j]*y[j]
    )

    @objective(m, Min,
        sum(FL.FixedCost[j]*y[j] for j in FL.Facilities)
        +sum(FL.Distance[(i,j)]*x[i,j] for i in FL.Customers, j in FL.Facilities)
        +sum(FL.PenaltyCost[i]*u[i] for i in FL.Customers)
        +sum(rho*z[j]*(sum(x[i,j] for i in FL.Customers)) for j in FL.Facilities)
    )
    optimize!(m)
    step = sum(z[j]*(sum(value(x[i,j]) for i in FL.Customers)) for j in FL.Facilities)
    return step
end

function record_scenario(FL::FacilityLocation, SP::JuMP.Model, scenario_list::Dict)
    z = JuMP.value.(SP[:z])
    z = Int.(round.(z)) # should be 0-1
    failedFacilities = Vector{Int}()
    for j in FL.Facilities
        if z[j] == 1
            push!(failedFacilities, j)
        end
    end
    in_list = true
    if !haskey(scenario_list, failedFacilities)
        in_list = false
        scenario_list[failedFacilities] = true
    end

    return in_list
end

