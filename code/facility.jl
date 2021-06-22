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

function init_master(FLP::FacilityLocation)
    m = initializeJuMPModel()

    @variable(m, y[FLP.Facilities], Bin)
    @variable(m, s >= 0)
    @objective(m, Min, s)

    return m
end

function update_master(FLP::FacilityLocation, MP::JuMP.Model, SP::JuMP.Model, master::MasterType, subproblem::SubproblemType)
    s = MP[:s]
    y = MP[:y]

    if master == CCG
        z = JuMP.value.(SP[:z])
        z = Float64.(Int.(round.(z))) # should be 0-1
        x = @variable(MP, [FLP.Customers, FLP.Facilities], lower_bound = 0)
        u = @variable(MP, [FLP.Customers], lower_bound = 0)

        @constraint(MP, [i in FLP.Customers],
            sum(x[i,j] for j in FLP.Facilities) + u[i] >= FLP.Demand[i]
        )
        @constraint(MP, [j in FLP.Facilities],
            sum(x[i,j] for i in FLP.Customers) <= FLP.Capacity[j]*y[j]*(1 - z[j])
        )
        @constraint(MP,
            s >= sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
                +sum(FLP.Distance[(i,j)]*x[i,j] for i in FLP.Customers, j in FLP.Facilities)
                +sum(FLP.PenaltyCost[i]*u[i] for i in FLP.Customers)
        )
    end

    if master == BD
        if subproblem == Penalty
            @warn("to do: optimality cut formuation not finalized for $(master)-$(subproblem)")
        end

        if subproblem == Penalty
            α = JuMP.value.(SP[:α])
            β = JuMP.value.(SP[:β])
            @constraint(MP,
                s >= sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
                    +sum(FLP.Demand[i]*α[i] for i in FLP.Customers)
                    -sum(FLP.Capacity[j]*y[j]*β[j] for j in FLP.Facilities)
            )
        end

        if subproblem ∈ [Penalty, LinearizedDual, IndicatorDual]
            α = JuMP.value.(SP[:α])
            β = JuMP.value.(SP[:β])
            z = JuMP.value.(SP[:z])
            z = Float64.(Int.(round.(z))) # should be 0-1
            @constraint(MP,
                s >= sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
                    +sum(FLP.Demand[i]*α[i] for i in FLP.Customers)
                    -sum(FLP.Capacity[j]*y[j]*β[j]*(1-z[j]) for j in FLP.Facilities)
            )
        end
    end
end

function build_sp_kkt(FLP::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1

    SP = initializeJuMPModel()
    @variable(SP, x[FLP.Customers, FLP.Facilities] >= 0)
    @variable(SP, u[FLP.Customers] >= 0)
    @variable(SP, α[FLP.Customers] >= 0)
    @variable(SP, β[FLP.Facilities] >= 0)
    @variable(SP, wx[FLP.Customers, FLP.Facilities], Bin)
    @variable(SP, wu[FLP.Customers], Bin)
    @variable(SP, wα[FLP.Customers], Bin)
    @variable(SP, wβ[FLP.Facilities], Bin)
    @variable(SP, z[FLP.Facilities], Bin)

    # objective
    @objective(SP, Max,
        +sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
        +sum(FLP.Distance[(i,j)]*x[i,j] for i in FLP.Customers, j in FLP.Facilities)
        +sum(FLP.PenaltyCost[i]*u[i] for i in FLP.Customers)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FLP.Facilities) <= FLP.budget)

    # primal feasibility
    @constraint(SP, [i in FLP.Customers],
        sum(x[i,j] for j in FLP.Facilities) + u[i] >= FLP.Demand[i]
    )
    @constraint(SP, [j in FLP.Facilities],
        sum(x[i,j] for i in FLP.Customers) <= FLP.Capacity[j]*y[j]*(1 - z[j])
    )
    # dual feasibility
    @constraint(SP, [i in FLP.Customers, j in FLP.Facilities],
        α[i] - β[j] <= FLP.Distance[(i,j)]
    )
    @constraint(SP, [i in FLP.Customers],
        α[i] <= FLP.PenaltyCost[i]
    )

    # complementarity - primal feas
    @constraint(SP, [i in FLP.Customers],
        sum(x[i,j] for j in FLP.Facilities) + u[i] <= FLP.Demand[i] + (FLP.Mα[i]*(1 - wα[i]))
    )
    @constraint(SP, [i in FLP.Customers],
        α[i] <= FLP.Mα[i]*wα[i]
    )
    @constraint(SP, [j in FLP.Facilities],
        sum(x[i,j] for i in FLP.Customers) >= ((FLP.Capacity[j]*y[j]*(1 - z[j])) + (FLP.Mβ[j]*(wβ[j] - 1)))
    )
    @constraint(SP, [j in FLP.Facilities],
        β[j] <= (FLP.Mβ[j]*wβ[j])
    )
    # complementarity - dual feas
    @constraint(SP, [i in FLP.Customers, j in FLP.Facilities],
        α[i] - β[j] >= (FLP.Distance[(i,j)] + (FLP.Mx[(i,j)]*(wx[i,j] - 1)))
    )
    @constraint(SP, [i in FLP.Customers, j in FLP.Facilities],
        x[i,j] <= (FLP.Mx[(i,j)]*wx[i,j])
    )
    @constraint(SP, [i in FLP.Customers],
        α[i] >= FLP.PenaltyCost[i] + (FLP.Mu[i]*(wu[i] - 1))
    )
    @constraint(SP, [i in FLP.Customers],
        u[i] <= FLP.Mu[i]*wu[i]
    )

    return SP
end

function build_sp_fixed_penalty(FLP::FacilityLocation, MP::JuMP.Model, rho::Float64)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    SP = initializeJuMPModel()
    @variable(SP, α[FLP.Customers] >= 0)
    @variable(SP, β[FLP.Facilities] >= 0)
    @variable(SP, z[FLP.Facilities], Bin)

    # objective
    @objective(SP, Max,
        +sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
        +sum(FLP.Demand[i]*α[i] for i in FLP.Customers)
        -sum(FLP.Capacity[j]*y[j]*β[j] for j in FLP.Facilities)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FLP.Facilities) <= FLP.budget)

    # dual feasibility
    @constraint(SP, [i in FLP.Customers, j in FLP.Facilities],
        α[i] - β[j] <= FLP.Distance[(i,j)] + (rho*z[j])
    )
    @constraint(SP, [i in FLP.Customers],
        α[i] <= FLP.PenaltyCost[i]
    )

    return SP
end

function build_sp_linearized_dual(FLP::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    SP = initializeJuMPModel()
    @variable(SP, α[FLP.Customers] >= 0)
    @variable(SP, β[FLP.Facilities] >= 0)
    @variable(SP, z[FLP.Facilities], Bin)
    @variable(SP, βtimes1minusz[FLP.Facilities] >= 0)

    # objective
    @objective(SP, Max,
        +sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
        +sum(FLP.Demand[i]*α[i] for i in FLP.Customers)
        -sum(FLP.Capacity[j]*y[j]*βtimes1minusz[j] for j in FLP.Facilities)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FLP.Facilities) <= FLP.budget)

    # dual feasibility
    @constraint(SP, [i in FLP.Customers, j in FLP.Facilities],
        α[i] - β[j] <= FLP.Distance[(i,j)]
    )
    @constraint(SP, [i in FLP.Customers],
        α[i] <= FLP.PenaltyCost[i]
    )

    # bilinearities
    @constraint(SP, [j in FLP.Facilities],
        βtimes1minusz[j] <= FLP.Mβ[j]*(1 - z[j])
    )
    @constraint(SP, [j in FLP.Facilities],
        βtimes1minusz[j] >= β[j] - (FLP.Mβ[j]*z[j])
    )
    @constraint(SP, [j in FLP.Facilities],
        βtimes1minusz[j] <= β[j]
    )

    return SP
end

function build_sp_indicator_dual(FLP::FacilityLocation, MP::JuMP.Model)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    SP = initializeJuMPModel()
    @variable(SP, α[FLP.Customers] >= 0)
    @variable(SP, β[FLP.Facilities] >= 0)
    @variable(SP, z[FLP.Facilities], Bin)
    @variable(SP, βtimes1minusz[FLP.Facilities] >= 0)

    # objective
    @objective(SP, Max,
        +sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
        +sum(FLP.Demand[i]*α[i] for i in FLP.Customers)
        -sum(FLP.Capacity[j]*y[j]*βtimes1minusz[j] for j in FLP.Facilities)
    )

    # uncertainty set
    @constraint(SP, sum(z[j] for j in FLP.Facilities) <= FLP.budget)

    # dual feasibility
    @constraint(SP, [i in FLP.Customers, j in FLP.Facilities],
        α[i] - β[j] <= FLP.Distance[(i,j)]
    )
    @constraint(SP, [i in FLP.Customers],
        α[i] <= FLP.PenaltyCost[i]
    )

    # indicators
    @constraint(SP, [j in FLP.Facilities],
        z[j] => {βtimes1minusz[j] == 0}
    )
    @constraint(SP, [j in FLP.Facilities],
        !z[j] => {βtimes1minusz[j] == β[j]}
    )

    return SP
end

function build_sp(FLP::FacilityLocation, MP::JuMP.Model, subproblem::SubproblemType, rho::Float64 = 1.0)
    if subproblem == KKT
        return build_sp_kkt(FLP, MP)
    end

    if subproblem == Penalty
        return build_sp_fixed_penalty(FLP, MP, rho)
    end

    if subproblem == LinearizedDual
        return build_sp_linearized_dual(FLP, MP)
    end

    if subproblem == IndicatorDual
        return build_sp_indicator_dual(FLP, MP)
    end
end

function solve_second_stage_problem_penalty(FLP::FacilityLocation, MP::JuMP.Model, SP::JuMP.Model, rho::Float64)
    y = JuMP.value.(MP[:y])
    y = Float64.(Int.(round.(y))) # should be 0-1
    z = JuMP.value.(SP[:z])
    z = Float64.(Int.(round.(z))) # should be 0-1
    
    m = initializeJuMPModel()
    @variable(m, x[FLP.Customers, FLP.Facilities] >= 0)
    @variable(m, u[FLP.Customers] >= 0)
    @constraint(m, [i in FLP.Customers],
        sum(x[i,j] for j in FLP.Facilities) + u[i] >= FLP.Demand[i]
    )
    @constraint(m, cap[j in FLP.Facilities],
        sum(x[i,j] for i in FLP.Customers) <= FLP.Capacity[j]*y[j]
    )

    @objective(m, Min,
        sum(FLP.FixedCost[j]*y[j] for j in FLP.Facilities)
        +sum(FLP.Distance[(i,j)]*x[i,j] for i in FLP.Customers, j in FLP.Facilities)
        +sum(FLP.PenaltyCost[i]*u[i] for i in FLP.Customers)
        +sum(rho*z[j]*(sum(x[i,j] for i in FLP.Customers)) for j in FLP.Facilities)
    )
    optimize!(m)
    step = sum(z[j]*(sum(value(x[i,j]) for i in FLP.Customers)) for j in FLP.Facilities)
    return step
end

function record_scenario(FLP::FacilityLocation, SP::JuMP.Model, scenario_list::Dict)
    z = JuMP.value.(SP[:z])
    z = Int.(round.(z)) # should be 0-1
    failedFacilities = Vector{Int}()
    for j in FLP.Facilities
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

