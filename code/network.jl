struct NetworkDesign <: Problem
    NumNodes::Int
    NumEdges::Int
    Nodes::Vector{String}
    Edges::Vector{String}
    Orig::Dict{String,String}
    Dest::Dict{String,String}
    PreInstalledCap::Dict{String,Float64}
    PerUnitCapCost::Dict{String,Float64}
    PerUnitTrCost::Dict{String,Float64}
    DisAllowCap::Vector{String}
    Demand::Dict{String,Float64}
    Multiplier::Float64
    budget::Int
    ObjScale::Float64
    CompleteRecourse::Bool

    function NetworkDesign(filename::String, Multiplier::Float64, budget::Int)
        data = readdlm(filename)
        pos = findall(x->x=="NODES",data); @assert length(pos)==1; (r,c) = Tuple(pos[1]); @assert c==1

        Nodes = Vector{String}()
        for rowIdx in r+1:size(data,1)
            if data[rowIdx,1] == ")"
                break
            end
            push!(Nodes, string(data[rowIdx,1]))
        end

        pos = findall(x->x=="LINKS",data); @assert length(pos)==1; (r,c) = Tuple(pos[1]); @assert c==1
        Edges = Vector{String}()
        Orig = Dict{String,String}()
        Dest = Dict{String,String}()
        PreInstalledCap = Dict{String,Float64}()
        PerUnitCapCost = Dict{String,Float64}()
        PerUnitTrCost = Dict{String,Float64}()
        DisAllowCap = Vector{String}()
        for rowIdx in r+1:size(data,1)
            if data[rowIdx,1] == ")"
                break
            end
            e = string(data[rowIdx,1])
            push!(Edges, e)
            Orig[e] = string(data[rowIdx,3])
            Dest[e] = string(data[rowIdx,4])
            PreInstalledCap[e] = data[rowIdx,6]
            if isa(data[rowIdx,11], Number) && isa(data[rowIdx,12], Number)
                ModuleCap = data[rowIdx,11]
                ModuleCost = data[rowIdx,12]
                PerUnitCapCost[e] = ModuleCost/ModuleCap
                PerUnitTrCost[e] = Multiplier*PerUnitCapCost[e]
            else
                PerUnitCapCost[e] = NaN
                PerUnitTrCost[e] = 0
                push!(DisAllowCap, e)
            end
        end

        pos = findall(x->x=="DEMANDS",data); @assert length(pos)==1; (r,c) = Tuple(pos[1]); @assert c==1
        Demand = Dict{String,Float64}()
        for i in Nodes
            Demand[i] = 0
        end
        for rowIdx in r+1:size(data,1)
            if data[rowIdx,1] == ")"
                break
            end
            src = string(data[rowIdx,3])
            tgt = string(data[rowIdx,4])
            demand = data[rowIdx,7]
            Demand[tgt] += demand
            Demand[src] -= demand
        end
        ObjScale = max(1.0, maximum(abs.(values(Demand))))
        for i in Nodes
            Demand[i] /= ObjScale
        end
        for e in Edges
            PreInstalledCap[e] /= ObjScale
        end
        NumNodes = length(Nodes)
        NumEdges = length(Edges)

        new(NumNodes,
            NumEdges,
            Nodes,
            Edges,
            Orig,
            Dest,
            PreInstalledCap,
            PerUnitCapCost,
            PerUnitTrCost,
            DisAllowCap,
            Demand,
            Multiplier,
            budget,
            ObjScale,
            false)
    end
end

function solve_deterministic_problem(ND::NetworkDesign)
    m = initializeJuMPModel()

    @variable(m, u[ND.Edges] >= 0)
    @variable(m, ff[ND.Edges] >= 0)
    @variable(m, fb[ND.Edges] >= 0)
    for e in ND.DisAllowCap
        fix(u[e], 0.0; force = true)
    end

    @objective(m, Min,
        sum((ND.PerUnitCapCost[e]*u[e]) + (ND.PerUnitTrCost[e]*(ff[e] + fb[e])) for e in ND.Edges if e ∉ ND.DisAllowCap)
    )

    @constraint(m, [e in ND.Edges],
        ff[e] <= u[e] + ND.PreInstalledCap[e]
    )
    @constraint(m, [e in ND.Edges],
        fb[e] <= u[e] + ND.PreInstalledCap[e]
    )
    @constraint(m, [i in ND.Nodes],
        sum(ff[e] - fb[e] for e in ND.Edges if ND.Dest[e]==i) + sum(fb[e]-ff[e] for e in ND.Edges if ND.Orig[e]==i) == ND.Demand[i]
    )

    optimize!(m)

    return ND.ObjScale * objective_value(m)
end

function init_master(ND::NetworkDesign)
    m = initializeJuMPModel()

    @variable(m, u[ND.Edges] >= 0)
    for e in ND.DisAllowCap
        fix(u[e], 0.0; force = true)
    end
    @variable(m, s >= 0)
    @objective(m, Min, s)

    return m
end

function update_master(ND::NetworkDesign, MP::JuMP.Model, SP::JuMP.Model, master::MasterType, subproblem::SubproblemType)
    s = MP[:s]
    u = MP[:u]

    if master == CCG
        z = JuMP.value.(SP[:z])
        z = Float64.(Int.(round.(z))) # should be 0-1
        ff = @variable(MP, [ND.Edges], lower_bound=0)
        fb = @variable(MP, [ND.Edges], lower_bound=0)

        @constraint(MP, [e in ND.Edges],
            ff[e] <= (u[e] + ND.PreInstalledCap[e])*(1 - z[e])
        )
        @constraint(MP, [e in ND.Edges],
            fb[e] <= (u[e] + ND.PreInstalledCap[e])*(1 - z[e])
        )
        @constraint(MP, [i in ND.Nodes],
            sum(ff[e] - fb[e] for e in ND.Edges if ND.Dest[e]==i) + sum(fb[e]-ff[e] for e in ND.Edges if ND.Orig[e]==i) == ND.Demand[i]
        )
        @constraint(MP,
            s >= sum((ND.PerUnitCapCost[e]*u[e]) + (ND.PerUnitTrCost[e]*(ff[e] + fb[e])) for e in ND.Edges if e ∉ ND.DisAllowCap)
        )
    end

    if master == BD
        @assert false
    end
end

function build_sp_indicator_dual(ND::NetworkDesign, MP::JuMP.Model, feasibility::Bool)
    u = JuMP.value.(MP[:u])

    SP = initializeJuMPModel()
    if feasibility
        @variable(SP, -1 <= μ[ND.Nodes] <= 1)
        @variable(SP, -1 <= vb[ND.Edges] <= 0)
        @variable(SP, -1 <= vf[ND.Edges] <= 0)
    else
        @variable(SP, μ[ND.Nodes])
        @variable(SP, vb[ND.Edges] <= 0)
        @variable(SP, vf[ND.Edges] <= 0)
    end
    @variable(SP, zb[ND.Edges] <= 0)
    @variable(SP, zf[ND.Edges] <= 0)
    @variable(SP, z[ND.Edges], Bin)

    # objective
    if feasibility
        @objective(SP, Max,
            +sum((ND.PreInstalledCap[e]+u[e])*(vb[e]-zb[e]) for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum((ND.PreInstalledCap[e]+u[e])*(vf[e]-zf[e]) for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum(ND.Demand[i]*μ[i] for i in ND.Nodes)
        )
    else
        @objective(SP, Max,
            +sum((ND.PerUnitCapCost[e]*u[e]) for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum((ND.PreInstalledCap[e]+u[e])*(vb[e]-zb[e]) for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum((ND.PreInstalledCap[e]+u[e])*(vf[e]-zf[e]) for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum(ND.Demand[i]*μ[i] for i in ND.Nodes)
        )
    end

    # uncertainty set
    @constraint(SP, sum(z[e] for e in ND.Edges) <= ND.budget)
    @constraint(SP, [i in ND.Nodes; abs(ND.ObjScale*ND.Demand[i]) > 1e-3],
        sum(z[e]-1 for e in ND.Edges if ND.Dest[e] == i) + sum(z[e]-1 for e in ND.Edges if ND.Orig[e] == i) <= - 1
    )


    # dual feasibility
    @constraint(SP, [e in ND.Edges],
        vb[e] + μ[ND.Orig[e]] - μ[ND.Dest[e]] <= (feasibility ? 0 : ND.PerUnitTrCost[e])
    )
    @constraint(SP, [e in ND.Edges],
        vf[e] - μ[ND.Orig[e]] + μ[ND.Dest[e]] <= (feasibility ? 0 : ND.PerUnitTrCost[e])
    )

    # indicators
    @constraint(SP, [e in ND.Edges],
        !z[e] => {zb[e] >= 0}
    )
    @constraint(SP, [e in ND.Edges],
        z[e] => {zb[e] == vb[e]}
    )
    @constraint(SP, [e in ND.Edges],
        !z[e] => {zf[e] >= 0}
    )
    @constraint(SP, [e in ND.Edges],
        z[e] => {zf[e] == vf[e]}
    )

    return SP
end

function build_sp_fixed_penalty(ND::NetworkDesign, MP::JuMP.Model, rho::Float64, feasibility::Bool)
    u = JuMP.value.(MP[:u])

    SP = initializeJuMPModel()
    if feasibility
        @variable(SP, -1 <= μ[ND.Nodes] <= 1)
        @variable(SP, -1 <= vb[ND.Edges] <= 0)
        @variable(SP, -1 <= vf[ND.Edges] <= 0)
    else
        @variable(SP, μ[ND.Nodes])
        @variable(SP, vb[ND.Edges] <= 0)
        @variable(SP, vf[ND.Edges] <= 0)
    end
    @variable(SP, z[ND.Edges], Bin)
    
    # objective
    if feasibility
        @objective(SP, Max,
            +sum((ND.PreInstalledCap[e]+u[e])*vb[e] for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum((ND.PreInstalledCap[e]+u[e])*vf[e] for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum(ND.Demand[i]*μ[i] for i in ND.Nodes)
        )
    else
        @objective(SP, Max,
            +sum((ND.PerUnitCapCost[e]*u[e]) for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum((ND.PreInstalledCap[e]+u[e])*vb[e] for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum((ND.PreInstalledCap[e]+u[e])*vf[e] for e in ND.Edges if e ∉ ND.DisAllowCap)
            +sum(ND.Demand[i]*μ[i] for i in ND.Nodes)
        )
    end

    # uncertainty set
    @constraint(SP, sum(z[e] for e in ND.Edges) <= ND.budget)
    @constraint(SP, [i in ND.Nodes; abs(ND.ObjScale*ND.Demand[i]) > 1e-3],
        sum(z[e]-1 for e in ND.Edges if ND.Dest[e] == i) + sum(z[e]-1 for e in ND.Edges if ND.Orig[e] == i) <= - 1
    )

    # dual feasibility
    @constraint(SP, [e in ND.Edges],
        vb[e] + μ[ND.Orig[e]] - μ[ND.Dest[e]] <= (feasibility ? 0 : ND.PerUnitTrCost[e]) + (rho*z[e])
    )
    @constraint(SP, [e in ND.Edges],
        vf[e] - μ[ND.Orig[e]] + μ[ND.Dest[e]] <= (feasibility ? 0 : ND.PerUnitTrCost[e]) + (rho*z[e])
    )

    return SP
end

function build_sp(ND::NetworkDesign, MP::JuMP.Model, subproblem::SubproblemType, rho::Float64 = 1.0, feasibility::Bool=true)
    if subproblem == LinearizedKKT
        @assert false
    end

    if subproblem == Penalty
        return build_sp_fixed_penalty(ND, MP, rho, feasibility)
    end

    if subproblem == LinearizedDual
        @assert false
    end

    if subproblem == IndicatorDual
        return build_sp_indicator_dual(ND, MP, feasibility)
    end
end

function solve_second_stage_problem_penalty(ND::NetworkDesign, MP::JuMP.Model, SP::JuMP.Model, rho::Float64)
    u = JuMP.value.(MP[:u])
    z = JuMP.value.(SP[:z])
    z = Float64.(Int.(round.(z))) # should be 0-1
    
    m = initializeJuMPModel()
    @variable(m, ff[ND.Edges] >= 0)
    @variable(m, fb[ND.Edges] >= 0)

    @objective(m, Min,
        +sum((ND.PerUnitCapCost[e]*u[e]) + (ND.PerUnitTrCost[e]*(ff[e] + fb[e])) for e in ND.Edges if e ∉ ND.DisAllowCap)
        +sum(rho*z[e]*(ff[e] + fb[e]) for e in ND.Edges)
    )

    @constraint(m, [e in ND.Edges],
        ff[e] <= u[e] + ND.PreInstalledCap[e]
    )
    @constraint(m, [e in ND.Edges],
        fb[e] <= u[e] + ND.PreInstalledCap[e]
    )
    @constraint(m, [i in ND.Nodes],
        sum(ff[e] - fb[e] for e in ND.Edges if ND.Dest[e]==i) + sum(fb[e]-ff[e] for e in ND.Edges if ND.Orig[e]==i) == ND.Demand[i]
    )

    optimize!(m)
    step = sum(z[e]*(value(ff[e]) + value(fb[e])) for e in ND.Edges)
    return step
end

function record_scenario(ND::NetworkDesign, SP::JuMP.Model, scenario_list::Dict)
    z = JuMP.value.(SP[:z])
    z = Int.(round.(z)) # should be 0-1
    failedEdges = Vector{String}()
    for e in ND.Edges
        if z[e] == 1
            push!(failedEdges, e)
        end
    end
    in_list = true
    if !haskey(scenario_list, failedEdges)
        in_list = false
        scenario_list[failedEdges] = true
    end

    return in_list
end