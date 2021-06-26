include("main.jl")

using ArgParse, JLD2, DelimitedFiles


function parse_commandline(args)
    s = ArgParseSettings()

    @add_arg_table s begin
        "--type"
            range_tester = x -> x in ["facility", "network", "rostering"]
        "--instance"
            arg_type = String
        "--budget"
            arg_type = Int
            range_tester = x -> x > 0
        "--subproblem"
            arg_type = String
            range_tester = x -> x in ["LinearizedKKT", "IndicatorKKT", "LinearizedDual", "IndicatorDual", "PenaltyDual"]
        "--seed"
            default = 1
            arg_type = Int
            range_tester = x -> x > 0
    end

    return parse_args(args, s; as_symbols = true)
end

function run(args)
    input = parse_commandline(args)

    instance = input[:instance]
    budget = input[:budget]
    seed = input[:seed]
    inst = instances(SubproblemType)
    syms = Symbol.(inst)
    smap = Dict(zip(syms, inst))
    subproblem = smap[Symbol(input[:subproblem])]

    problem = nothing
    tiny = nothing
    time_limit = 3600.0
    if input[:type] == "facility"
        filename = joinpath(dirname(@__FILE__), "data/CFLP/$(instance).txt")
        problem = FacilityLocation(filename, budget)
        tiny = FacilityLocation(joinpath(dirname(@__FILE__), "data/CFLP/Cap_F10_C10.txt"), 1)
    elseif input[:type] == "network"
        filename = joinpath(dirname(@__FILE__), "data/SNDLIB/$(instance).txt")
        problem = NetworkDesign(filename, 1.0, budget)
        tiny = NetworkDesign(joinpath(dirname(@__FILE__), "data/SNDLIB/dfn-bwin.txt"), 1.0, 1)
    elseif input[:type] == "rostering"
        scale = parse(Int, instance)
        problem = Rostering(budget, scale, scale, seed)
        tiny = Rostering(1, 1, 1, 1)
    end

    outdir = joinpath(dirname(@__FILE__), "output", input[:type])
    outfile = joinpath(outdir, "$(instance)_s=$(seed)_b=$(budget)_sp=$(input[:subproblem]).jld2")
    if !isdir(outdir)
        mkpath(outdir)
    end

    # first pre-compile
    run_ccg(tiny, subproblem, time_limit)

    # main run
    reset_timer!()
    iter, lb, ub, tottime, iter_inner = run_ccg(problem, subproblem, time_limit)

    f = jldopen(outfile, "w")
    f["timer"] = TimerOutputs.get_defaulttimer()
    f["iter"] = iter
    f["lb"] = lb
    f["ub"] = ub
    f["tottime"] = tottime
    f["iter_inner"] = iter_inner
    close(f)
end

function batch_gen()
    basecmd = "julia --project run.jl"
    runs = Vector{String}()

    # facility
    files = readdir(joinpath(dirname(@__FILE__), "data/CFLP"))
    for f in files
        instance = splitext(f)[1]
        for budget in 1:4, subproblem in ["LinearizedKKT", "PenaltyDual"]
            cmd = "$basecmd --type=facility --instance=$instance --budget=$budget --subproblem=$subproblem"
            push!(runs, cmd)
        end
    end

    # network
    for instance in ["dfn-bwin", "dfn-gwin", "di-yuan"]
        for budget in 1:16, subproblem in ["IndicatorDual", "PenaltyDual"]
            cmd = "$basecmd --type=network --instance=$instance --budget=$budget --subproblem=$subproblem"
            push!(runs, cmd)
        end
    end

    # rostering
    for instance in [1, 2], seed in 1:10
        for budget in 3:3:18, subproblem in ["LinearizedKKT", "IndicatorKKT", "PenaltyDual"]
            cmd = "$basecmd --type=rostering --instance=$instance --seed=$seed --budget=$budget --subproblem=$subproblem"
            push!(runs, cmd)
        end
    end

    writedlm("batch_runs.txt", runs)
end

run(ARGS)



    
