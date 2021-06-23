abstract type Problem end

@enum(SubproblemType,
    LinearizedKKT,
    IndicatorKKT,
    LinearizedDual,
    IndicatorDual,
    Penalty,
)
@enum(MasterType,
    CCG,
    BD,
)
