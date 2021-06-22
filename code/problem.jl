abstract type Problem end

@enum(SubproblemType,
    LinearizedKKT,
    LinearizedDual,
    IndicatorDual,
    Penalty,
)
@enum(MasterType,
    CCG,
    BD,
)
