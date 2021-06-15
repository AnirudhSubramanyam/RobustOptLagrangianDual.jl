abstract type Problem end

@enum(SubproblemType,
    KKT,
    Penalty,
    LinearizedDual,
    IndicatorDual,
)
@enum(MasterType,
    CCG,
    BD,
)
