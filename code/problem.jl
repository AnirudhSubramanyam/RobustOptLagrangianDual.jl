abstract type Problem end

@enum(SubproblemType,
    LinearizedKKT,
    IndicatorKKT,
    LinearizedDual,
    IndicatorDual,
    PenaltyDual,
)
@enum(MasterType,
    CCG,
    Benders,
)
