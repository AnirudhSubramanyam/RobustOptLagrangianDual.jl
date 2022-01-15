# RobustOptLagrangianDual.jl

[RobustOptLagrangianDual.jl](https://github.com/AnirudhSubramanyam/RobustOptLagrangianDual.jl) is a Julia implementation of a Lagrangian dual algorithm for solving two-stage robust optimization problems. The details of the algorithm and problem formulation can be found in the following paper:
```
@article{subramanyam2021lagrangian,
  title={A Lagrangian Dual Method for Two-Stage Robust Optimization with Binary Uncertainties},
  author={Subramanyam, Anirudh},
  journal={arXiv preprint arXiv:2112.13138},
  year={2021}
}
```

## Installation
The package can be installed via:
```shell
$ git clone https://github.com/AnirudhSubramanyam/RobustOptLagrangianDual.jl
$ cd RobustOptLagrangianDual.jl
$ julia --project
```

Test the installation by running:
```shell
$ julia --project test/test.jl
```

By default, the package uses the open-source solver [Cbc](https://github.com/jump-dev/Cbc.jl).
One can optionally use one of the following solvers as well, which must be installed separately:
* [CPLEX](https://github.com/jump-dev/CPLEX.jl)
* [Gurobi](https://github.com/jump-dev/Gurobi.jl)
* [Mosek](https://github.com/MOSEK/Mosek.jl)

See also [Solvers](@ref).

## Documentation and Usage
```@contents
Pages = [
    "algorithms.md",
    "solvers.md",
    "internals.md",
    "examples.md",
]
Depth = 1
```

Example problem implementations can be found in `src/network.jl`, `src/facility.jl` and `src/rostering.jl`.

## Funding

This research was supported by the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research, under Contract DE-AC02-06CH11357.
