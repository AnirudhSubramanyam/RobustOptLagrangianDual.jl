# RobustOptLagrangianDual.jl

[RobustOptLagrangianDual.jl](https://github.com/AnirudhSubramanyam/RobustOpt_LagrangianDual) is a Julia implementation of a Lagrangian dual algorithm for solving two-stage robust optimization problems. The details of the algorithm and problem formulation can be found in the following paper:
```
@article{subramanyam2021lagrangian,
  title={A Lagrangian Dual Method for Two-Stage Robust Optimization with Binary Uncertainties},
  author={Subramanyam, Anirudh},
  journal={arXiv preprint arXiv:2112.13138},
  year={2021}
}
```

## Installation
The `RobustOptLagrangianDual.jl` package requires at least one of the following solvers to be installed:
* [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)
* [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)
* [Mosek.jl](https://github.com/MOSEK/Mosek.jl)

After installing one of the above, run:
```shell
$ git clone https://github.com/AnirudhSubramanyam/RobustOpt_LagrangianDual.git
$ cd RobustOpt_LagrangianDual
```
Depending on the solver you want to use, set the `SOLVER` variable in `src/RobustOptLagrangianDual.jl`.

Test the installation by running:
```shell
$ julia --project test/test.jl
```

## Documentation and Usage
```@contents
Pages = [
    "algorithms.md",
    "internals.md",
    "examples.md",
]
Depth = 1
```

Example problem implementations can be found in `src/network.jl`, `src/facility.jl` and `src/rostering.jl`.

## Funding

This research was supported by the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research, under Contract DE-AC02-06CH11357.
