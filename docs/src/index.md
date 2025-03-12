# RobustOptLagrangianDual.jl

[RobustOptLagrangianDual.jl](https://github.com/AnirudhSubramanyam/RobustOptLagrangianDual.jl) is a Julia implementation of a Lagrangian dual algorithm for solving two-stage robust optimization problems. The details of the algorithm and problem formulation can be found in the following papers:
```
@article{subramanyam2022lagrangian,
  title={A Lagrangian dual method for two-stage robust optimization with binary uncertainties},
  author={Subramanyam, Anirudh},
  journal={Optimization and Engineering},
  volume={23},
  number={4},
  pages={1831--1871},
  year={2022},
  publisher={Springer}
}
@article{lefebvre2024correction,
  title={Correction to: A Lagrangian dual method for two-stage robust optimization with binary uncertainties},
  author={Lefebvre, Henri and Subramanyam, Anirudh},
  journal={arXiv preprint arXiv:2411.04307},
  year={2024}
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

This research was supported in part by the U.S. Department of Energy, Office of Science, Advanced Scientific Computing Research, under Contract DE-AC02-06CH11357, and in part by the U.S. National Science Foundation under grant DMS-2229408.
