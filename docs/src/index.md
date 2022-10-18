# Synthetic Power Grids 

A package to generate synthetic power grid dynamics that is build on top of [PowerDynamics.jl](https://github.com/JuliaEnergy/PowerDynamics.jl). 

## Installation
You can directly add the package from our github / gitlab repository.
```@julia
] add https://gitlab.pik-potsdam.de/buettner/syntheticpowergrids
```

## Overview
A main feature of the package is the struct `PGGeneration` that contains all relevant data about the synthetic grid, such as the dynamics on the nodes. Several fields in the struct have predefined values, such as the per unit system, but can be changed by the user. The second main feature is the `generate_powergrid_dynamics` function which generates a `PowerGrid` and its operation point using the information given in the `PGGeneration` struct. A number of validators are run internally to ensure the validity of the generated synthetic dynamics.

As simple example to generate a synthetic grid would look like this:

```@julia
parameters = Dict(:η => 3 * 10^(-3), :α => 5.0, :κ => π/2, :Ω => 0)
dynamics = [(1.0, get_dVOCapprox, parameters)]
pg_struct = PGGeneration(num_nodes = 10, nodal_dynamics = dynamics)
pg, op, pg_struct_new, rejections = generate_powergrid_dynamics(pg_struct)

```
The default of the grid is to generate a grid topology using [SyntheticNetworks.jl](https://github.com/PIK-ICoNe/SyntheticNetworks.jl) but it is also possible to use your own predefined topology.

Furthermore the default active power set points are sampled from a bimodal distribution, as explained in the paper, and the reactive power set points are calculated such that the voltage in the grid is close to 1 p.u., but here it is also possible to use predefined set-points.