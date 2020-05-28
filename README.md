# LatticeRules

| **Documentation** | **Build Status** | **Coverage** |
|-------------------|------------------|--------------|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PieterjanRobbe.github.io/LatticeRules.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PieterjanRobbe.github.io/LatticeRules.jl/dev) | [![Build Status](https://github.com/PieterjanRobbe/LatticeRules.jl/workflows/CI/badge.svg)](https://github.com/PieterjanRobbe/LatticeRules.jl/actions) [![Build Status](https://travis-ci.com/PieterjanRobbe/LatticeRules.jl.svg?branch=master)](https://travis-ci.com/PieterjanRobbe/LatticeRules.jl) [![Build Status](https://ci.appveyor.com/api/projects/status/github/PieterjanRobbe/LatticeRules.jl?svg=true)](https://ci.appveyor.com/project/PieterjanRobbe/LatticeRules-jl) | [![Coverage](https://codecov.io/gh/PieterjanRobbe/LatticeRules.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PieterjanRobbe/LatticeRules.jl) [![Coverage](https://coveralls.io/repos/github/PieterjanRobbe/LatticeRules.jl/badge.svg?branch=master)](https://coveralls.io/github/PieterjanRobbe/LatticeRules.jl?branch=master) |

This module provides an implementation of rank-1 lattice rules. Lattice rules generate "quasi-random" sequences of points in `d` dimensions which are equally distributed over the `d`-dimensional unit cube [0,1]<sup>d</sup>.

## Usage

To initialize a rank-1 lattice rule `lattice_rule` in `d` dimensions, use
```julia
using LatticeRules
my_lattice = LatticeRule(d)
``` 

Then
```julia
getpoint(my_lattice, 0)
```
or
```julia
my_lattice[0]
```
returns the first point of the lattice, and 
```julia
my_lattice[k]
```
returns the `k`th point of the lattice. 

For a `d`-dimensional function `f`,
```julia
f.(my_lattice[1:N])
```
gives an approximation for the integral of `f` using `N` lattice points.

Providing your own generating vector `z` is possible with
```julia
my_other_lattice = LatticeRule(z, d, n)
```
where `n` is the maximum number of points in the lattice.

In practice, it is more useful to have a shifted rank-1 lattice rule
```julia
my_shifted_lattice = ShiftedLatticeRule(d)
```
to obtain an error estimate in the same way as in the Monte Carlo method.

An existing lattice rule can be turned into a randomly shifted lattice rule using
```julia
my_other_shifted_lattice = ShiftedLatticeRule(my_lattice)
```
or
```julia
shift = rand(ndims(my_lattice)) 
my_final_shifted_lattice = ShiftedLatticeRule(my_lattice, shift)
```
optionally providing the random shift vector `shift`.

More extensive documentation can be found [here](https://PieterjanRobbe.github.io/LatticeRules.jl/dev).

## Example

A classic toy example to illustrate the Monte Carlo method is to approximate the value of &pi; by throwing random darts on a square board. Suppose we draw a circle on the board with a diameter equal to the length of one side of the square. Then, the ratio of the area of the circle to the area of the square is &pi;/4. If we now repeatedly throw darts at random on the board, the ratio of the number of darts that landed inside the circle and the total number of darts, multiplied by 4, is an approximation to &pi;.

First, generate a lattice rule in two dimensions.
```julia
using LatticeRules, Statistics
my_lattice = LatticeRule(2)
```

The function `inside_circle` checks if a dart is inside the circle:
```julia
inside_circle(x) = x[1]*x[1] + x[2]*x[2] < 1
```

Our approximation for the value of &pi; is
```
Q = 4 * mean(inside_circle.(collect(my_lattice)))
``` 
with `Q = 3.1416015625`.

## See also

- [The "Magic Point Shop" of QMC point generators and generating vectors](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/) by D. Nuyens
- [Lattice rule generating vectors](https://web.maths.unsw.edu.au/~fkuo/lattice/index.html) by F. Y. Kuo.
