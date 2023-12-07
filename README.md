# SpectraShooter.jl
Parameter searching for possible spectra(a collection of eigenvalues like masses) of a Boundary Value Problem(BVP) via the shooting method. Originally designed for analysis of the radion stabilization problem in extra-dimensions, see [cite]
## Features
- [x] Shooting method for solving Goldberger-Wise EoMs with the given boundary conditions
- [x] Extract the mass spectrum w.r.t to the given parameter(including the backreaction parameter $l^2$ and the quadratic superpotential parameter $\gamma^2$)
- [x] Optimized the shooting method by using the lazy matrix method, which reduces most of the unnecessary calls of the ODE solver
## Usage
### Pre-requisites
- `DifferentialEquations.jl`
### Installing
- Clone this repository and switch to the directory
- Start Julia and run the following commands
```jldoctest
julia> include("LazyMatrices.jl", "TurtleSearch.jl", "eom.jl")
julia> using .LazyMatrices, .TurtleSearch, .GoldbergerWiseEom
```

### Running

## (on dev) Requirements List

### Workflows

### Data structure

### Data handling logic
### Functions and their behaviors
## Trouble Shooting