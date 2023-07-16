# SpectraShooter.jl
Parameter searching for possible spectra(a collection of eigenvalues like masses) of a Boundary Value Problem(BVP) via the shooting method. Originally designed for analysis of the radion stablization problem in extra-dimensions, see [cite]
## Features
- [ ] Visualization of the spectra with respect to the user-defined parameter
- [ ] An easy to use interface for shooting method with a given rigion
- [ ] edge growth algorithm for the searching of the spectrab 
<!-- - optimize the resolusion of the spectra
- tanh activation
- finetuning on interested region
- abstract matrix output
- an api for root searching like `fzero(aiming(config), (xmin, xmax))` -->
## Usage
### Pre-requisites
- `DifferentialEquations.jl`
- 
### Installing
### Running
```julia
fzero(shooter(...), (xmin, xmax))
```
## (on dev) Requirements List

### Workflows
- setup the problem for shooting method
    <!-- - define the bulk equation as `ODEProblem` via `DifferentialEquations.jl` -->
    - define the bulk equations in the following form:
        ```julia
        function bulk!(du, u, p::Params, t)
            m = p.eigval # eigenvalue parameter
            l = p.freeval # free parameter of interest
            p1, p2,... = p.consts # other parameters
            du = ... # bulk equation
            ...
        end
        ```
        <!-- - `PDEProblem` generalization? -->
    - define the boundary conditions at both sides in the following form:
        ```
        bc1(u1, p::Params) = ... return du1
        bc2(u2, p::Params) = ... return du2
        ``` 
        - `bc1` and `bc2` are the boundary conditions at the left and right side of the bulk, respectively
        - `du1` and `du2` are the derivatives of the solution the bulk equation at the left and right side of the bulk, respectively
        - `u1` and `u2` are the values of the solution at the left and right side of the bulk, respectively
        - `p` is the collection of parameters of the bulk equation
        - currently only support linear boundary conditions like $a u + b u' = 0$, and two `bc`s must be independent(i.e. `bc1` does not depend `u2`s and `bc2` does not depends `u1`s)
    - identify the eigenvalue parameter and the free parameter of interest
    - wrap the bulk equation, boundary conditions and parameters into a evaluation function `shooter`:
        ```julia
        function shooter(i::Int, j::Int)
            
        end
        ```
- initialize the seeds of spectra(intervals containing at least one eigenvalue) by scanning from a specific initial parameter value
- let the seeds grow along the parameter space
    - iteration algorithm: wall following
    - termination condition:scanned or out of the region
- extract spectra lines
- find a specific solution with given parameters
- normalize the solution from the given inner product
### Data structure
- [ ] `Spectra` type: include `Params`, `prob` and `bc`s
- [ ] `Params` type: identify the eigenvalue parameter, the free parameter of interest
    ```julia
    struct Params
        consts::Tuple
        loc::CartesianIndex
        logeigvals::StepRangeLen
        logfreevals::StepRangeLen
    end
    Base.getproperty(p::Params, ::Val{:eigval}) = exp10(p.logeigvals[p.loc[1]])
    Base.getproperty(p::Params, ::Val{:freeval}) =  exp10(p.logfreevals[p.loc[2]])
    ```
### Data handling logic
### Functions and their behaviors
- `Spectra(bulk!, p::Params, bc1, bc2)::Spectra` 
- `solve!(s::Spectra)`
- `monitor(s::Spectra, p::Params)::Tuple{Spectra, ODESolution}`
## Trouble Shooting