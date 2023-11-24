module GoldbergerWiseEoM
export errBCwithφ, M_IR, k, u, ϕP, yₘ, solveODE, getφ
using DifferentialEquations

# parameters
const yₘ = 1e0π #overall normalization s.t. y_m * Lambda = pi
const u  = 1.0e-1 #   log(ϕT / ϕP)/yₘ, this parameter should be fine-tuned to satisfy k*ym ~ O(50), but this causes instability by inrtoducing such a big hierarchy in numerical computation. Therefore here k*ym is set at O(10)
const ϕP = 1.e-1 # The scalar field value at Plank-brane
const ϕT = exp(-u * yₘ)*ϕP
const k = 37u #pp13 below eq(6.6)
M_IR = exp(-k*yₘ) #IR brane scale;(with M_Pl=1) note that k*ym ~ O(50) is required to get the correct IR scale, but for numerical stability, k*ym is set at O(10)
# l²= 1e-3 #kappa^2 * phiP^2 / 2
# γ² = 1e9 at large gamma limit

#static profile
@inline ϕ0(  y::Number)                                    = ϕP * (ϕT/ϕP)^(y/yₘ) #ϕ0' = -u ϕ0
@inline A(   y::Number, l²::Number, k::Number, γ²::Number) = k * y + l²/6 * (ϕT/ϕP)^(2y/yₘ)
@inline A′(  y::Number, l²::Number, k::Number, γ²::Number) = k     + l²/6 * (ϕT/ϕP)^(2y/yₘ)* (-2u)
@inline A′′( y::Number, l²::Number, k::Number, γ²::Number) =         l²/6 * (ϕT/ϕP)^(2y/yₘ)* 4u^2
@inline V′(  y::Number, l²::Number, k::Number, γ²::Number) = u*(4k + u)*ϕ0(y) - 2/3*u^2 * 2l²/ϕP^2*ϕ0(y)^3 #κ^2 = 2l²/ϕP^2
@inline V′′( y::Number, l²::Number, k::Number, γ²::Number) = u*(4k + u -2u*2l²*(ϕT/ϕP)^(2y/yₘ))
@inline λP′( φ::Number, l²::Number, k::Number, γ²::Number) = -2u * ϕP  + 2γ² * φ
# @inline λP′( φ::Number, l²::Number, k::Number, γ²::Number) = -2u * ϕP * φ + 2γ² * φ
@inline λT′( φ::Number, l²::Number, k::Number, γ²::Number) =  2u * ϕT  + 2γ² * φ
# @inline λT′( φ::Number, l²::Number, k::Number, γ²::Number) =  2u * ϕT * φ + 2γ² * φ
@inline λP′′(φ::Number, l²::Number, k::Number, γ²::Number) = 2γ²
@inline λT′′(φ::Number, l²::Number, k::Number, γ²::Number) = 2γ²

# EoM of φ, obtained from the zero solution of the Einstein eq #(3.12)
function getφ(sol_of_F::ODESolution, params) 
    mF, l², γ²= params
    F(y::Number)  = sol_of_F(y)[2]
    F′(y::Number) = sol_of_F(y)[1]
    φ(y::Number)  = 3ϕP^2/(2u*l²*ϕ0(y)) * (F′(y) - 2A′(y, l², k, γ²)*F(y))  #(3.12)
    return φ
end

#     # φ(y::Number)  = -6ϕP^2/(u*l²*ϕ0(y)) * (F′(y) - 2A′(y, l², k, γ²)*F(y))  #(3.12)


#BCs
# at infinite gamma limit:
@inline dFP_frozen(φ, F, l², k, γ²) = 2A′(0 , l², k, γ²)*F
@inline dFT_frozen(φ, F, l², k, γ²) = 2A′(yₘ, l², k, γ²)*F
# non-perturbative gamma:
@inline dφP(φ, F, l², k, γ²) = 0.5λP′′(φ, l², k, γ²) * φ - 2u * ϕP * F  #(3.14)
    # @inline dφP(φ, F, l², k, γ²) = 0.5( λP′′(φ, l², k, γ²) * φ + 2λP′(0, l², k, γ²) * F ) #(3.14)(wrong sign)
    # @inline dφP(φ, F, l², k, γ²) = 0.5( λP′′(φ, l², k, γ²) * φ + 2λP′(φ, l², k, γ²) * F ) #(3.14)
@inline dφT(φ, F, l², k, γ²) =-0.5λT′′(φ, l², k, γ²) * φ - 2u * ϕT * F  #(3.14)
    # @inline dφT(φ, F, l², k, γ²) =-0.5( λT′′(φ, l², k, γ²) * φ + 2λT′(0, l², k, γ²) * F ) #(3.14)(wrong sign)
    # @inline dφT(φ, F, l², k, γ²) =-0.5( λT′′(φ, l², k, γ²) * φ + 2λT′(φ, l², k, γ²) * F ) #(3.14)
@inline dFP(φ, F, l², k, γ²) = 2A′(0 , l², k, γ²)*F - u*l²*ϕ0(0)/ϕP^2/6 * φ
@inline dFT(φ, F, l², k, γ²) = 2A′(yₘ, l², k, γ²)*F - u*l²*ϕ0(yₘ)/ϕP^2/6 * φ

# EoM of F(perturbation on redshift factor A in metric)
function radionSpectrum_secondOrder!(ddf, df, f, params, y)
    F  = f[1]
    F′ = df[1]
    mF, l², γ²= params
    dF′ = 2A′(y, l², k, γ²)*F′ + 4A′′(y, l², k, γ²)*F - 2u*F′ + 4u*A′(y, l², k, γ²)*F - mF^2 * exp(2A(y, l², k, γ²))*F #(3.17)
    ddf[1]=dF′
end

W(ϕ, κ, k, u) = 6k/κ^2  - u*ϕ^2
W′(ϕ, κ, k, u) =        - 2u*ϕ
V(W, ϕ, κ, k, u) = 1/8 * (W′(ϕ, κ, k, u))^2 - κ^2 / 6 * W(ϕ, κ, k, u)^2
    
# function errBCwithφ(FP, params; φP = 0)
#     mF, l², γ²= params
#     yspan = (0.0,yₘ)
    
#     F′P= dFP(φP, FP, l², k, γ²)
#     prob = SecondOrderODEProblem(radionSpectrum_secondOrder!,[F′P], [FP],yspan, params)
#     Fsol = solve(prob, ImplicitEuler())
#     # Fsol = solve(prob, Tsit5())
#     F′T, FT = Fsol(yₘ)
#     φ = getφ(Fsol, params)
#     φT = φ(yₘ)
#     # φ(0)=φP from solving F(y) with BC
#     ys = range(yₘ*(1-1e-6) , yₘ, 2)
#     φ′T = (diff(φ.(ys))/diff(ys))[1]
#     # dφT(φT, FT, l², k, γ²)#(3.14)
#     ΔφT =  (dφT(φT, FT, l², k, γ²) + φ′T)/sqrt(   ( λT′′(0, l², k, γ²)/2)^2+  λT′(0, l², k, γ² )^2+1)  
#     return ΔφT #distance from (F'T, FT) to the TeV brane BC line in phase diagram
# end
function solveODE(FP, φP, params)
    mF, l², γ²= params
    yspan = (0.0,yₘ)
    F′P= dFP(φP, FP, l², k, γ²)
    prob = SecondOrderODEProblem(radionSpectrum_secondOrder!,[F′P], [FP],yspan, params)
    # return solve(prob, Tsit5())
    return solve(prob, ImplicitEuler())
    # return solve(prob, Rosenbrock23())
end

function calculateΔφT(Fsol, params)
    mF, l², γ²= params
    F′T, FT = Fsol(yₘ)
    φ = getφ(Fsol, params)
    φT = φ(yₘ)
    ys = range(yₘ*(1-1e-6) , yₘ, 2)
    φ′T = (diff(φ.(ys))/diff(ys))[1]
    # return (dφT(φT, FT, l², k, γ²) + φ′T)/sqrt((λT′′(0, l², k, γ²)/2)^2+ λT′(φT, l², k, γ² )^2+1)
    return (dφT(φT, FT, l², k, γ²) + φ′T)/sqrt((λT′′(0, l², k, γ²)/2)^2+ λT′(0, l², k, γ² )^2+1)
end

function errBCwithφ(FP, params; φP = 0)
    Fsol = solveODE(FP, φP, params)
    return calculateΔφT(Fsol, params)
end
end