export SaturationFunction, Atan, Tanh, Algebraic

"""
    SaturationFunction

Abstract type for saturation functions. A saturation function is a smooth, bounded, and odd function.
It is used to limit the output of a controller to a certain range.

The subtypes of `SaturationFunction` are defined using two parameters: the derivative at zero, `k_0`, and the limit of the function, `A_max`.
"""
abstract type SaturationFunction <: Function end

"""
    Atan(k_0, A_max)

Arctangent saturation function.

The function is given by
```math
f(x) = A \\tan^{-1}(kx)
```
where the constants ``A`` and ``k`` are chosen to satisfy the conditions on the derivative at zero and the saturation level.
"""
struct Atan <: SaturationFunction
    "Near-zero gain"
    k_0::Real
    "Saturation level"
    A_max::Real
end

"""
    Tanh(k_0, A_max)

Hyperbolic tangent saturation function.

The function is given by
```math
f(x) = A_{\\max} \\tanh(kx)
```
where the constant ``k`` is chosen to satisfy the condition on the derivative at zero.
"""
struct Tanh <: SaturationFunction
    "Near-zero gain"
    k_0::Real
    "Saturation level"
    A_max::Real
end

"""
    Algebraic(k_0, A_max)

Algebraic saturation function.

The function is given by
```math
f(x) = \\frac{kx}{\\sqrt{\\Delta^2 + x^2}}
```
where the constants ``k`` and ``\\Delta`` are chosen to satisfy the conditions on the derivative at zero and the saturation level.
"""
struct Algebraic <: SaturationFunction
    "Near-zero gain"
    k_0::Real
    "Saturation level"
    A_max::Real
end

function (sat::Atan)(x::Real)
    A = sat.A_max * 2 / π
    k = sat.k_0 / A    
    return A * atan(k * x)
end

function (sat::Tanh)(x::Real)
    k = sat.k_0 / sat.A_max
    return sat.A_max * tanh(k * x)
end

function (sat::Algebraic)(x::Real)
    k = sat.A_max
    Δ = sat.A_max / sat.k_0
    return k * x / sqrt(Δ^2 + x^2)
end
