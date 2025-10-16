export SignFunction, Sign, HighGainSaturation

"""
    SignFunction

Abstract base type for sign-like functions used in sliding mode control.
"""
abstract type SignFunction <: Function end

"""
    Sign <: SignFunction

Standard sign function implementation.
"""
struct Sign <: SignFunction
end

function (s::Sign)(x::Real)
    return sign(x)
end

"""
    HighGainSaturation(sat::Type{<:SaturationFunction}, k::Real) <: SignFunction

Smooth approximation of the sign function using a high-gain saturation function.
"""
struct HighGainSaturation <: SignFunction
    "Saturation function"
    sat::Type{<:SaturationFunction}
    "Gain"
    k::Real
end

function (s::HighGainSaturation)(x::Real)
    fun = s.sat(s.k, 1)
    return fun(x)
end
