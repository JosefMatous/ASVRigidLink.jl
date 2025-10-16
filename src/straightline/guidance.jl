export StraightLineParameters
export SpeedMode, SurgeVelocity, SpeedOverGround
export LOSGuidance, guidance_law

"""
    StraightLineParameters

Parameters of a straight-line path.

# Fields
$(FIELDS)
"""
struct StraightLineParameters
    "Starting point"
    p0::Rn
    "Orientation"
    ψ::Real

    function StraightLineParameters(p0::Rn=[0,0], ψ::Real=0)
        @assert length(p0) == 2 "Starting point must be 2D."
        return new(p0, ψ)
    end
end

abstract type SpeedMode end

struct SurgeVelocity <: SpeedMode end

function _surge_reference(x::Rn, U::Real, mode::SurgeVelocity)
    return U
end

struct SpeedOverGround <: SpeedMode
    "Minimum speed"
    u_min::Real
end

function _surge_reference(x::Rn, U::Real, mode::SpeedOverGround)
    # Unpack state
    x_dot = x[5]
    y_dot = x[6]
    ψ = x[3]
    #u = x_dot * cos(ψ) + y_dot * sin(ψ)
    v = -x_dot * sin(ψ) + y_dot * cos(ψ)
    u_sq = U^2 - v^2
    u_ref = u_sq > mode.u_min^2 ? sqrt(u_sq) : mode.u_min
    return u_ref
end

struct LOSGuidance
    "Reference speed"
    U::Real
    "Look-ahead distance"
    Δ::Real
    "Adaptation position gain"
    γ_y::Real
    "Adaptation velocity gain"
    γ_ydot::Real
    "Speed mode"
    speed_mode::SpeedMode

    function LOSGuidance(U::Real, Δ::Real, γ_y::Real, γ_ydot::Real, speed_mode::SpeedMode=SurgeVelocity())
        @assert U > 0 "Reference speed must be positive."
        @assert Δ > 0 "Look-ahead distance must be positive."
        @assert γ_y >= 0 "Adaptation position gain must be non-negative."
        @assert γ_ydot >= 0 "Adaptation velocity gain must be non-negative."
        return new(U, Δ, γ_y, γ_ydot, speed_mode)
    end
end

function guidance_law(x::Rn, y_add::Real, L::Real, path::StraightLineParameters, guidance::LOSGuidance)
    # Unpack state
    p_asv = x[1:2]
    #ψ = x[3]
    θ = x[4]
    v_asv = x[5:6]
    #r = x[7]
    θ_dot = x[8]

    # Unpack path
    p_path = path.p0
    ψ_path = path.ψ
    #Γ_path = [cos(ψ_path), sin(ψ_path)]
    dΓ_path = [-sin(ψ_path), cos(ψ_path)]

    # Serret-Frenet frame
    y = dΓ_path ⋅ (p_asv - p_path)
    p_tow = p_asv + L * [cos(θ), sin(θ)]
    y_tow = dΓ_path ⋅ (p_tow - p_path)

    y_dot = dΓ_path ⋅ v_asv
    v_tow = v_asv + L * θ_dot * [-sin(θ), cos(θ)]
    y_tow_dot = dΓ_path ⋅ v_tow

    # Unpack guidance
    U = guidance.U
    Δ = guidance.Δ
    γ_y = guidance.γ_y
    γ_ydot = guidance.γ_ydot
    speed_mode = guidance.speed_mode

    # LOS guidance law
    ψ_ref = ψ_path - atan((y + y_add) / Δ) # desired yaw angle
    ∂atan = ForwardDiff.derivative(atan, (y + y_add) / Δ)

    # Adaptation law
    y_add_dot = γ_y * y_tow + γ_ydot * y_tow_dot

    # Reference signals
    r_ref = -∂atan * (y_dot + y_add_dot) / Δ # desired yaw rate
    u_ref = _surge_reference(x, U, speed_mode)

    return LowLevelReference(u_ref, ψ_ref, r_ref), y_add_dot
end
