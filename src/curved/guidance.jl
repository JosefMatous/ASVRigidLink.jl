export VirtualPointGuidance, VirtualPointLOS, virtual_point_guidance_law

"""
Abstract type for virtual point guidance laws.

Let `p_0` be the position of the ASV, `θ` the angle of the cable, and `ε` a constant between 0 and 1. The virtual point is located at
```math
p_ε = p_0 + ε L [cos(θ), sin(θ)]
```
where `L` is the cable length. The guidance law computes a desired velocity for the virtual point, as well as the path parameter rate.
"""
abstract type VirtualPointGuidance end

get_num_states(guidance::VirtualPointGuidance) = error("get_num_states not implemented for type $(typeof(guidance))")

"""
LOS guidance law for the virtual point.

# Fields
$(FIELDS)
"""
struct VirtualPointLOS <: VirtualPointGuidance
    "Reference speed"
    U::Real
    "Look-ahead distance"
    Δ::Real
end

get_num_states(guidance::VirtualPointLOS) = 0

"""
Compute the desired velocity for the virtual point and the path parameter rate.

    v_ref, s_dot, ζ_dot = virtual_point_guidance_law(
        x,
        s,
        ζ,
        path_fcn,
        guidance,
        ε,
        model,
        V_c
    )

Returns the desired velocity `v_ref` for the virtual point, the path parameter rate `s_dot`, and the derivative of the guidance internal state `ζ_dot`.

# Arguments
- `x::Rn`: Current system state (cf., `asv_ode`).
- `s::Real`: Current path parameter.
- `ζ::Rn`: Current guidance internal state (empty vector if the guidance law has no internal states).
- `path_fcn::Function`: Path function mapping path parameter to position.
- `guidance::VirtualPointGuidance`: Virtual point guidance parameters.
- `ε::Real`: Virtual point parameter.
- `model::ASVTowingModel`: ASV and cable model parameters.
- `V_c::Rn`: Ocean current (typically unused).
"""
function virtual_point_guidance_law(
    x::Rn,
    s::Real,
    ζ::Rn,
    path_fcn::Function,
    guidance::VirtualPointLOS,
    ε::Real,
    model::ASVTowingModel,
    V_c::Rn
)

    # Unpack state
    p_asv = x[1:2]
    #ψ = x[3]
    θ = x[4]
    #v_asv = x[5:6]
    #r = x[7]
    #θ_dot = x[8]

    # Virtual output
    L = model.L
    Γ_cable = [cos(θ), sin(θ)]
    #dΓ_cable = [-sin(θ), cos(θ)]
    p_ε = p_asv + ε * L * Γ_cable
    #v_ε = v_asv + ε * L * dΓ_cable * θ_dot

    # LOS guidance
    Δ = guidance.Δ
    p_path = path_fcn(s)
    ψ_path = path_angle(s, path_fcn)
    R_path = [cos(ψ_path) -sin(ψ_path);
              sin(ψ_path)  cos(ψ_path)]
    ∇_path = ForwardDiff.derivative(path_fcn, s)
    path_error = R_path' * (p_ε - p_path)
    e_x, e_y = path_error
    D_x = sqrt(e_x^2 + Δ^2)
    D_y = sqrt(e_y^2 + Δ^2)
    v_LOS = guidance.U * R_path * [Δ, -e_y] ./ D_y
    s_dot = guidance.U * (Δ / D_y + e_x / D_x) / norm(∇_path)

    T = eltype(ζ)

    return v_LOS, s_dot, zeros(T, length(ζ))
end
