export VirtualPointController, virtual_point_control_law, VirtualPointLinearizingController

"""
Abstract type for virtual point controllers.

Let `p_0` be the position of the ASV, `θ` the angle of the cable, and `ε` a constant between 0 and 1. The virtual point is located at
```math
p_ε = p_0 + ε L [cos(θ), sin(θ)]
```
where `L` is the cable length. The controller computes the control input `τ` for the ASV to track a desired velocity for the virtual point.
"""
abstract type VirtualPointController end

get_num_states(controller::VirtualPointController) = error("get_num_states not implemented for type $(typeof(controller))")

"""
Feedback linearizing controller for the virtual point.

Let `e_v` be the velocity error of the virtual point. The control input is calculated such that the derivative of `e_v` is given by:
```math
ė_v = -K_v e_v
```

# Fields
$(FIELDS)
"""
struct VirtualPointLinearizingController <: VirtualPointController
    "Velocity gain"
    K_v::Real
end

get_num_states(::VirtualPointLinearizingController) = 0

"""
Compute the control input for the ASV to track the desired virtual point velocity.

    τ, x_i_dot = virtual_point_control_law(
        x,
        x_i,
        v_ref,
        v_ref_dot,
        controller,
        ε,
        model,
        V_c
    )

Returns the control input `τ` for the ASV and the derivative of the controller internal state `x_i_dot` (empty vector if the controller has no internal states).

# Arguments
- `x::Rn`: Current system state (cf., `asv_ode`).
- `x_i::Rn`: Current controller state.
- `v_ref::Rn`: Desired virtual point velocity.
- `v_ref_dot::Rn`: Desired virtual point acceleration.
- `controller::VirtualPointController`: Virtual point controller parameters.
- `ε::Real`: Virtual point parameter.
- `model::ASVTowingModel`: ASV and cable model parameters.
- `V_c::Rn`: Ocean current.
"""
function virtual_point_control_law(
    x::Rn,
    x_i::Rn,
    v_ref::Rn,
    v_ref_dot::Rn,
    controller::VirtualPointLinearizingController,
    ε::Real,
    model::ASVTowingModel,
    V_c::Rn
)
    
    # Unpack state
    ψ = x[3]
    θ = x[4]
    v_ASV = x[5:6]
    #ψ_dot = x[7]
    θ_dot = x[8]

    Γ_ψ = [cos(ψ), sin(ψ)]
    Γ_θ = [cos(θ), sin(θ)]
    dΓ_ψ = [-sin(ψ), cos(ψ)]
    dΓ_θ = [-sin(θ), cos(θ)]

    # Virtual point
    L = model.L
    v = v_ASV + ε * L * θ_dot * dΓ_θ
    v_err = v_ref - v

    # Zero-input dynamics
    dx0 = asv_ode(x, zeros(2), V_c, model)
    v_ASV_dot0 = dx0[5:6]
    θ_ddot0 = dx0[8]
    v_dot0 = v_ASV_dot0 + ε * L * (θ_ddot0 * dΓ_θ - θ_dot^2 * Γ_θ)

    # Input matrix
    # Calculate A such that v_dot = v_dot0 + A * u
    #  A = C * inv(M) * B, where
    # C is the Jacobian of v w.r.t. generalized velocities
    C = [I(2) zeros(2) ε * L * dΓ_θ]
    # M is the mass matrix of the system
    m_x = model.m_x
    m_y = model.m_y
    I_z = model.I_z
    x_G = model.x_G
    M_ASV = [m_x 0 0; 0 m_y m_y*x_G; 0 m_y*x_G I_z+m_y*x_G^2]
    J_ASV = [Γ_ψ dΓ_ψ zeros(2); 0 0 1]
    M = zeros(4,4)
    M[1:3,1:3] = J_ASV * M_ASV * J_ASV'
    m_T = model.m_T
    M[1, 1] += m_T
    M[2, 2] += m_T
    M[1:2, 4] = m_T * L * dΓ_θ
    M[4, 1:2] = M[1:2, 4]'
    M[4,4] = m_T * L^2
    # B is the input matrix of the ASV
    B = [Γ_ψ zeros(2); 0 1; 0 0]

    A = C * inv(M) * B

    # Control law
    #  v_dot = v_dot0 + A * u = K_v * v_err + v_ref_dot
    Au = controller.K_v * v_err + v_ref_dot - v_dot0
    u = A \ Au

    T = eltype(x_i)
    return u, zeros(T, length(x_i)) # no integral action
end
