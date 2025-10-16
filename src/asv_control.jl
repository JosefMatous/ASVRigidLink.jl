export LowLevelController, LowLevelPID, LowLevelST, LowLevelReference, control_law
export get_num_states

abstract type LowLevelController end

struct LowLevelPID <: LowLevelController
    "Surge proportional gain"
    Kp_u::Real
    "Surge integral gain"
    Ki_u::Real
    "Yaw proportional gain"
    Kp_r::Real
    "Yaw integral gain"
    Ki_r::Real
    "Yaw derivative gain"
    Kd_r::Real
end

get_num_states(controller::LowLevelPID) = 2

struct LowLevelST <: LowLevelController
    "Surge proportional gain"
    Kp_u::Real
    "Surge integral gain"
    Ki_u::Real
    "Yaw proportional gain"
    Kp_r::Real
    "Yaw integral gain"
    Ki_r::Real
    "Yaw derivative gain"
    Kd_r::Real
    "Sign function"
    sgn::SignFunction
end

get_num_states(controller::LowLevelST) = 2

struct LowLevelReference
    "Surge reference speed"
    u_ref::Real
    "Yaw reference angle"
    ψ_ref::Real
    "Yaw reference angular velocity"
    r_ref::Real
end

function control_law(x::Rn, x_i::Rn, ref::LowLevelReference, controller::LowLevelPID)
    # Unpack state
    x_dot = x[5]
    y_dot = x[6]
    r = x[7]
    ψ = x[3]
    u = x_dot * cos(ψ) + y_dot * sin(ψ)

    # Unpack integral state
    u_i = x_i[1]
    r_i = x_i[2]

    # Unpack reference
    u_ref = ref.u_ref
    ψ_ref = ref.ψ_ref
    r_ref = ref.r_ref

    # Control law
    ψ_err = mod(ψ_ref - ψ + π, 2π) - π # wrap to [-π, π]
    τ_u = controller.Kp_u * (u_ref - u) + u_i
    τ_r = controller.Kp_r * ψ_err + controller.Kd_r * (r_ref - r) + r_i

    # Integral state dynamics
    u_i_dot = controller.Ki_u * (u_ref - u)
    r_i_dot = controller.Ki_r * ψ_err

    return [τ_u, τ_r], [u_i_dot, r_i_dot]
end

spow(x::Real, a::Real) = sign(x) * abs(x)^a

function control_law(x::Rn, x_i::Rn, ref::LowLevelReference, controller::LowLevelST)
    # Unpack state
    x_dot = x[5]
    y_dot = x[6]
    r = x[7]
    ψ = x[3]
    u = x_dot * cos(ψ) + y_dot * sin(ψ)

    # Unpack integral state
    u_i = x_i[1]
    r_i = x_i[2]

    # Unpack reference
    u_ref = ref.u_ref
    ψ_ref = ref.ψ_ref
    r_ref = ref.r_ref

    # Control law
    τ_u = controller.Kp_u * spow(u_ref - u, 1/2) + u_i
    k1 = controller.Kp_r / controller.Kd_r
    ϕ_r = k1 * spow(ψ_ref - ψ, 2/3) + r_ref - r # sliding variable
    τ_r = controller.Kd_r * spow(ϕ_r, 1/2) + r_i

    # Integral state dynamics
    sgn = controller.sgn
    u_i_dot = controller.Ki_u * sgn(u_ref - u)
    k3 = controller.Ki_r / k1
    r_i_dot = k3 * sgn(ϕ_r)

    return [τ_u, τ_r], [u_i_dot, r_i_dot]
end
