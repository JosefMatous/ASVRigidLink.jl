export ASVTowingModel, asv_ode

struct ASVTowingModel
    "Vehicle inertia in x-axis"
    m_x::Real
    "Vehicle inertia in y-axis"
    m_y::Real
    "Vehicle inertia in yaw"
    I_z::Real
    "Distance from CG to towing point"
    x_G::Real
    "Hydrodynamic damping"
    D::Rmxn
    "Cable length"
    L::Real
    "Payload mass"
    m_T::Real
    "Payload damping"
    b_T::Real
end

function asv_lagrangian(X::Rn, model::ASVTowingModel)
    # Unpack generalized coordinates and velocities
    _, _, ψ, θ, x_dot, y_dot, ψ_dot, θ_dot = X

    # ASV body-fixed velocities
    J_ASV = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
    v_ASV = J_ASV' * [x_dot; y_dot; ψ_dot]

    # Kinetic energy of ASV
    # ASV inertia matrix
    m_x = model.m_x
    m_y = model.m_y
    I_z = model.I_z
    x_G = model.x_G
    M_ASV = [m_x 0 0; 0 m_y m_y*x_G; 0 m_y*x_G I_z]
    T_ASV = 0.5 * v_ASV' * M_ASV * v_ASV

    # Payload velocity
    v_T = [x_dot - θ_dot * model.L * sin(θ), y_dot + θ_dot * model.L * cos(θ)]

    # Kinetic energy of payload
    m_T = model.m_T
    T_T = 0.5 * m_T * (v_T[1]^2 + v_T[2]^2)

    # Total kinetic energy
    T = T_ASV + T_T

    return T    
end

function asv_ode(x::Rn, u::Rn, V_c::Rn, model::ASVTowingModel)
    # Unpack state
    ψ = x[3]
    θ = x[4]
    v_ASV = x[5:6]
    ψ_dot = x[7]
    θ_dot = x[8]

    # ASV body-fixed velocities
    J_ASV = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]
    ν_ASV_rel = J_ASV' * [v_ASV - V_c; ψ_dot] # relative to ocean current

    # Payload velocity
    v_T = v_ASV + θ_dot * model.L * [-sin(θ), cos(θ)]

    # Mass and Coriolis matrices
    L_fcn = _x -> asv_lagrangian(_x, model)
    res = DiffResults.HessianResult(x)
    ForwardDiff.hessian!(res, L_fcn, x)

    H = DiffResults.hessian(res)
    ∇L = DiffResults.jacobian(res)
    M = H[5:8, 5:8]
    b = ∇L[1:4] - H[5:8, 1:4] * x[5:8]

    # Damping forces
    # ASV hydrodynamic damping
    Q_ASV = [-J_ASV * model.D * ν_ASV_rel; 0]
    # Payload damping
    F_T = -model.b_T * (v_T - V_c)
    Jac_T = [I(2); zeros(1,2); model.L * [-sin(θ) cos(θ)]] # Jacobian of payload position w.r.t. generalized coordinates
    Q_T = Jac_T * F_T

    # Control inputs
    τ_u, τ_r = u
    Q_u = [J_ASV * [τ_u, 0, τ_r]; 0]

    # ODEs
    dx = similar(x)
    dx[1:4] = x[5:8]
    dx[5:8] = M \ (b + Q_ASV + Q_T + Q_u)
    return dx
end
