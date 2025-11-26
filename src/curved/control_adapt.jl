export ASVDamping, DIAGONAL, SWAY_YAW, FULL, detect_damping
export VirtualPointAdaptiveController

@enum ASVDamping DIAGONAL SWAY_YAW FULL

function detect_damping(model::ASVTowingModel, rel_tol::Real=1e-6)
    D = model.D
    abs_tol = rel_tol * tr(D) / 3

    ind_offdiag = (
        CartesianIndex(1,2),
        CartesianIndex(2,1),
        CartesianIndex(1,3),
        CartesianIndex(3,1),
        CartesianIndex(2,3),
        CartesianIndex(3,2),
    )

    ind_not_sway_yaw = ind_offdiag[1:4]
        
    if all([abs(D[i]) < abs_tol for i in ind_offdiag])
        return DIAGONAL
    elseif all([abs(D[i]) < abs_tol for i in ind_not_sway_yaw])
        return SWAY_YAW
    else
        return FULL
    end
end

include("control_adapt_regressor.jl")

struct VirtualPointAdaptiveController <: VirtualPointController
    "Velocity gain"
    k_v::Real
    "Adaptation gain"
    γ::Union{Real, Rmxn}
    "Structure of the damping matrix"
    damping::ASVDamping
end

function get_num_states(c::VirtualPointAdaptiveController)
    dmp = c.damping
    if dmp == DIAGONAL
        N_d = 4
        N_x = 3
    elseif dmp == SWAY_YAW
        N_d = 6
        N_x = 4
    elseif dmp == FULL
        N_d = 10
        N_x = 7
    else
        error("Unknown damping structure")
    end
    return N_d + 2N_x
end

function virtual_point_control_law(
    x::Rn,
    x_i::Rn,
    v_ref::Rn,
    v_ref_dot::Rn,
    controller::VirtualPointAdaptiveController,
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
    dΓ_θ = [-sin(θ), cos(θ)]

    # Virtual point
    L = model.L
    v = v_ASV + ε * L * θ_dot * dΓ_θ
    v_err = v_ref - v

    # Regressor matrix
    Y = regressor_matrix(x, model, controller.damping)

    # Matrices
    M, b = asv_mass_coriolis(x, model)
    Minv = inv(M)
    C = [I(2) zeros(2) ε * L * dΓ_θ]
    B = [Γ_ψ zeros(2); 0 1; 0 0]

    A = C * Minv * B
    
    # Control law (Coriolis cancellation + feedforward)
    k_v = controller.k_v
    v_dot_ff = v_ref_dot  +  k_v * v_err  +  ε * L * θ_dot^2 * Γ_θ
    u = A \ (v_dot_ff  +  C * Minv * (Y * x_i - b))

    # Parameter adaptation law
    γ = controller.γ
    x_i_dot = γ * (C * Minv * Y)' * v_err

    return u, x_i_dot
end
