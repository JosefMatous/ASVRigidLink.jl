export VirtualPointSimulationState
export VirtualPointSimulationParameters

struct VirtualPointSimulationState
    "ASV position"
    p_asv::Rn
    "ASV orientation"
    ψ::Real
    "Cable angle"
    θ::Real
    "ASV velocity"
    v_asv::Rn
    "ASV angular velocity"
    r::Real
    "Cable angular velocity"
    θ_dot::Real
    "Path parameter"
    s::Real
    "Guidance internal state"
    ζ::Rn
    "Low-level controller state"
    x_i::Rn
end

function vectorize(state::VirtualPointSimulationState)
    return [state.p_asv; state.ψ; state.θ; state.v_asv; state.r; state.θ_dot; state.s; state.ζ; state.x_i]
end

struct VirtualPointSimulationParameters
    "ASV and cable parameters"
    model::ASVTowingModel
    "Virtual point"
    ε::Real
    "Path function"
    path_fcn::Function
    "Guidance law"
    guidance::VirtualPointGuidance
    "Low-level controller"
    controller::VirtualPointController
    "Ocean current"
    V_c::Rn
end

function closed_loop_ode!(dx::Rn, x::Rn, p::VirtualPointSimulationParameters, t::Real)
    # Unpack parameters
    model = p.model
    ε = p.ε
    path_fcn = p.path_fcn
    guidance = p.guidance
    controller = p.controller
    V_c = p.V_c

    # Unpack state
    x_system = x[1:8]
    s = x[9]
    N_ζ = get_num_states(guidance)
    ζ = x[(1:N_ζ) .+ 9]
    N_i = get_num_states(controller)
    x_i = x[(1:N_i) .+ (9 + N_ζ)]

    # Guidance law
    v_ref, s_dot, ζ_dot = virtual_point_guidance_law(
        x_system,
        s,
        ζ,
        path_fcn,
        guidance,
        ε,
        model,
        V_c
    )
    v_ref_fcn = (_ξ) -> virtual_point_guidance_law(
        _ξ[1:8],
        _ξ[9],
        _ξ[(1:N_ζ) .+ 9],
        path_fcn,
        guidance,
        ε,
        model,
        V_c
    )[1]
    ξ = [x_system; s; ζ]
    ξ_dot = [asv_ode(x_system, zeros(2), V_c, model); s_dot; ζ_dot]
    v_ref_dot = ForwardDiff.jacobian(v_ref_fcn, ξ) * ξ_dot

    # Control law
    u, x_i_dot = virtual_point_control_law(
        x_system,
        x_i,
        v_ref,
        v_ref_dot,
        controller,
        ε,
        model,
        V_c
    )

    # System dynamics
    dx[1:8] = asv_ode(x_system, u, V_c, model)
    dx[9] = s_dot
    dx[(1:N_ζ) .+ 9] = ζ_dot
    dx[(1:N_i) .+ (9 + N_ζ)] = x_i_dot

    nothing
end
