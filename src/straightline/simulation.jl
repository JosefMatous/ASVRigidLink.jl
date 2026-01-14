export SimulationState, vectorize
export SimulationParameters, closed_loop_ode!

"""
Current state of the straight-line towing simulation.

# Fields
$(FIELDS)
"""
struct SimulationState
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
    "LOS integrator"
    y_add::Real
    "Low-level controller integrator"
    x_i::Rn
end

"""
Convert a `SimulationState` to a vector.
"""
function vectorize(state::SimulationState)
    return [state.p_asv; state.ψ; state.θ; state.v_asv; state.r; state.θ_dot; state.y_add; state.x_i]
end

"""
Create a `SimulationState` from a vector.
"""
function SimulationState(x::Rn)
    return SimulationState(
        x[1:2],
        x[3],
        x[4],
        x[5:6],
        x[7],
        x[8],
        x[9],
        x[10:end]
    )
end
Base.show(io::IO, state::SimulationState) = print(io, "p_asv=$(state.p_asv), ψ=$(state.ψ), θ=$(state.θ), v_asv=$(state.v_asv), r=$(state.r), θ_dot=$(state.θ_dot), y_add=$(state.y_add), x_i=$(state.x_i)")

"""
Simulation parameters for the straight-line towing simulation.

# Fields
$(FIELDS)
"""
struct SimulationParameters
    "ASV and cable parameters"
    model::ASVTowingModel
    "Path"
    path::StraightLineParameters
    "Guidance law"
    guidance::StraightLineGuidance
    "Low-level controller"
    controller::LowLevelController
    "Ocean current"
    V_c::Rn
end

"""
Closed-loop ODE for the straight-line towing simulation.

    closed_loop_ode!(dx, x, p, t)

# Arguments
- `dx::Rn`: Derivative of the state.
- `x::Rn`: Current system state.
- `p::SimulationParameters`: Simulation parameters.
- `t::Real`: Current time.
"""
function closed_loop_ode!(dx::Rn, x::Rn, p::SimulationParameters, ::Real)
    # Unpack parameters
    model = p.model
    path = p.path
    guidance = p.guidance
    controller = p.controller

    # Unpack state
    x_system = x[1:8]
    y_add = x[9]
    x_i = x[10:end]

    # Guidance law
    ref, y_add_dot = guidance_law(x_system, y_add, model.L, path, guidance)
    dx[9] = y_add_dot

    # Control law
    τ, x_i_dot = control_law(x_system, x_i, ref, controller)
    dx[10:end] = x_i_dot

    # System dynamics
    x_system_dot = asv_ode(x_system, τ, p.V_c, model)
    dx[1:8] = x_system_dot
    nothing
end
