using Revise
using ASVRigidLink
using LinearAlgebra
import ForwardDiff

# Parameters
m_x = 1.5e3
m_y = 1.2e3
I_z = m_x*25/6
T = [5.0, 2.0, 1.0] # Damping time constants
D = Diagonal([m_x, m_y, I_z] ./ T)
L = 100.0
x_G = -2.0
m_T = 250.0
b_T = 1.0 * m_T
mdl = ASVTowingModel(m_x, m_y, I_z, x_G, D, L, m_T, b_T)

V_c = [-0.5, 0.5]

# Controller
path = StraightLineParameters([0.0, 0.0], 0.0)
k = 0.5 # PID pole
Ki_r = k^3 * I_z
Kp_r = 3*k^2 * I_z
Kd_r = 3*k * I_z
Kp_u = 2*k * m_x
Ki_u = k^2 * m_x
pid = LowLevelPID(Kp_u, Ki_u, Kp_r, Ki_r, Kd_r)

# Guidance algorithms
struct GenericLOS <: StraightLineGuidance
    los_fcn::Function
    parameters::Any
    speed_mode::SpeedMode
end

function ASVRigidLink.guidance_law(x::Rn, y_add::Real, L::Real, path::StraightLineParameters, guidance::GenericLOS)
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

    ψ_add, r_ref, y_add_dot = guidance.los_fcn(
        y, y_dot, y_tow, y_tow_dot, y_add, guidance.parameters
    )

    U = reference_speed(guidance.parameters)
    u_ref = ASVRigidLink._surge_reference(x, U, guidance.speed_mode)

    return LowLevelReference(u_ref, ψ_path + ψ_add, r_ref), y_add_dot
end

# LOS + Linear integrator
struct LOSLinearParams
    U::Real
    Δ::Real
    k_y::Real
    k_ydot::Real
    on_payload::Bool
end

reference_speed(params::LOSLinearParams) = params.U

function los_fcn(y_asv::Real, y_asv_dot::Real, y_tow::Real, y_tow_dot::Real, y_add::Real, params::LOSLinearParams)
    Δ = params.Δ
    k_y = params.k_y
    k_ydot = params.k_ydot
    on_payload = params.on_payload

    y = on_payload ? y_tow : y_asv
    y_dot = on_payload ? y_tow_dot : y_asv_dot

    ψ_add = -atan((y + y_add) / Δ)
    ∂atan = ForwardDiff.derivative(atan, (y + y_add) / Δ)
    
    y_add_dot = k_y * y_tow + k_ydot * y_tow_dot

    r_ref = -∂atan * (y_dot + y_add_dot) / Δ

    return ψ_add, r_ref, y_add_dot
end

# ILOS
struct ILOSParams
    U::Real
    Δ::Real
    k_y::Real
    k_ydot::Real
    on_payload::Bool
end

reference_speed(params::ILOSParams) = params.U

function los_fcn(y_asv::Real, y_asv_dot::Real, y_tow::Real, y_tow_dot::Real, y_add::Real, params::ILOSParams)
    Δ = params.Δ
    k_y = params.k_y
    k_ydot = params.k_ydot
    on_payload = params.on_payload

    y = on_payload ? y_tow : y_asv
    y_dot = on_payload ? y_tow_dot : y_asv_dot

    ψ_add = -atan((y + k_y * y_add) / Δ)
    ∂atan = ForwardDiff.derivative(atan, (y + k_y * y_add) / Δ)

    y_add_dot = Δ * y_tow / (Δ^2 + (y_tow + k_y * y_add)^2) + k_ydot * y_tow_dot
    r_ref = -∂atan * (y_dot + y_add_dot) / Δ

    return ψ_add, r_ref, y_add_dot
end

# Adaptive ILOS
struct AILOSParams
    U::Real
    Δ::Real
    k_y::Real
    k_ydot::Real
    on_payload::Bool
end

reference_speed(params::AILOSParams) = params.U

function los_fcn(y_asv::Real, y_asv_dot::Real, y_tow::Real, y_tow_dot::Real, β::Real, params::AILOSParams)
    U = params.U
    Δ = params.Δ
    k_y = params.k_y
    k_ydot = params.k_ydot
    on_payload = params.on_payload

    y = on_payload ? y_tow : y_asv
    y_dot = on_payload ? y_tow_dot : y_asv_dot

    ψ_add = -atan(y / Δ + β)
    ∂atan = ForwardDiff.derivative(atan, y / Δ + β)

    β_dot = k_y * U * Δ * y_tow / (Δ^2 + (y_tow + β * Δ)^2) + k_ydot * y_tow_dot
    r_ref = -∂atan * (y_dot + Δ * β_dot) / Δ

    return ψ_add, r_ref, β_dot
end

# Adaptive LOS
struct ALOSParams
    U::Real
    Δ::Real
    k_y::Real
    k_ydot::Real
    on_payload::Bool
end

reference_speed(params::ALOSParams) = params.U

function los_fcn(y_asv::Real, y_asv_dot::Real, y_tow::Real, y_tow_dot::Real, β::Real, params::ALOSParams)
    Δ = params.Δ
    k_y = params.k_y
    k_ydot = params.k_ydot
    on_payload = params.on_payload

    y = on_payload ? y_tow : y_asv
    y_dot = on_payload ? y_tow_dot : y_asv_dot

    ψ_add = -atan(y / Δ) - β
    ∂atan = ForwardDiff.derivative(atan, y / Δ)

    β_dot = k_y * Δ / sqrt(Δ^2 + y_tow^2) * y_tow + k_ydot * y_tow_dot
    r_ref = -∂atan * (y_dot + Δ * β_dot) / Δ

    return ψ_add, r_ref, β_dot
end

guidance = GenericLOS(los_fcn, LOSLinearParams(3.0, 50.0, 0.1, 1.0, false), SurgeVelocity())
# Initial conditions
p_asv_0 = [0, 25.0]
ψ_0 = 0.0
θ_0 = π
v_asv_0 = [1.0, 0.0]
r_0 = 0.0
θ_dot_0 = 0.0
init_state = SimulationState(p_asv_0, ψ_0, θ_0, v_asv_0, r_0, θ_dot_0, 0.0, zeros(2))

# Convergence analysis
using DifferentialEquations
function convergence_time(los_params; T_stop=1e4, y_tol=0.1, vy_tol=0.1, reltol=1e-6, abstol=1e-6, solver=Tsit5())
    guidance = GenericLOS(los_fcn, los_params, SurgeVelocity())
    sim_params = SimulationParameters(mdl, path, guidance, pid, V_c)
    prob = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim_params)
    sol = solve(prob, solver, reltol=reltol, abstol=abstol)

    y0 = [state[2] for state in sol.u]
    θ = [state[4] for state in sol.u]
    y = y0 + mdl.L * sin.(θ)

    vy0 = [state[6] for state in sol.u]
    θ_dot = [state[8] for state in sol.u]
    vy = vy0 + mdl.L * θ_dot .* cos.(θ)

    ind_y = findlast(abs.(y) .> y_tol)
    if isnothing(ind_y)
        ind_y = 0
    end
    ind_vy = findlast(abs.(vy) .> vy_tol)
    if isnothing(ind_vy)
        ind_vy = 0
    end
    ind = max(ind_y, ind_vy) + 1
    if ind > length(sol.t)
        return T_stop
    else
        return sol.t[ind]
    end
end

# Optimization of LOS parameters
using NLopt
function optimize_los_params(p0::Vector, param_generator::Function; algorithm=:LN_BOBYQA, lower_bounds=nothing, upper_bounds=nothing, T_stop=1e3, maxevals=1000, args::Tuple=())
    opt = Opt(algorithm, length(p0))
    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds
    opt.min_objective = (p, _) -> begin
        los_params = param_generator(p, args...)
        t_conv = convergence_time(los_params; T_stop=T_stop)
        return t_conv
    end
    opt.maxeval = maxevals
    (minf, minx, ret) = optimize(opt, p0)
    los_params_opt = param_generator(minx, args...)
    return los_params_opt, minf, ret
end

U = 3.0
function generate_LOSParams(p::Vector, ptype::Type, on_payload::Bool)
    return ptype(U, p[1], p[2], p[3], on_payload)
end

# Optimize LOS parameters
# Linear LOS
p0 = [50.0, 0.1, 1.0]
lower_bounds = [10.0, 0.0, 0.0]
upper_bounds = [100.0, 10.0, 10.0]
popt, minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(LOSLinearParams, false), algorithm=:GN_ESCH)
p0 = [popt.Δ, popt.k_y, popt.k_ydot]
linlos_params_opt, linlos_minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(LOSLinearParams, false), algorithm=:LN_BOBYQA)

# ILOS
p0 = [50.0, 1.0, 1.0]
popt, minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(ILOSParams, false), algorithm=:GN_ESCH)
p0 = [popt.Δ, popt.k_y, popt.k_ydot]
ilos_params_opt, ilos_minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(ILOSParams, false), algorithm=:LN_BOBYQA)

# Adaptive ILOS
p0 = [50.0, 1.0, 1.0]
popt, minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(AILOSParams, false), algorithm=:GN_ESCH)
p0 = [popt.Δ, popt.k_y, popt.k_ydot]
ailos_params_opt, ailos_minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(AILOSParams, false), algorithm=:LN_BOBYQA)

# Adaptive LOS
p0 = [50.0, 1.0, 1.0]
popt, minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(ALOSParams, false), algorithm=:GN_ESCH)
p0 = [popt.Δ, popt.k_y, popt.k_ydot]
alos_params_opt, alos_minf, ret = optimize_los_params(p0, generate_LOSParams; lower_bounds=lower_bounds, upper_bounds=upper_bounds, args=(ALOSParams, false), algorithm=:LN_BOBYQA)

# Save optimized parameters
fd = open(joinpath(@__DIR__, "optimized_los_params.jl"), "w")
println(fd, "# Optimized LOS parameters")
println(fd, "linlos_params_opt = $(linlos_params_opt)")
println(fd, "ilos_params_opt = $(ilos_params_opt)")
println(fd, "ailos_params_opt = $(ailos_params_opt)")
println(fd, "alos_params_opt = $(alos_params_opt)")
close(fd)

# Errors
function get_plot_values(sol)
    y0 = [state[2] for state in sol.u]
    θ = [state[4] for state in sol.u]
    y = y0 + mdl.L * sin.(θ)

    vx0 = [state[5] for state in sol.u]
    θ_dot = [state[8] for state in sol.u]
    vx = vx0 - mdl.L * θ_dot .* sin.(θ)

    return y, vx, θ, θ_dot
end

using Plots, ColorSchemes

function simulate_and_plot(los_params, p1, p2, p3, p4; label="", color=:blue, T_stop=200.0)
    guidance = GenericLOS(los_fcn, los_params, SurgeVelocity())
    sim_params = SimulationParameters(mdl, path, guidance, pid, V_c)
    prob = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim_params)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

    y, vx, θ, θ_dot = get_plot_values(sol)

    t = sol.t

    plot!(p1, t, y, label=label, color=color)
    plot!(p2, t, vx, label=label, color=color)
    plot!(p3, t, θ, label=label, color=color)
    plot!(p4, t, θ_dot, label=label, color=color)
    return p1, p2, p3, p4
end

scheme = ColorSchemes.Dark2_8
p1 = plot(;xlabel="Time (s)", ylabel="Lateral position (m)", title="Lateral position")
p2 = plot(;xlabel="Time (s)", ylabel="Lateral velocity (m/s)", title="Lateral velocity")
p3 = plot(;xlabel="Time (s)", ylabel="Tow angle (rad)", title="Tow angle")
p4 = plot(;xlabel="Time (s)", ylabel="Tow angular velocity (rad/s)", title="Tow angular velocity")

simulate_and_plot(linlos_params_opt, p1, p2, p3, p4; label="Linear LOS", color=scheme.colors[1])
simulate_and_plot(ilos_params_opt, p1, p2, p3, p4; label="ILOS", color=scheme.colors[2])
simulate_and_plot(ailos_params_opt, p1, p2, p3, p4; label="Adaptive ILOS", color=scheme.colors[3])
simulate_and_plot(alos_params_opt, p1, p2, p3, p4; label="Adaptive LOS", color=scheme.colors[4])
plot(p1, p2, p3, p4; layout=(2,2), size=(800, 600))
