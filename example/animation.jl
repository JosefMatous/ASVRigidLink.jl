using Revise
using ASVRigidLink
using LinearAlgebra

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

# Controller
path = StraightLineParameters([0.0, 0.0], 0.0)
U = 3.0
Δ = 50.0
k_i = 0.1
Δ_a = 10.0
k_y = 0.1
k_ydot = 1.0
speed_mode = SurgeVelocity()
guidance = LOSGuidance(U, Δ, k_y, k_ydot, speed_mode)

k = 0.5 # PID pole
Ki_r = k^3 * I_z
Kp_r = 3*k^2 * I_z
Kd_r = 3*k * I_z
Kp_u = 2*k * m_x
Ki_u = k^2 * m_x
ctrl = LowLevelPID(Kp_u, Ki_u, Kp_r, Ki_r, Kd_r)

V_c = [-0.5, 0.5]

# Initial conditions
p_asv_0 = [0, 25.0]
ψ_0 = 0.0
θ_0 = π
v_asv_0 = [1.0, 0.0]
r_0 = 0.0
θ_dot_0 = 0.0
init_state = SimulationState(p_asv_0, ψ_0, θ_0, v_asv_0, r_0, θ_dot_0, 0.0, zeros(2))

# Simulate
sim_params = SimulationParameters(mdl, path, guidance, ctrl, V_c)
T_stop = 250.0
using DifferentialEquations
prob = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim_params)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# Plot
using Plots, ColorSchemes, Printf
asv_body_length = 10.0
asv_bow_length = 3.0
asv_width = 4.0
dims = SystemDimensions(asv_body_length, asv_bow_length, asv_width, mdl)
scheme = ColorSchemes.Dark2_8
asv_param = (fill=scheme.colors[1], line=:black, linewidth=1, alpha=0.8, label="ASV")
cable_param = (line=scheme.colors[2], linewidth=1.5, label="Cable")
payload_param = (marker=:circle, markersize=6, markercolor=scheme.colors[3], line=:black, label="Payload")
function plot_frame(t::Real)
    title_string = @sprintf("Trajectory (t = %.2f s)", t)
    plt = plot(; title=title_string, xlabel="East [m]", ylabel="North [m]", aspect_ratio=1, legend=:outertopright)
    plot_system!(plt, sol(t), dims; asv_params=asv_param, cable_params=cable_param, payload_params=payload_param)
    return plt
end

fps = 20
spd = 5.0
anim = @animate for t in 0:spd/fps:T_stop
    plot_frame(t)
end
gif(anim, "asv_towing_simulation.gif", fps=fps)
