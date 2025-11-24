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
pid = LowLevelPID(Kp_u, Ki_u, Kp_r, Ki_r, Kd_r)
sgn = HighGainSaturation(Tanh, 1e3)
smc = LowLevelST(Kp_u, Ki_u, Kp_r, Ki_r, Kd_r, sgn)

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
sim_params_pid = SimulationParameters(mdl, path, guidance, pid, V_c)
sim_params_smc = SimulationParameters(mdl, path, guidance, smc, V_c)
T_stop = 250.0
using DifferentialEquations
prob_pid = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim_params_pid)
prob_smc = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim_params_smc)
sol_pid = solve(prob_pid, Tsit5(), reltol=1e-6, abstol=1e-6)
sol_smc = solve(prob_smc, Tsit5(), reltol=1e-6, abstol=1e-6)

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

y_pid, vx_pid, θ_pid, θ_dot_pid = get_plot_values(sol_pid)
y_smc, vx_smc, θ_smc, θ_dot_smc = get_plot_values(sol_smc)
using Plots, ColorSchemes

scheme = ColorSchemes.Dark2_8

plt_y = plot(sol_pid.t, y_pid, color=scheme.colors[1], label="PID", xlabel="Time (s)", ylabel="Error (m)", title="Cross-track Error")
plot!(plt_y, sol_smc.t, y_smc, color=scheme.colors[2], label="ST-SMC")

plt_vx = plot(sol_pid.t, vx_pid, color=scheme.colors[1], label="PID", xlabel="Time (s)", ylabel="Velocity (m/s)", title="Payload x-velocity")
plot!(plt_vx, sol_smc.t, vx_smc, color=scheme.colors[2], label="ST-SMC")
plot!(plt_vx, [0, T_stop], [U, U], linestyle=:dash, color=:black, label="Reference Speed")

plt_θ = plot(sol_pid.t, θ_pid, color=scheme.colors[1], label="PID", xlabel="Time (s)", ylabel="Angle (rad)", title="Payload Angle")
plot!(plt_θ, sol_smc.t, θ_smc, color=scheme.colors[2], label="ST-SMC")
plt_θdot = plot(sol_pid.t, θ_dot_pid, color=scheme.colors[1], label="PID", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", title="Payload Angular Velocity")
plot!(plt_θdot, sol_smc.t, θ_dot_smc, color=scheme.colors[2], label="ST-SMC")

plot(plt_y, plt_vx, plt_θ, plt_θdot; layout=(2,2), size=(800,600))
