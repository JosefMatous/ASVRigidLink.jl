using Revise
using ASVRigidLink
using LinearAlgebra

# Parameters
m_x = 1.0e3
m_y = m_x
I_z = m_x*25/6
T = [5.0, 2.0, 1.0] # Damping time constants
D = Diagonal([m_x, m_y, I_z] ./ T)
L = 100.0
x_G = -2.0
m_T = 250.0
b_T = 1.0 * m_T
mdl = ASVTowingModel(m_x, m_y, I_z, x_G, D, L, m_T, b_T)

# Path
R_path = 100.0
path_fcn(s) = R_path * [cos(s), sin(s)]
# Controller
U_ref = 3.0
Δ = 50.0
ε = 0.5
guidance = VirtualPointLOS(U_ref, Δ)

k_v = 0.5
ctrl = VirtualPointLinearizingController(k_v)

V_c = [-0.5, 0.5]

# Initial conditions
p_asv_0 = path_fcn(0) + [50.0, 50.0]
ψ_0 = π/2
θ_0 = -π/2
v_asv_0 = [1.0, 0.0]
r_0 = 0.0
θ_dot_0 = 0.0
ζ0 = zeros(get_num_states(guidance)) # guidance states
# Initial controller states
x_i0 = zeros(get_num_states(ctrl)) # controller states
init_state = VirtualPointSimulationState(p_asv_0, ψ_0, θ_0, v_asv_0, r_0, θ_dot_0, 0.0, ζ0, x_i0)

# Simulate
sim_params = VirtualPointSimulationParameters(mdl, ε, path_fcn, guidance, ctrl, V_c)
T_stop = 400.0
using DifferentialEquations
prob = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim_params)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# Visualize
using Plots

x_asv = [u[1] for u in sol.u]
y_asv = [u[2] for u in sol.u]
s_sol = [u[9] for u in sol.u]
x_path = [path_fcn(s)[1] for s in s_sol]
y_path = [path_fcn(s)[2] for s in s_sol]
θ = [u[4] for u in sol.u]
x_tow = x_asv .+ mdl.L * cos.(θ)
y_tow = y_asv .+ mdl.L * sin.(θ)
x_ε = x_asv .+ mdl.L * ε * cos.(θ)
y_ε = y_asv .+ mdl.L * ε * sin.(θ)
plt1 = plot(y_path, x_path, label="Path", color=:black, lw=2, ls=:dash, aspect_ratio=:equal)
plot!(y_asv, x_asv, label="ASV")
plot!(y_tow, x_tow, label="Payload")
plot!(y_ε, x_ε, label="Virtual output")
xlabel!("East [m]")
ylabel!("North [m]")
title!("Trajectory")

# Path-following error
ψ_path = [ASVRigidLink.path_angle(s, path_fcn) for s in s_sol]
cψ = cos.(ψ_path)
sψ = sin.(ψ_path)
e_x = cψ .* (x_ε - x_path) + sψ .* (y_ε - y_path)
e_y = -sψ .* (x_ε - x_path) + cψ .* (y_ε - y_path)
plt2 = plot(sol.t, e_x, label="e_x", lw=2)
plot!(sol.t, e_y, label="e_y", lw=2)
xlabel!("Time [s]")
ylabel!("Path-following error [m]")
title!("Path-following error")

plot(plt1, plt2, layout=(2,1), size=(600,800))
