using Revise
using ASVRigidLink
using LinearAlgebra
using DifferentialEquations
using ForwardDiff

# Parameters
m_x = 1.5e3
m_y = 1.5e3
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
#speed_mode = SpeedOverGround(1.0)
speed_mode = SurgeVelocity()
guidance = LOSGuidance(U, Δ, k_y, k_ydot, speed_mode)

k = 0.5 # PID pole
Ki_r = k^3 * I_z
Kp_r = 3*k^2 * I_z
Kd_r = 3*k * I_z
Kp_u = 2*k * m_x
Ki_u = k^2 * m_x
ctrl = LowLevelPID(Kp_u, Ki_u, Kp_r, Ki_r, Kd_r)

V_c = [-0.5, 1.0]
sim_params = SimulationParameters(mdl, path, guidance, ctrl, V_c)

function steady_state_analytical(sim::SimulationParameters; N_iter::Integer=10, tol::Real=1e-6)
    model = sim.model
    U = sim.guidance.U
    D = model.D
    V_c = sim.V_c
    b_T = model.b_T

    if isa(sim.guidance.speed_mode, SpeedOverGround)
        U_ss = U

        # Damping force balance
        # ASV damping torque balance
        a_x = D[2,1] * (U_ss - V_c[1]) - (D[2,2] + b_T) * V_c[2]
        a_y = -(D[2,2] + b_T) * (U_ss - V_c[1]) - D[2,1] * V_c[2]
        # a_x * cos(ψ) + a_y * sin(ψ) = 0
        ψ = atan(a_x, -a_y)
    else
        # Yaw angle; let t = tan(ψ)
        f_t(t) = t * (U*sqrt(1 + t^2) - V_c[1]) + V_c[2] # f_t(t) = 0
        t_sol = - V_c[2] / (U - V_c[1])
        for i in 1:N_iter
            t_sol = t_sol - f_t(t_sol) / ForwardDiff.derivative(f_t, t_sol)
            if abs(f_t(t_sol)) < tol
                @debug "Converged in $i iterations"
                break
            end
        end
        ψ = atan(t_sol)
        U_ss = U * sqrt(1 + t_sol^2)

    end
    v_asv = [U_ss, 0.0]
    # Payload damping
    # [-sin(θ), cos(θ)] ⋅ (v_ASV - V_c) = 0
    θ = atan(V_c[2], V_c[1] - U_ss)

    # Control input
    F_u = [(D[1,1] + b_T) * cos(ψ) - D[1,2] * sin(ψ), (D[1,1] + b_T) * sin(ψ) + D[1,2] * cos(ψ)] ⋅ (v_asv - V_c)

    return ψ, θ, v_asv, F_u
end

function steady_state_simulation(sim::SimulationParameters; T_stop::Real=500.0)
    # Initial conditions
    p_asv_0 = [0, 0.0]
    ψ_0 = 0.0
    θ_0 = π
    v_asv_0 = [3.0, 0.0]
    r_0 = 0.0
    θ_dot_0 = 0.0
    init_state = SimulationState(p_asv_0, ψ_0, θ_0, v_asv_0, r_0, θ_dot_0, 0.0, zeros(2))

    # Simulate
    prob = ODEProblem(closed_loop_ode!, vectorize(init_state), (0.0, T_stop), sim)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, save_everystep=false, save_start=false)
    x_ss = sol.u[end]

    ψ = x_ss[3]
    θ = x_ss[4]
    v_asv = x_ss[5:6]
    return ψ, θ, v_asv
end

ψ, θ, v_asv, F_u = steady_state_analytical(sim_params)
x_ASV = [0.0, -L*sin(θ), ψ, θ, v_asv[1], v_asv[2], 0.0, 0.0]
dx_ASV = asv_ode(x_ASV, [F_u, 0.0], V_c, sim_params.model)
