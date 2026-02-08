module Controller

using LinearAlgebra
using BlockDiagonals
using ForwardDiff
using ControlSystems
import Convex as cvx
using ECOS
import DifferentialEquations: solve as DEsolve
import DifferentialEquations: ODEProblem, Tsit5, PeriodicCallback, CallbackSet

using ..QuadSys
using ..QuaternionUtils
using ..RefTrajLib
using ..SimulateCoverge

# Compute Gain Matrix for LQR
A̅, B̅ = QuadSys.get_linearized_model()
Q, R = QuadSys.cost_matrices()
const K = dlqr(A̅, B̅, Q, R)

uhover = QuadSys.uhover
xeq = QuadSys.x0
req = xeq[1:3]
qeq= xeq[4:7]
q̅eq = QuaternionUtils.qtorp(qeq)
veq = xeq[8:10]
ωeq = xeq[11:13]
x̅eq = [req; q̅eq; veq ; ωeq]

Nx = QuadSys.Nx
Nx̅ = QuadSys.Nx̅
Nu = QuadSys.Nu
umin = QuadSys.umin
umax = QuadSys.umax

struct NomTrajectory{F}
    ts::Vector{F}
    xs::Vector{Vector{F}}
    us::Vector{Vector{F}}
end

NomTrajectory() = NomTrajectory(Float64[],Vector{Float64}[],Vector{Float64}[])

function get_nominal_mpc_π(x0, DI_ergodic_traj; t0 = 0.0, N_mpc = 40, T_H_sec = 2.0, TRAJ = "NM")
    r0 = x0[1:3]
    q0 = x0[4:7]
    q̅0 = QuaternionUtils.qtorp(q0)
    v0 = x0[8:10]
    ω0 = x0[11:13]
    x̅0 = [r0; q̅0; v0; ω0]

    dt = T_H_sec / N_mpc

    if TRAJ == "NM"
        X_ref, U_ref = SimulateCoverge.generate_reference_trajectory(t0, N_mpc, dt, DI_ergodic_traj)
    else
        X_ref, U_ref = RefTrajLib.get_desired_trajectory_point(t0, N_mpc, dt)
    end

    X_traj = [cvx.Variable(Nx̅) for i=1:N_mpc]
    U_traj = [cvx.Variable(Nu) for i=1:(N_mpc-1)];

    cost = 0.0
    for i = 1:N_mpc
        cost += 0.5*cvx.quadform(X_traj[i] - X_ref[i], Q)
    end
    for i = 1:(N_mpc - 1)
        cost += 0.5*cvx.quadform(U_traj[i] - U_ref[i], R)
    end

    prob = cvx.minimize(cost)
    prob.constraints += X_traj[1] == x̅0

    xrk4 = QuadSys.quad_dynamics_rk4(xeq, uhover, dt)
    x̅rk4 = [xrk4[1:3]; QuaternionUtils.qtorp(xrk4[4:7]); xrk4[8:10]; xrk4[11:13]]
    for i = 1:(N_mpc - 1)
        prob.constraints += X_traj[i+1] == x̅rk4 + A̅*(X_traj[i] - x̅eq) + B̅*(U_traj[i] - uhover)
    end

    for i = 1:(N_mpc - 1)
        prob.constraints += U_traj[i] <= umax
        prob.constraints += U_traj[i] >= umin
    end
    cvx.solve!(prob, ECOS.Optimizer; silent_solver = true)

    X_traj_sol = cvx.evaluate.(X_traj)
    U_traj_sol = cvx.evaluate.(U_traj)

    return NomTrajectory(
        [t0 + (i-1) * dt for i=1:N_mpc],
        X_traj_sol,
        U_traj_sol
    )
end

function interpolate_traj(t, traj::NomTrajectory)
    @assert t >= traj.ts[1]

    if t >= traj.ts[end]
        h = t - traj.ts[end]
        A̅, B̅ = QuadSys.get_linearized_model(h)
        di_state_ref = A̅ * traj.xs[end] + B̅ * traj.us[end]
        di_u_ref = traj.us[end]
    else
        N = length(traj.ts)
        ind = findfirst(traj.ts .> t) - 1
        @assert !isnothing(ind)
        h = t - traj.ts[ind]
        A̅, B̅ = QuadSys.get_linearized_model(h)
        di_state_ref = A̅ * traj.xs[ind] + B̅ * traj.us[ind]
        di_u_ref = traj.us[ind]
    end

    return di_state_ref, di_u_ref
end

function LQR_π_tracker(t, x, traj::NomTrajectory)
    x_ref, u_ref = interpolate_traj(t, traj)

    r_ref = x_ref[1:3]
    p_ref = x_ref[4:6]
    q_ref = QuaternionUtils.rptoq(p_ref)
    v_ref = x_ref[7:9]
    ω_ref = x_ref[10:12]

    r = x[1:3]
    q = x[4:7]
    v = x[8:10]
    ω = x[11:13]
    ϕ = QuaternionUtils.qtorp(QuaternionUtils.L(q_ref)'*q)
    Δx̅ = [r-r_ref; ϕ; v-v_ref; ω-ω_ref]
    u = u_ref - K*Δx̅
    return u
end

end
