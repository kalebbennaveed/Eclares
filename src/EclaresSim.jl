module EclaresSim

using DifferentialEquations
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
using ..Controller
using ..OceanographicData
using ..SimulateCoverge

NomTrajectory = Controller.NomTrajectory
DI_erg_Traj = SimulateCoverge.DI_erg_Traj

function reached_charging_station(comb_state, flight_mode)
    robot_pos = comb_state[1:3]
    if (norm(robot_pos) <= 0.115) && (flight_mode == "Landing")
        return true
    else
        return false
    end
end

function recompute_u_track_cand!(integrator)
    t = integrator.t
    comb_state = integrator.u
    quad_state = comb_state[1:13]
    u = Controller.LQR_π_tracker(t, quad_state, integrator.p[2])
    integrator.p[1] = u
end

function simulate_u_track_cand!(integrator)
    t = integrator.t
    comb_state = integrator.u
    quad_state = comb_state[1:13]
    u = Controller.LQR_π_tracker(t, quad_state, integrator.p[2])
    integrator.p[1] = u
end

function simulate_tracking_candidate(t0, initial_state, cand_traj::NomTrajectory)
    tspan = (t0, cand_traj.ts[end])
    params = [[0,0.,0,0], cand_traj]
    prob = ODEProblem(QuadSys.sim_closed_loop_comb!, initial_state, tspan, params)
    cb_u = PeriodicCallback(simulate_u_track_cand!, 0.05; initial_affect=true)
    sol = solve(prob, Tsit5(), callback=cb_u)
    return sol
end

function get_committed_traj!(integrator)
    if integrator.p[6] != "Mission"
        return
    end

    t_k = integrator.t
    state_k = integrator.u
    state_quad = state_k[1:13]
    traj_nom_k = integrator.p[3]

    t_kN = traj_nom_k.ts[end]
    state_kN_red = traj_nom_k.xs[end]
    state_kN = [state_kN_red[1:3]; QuaternionUtils.rptoq(state_kN_red[4:6]); state_kN_red[7:9]; state_kN_red[10:12]]

    comp_time_b2b = @elapsed traj_b2b_kn = Controller.get_nominal_mpc_π(state_kN, integrator.p[7]; t0 = t_kN, T_H_sec = 10.0, N_mpc = 200, TRAJ = "B2B")
    push!(integrator.p[5], comp_time_b2b)

    traj_cand = deepcopy(traj_nom_k)
    pop!(traj_cand.xs)
    pop!(traj_cand.ts)

    for i = 1:length(traj_b2b_kn.ts)
        push!(traj_cand.ts, traj_b2b_kn.ts[i])
        push!(traj_cand.xs, traj_b2b_kn.xs[i])
    end

    for i = 1:length(traj_b2b_kn.us)
        push!(traj_cand.us, traj_b2b_kn.us[i])
    end

    traj_cand_sol_k = simulate_tracking_candidate(t_k, state_k, traj_cand)
    combined_state_sol = traj_cand_sol_k.u
    SOC_at_END = combined_state_sol[end][14]

    if SOC_at_END > 0.
        integrator.p[2] = traj_cand
    else
        integrator.p[6] = "Landing"
        new_current_map = SimulateCoverge.update_current_map(t_k, integrator.p[7], integrator.p[8], integrator.p[9])
        if new_current_map !== nothing
            integrator.p[9] = new_current_map
        end
    end
end

function recompute_committed_timed!(integrator)
    if integrator.p[6] != "Mission"
        return
    end
    t_k = integrator.t
    state_k = integrator.u
    state_quad = state_k[1:13]
    comp_time_nom = @elapsed traj_nom_k = Controller.get_nominal_mpc_π(state_quad, integrator.p[7]; t0 = t_k)
    integrator.p[3] = traj_nom_k
    get_committed_traj!(integrator)
    return
end

function recompute_ergodic_timed!(integrator)
    if integrator.p[6] != "Mission"
        return
    end

    if integrator.p[8] > 0
        DI_ergodic_traj, new_ϕ_map = SimulateCoverge.generate_DI_ergodic_trajectory(integrator.p[9], integrator.p[7])
        integrator.p[8] += 1
        integrator.p[7] = DI_ergodic_traj
        integrator.p[9] = new_ϕ_map
    else
        DI_ergodic_traj, new_ϕ_map = SimulateCoverge.generate_initial_DI_ergodic_trajectory(integrator.p[9])
        integrator.p[8] += 1
        integrator.p[7] = DI_ergodic_traj
        integrator.p[9] = new_ϕ_map
    end
    return
end

function condition!(u, t, integrator)
    at_charging = reached_charging_station(u, integrator.p[6])
    return at_charging == true
end

function simulate(initial_state, t0, t_max)
    tspan = (t0, t_max)

    u0 = [0., 0, 0, 0]
    com_traj = NomTrajectory()
    nom_traj = NomTrajectory()
    DI_ergodic_traj = DI_erg_Traj()
    erg_iter = 0
    flight_mode = "Mission"
    ϕ_current = (OceanographicData.get_target_distribution(distribution_type = "current")')

    params = [u0, com_traj, nom_traj, Float64[], Float64[], flight_mode, DI_ergodic_traj, erg_iter, ϕ_current]

    prob = ODEProblem(QuadSys.closed_loop_comb!, initial_state, tspan, params)

    cb_traj_erg = PeriodicCallback(recompute_ergodic_timed!, 30.0; initial_affect=true)
    cb_traj_com = PeriodicCallback(recompute_committed_timed!, 1.0; initial_affect=true)
    cb_u = PeriodicCallback(recompute_u_track_cand!, 0.05; initial_affect=true, save_positions=(false, false))
    affect!(integrator) = terminate!(integrator)
    cb_terminate = DiscreteCallback(condition!, affect!)

    cbs = CallbackSet(cb_traj_erg, cb_traj_com, cb_u, cb_terminate)

    sol = DEsolve(prob, Tsit5(), callback=cbs)
    return sol
end

end
