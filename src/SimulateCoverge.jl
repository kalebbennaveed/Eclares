module SimulateCoverge

using LinearAlgebra
using Plots
using NetCDF
using Distributions
using Contour
using Dates
using StatsBase
using DifferentialEquations
using ForwardDiff

using ..DomainMod
using ..MathUtils
using ..ErgodicManager
using ..TrajectoryManagerLib
using ..OceanographicData
using ..TrajectoryGeneration
using ..QuadSys

ErgodicInstance = ErgodicManager.ErgodicInstance

struct DI_erg_Traj{F}
    ts::Vector{F}
    xs::Vector{Vector{F}}
    us::Vector{Vector{F}}
end

DI_erg_Traj() = DI_erg_Traj(Float64[],Vector{Float64}[],Vector{Float64}[])

function DI_Dynamics(x, u)
    A = [0. 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0]
    B = [0 0;0 0;1 0; 0 1]
    dx = A*x + B*u
end

function DI_dynamics_rk4(x,u,h)
    f1 = DI_Dynamics(x, u)
    f2 = DI_Dynamics(x + 0.5*h*f1, u)
    f3 = DI_Dynamics(x + 0.5*h*f2, u)
    f4 = DI_Dynamics(x + h*f3, u)
    xn = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
    return xn
end

const C = 1.0
const R = 0.5
const Q = 0.0

clarity(q) = (C^2 / R) * (1-q)^2
clarity_with_decay(q, Q) = (C^2 / R) * (1-q)^2 - Q * q^2
clarity_only_decay(q, Q) = - Q * q^2
Clarity(u,p,t) = (C^2 / R) * (1-u)^2
Clarity_decay(u,p,t) = (C^2 / R) * (1-u)^2 - Q * u^2

function clarity_dynamics_rk4(q,h; fn::Function = clarity)
    f1 = fn(q)
    f2 = fn(q + 0.5*h*f1)
    f3 = fn(q + 0.5*h*f2)
    f4 = fn(q + h*f3)
    qn = q + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
    return qn
end

function clarity_dynamics_rk4_decay(q,h,Q_val; fn::Function = clarity_with_decay)
    f1 = fn(q, Q_val)
    f2 = fn(q + 0.5*h*f1, Q_val)
    f3 = fn(q + 0.5*h*f2, Q_val)
    f4 = fn(q + h*f3, Q_val)
    qn = q + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
    return qn
end

function update_clarity(em::ErgodicInstance, xds, ϕ_curr, h)
    ϕ_new = 1.0 * ϕ_curr
    Lcellx = em.domain.cell_lengths[1]
    Lcelly = em.domain.cell_lengths[2]
    @assert Lcellx == Lcelly

    for i in 1:length(xds)
        x_pos_i = xds[i][1:2]
        ϕ_xi = ceil(Int, (x_pos_i[1]/Lcellx) )
        ϕ_yi = ceil(Int, (x_pos_i[2]/Lcelly) )

        if ϕ_yi == 1.0
            row_index = ϕ_yi:ϕ_yi+1
        elseif ϕ_yi == size(ϕ_new, 1)
            row_index = ϕ_yi - 1:ϕ_yi
        else
            row_index = ϕ_yi-1:ϕ_yi+1
        end

        if ϕ_xi == 1.0
            col_index = ϕ_xi:ϕ_xi+1
        elseif ϕ_xi == size(ϕ_new, 2)
            col_index = ϕ_xi - 1:ϕ_xi
        else
            col_index = ϕ_xi-1:ϕ_xi+1
        end

        ϕ_sensed = ϕ_new[row_index, col_index]
        for cell in CartesianIndices(ϕ_sensed)
            ϕ_sensed[cell] = clarity_dynamics_rk4(ϕ_sensed[cell],h)
        end
        ϕ_new[row_index, col_index] = ϕ_sensed
    end
    return ϕ_new
end

function update_clarity_decay(em::ErgodicInstance, xds, ϕ_curr, h)
    Q_mat = (OceanographicData.get_Q_matrix()')
    ϕ_new = 1.0 * ϕ_curr
    Lcellx = em.domain.cell_lengths[1]
    Lcelly = em.domain.cell_lengths[2]
    @assert Lcellx == Lcelly

    for i in 1:length(xds)
        x_pos_i = @view xds[i][1:2]
        ϕ_xi = ceil(Int, (x_pos_i[1]/Lcellx) )
        ϕ_yi = ceil(Int, (x_pos_i[2]/Lcelly) )

        if ϕ_yi == 1.0
            row_index = ϕ_yi:ϕ_yi+1
        elseif ϕ_yi == size(ϕ_new, 1)
            row_index = ϕ_yi - 1:ϕ_yi
        else
            row_index = ϕ_yi-1:ϕ_yi+1
        end

        if ϕ_xi == 1.0
            col_index = ϕ_xi:ϕ_xi+1
        elseif ϕ_xi == size(ϕ_new, 2)
            col_index = ϕ_xi - 1:ϕ_xi
        else
            col_index = ϕ_xi-1:ϕ_xi+1
        end

        ϕ_sensed_idx = CartesianIndices((row_index, col_index))
        for cell in CartesianIndices(ϕ_new)
            if cell in ϕ_sensed_idx
                Q_val = 0.0015
                ϕ_new[cell] = clarity_dynamics_rk4_decay(ϕ_new[cell],h,Q_val; fn = clarity_with_decay)
            else
                Q_val = 0.0015
                ϕ_new[cell] = clarity_dynamics_rk4_decay(ϕ_new[cell],h,Q_val; fn = clarity_only_decay)
            end
        end
    end
    return ϕ_new
end

function remove_NaNs(mat)
    for cell in CartesianIndices(mat)
        if isnan(mat[cell])
            mat[cell] = 0
        else
            mat[cell] = max(0, mat[cell])
        end
    end
    return mat
end

function Forward_Simulate_Clarity(current_clarity, target_clarity; fn::Function = Clarity)
    tspan = (0.0, 10.0)
    u0 = current_clarity
    prob = ODEProblem(fn, u0, tspan)
    condition(u, t, integrator) = u ≥ target_clarity
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)
    sol = solve(prob, Tsit5(), callback = cb);
    return sol.t[end]
end

function get_target_clarity(target_clarity; Q_val = 0.0015)
    k = C / sqrt(Q_val * R)
    q∞ = k / (k + 1)
    target_clarity = min(target_clarity, q∞)
    return target_clarity
end

function get_target_distribution(ϕ_init_target, ϕ_current; case = "spatiostatic")
    ϕ_time_to_spent = 1.0 * ϕ_current
    for cell in CartesianIndices(ϕ_current)
        target_clarity = ϕ_init_target[cell]
        if case == "spatiotemporal"
            target_clarity = get_target_clarity(target_clarity)
        end
        current_clarity = ϕ_current[cell]
        if target_clarity > current_clarity
            ϕ_time_to_spent[cell] = Forward_Simulate_Clarity(current_clarity, target_clarity; fn = Clarity_decay)
        else
            ϕ_time_to_spent[cell] = 0.0
        end
    end
    ϕ_target_distribution = 1.0 * ϕ_time_to_spent
    for cell in CartesianIndices(ϕ_target_distribution)
        ϕ_target_distribution[cell] = ϕ_target_distribution[cell]/sum(ϕ_target_distribution)
    end
    return ϕ_time_to_spent
end

function generate_initial_DI_ergodic_trajectory(ϕ_current; N::Int64 = 150, h::Float64 = 0.2)
    x0 = [0.7, 0.6, 0.0, 0.0]
    t0 = 0.0
    ϕ_final_target = OceanographicData.get_target_distribution(distribution_type = "init_target")
    ϕ_time_difference = get_target_distribution(ϕ_final_target, ϕ_current')
    ems = ErgodicManager.ErgodicManagerR2(ϕ_time_difference)
    tm = TrajectoryManagerLib.TrajectoryManager(x0, h, N; i = TrajectoryManagerLib.ConstantInitializer([0.0,0.0]) )
    xds, uds = TrajectoryGeneration.pto_linear(ems, tm)
    ϕ_current_new = update_clarity(ems, xds, ϕ_current, h)
    DI_ergodic_traj = DI_erg_Traj(
        [t0 + (i-1) * h for i=1:length(xds)],
        xds,
        uds
        )
    return DI_ergodic_traj, ϕ_current_new
end

function generate_DI_ergodic_trajectory(ϕ_current, past_DI_erg_traj ; N::Int64 = 150, h::Float64 = 0.2)
    init_time = past_DI_erg_traj.ts[end]
    x0 = past_DI_erg_traj.xs[end]
    ϕ_final_target = OceanographicData.get_target_distribution(distribution_type = "init_target")
    ϕ_time_difference = get_target_distribution(ϕ_final_target, ϕ_current')
    ems = ErgodicManager.ErgodicManagerR2(ϕ_time_difference)
    tm = TrajectoryManagerLib.TrajectoryManager(x0, h, N; i = TrajectoryManagerLib.ConstantInitializer([0.0,0.0]) )
    xds, uds = TrajectoryGeneration.pto_linear(ems, tm)
    ϕ_current_new = update_clarity(ems, xds, ϕ_current, h)
    ts_new = [init_time + (i-1) * h for i=1:length(xds)]
    extented_time = vcat(past_DI_erg_traj.ts, ts_new[2:end])
    extended_DI_xtraj = vcat(past_DI_erg_traj.xs, xds[2:end])
    extended_DI_utraj = vcat(past_DI_erg_traj.us, uds[2:end])
    DI_ergodic_traj = DI_erg_Traj(
        extented_time,
        extended_DI_xtraj,
        extended_DI_utraj
        )
    return DI_ergodic_traj, ϕ_current_new
end

function generate_reference_trajectory(t0, N, dt, DI_ergodic_traj)
    DI_reference_traj = Vector{Vector{Float64}}(undef, N+1)

    for i = 1:(N+1)
        t = t0 + (i-1)*dt
        if t >= DI_ergodic_traj.ts[end]
            h = t - DI_ergodic_traj.ts[end]
            xn = DI_dynamics_rk4(DI_ergodic_traj.xs[end], DI_ergodic_traj.us[end], h)
            DI_reference_traj[i] = xn
        else
            ind = findfirst(DI_ergodic_traj.ts .> t) - 1
            h = t - DI_ergodic_traj.ts[ind]
            xn = DI_dynamics_rk4(DI_ergodic_traj.xs[ind], DI_ergodic_traj.us[ind], h)
            DI_reference_traj[i] = xn
        end
    end

    DI_reference_traj_T = [20.0 .* x for x in DI_reference_traj]
    pos_states = [state[1:2] for state in DI_reference_traj_T]
    pos_states_T = [[pos_vec[1] - 10.0, pos_vec[2] - 10.0] for pos_vec in pos_states]

    Quad_Xerg = [zeros(12) for i=1:length(DI_reference_traj)]
    for i = 1:(N+1)
        Quad_Xerg[i][1:2] = pos_states_T[i]
        Quad_Xerg[i][3] = 4.0
    end

    for i = 1:N
        Quad_Xerg[i][7:9] = (Quad_Xerg[i+1][1:3] - Quad_Xerg[i][1:3])/dt
    end

    Quad_Uerg = [(QuadSys.g*QuadSys.m/4)*ones(4) for i = 1:(length(DI_reference_traj)-1)]
    return Quad_Xerg[1:N], Quad_Uerg[1:N-1]
end

function update_current_map(tk, DI_ergodic_traj, ergodic_iterations, ergodic_)
    # Placeholder - can be extended for map updates
    return nothing
end

end
