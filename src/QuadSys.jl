module QuadSys

using LinearAlgebra
using BlockDiagonals
using ForwardDiff
using ControlSystems

using ..QuaternionUtils

# Quad Parameters
m = 0.5
ℓ = 0.1750
J = Diagonal([0.0023, 0.0023, 0.004])
g = 9.81
kt=1.0
km=0.0245

# Quad Hover Condition (Equilibrium point)
uhover = (m*g/4)*ones(4)

r0 = [0.0; 0; 1.0]
q0 = [1.0; 0; 0; 0]
v0 = zeros(3)
ω0 = zeros(3)
x0 = [r0; q0; v0; ω0];

#Thrust limits
umin = [0.2*m*g; 0.2*m*g; 0.2*m*g; 0.2*m*g]
umax = [0.6*m*g; 0.6*m*g; 0.6*m*g; 0.6*m*g]

# ================================================
h = 0.05 #20 Hz

Nx = 13     # number of states (quaternion)
Nx̅ = 12     # number of states (linearized)
Nu = 4     # number of controls
# ================================================

function quad_dynamics(x,u)
    r = x[1:3]
    q = x[4:7]/norm(x[4:7]) #normalize q just to be careful
    v = x[8:10]
    ω = x[11:13]
    Q = QuaternionUtils.qtoQ(q)

    ṙ = Q*v
    q̇ = 0.5*QuaternionUtils.L(q)*QuaternionUtils.H*ω

    v̇ = Q'*[0; 0; -g] + (1/m)*[zeros(2,4); kt*ones(1,4)]*u - QuaternionUtils.hat(ω)*v

    ω̇ = J\(-QuaternionUtils.hat(ω)*J*ω + [0 ℓ*kt 0 -ℓ*kt; -ℓ*kt 0 ℓ*kt 0; km -km km -km]*u)

    return [ṙ; q̇; v̇; ω̇]
end


function dynamics!(dx, x, u)
    r = x[1:3]
    q = x[4:7]/norm(x[4:7]) #normalize q just to be careful
    v = x[8:10]
    ω = x[11:13]
    Q = QuaternionUtils.qtoQ(q)

    dx[1:3] = Q*v
    dx[4:7]= 0.5*QuaternionUtils.L(q)*QuaternionUtils.H*ω

    dx[8:10] = Q'*[0; 0; -g] + (1/m)*[zeros(2,4); kt*ones(1,4)]*u - QuaternionUtils.hat(ω)*v

    dx[11:13] = J\(-QuaternionUtils.hat(ω)*J*ω + [0 ℓ*kt 0 -ℓ*kt; -ℓ*kt 0 ℓ*kt 0; km -km km -km]*u)

    return
end


function combined_dynamics!(dx, x, u; flight_mode = "Mission")
    r = x[1:3]
    q = x[4:7]/norm(x[4:7]) #normalize q just to be careful
    v = x[8:10]
    ω = x[11:13]
    Q = QuaternionUtils.qtoQ(q)

    dx[1:3] = Q*v
    dx[4:7]= 0.5*QuaternionUtils.L(q)*QuaternionUtils.H*ω

    dx[8:10] = Q'*[0; 0; -g] + (1/m)*[zeros(2,4); kt*ones(1,4)]*u - QuaternionUtils.hat(ω)*v

    dx[11:13] = J\(-QuaternionUtils.hat(ω)*J*ω + [0 ℓ*kt 0 -ℓ*kt; -ℓ*kt 0 ℓ*kt 0; km -km km -km]*u)

    # Battery dynamics
    if flight_mode == "Battery_Replacement"
        dx[14] = 0.0
    else
        dx[14] = -0.1*norm(u, 2)^2
    end
    return
end

#RK4 integration with zero-order hold on u
function quad_dynamics_rk4(x,u,h)
    f1 = quad_dynamics(x, u)
    f2 = quad_dynamics(x + 0.5*h*f1, u)
    f3 = quad_dynamics(x + 0.5*h*f2, u)
    f4 = quad_dynamics(x + h*f3, u)
    xn = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)

    #re-normalize quaternion
    xn[4:7] .= xn[4:7]/norm(xn[4:7])
    return xn
end

function get_linearized_model(h = 0.1)
    A = ForwardDiff.jacobian(x->quad_dynamics_rk4(x,uhover,h), x0)
    B = ForwardDiff.jacobian(u->quad_dynamics_rk4(x0, u, h), uhover)

    A̅ = Array(QuaternionUtils.E(q0)'*A*QuaternionUtils.E(q0))
    B̅ = Array(QuaternionUtils.E(q0)'*B)

    return A̅, B̅
end

function cost_matrices()
    Q = Array(3.0*I(Nx̅ ))
    R = Array(.1*I(Nu));
    return Q, R
end

function closed_loop!(dx, x, params, t)
    dynamics!(dx, x, params[1])
    return
end

function closed_loop_comb!(dx, x, params, t)
    combined_dynamics!(dx, x, params[1]; flight_mode = params[6])
    return
end

function sim_closed_loop_comb!(dx, x, params, t)
    combined_dynamics!(dx, x, params[1])
    return
end

end
