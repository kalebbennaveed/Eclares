module TrajectoryGeneration

using LinearAlgebra

using ..DomainMod
using ..MathUtils
using ..ErgodicManager
using ..TrajectoryManagerLib

F = Float64
VF = Vector{Float64}
MF = Matrix{Float64}
VVF = Vector{VF}
VVVF = Vector{VVF}
VMF = Vector{MF}

DomainInstance = DomainMod.DomainInstance
ErgodicInstance = ErgodicManager.ErgodicInstance
TrajectoryInstanse = TrajectoryManagerLib.TrajectoryInstanse
ArmijoLineSearch = TrajectoryManagerLib.ArmijoLineSearch

function LQ2(A::MF, B::MF, Q::MF, R::MF, N::Int)
    P = [zeros(4,4) for i=1:N+1]
    G = [zeros(4,4) for i=1:N]
    Ginv = [zeros(4,4) for i=1:N]
    K = [zeros(4,4) for i=1:N]
    P[N+1] = Q

    for n = N:-1:1
        G[n] = R + (B' * P[n+1] * B)
        Ginv[n] = inv(G[n])
        K[n] = Ginv[n] * B' * P[n+1] * A
        Kp = K[n]'
        P[n] = Q + (A' * P[n+1] * A) - (Kp * G[n] * K[n])
    end
    return K, Ginv
end

function compute_ans(em::ErgodicInstance, xd::VVF, tm::TrajectoryInstanse, n::Int, ck::MF)
    x = xd[n + 1][1]
    y = xd[n + 1][2]
    Lx = em.domain.lengths[1]
    Ly = em.domain.lengths[2]
    an_x = 0.0
    an_y = 0.0
    xm = ErgodicManager.x_min(em)
    ym = ErgodicManager.y_min(em)

    for k1 = 0:em.K
        for k2 = 0:em.K
            hₖ = em.hₖ[k1+1,k2+1]
            dFk_dxn1 = -k1*pi*sin(k1*pi*(x-xm)/Lx)*cos(k2*pi*(y-ym)/Ly) / (hₖ*Lx)
            dFk_dxn2 = -k2*pi*cos(k1*pi*(x-xm)/Lx)*sin(k2*pi*(y-ym)/Ly) / (hₖ*Ly)
            c = em.λ[k1+1,k2+1] * (ck[k1+1,k2+1] - em.ϕₖ[k1+1,k2+1])
            an_x += c*dFk_dxn1
            an_y += c*dFk_dxn2
        end
    end
    an_x *= 2.0/(tm.N+1)
    an_y *= 2.0/(tm.N+1)
    return an_x, an_y
end

function gradients!(ad::MF, bd::MF, em::ErgodicInstance, tm::TrajectoryInstanse, xd::VVF, ud::VVF)
    ck = ErgodicManager.decompose(em, xd)
    gradients!(ad, bd, em, tm, xd, ud, ck)
end

function gradients!(ad::MF, bd::MF, em::ErgodicInstance, tm::TrajectoryInstanse, xd::VVF, ud::VVF, ck)
    ni =  1
    for n = 0:(tm.N-1)
        an = compute_ans(em, xd, tm, n, ck)
        for i = 1:length(an)
            ad[i,ni] = an[i]
        end
        if tm.barrier_cost > 0.0
            xnx = xd[n+1][1]
            xny = xd[n+1][2]
            xmax = ErgodicManager.x_max(em)
            xmin = ErgodicManager.x_min(em)
            ymax = ErgodicManager.y_max(em)
            ymin = ErgodicManager.y_min(em)
            if (xnx > xmax)
                ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmax)
            elseif xnx < xmin
                ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmin)
            end
            if xny > ymax
                ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymax)
            elseif xny < ymin
                ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymin)
            end
        end
        bd[:,ni] = tm.h * tm.R * ud[ni]
        ni += 1
    end

    n = tm.N
    an = compute_ans(em, xd, tm, n, ck)
    for i = 1:length(an)
        ad[i,ni] = an[i]
    end
    if tm.barrier_cost > 0.0
        xnx = xd[n+1][1]
        xny = xd[n+1][2]
        xmax = ErgodicManager.x_max(em)
        xmin = ErgodicManager.x_min(em)
        ymax = ErgodicManager.y_max(em)
        ymin = ErgodicManager.y_min(em)
        if (xnx > xmax)
            ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmax)
        elseif xnx < xmin
            ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmin)
        end
        if xny > ymax
            ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymax)
        elseif xny < ymin
            ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymin)
        end
    end
end

function LQ(K::VMF, Ginv::VMF, A::MF, B::MF, a::MF, b::MF, Q::MF, R::MF, N::Int)
    r = [zeros(4) for i=1:N+1]
    r[N+1] = 0.5*a[:,N+1]
    C = [zeros(4) for i=1:N]

    for n = N:-1:1
        Kp = K[n]'
        r[n] = (A'-Kp*B')*r[n+1] + .5*a[:,n] - .5*Kp*b[:,n]
        C[n] = Ginv[n] * (B'*r[n+1] + .5*b[:,n])
    end
    return C
end

function apply_LQ_gains(A::MF, B::MF, K::VMF, C::VVF)
    N = length(K)
    z = [ zeros(size(B,1) ) ]
    v = Array{Vector{Float64}}(undef, 0)
    for n = 1:N
        push!(v, -K[n]*z[n] - C[n])
        push!(z, A*z[n] + B*v[n])
    end
    return z, v
end

function ergodic_score(em::ErgodicInstance, traj::VVF)
    ck = ErgodicManager.decompose(em, traj)
    return ergodic_score(em, ck)
end

function ergodic_score(em::ErgodicInstance, ck)
    val = 0.0
    for (i, L) in enumerate(em.λ)
        d = em.ϕₖ[i] - ck[i]
        val += L * d * d
    end
    return val
end

function control_score(ud::VVF, R::MF, h::F)
    cs = 0.0
    num_u = length(ud[1])
    for ui in ud
        for j = 1:num_u
           cs += R[j,j] * ui[j] * ui[j]
        end
    end
    return 0.5 * h * cs
end

function barrier_score(em::ErgodicInstance, xd::VVF, c::F)
    if c == 0.0; return 0.0; end
    bs = 0.0
    xmax = ErgodicManager.x_max(em)
    ymax = ErgodicManager.y_max(em)
    xmin = ErgodicManager.x_min(em)
    ymin = ErgodicManager.y_min(em)
    for xi in xd
        if (xi[1] > xmax)
            dx = xi[1] - xmax
            bs += c * dx * dx
        elseif (xi[1] < xmin)
            dx = xi[1] - xmin
            bs += c * dx * dx
        end
        if (xi[2] > ymax)
            dy = xi[2] - ymax
            bs += c * dy * dy
        elseif (xi[2] < ymin)
            dy = xi[2] - ymin
            bs += c * dy * dy
        end
    end
    return bs
end

function total_score(em::ErgodicInstance, tm::TrajectoryInstanse, xd::VVF, ud::VVF)
    es = tm.q * ergodic_score(em, xd)
    cs = control_score(ud, tm.R, tm.h)
    bs = barrier_score(em, xd, tm.barrier_cost)
    return es + cs + bs
end

function all_scores(em::ErgodicInstance, tm::TrajectoryInstanse, xd::VVF, ud::VVF)
    es = tm.q * ergodic_score(em, xd)
    cs = control_score(ud, tm.R, tm.h)
    ts = es + cs
    return es, cs, ts
end

function project2(em::ErgodicInstance, tm::TrajectoryInstanse, K::VMF, xd::VVF, ud::VVF, zd::VVF, vd::VVF, step_size::F)
    xdn = VVF(undef, 0)
    udn = VVF(undef, 0)
    for n = 1:tm.N
        push!(udn, ud[n] + step_size*vd[n])
        push!(xdn, xd[n] + step_size*zd[n])
    end
    push!(xdn, xd[tm.N+1] + step_size*zd[tm.N+1])
    return xdn, udn
end

function get_step_size2(als::ArmijoLineSearch, em::ErgodicInstance, tm::TrajectoryInstanse, xd::VVF, ud::VVF, zd::VVF, vd::VVF, ad::MF, bd::MF, K::VMF, i::Int)
    tau = 0.5
    step_size = als.initial_step
    m = MathUtils.directional_derivative(ad, bd, zd, vd)
    f_x = total_score(em, tm, xd, ud)
    xdn, udn = project2(em, tm, K, xd, ud, zd, vd, step_size)
    armijo_index = 0
    while (total_score(em, tm, xdn, udn) > f_x + step_size*als.c*m) && (armijo_index < als.max_iters)
        step_size *= tau
        xdn, udn = project2(em, tm, K, xd, ud, zd, vd, step_size)
        armijo_index += 1
    end
    return step_size
end

function check_convergence(es::F, es_crit::F, i::Int, max_iters::Int, dd::F, dd_crit::F, es_count::Int)
    not_finished = true
    if es < es_crit; not_finished = false; end
    if i > max_iters; not_finished = false; end
    if abs(dd) < dd_crit; not_finished = false; end
    if es_count > 50; not_finished = false; end
    return not_finished
end

function pto_linear(em::ErgodicInstance, tm::TrajectoryInstanse; max_iters::Int =10000, es_crit::F = 0.0, dd_crit::F = 1e-6)
    xd0, ud0 = TrajectoryManagerLib.initialize(tm.initializer, em, tm)
    pto_linear(em, tm, xd0, ud0; max_iters=max_iters, es_crit=es_crit, dd_crit = dd_crit)
end

function pto_linear(em::ErgodicInstance, tm::TrajectoryInstanse, xd0::VVF, ud0::VVF ; max_iters::Int=10000, es_crit::F=0.0, dd_crit::F=1e-6)
    xd = deepcopy(xd0)
    ud = deepcopy(ud0)
    N = tm.N
    ad = zeros(tm.dynamics.n, N+1)
    bd = zeros(tm.dynamics.m, N)
    i = 1
    not_finished = true
    es = 0.; cs = 0.; ts = 0.; dd = 0.; step_size = 0.
    es_prev = 0.0
    es_count = 1
    A = tm.dynamics.A
    B = tm.dynamics.B
    K, Ginv = LQ2(A, B, tm.Qn, tm.Rn, tm.N)

    while not_finished
        gradients!(ad, bd, em, tm, xd, ud)
        C = LQ(K, Ginv, A, B, ad, bd, tm.Qn, tm.Rn, tm.N)
        zd, vd = apply_LQ_gains(A, B, K, C)
        step_size = get_step_size2(tm.descender, em, tm, xd, ud, zd, vd, ad, bd, K, i)
        xd, ud = project2(em, tm, K, xd, ud, zd, vd, step_size)
        es, cs, ts = all_scores(em, tm, xd, ud)
        dd = MathUtils.directional_derivative(ad, bd, zd, vd)
        i += 1
        es_count = abs(es - es_prev) < 1e-7 ? es_count + 1 : 0
        es_prev = es
        not_finished = check_convergence(es,es_crit,i,max_iters,dd,dd_crit, es_count)
    end
    return xd, ud
end

end
