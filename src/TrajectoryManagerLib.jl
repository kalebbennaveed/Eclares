module TrajectoryManagerLib

using LinearAlgebra

using ..DomainMod
using ..MathUtils
using ..ErgodicManager

F = Float64
VF = Vector{Float64}
MF = Matrix{Float64}
VVF = Vector{VF}
VVVF = Vector{VVF}
VMF = Vector{MF}

DomainInstance = DomainMod.DomainInstance
ErgodicInstance = ErgodicManager.ErgodicInstance

abstract type Dynamics end
abstract type Initializer end
abstract type Descender end
abstract type IntegrationScheme end

Base.@kwdef mutable struct TrajectoryInstanse
    N::Int
    h::F
    x0::VF
    q::F
    R::MF
    Qn::MF
    Rn::MF
    barrier_cost::F
    initializer::Initializer
    descender::Descender
    dynamics::Dynamics
    int_scheme::IntegrationScheme
end

struct LinearDynamics <: Dynamics
    n::Int
    m::Int
    A::MF
    B::MF

    function LinearDynamics(A, B)
        n,m = size(B)
        return new(n,m,deepcopy(A),deepcopy(B) )
    end
end

function DoubleIntegrator(n::Int, h::F)
    A = MF(I(2*n))
    for i = 1:n
        A[i, n+i] = h
    end

    B = zeros(2*n, n)
    for i = 1:n
        B[n+i, i] = h
    end

    return LinearDynamics(A,B)
end

struct ConstantInitializer <: Initializer
    action::Vector{Float64}
end

struct RandomInitializer <: Initializer
end

function initialize(ci::ConstantInitializer, em::ErgodicInstance, tm::TrajectoryInstanse)
    xd = [[tm.x0[1] + 0.0 * i / (tm.N+1), tm.x0[2] + 0.0 * i / (tm.N+1), (0.0)/((tm.N + 1) * tm.h),  (0.0)/((tm.N + 1) * tm.h)] for i=1:tm.N+1]
    ud = [[0., 0.] for i=1:tm.N]
    xd[1] = deepcopy(tm.x0)
    ud[1] = deepcopy(ci.action)
    for i = 1:(tm.N-1)
        xd[i+1] = integrate(tm, xd[i], ud[i])
        ud[i+1] = deepcopy(ci.action)
    end
    xd[tm.N+1] = integrate(tm, xd[tm.N], ud[tm.N])
    return xd, ud
end

function initialize(ri::RandomInitializer, em::ErgodicInstance, tm::TrajectoryInstanse)
    initialize(ConstantInitializer([0.0, 0.0]), em, tm)
end

mutable struct ArmijoLineSearch <: Descender
    initial_step::F
    c::F
    max_iters::F

    function ArmijoLineSearch(initial_step::Real, c::Real, mi::Real)
        return new(float(initial_step), float(c), float(mi) )
    end
    function ArmijoLineSearch(initial_step::Real, c::Real)
        return new(float(initial_step), float(c), 50.)
    end
    ArmijoLineSearch() = ArmijoLineSearch(1, 0.01, 50.)
end

struct ForwardEuler <: IntegrationScheme end

function integrate(tm::TrajectoryInstanse, x::VF, u::VF)
    integrate(tm.int_scheme, tm.dynamics, x, u, tm.h)
end

function integrate(::ForwardEuler, d::Dynamics, x::VF, u::VF, h::Float64)
    forward_euler(d, x, u, h)
end

function forward_euler(ld::LinearDynamics, x::VF, u::VF, h::Float64)
    return ld.A*x + ld.B*u
end

function TrajectoryManager(x0::VF, h::Real, N::Int; i::Initializer=RandomInitializer())
    tm = TrajectoryInstanse(
            N = N,
            h = h,
            x0 = x0,
            q = 100.0,
            Qn = MF(I(4)),
            R = 0.1*h * MF(I(2)),
            Rn = 0.1*h * MF(I(2)),
            barrier_cost = 38.0,
            initializer = i,
            descender = ArmijoLineSearch(),
            dynamics = DoubleIntegrator(2, h),
            int_scheme = ForwardEuler()
        )
    return tm
end

end
