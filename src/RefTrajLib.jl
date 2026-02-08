module RefTrajLib

using ComponentArrays
using Parameters: @unpack
using StaticArrays

using ..QuadSys

m = QuadSys.m
g = QuadSys.g

# ****** Uref is given is reduced state here with [r, p, v, ω] ******

trajectory_params = ComponentArray(
    A = [5,5, 0.7],
    ω = [5/8,4/8,6/8],
    ψ = [π/2, 0, 0],
    off = [0.,0,5.]
)

function lissajous(t, params, D::Integer = 0)
    @unpack A, ω, ψ, off = params
    return lissajous(t, A, ω, ψ, off, D)
end

function lissajous(t, A, ω, ψ, off, D = 0)
    N = length(A)
    if D == 0
        return SVector{N}(off[i] + A[i] * sin(ω[i] * t + ψ[i]) for i = 1:N)
    else
        return SVector{N}(A[i] * ω[i]^D * sin(ω[i] * t  + π * D /2 + ψ[i]) for i = 1:N)
    end
end

function get_desired_trajectory(t0, traj_type, N, dt)
    if (traj_type == 1)
        Xref = [zeros(12) for i=1:N]
        for i = 1:(N)
            ti = t0 + (i-1)*dt
            Xref[i][1:3] = lissajous(ti, trajectory_params)
            Xref[i][7:9] = lissajous(ti, trajectory_params, 1)
        end
        Uref = [(g*m/4)*ones(4) for i = 1:(N-1)]
        return Xref, Uref
    end
end

function get_desired_trajectory_point(t0, N, dt)
    Xref = [zeros(12) for i=1:N]
    for i = 1:N
        ti = t0 + (i-1)*dt
        Xref[i][1:3] = [0.0; 0.0; 0.11]
    end
    Uref = [(g*m/4)*ones(4) for i = 1:(N-1)]
    return Xref, Uref
end

end
