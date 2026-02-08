module MathUtils

using LinearAlgebra

F = Float64
VF = Vector{Float64}
MF = Matrix{Float64}
VVF = Vector{VF}
VMF = Vector{MF}

function fast_pdf(x::VF, m::VF, inv_S::MF, den::Float64)
    xm = x - m
    return exp(-0.5 * dot(xm, inv_S*xm)) / den
end

function normalize!(mat::Matrix{Float64}, dA::Float64, dom_max_size::Float64)
    num_x, num_y = size(mat)
    c = 1.0 / (sum(mat) * dA)
    for xi = 1:num_x
        for yi = 1:num_y
            mat[xi,yi] *= c
        end
    end
end

function directional_derivative(g_f1::MF, g_f2::MF, u1::VVF, u2::VVF)
    N = size(g_f2, 2)
    dd = 0.0
    for i = 1:N
        dd += dot(g_f1[:,i], u1[i]) + dot(g_f2[:,i], u2[i])
    end
    dd += dot(g_f1[:,N+1], u1[N+1])
    return dd
end

end
