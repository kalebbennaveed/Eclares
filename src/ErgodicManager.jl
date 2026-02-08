module ErgodicManager

using LinearAlgebra

using ..DomainMod
using ..MathUtils
using ..OceanographicData

F = Float64
VF = Vector{Float64}
MF = Matrix{Float64}
VVF = Vector{VF}
VVVF = Vector{VVF}
VMF = Vector{MF}

DomainInstance = DomainMod.DomainInstance

Base.@kwdef mutable struct ErgodicInstance
    domain::DomainInstance
    K::Int64
    ϕ::Matrix{F}
    ϕₖ::Matrix{F}
    λ::Matrix{F}
    hₖ::Matrix{F}
    kpixl::Matrix{F}
    kpiyl::Matrix{F}
end

x_min(em::ErgodicInstance) = em.domain.mins[1]
y_min(em::ErgodicInstance) = em.domain.mins[2]
z_min(em::ErgodicInstance) = em.domain.mins[3]
x_max(em::ErgodicInstance) = em.domain.maxes[1]
y_max(em::ErgodicInstance) = em.domain.maxes[2]
z_max(em::ErgodicInstance) = em.domain.maxes[3]
x_size(em::ErgodicInstance) = em.domain.cell_lengths[1]
y_size(em::ErgodicInstance) = em.domain.cell_lengths[2]
z_size(em::ErgodicInstance) = em.domain.cell_lengths[3]
x_cells(em::ErgodicInstance) = em.domain.cells[1]
y_cells(em::ErgodicInstance) = em.domain.cells[2]
z_cells(em::ErgodicInstance) = em.domain.cells[3]

function Lambda!(em::ErgodicInstance)
    for k1 = 0:em.K, k2 = 0:em.K
        den = (1.0 + k1*k1 + k2*k2) ^ 1.5
        em.λ[k1+1, k2+1] = 1.0 / den
    end
end

function kpixl!(em::ErgodicInstance)
    Lx = em.domain.lengths[1]
    xmin = x_min(em)
    for xi = 1:x_cells(em)
        x = xmin + (xi-0.5)*x_size(em)
        for k = 0:em.K
            em.kpixl[k+1,xi] = cos(k*pi*(x-xmin) / Lx)
        end
    end

    Ly = em.domain.lengths[2]
    ymin = y_min(em)
    for yi = 1:y_cells(em)
        y = ymin + (yi-0.5)*y_size(em)
        for k = 0:em.K
            em.kpiyl[k+1,yi] = cos(k*pi*(y-ymin) / Ly)
        end
    end
end

function hk_ij(em::ErgodicInstance, k1::Int, k2::Int)
    val = 0.0
    for xi = 1:x_cells(em)
        cx = em.kpixl[k1+1,xi]
        cx2 = cx * cx
        for yi = 1:y_cells(em)
            cy = em.kpiyl[k2+1,yi]
            val += cx2 * cy * cy * em.domain.cell_size
        end
    end
    return sqrt(val)
end

function hk!(em::ErgodicInstance)
    for k1 = 0:em.K, k2 = 0:em.K
        em.hₖ[k1+1,k2+1] = hk_ij(em, k1, k2)
    end
end

function phi_ij(em::ErgodicInstance, k1::Int, k2::Int, d::MF)
    val = 0.0
    for xi = 1:x_cells(em)
        cx = em.kpixl[k1+1,xi]
        for yi = 1:y_cells(em)
            cy = em.kpiyl[k2+1,yi]
            val += d[xi,yi] * cx * cy * em.domain.cell_size
        end
    end
    return val / em.hₖ[k1+1,k2+1]
end

function decompose!(em::ErgodicInstance, d::MF)
    for k1 = 0:em.K, k2 = 0:em.K
        em.ϕₖ[k1+1,k2+1] = phi_ij(em, k1, k2, d)
    end
    em.ϕ = d
end

function decompose(em::ErgodicInstance, traj::VVF)
    N = length(traj)-1
    ck = zeros(em.K+1, em.K+1)
    Lx = em.domain.lengths[1]
    Ly = em.domain.lengths[2]
    xmin = x_min(em)
    ymin = y_min(em)

    for k1 = 0:em.K
        kpiL1 = k1 * π / Lx
        for k2 = 0:em.K
            kpiL2 = k2 * π / Ly
            hₖ = em.hₖ[k1+1, k2+1]
            fk_sum = 0.0
            for n = 0:N
                xn = traj[n + 1]
                c1 = cos(kpiL1 * (xn[1]-xmin))
                c2 = cos(kpiL2 * (xn[2]-ymin))
                fk_sum += c1*c2
            end
            ck[k1+1, k2+1] = fk_sum / (hₖ * (N+1))
        end
    end
    return ck
end

function ErgodicManagerR2(d::DomainInstance, ϕ::MF; K::Int=5)
    em = ErgodicInstance(
          domain = deepcopy(d),
          K = K,
          ϕ = deepcopy(ϕ),
          ϕₖ = zeros(K+1, K+1),
          λ = zeros(K+1, K+1),
          hₖ = zeros(K+1, K+1),
          kpixl = zeros(K+1, d.cells[1]),
          kpiyl = zeros(K+1, d.cells[2])
         )
    Lambda!(em)
    kpixl!(em)
    hk!(em)
    decompose!(em, em.ϕ)
    return em
end

function ErgodicManagerR2(ϕ::MF)
    K = 5
    d = DomainMod.Domain([0., 0.], [1.0, 1.0], collect(size(ϕ)))
    ϕn = convert(Matrix{Float64}, ϕ)
    MathUtils.normalize!(ϕn, d.cell_size, d.maxes[1])
    return ErgodicManagerR2(d, ϕn, K = K)
end

end
