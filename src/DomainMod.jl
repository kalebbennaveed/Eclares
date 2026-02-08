module DomainMod

using LinearAlgebra

using ..MathUtils

F = Float64
VF = Vector{Float64}
MF = Matrix{Float64}
VVF = Vector{VF}
VMF = Vector{MF}

Base.@kwdef mutable struct DomainInstance
    mins::Vector{F}
    maxes::Vector{F}
    cells::Vector{Int64}
    lengths::Vector{F}
    cell_lengths::Vector{F}
    num_dims::Int
    cell_size::F
end

x_min(domain::DomainInstance) = domain.mins[1]
y_min(domain::DomainInstance) = domain.mins[2]
x_max(domain::DomainInstance) = domain.maxes[1]
y_max(domain::DomainInstance) = domain.maxes[2]
x_size(domain::DomainInstance) = domain.cell_lengths[1]
y_size(domain::DomainInstance) = domain.cell_lengths[2]
x_cells(domain::DomainInstance) = domain.cells[1]
y_cells(domain::DomainInstance) = domain.cells[2]

function Domain(mins::Vector{F}, maxes::Vector{F}, cells::Vector{Int})
    num_dims = length(mins)
    d = DomainInstance(
        mins = deepcopy(mins),
        maxes = deepcopy(maxes),
        cells = deepcopy(cells),
        lengths = zeros(num_dims),
        num_dims = num_dims,
        cell_lengths = zeros(num_dims),
        cell_size = 1.0
    )

    for i = 1:num_dims
        d.lengths[i] = d.maxes[i] - d.mins[i]
        d.cell_lengths[i] = d.lengths[i] / d.cells[i]
        d.cell_size *= d.cell_lengths[i]
    end

    return d
end

end
