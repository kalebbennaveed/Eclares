module OceanographicData

using Plots
using NetCDF
using Distributions
using Contour
using Dates
using StatsBase

"""
    prob_in_range(μ, σ, l, h)

returns the probability that l <= x <= h for x ∼ N(μ, σ²)
"""
function prob_in_range(μ, σ, l, h)
    if isnan(μ) || isnan(σ) || isnan(l) || isnan(h)
        return NaN
    end
    g = Normal(μ, σ)
    return cdf(g, h) - cdf(g, l)
end

function clarity(σ)
    return 1/(1+σ^2)
end

"""
    indices_in_range(v, l, h)
given a sorted vector v, it returns the index range of v that lies in [l, h]
"""
function indices_in_range(v, l, h)
    low_index = findfirst(v .>= l)
    high_index = findlast(v .<= h)
    return low_index:high_index
end

"""
    rectangle!(xl, xh,  yl, yh; kwargs...)
plots a rectangle with corners (xl, yl) and (xh, yh)
"""
function rectangle!(xl, xh,  yl, yh; kwargs...)
    plot!([xl, xh, xh, xl, xl], [yl, yl, yh, yh, yl]; kwargs...)
end

"""
    (xs, ys) = get_contour(x,y,z,l)
returns the coordinates of the contour lines of (x, y, z) at level l
"""
function get_contour(x, y, z, l)
    conts = contours(x, y, z, [l])
    coords_x = Float64[]
    coords_y = Float64[]
    for l in lines(levels(conts)[1])
        xs, ys = coordinates(l)
        append!(coords_x, xs)
        push!(coords_x, NaN)
        append!(coords_y, ys)
        push!(coords_y, NaN)
    end
    return coords_x, coords_y
end

function get_dates(filename)
    t =  ncread(filename, "MT")
    dates = Date(1900,12,31) .+ Day.(t)
    return dates
end

function get_lat_lon(filename)
    lat = ncread(filename, "Latitude")
    lon = ncread(filename, "Longitude")
    return lat, lon
end

function get_salinity(filename)
    salinity_fill_value = ncgetatt(filename, "salinity", "_FillValue")
    salinity = ncread(filename, "salinity")[:,:,1,:]
    salinity[salinity .== salinity_fill_value] .= NaN
    return salinity
end

function make_domain_square!(ϕ)
    rows = size(ϕ,1)
    cols = size(ϕ,2)
    diff = abs(rows - cols)
    if rows != cols
        if rows > cols
            ϕ = [ϕ zeros(rows, diff)]
        else
            ϕ = [ϕ; zeros(diff, cols)]
        end
    end

    # Pad the domain by 10 on each side
    ϕ = [zeros(size(ϕ,1), 10) ϕ zeros(size(ϕ,1), 10)]
    ϕ = [zeros(10, size(ϕ,2) ); ϕ; zeros(10, size(ϕ,2) )]

    return ϕ
end

# Data file path - relative to project root
const DATA_DIR = joinpath(@__DIR__, "..", "data")
const DEFAULT_FILENAME = joinpath(DATA_DIR, "expt_32_2019.nc4")

# Load data on module initialization (lazy - only when needed)
function _load_oceanographic_data()
    filename = isfile(DEFAULT_FILENAME) ? DEFAULT_FILENAME : "expt_32_2019.nc4"
    if !isfile(filename)
        error("Oceanographic data file not found. Place expt_32_2019.nc4 in data/ or project root.")
    end
    dates = get_dates(filename)
    lat, lon = get_lat_lon(filename)
    salinity = get_salinity(filename)
    mean_salinity = mean(salinity, dims=3)[:,:,1]
    stddev_salinity = std(salinity, dims=3)[:,:,1]
    clarity_salinity = clarity.(stddev_salinity)
    prob_salinity_is_32 = similar(mean_salinity)
    for index in CartesianIndices(mean_salinity)
        prob_salinity_is_32[index] = prob_in_range(mean_salinity[index], stddev_salinity[index], 31, 33)
    end
    return filename, dates, lat, lon, salinity, mean_salinity, stddev_salinity, clarity_salinity, prob_salinity_is_32
end

# Lazy-loaded globals
const _data_cache = Ref{Union{Nothing,Tuple}}(nothing)

function _get_data()
    if _data_cache[] === nothing
        _data_cache[] = _load_oceanographic_data()
    end
    return _data_cache[]
end

function get_Q_matrix()
    _, _, _, _, _, _, stddev_salinity, _, _ = _get_data()

    roi1 = ((-90., -87.),  (28.5, 30.5))
    roi2 = ((-98., -90.), (25., 30.))
    roi = roi1

    lon_inds = indices_in_range(_get_data()[4], roi[1][1], roi[1][2])
    lat_inds = indices_in_range(_get_data()[3], roi[2][1], roi[2][2])

    Q = (stddev_salinity)[lon_inds, lat_inds]

    for cell in CartesianIndices(Q)
        if isnan(Q[cell])
            Q[cell] = 0
        else
            Q[cell] = max(0, Q[cell])
        end
    end

    Q = make_domain_square!(Q)
    return Q
end

function get_target_distribution(; distribution_type = "final_target")
    _, _, lat, lon, _, mean_salinity, stddev_salinity, clarity_salinity, prob_salinity_is_32 = _get_data()

    roi1 = ((-90., -87.),  (28.5, 30.5))
    roi2 = ((-98., -90.), (25., 30.))
    roi = roi1

    lon_inds = indices_in_range(lon, roi[1][1], roi[1][2])
    lat_inds = indices_in_range(lat, roi[2][1], roi[2][2])
    target_clarity = 1 .- (1 .-  prob_salinity_is_32)/2;

    if distribution_type == "init_target"
        ϕ = (target_clarity )[lon_inds, lat_inds]
    elseif distribution_type == "current"
        ϕ = (clarity_salinity)[lon_inds, lat_inds]
    elseif distribution_type == "final_target"
        ϕ = (target_clarity - clarity_salinity)[lon_inds, lat_inds]
    else
        error("Unknown distribution_type: $distribution_type")
    end

    for cell in CartesianIndices(ϕ)
        if isnan(ϕ[cell])
            ϕ[cell] = 0
        else
            ϕ[cell] = max(0, ϕ[cell])
        end
    end

    ϕ = make_domain_square!(ϕ)
    return ϕ
end

function get_uniform_distribution()
    ϕ = get_target_distribution(;distribution_type = "init_target")
    for cell in CartesianIndices(ϕ)
        if ϕ[cell] > 0.
            ϕ[cell] = 1.
        else
            ϕ[cell] = 0.
        end
    end
    return ϕ
end

end
