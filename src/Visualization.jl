module Visualization

import MeshCat as mc
using TrajOptPlots
using GeometryBasics
using LinearAlgebra
using StaticArrays
using CoordinateTransformations
using ColorTypes
using ColorSchemes
using Colors

using ..QuaternionUtils
using ..GeneralUtils

function vis_traj!(vis, name, X; R = 0.4, color = mc.RGBA(1.0, 0.8, 0.0, 5.0))
    for i = 1:length(X)
        a = X[i][1:3]
        sph = mc.HyperSphere(mc.Point(a...), R)
        mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color))
    end
end

function animate_system!(vis, Xsim, dt)
    utils_dir = joinpath(@__DIR__, "..", "utils")
    quad_obj_path = joinpath(utils_dir, "quadrotor_scaled.obj")
    if isfile(quad_obj_path)
        quad_obj = mc.MeshFileGeometry(quad_obj_path)
    else
        quad_obj = mc.HyperSphere(mc.Point(0, 0, 0), 0.2)
    end
    mc.setobject!(vis[:quad], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

    cone_radius = 0.70
    cone = mc.Cone(mc.Point(0, 0, -0.2), mc.Point(0, 0, 0.2), cone_radius)
    mc.setobject!(vis[:cone], cone, mc.MeshPhongMaterial(color= mc.RGBA(1.0, 0.7431372549, 0., 0.4)))

    anim = mc.Animation(floor(Int,1/dt))

    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            r_q = Xsim[k][1:3]
            q_q = Xsim[k][4:7]
            dcm_q = QuaternionUtils.qtoQ(q_q)
            cone_pos = deepcopy(r_q)
            cone_pos[3] -= 0.5 * cone_pos[3]
            mc.settransform!(vis[:quad], mc.compose(mc.Translation(r_q), mc.LinearMap(dcm_q)))
            cone_scaling_diag = SVector(1.0, 1.0, 5*cone_pos[3])
            cone_transformation = mc.compose(mc.Translation(cone_pos), mc.LinearMap(Diagonal(cone_scaling_diag)))
            mc.settransform!(vis[:cone], cone_transformation)
            mc.settransform!(vis["/Cameras/default"], Translation(r_q[1], r_q[2], 1))
        end
    end

    mc.setanimation!(vis, anim)
end

function create_cone(vis)
    radius = 0.5
    height = 1.0
    cone = mc.Cone(mc.Point(0, 0, -height / 2), mc.Point(0, 0, height / 2), radius)
    material = mc.MeshPhongMaterial(color= mc.RGBA(1.0, 0.6470, 0.0, 0.4))
    mc.setobject!(vis[:cone], cone, material)
end

end
