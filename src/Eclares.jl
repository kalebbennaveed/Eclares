"""
    Eclares

Clarity-driven Ergodic Search for autonomous coverage with quadrotor systems.
Combines ergodic trajectory generation, MPC control, and oceanographic data processing.
"""
module Eclares

# Include submodules in dependency order
include("QuaternionUtils.jl")
include("GeneralUtils.jl")
include("QuadSys.jl")
include("RefTrajLib.jl")
include("MathUtils.jl")
include("DomainMod.jl")
include("OceanographicData.jl")
include("ErgodicManager.jl")
include("TrajectoryManagerLib.jl")
include("TrajectoryGeneration.jl")
include("SimulateCoverge.jl")
include("Controller.jl")
include("Visualization.jl")

# Main Eclares simulation module (separate file to avoid circular dependency)
include("EclaresSim.jl")

# Export main types and functions
export QuaternionUtils, GeneralUtils, QuadSys, RefTrajLib
export MathUtils, DomainMod, OceanographicData, ErgodicManager
export TrajectoryManagerLib, TrajectoryGeneration, SimulateCoverge
export Controller, Visualization, EclaresSim

end
