# Eclares

**Clarity-driven Ergodic Search** for autonomous oceanographic coverage with quadrotor systems. This package combines ergodic trajectory generation, model predictive control (MPC), and oceanographic data processing for persistent marine monitoring missions.

For implementation details on ergodic trajectory generation and back-to-base (B2B) trajectory planning, see the Appendix of *Eclares_Extended_Version*. The paper has been published in IEEE—see the [Citation](#citation) section below for the link.

## Requirements

- **Julia** ≥ 1.9
- See [Project.toml](Project.toml) for package dependencies

### Key Dependencies

| Package | Purpose |
|---------|---------|
| DifferentialEquations | ODE integration, callbacks |
| ControlSystems | LQR controller design |
| Convex + ECOS | MPC optimization |
| NetCDF | Oceanographic data loading |
| MeshCat + TrajOptPlots | 3D visualization |
| ForwardDiff | Automatic differentiation |
| Plots, Distributions | Plotting and statistics |

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/Eclares.git
   cd Eclares
   ```

2. Start Julia and activate the project:
   ```julia
   import Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. Place the oceanographic data file `expt_32_2019.nc4` in the `data/` folder (or project root).

## Project Structure

```
Eclares/
├── src/
│   ├── Eclares.jl           # Main package entry
│   ├── QuaternionUtils.jl   # Quaternion math
│   ├── QuadSys.jl           # Quadrotor dynamics, linearization
│   ├── RefTrajLib.jl        # Reference trajectories (Lissajous)
│   ├── OceanographicData.jl # NetCDF salinity/clarity processing
│   ├── MathUtils.jl         # Math utilities
│   ├── DomainMod.jl         # Spatial domain representation
│   ├── ErgodicManager.jl    # Ergodic coverage optimization
│   ├── TrajectoryManagerLib.jl # Trajectory optimization
│   ├── TrajectoryGeneration.jl # Ergodic trajectory generation
│   ├── SimulateCoverge.jl   # Clarity dynamics, DI ergodic sim
│   ├── Controller.jl        # MPC + LQR controller
│   ├── EclaresSim.jl        # Main simulation loop
│   └── Visualization.jl     # MeshCat animation
├── notebooks/
│   └── Eclares_Full_Code_Working.ipynb
├── data/
│   └── expt_32_2019.nc4     # Oceanographic data (NetCDF)
├── utils/
│   └── quadrotor_scaled.obj # Quadrotor mesh for visualization
├── Project.toml
└── README.md
```

## Usage

### From Notebook

```julia
# Activate project
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Eclares

# Initial state: [position, quaternion, velocity, angular_velocity, state_of_charge]
r0 = [0.0; 0.0; 0.0]
q0 = Eclares.QuaternionUtils.L([1.0; 0.0; 0; 0]) * Eclares.QuaternionUtils.rptoq([0.; 0; 0])
x0 = [r0; q0; zeros(3); zeros(3); 40.0]

# Run simulation
sol = Eclares.EclaresSim.simulate(x0, 0.0, 60.0)
```

### From Julia REPL

```julia
using Pkg
Pkg.activate("/path/to/Eclares")
Pkg.instantiate()

using Eclares

x0 = [0; 0; 0; 1; 0; 0; 0; zeros(6); 40.0]
sol = Eclares.EclaresSim.simulate(x0, 0.0, 60.0)
```

## Data Requirements

The `expt_32_2019.nc4` file (HYCOM Gulf of Mexico salinity data) is required for oceanographic clarity distributions. Place it in:

- `data/expt_32_2019.nc4`, or  
- The project root directory

## License

See [LICENSE](LICENSE) for details.

## Citation

The paper has been published in IEEE ICRA 2024. If you use this code, please cite:

- **[Eclares: Energy-Aware Clarity-Driven Ergodic Search](https://ieeexplore.ieee.org/document/10611286/)**  
  *Kaleb Ben Naveed, Devansh Agrawal, Christopher Vermillion, Dimitra Panagou*, 2024 IEEE International Conference on Robotics and Automation (ICRA), pp. 14326-14332, 2024.

```bibtex
@INPROCEEDINGS{10611286,
  author={Naveed, Kaleb Ben and Agrawal, Devansh and Vermillion, Christopher and Panagou, Dimitra},
  booktitle={2024 IEEE International Conference on Robotics and Automation (ICRA)},
  title={Eclares: Energy-Aware Clarity-Driven Ergodic Search},
  year={2024},
  pages={14326-14332},
  doi={10.1109/ICRA57147.2024.10611286}
}
```
