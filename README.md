# Double Pendulum

This project simulates the dynamics of a double pendulum using the Runge-Kutta 4 (RK4) integration method and optimizes the system parameters (masses, initial velocities) to match experimental tracking data.

## Files Structure

- **`double_pendulum.jl`**: Core simulation logic. Contains the equations of motion, RK4 integrator, and visualization functions to generate MP4 animations and position plots.
- **`optimisation_pendule.jl`**: Optimization script. Uses the `Optim.jl` package (Particle Swarm Optimization followed by Nelder-Mead) to find the best-fitting parameters (`m1`, `m2`, `theta1_dot`, `theta2_dot`) that minimize the error (RMSE) between the simulation and real-world data.
- **`*.csv`**: Tracked position data files (e.g., `mass_a_200.csv`, `mass_b_200.csv`) used as ground truth for optimization. Files created by using [Tracker Online](https://opensourcephysics.github.io/tracker-online/index_ff.html) and selecting each position by hand.

## Dependencies

### Julia Packages
The project is written in Julia. You will need the following packages:
```julia
using DataFrames, CSV, Plots, LinearAlgebra, Statistics, Optim, LineSearches
```
Install them via the Julia REPL:
```julia
] add DataFrames CSV Plots Optim LineSearches
```

### System Tools
- **`make`**: To automate the execution workflow.
- **`ImageMagick`**: Required for the `montage` command used in the Makefile to combine plots.

## Usage

To run the parameter optimization and generate the comparison visualization:

```bash
make
```

This script will:
1.  Launch `optimisation_pendule.jl`
2.  Which will load the tracking data (`mass_a_200.csv`, `mass_b_200.csv`).
3.  Run a global optimization (Particle Swarm) followed by a local refinement (Nelder-Mead) to estimate masses and initial angular velocities.
4.  Output the best parameters found.
5.  Generate an animation (`double_pendulum.mp4`) and plots (`masse_a_x.png`, `masse_b_x.png`, etc.) comparing the simulation with the experimental data.
6.  Combine all simulation data plots (`masse_*_*.png`) into one single `resultat_combine.png` file.

Or else you can run :
```bash
julia double_pendulum.jl
```

## Methodology

- **Dynamics**: Solved using Lagrangian mechanics and Cramer's rule for the system of linear equations for accelerations.
- **Integration**: Fixed-step Runge-Kutta 4th order (RK4).
- **Optimization**: The cost function calculates the Root Mean Square Error (RMSE) between the simulated coordinates and the tracked coordinates from the CSV files.
