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
using DataFrames, CSV, Plots, LinearAlgebra, Statistics, Optim
```
Install them via the Julia REPL:
```julia
] add DataFrames CSV Plots Optim
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

### Lagrangian Mechanics
Taken the Lagrangian equations from [Wikipedia](https://fr.wikipedia.org/wiki/Pendule_double) and applied the formulas
$L=T-V$
$T={\frac {1}{2}}m_{1}l_{1}^{2}{\dot {\theta }}_{1}^{2}+{\frac {1}{2}}m_{2}[l_{1}^{2}{\dot {\theta }}_{1}^{2}+l_{2}^{2}{\dot {\theta }}_{2}^{2}+2l_{1}l_{2}{\dot {\theta }}_{1}{\dot {\theta }}_{2}\cos(\theta _{1}-\theta _{2})]$
$ V=-(m_{1}+m_{2})gl_{1}\cos(\theta _{1})-m_{2}gl_{2}\cos(\theta _{2})$

Obtained equations:
$$
{\begin{array}{l}(m_{1}+m_{2})l_{1}{\ddot {\theta }}_{1}+m_{2}l_{2}{\ddot {\theta }}_{2}\cos(\theta _{1}-\theta _{2})+m_{2}l_{2}{\dot {\theta }}_{2}^{2}\sin(\theta _{1}-\theta _{2})+(m_{1}+m_{2})g\sin(\theta _{1})=0\\l_{1}{\ddot {\theta }}_{1}\cos(\theta _{1}-\theta _{2})+l_{2}{\ddot {\theta }}_{2}-l_{1}{\dot {\theta }}_{1}^{2}\sin(\theta _{1}-\theta _{2})+g\sin(\theta _{2})=0\end{array}}
$$

Final equations (to solve for $\ddot{\theta}_1$ and $\ddot{\theta}_2$):

$M_{11} = (m_1 + m_2) * l_1$
$M_{12} = m_2 * l_2 * cos(\theta_1 - \theta_2)$
$R_1  = -m_2 * l_2 * \dot{\theta}_2^2 * sin(\theta_1 - \theta_2) - (m_1 + m_2) * g * sin(\theta_1)$

$M_{21} = l_1 * cos(\theta_1 - \theta_2)$
$M_{22} = l_2$
$R_2  = l_1 * \dot{\theta_1}^2 * sin(\theta_1 - \theta_2) - g * sin(\theta_2)$

### Runge Kutta
The system state is defined by the vector $u = [\theta_1, \theta_2, \dot{\theta}_1, \dot{\theta}_2]$. To evolve this state over time, we use the 4th-order Runge-Kutta method (RK4) with a fixed time step $\Delta t$:

$k_1 = f(u_n)$
$k_2 = f(u_n + \frac{\Delta t}{2} k_1)$
$k_3 = f(u_n + \frac{\Delta t}{2} k_2)$
$k_4 = f(u_n + \Delta t k_3)$

$u_{n+1} = u_n + \frac{\Delta t}{6} (k_1 + 2k_2 + 2k_3 + k_4)$

Where $f(u)$ is the function calculating the derivatives $[\dot{\theta}_1, \dot{\theta}_2, \ddot{\theta}_1, \ddot{\theta}_2]$.

### Optimization Strategy

To find the physical parameters that best match the experimental data, we use a two-step optimization approach using the `Optim.jl` package.

1.  **Global Search (Particle Swarm Optimization)**:
    -   We start with ParticleSwarm to explore the parameter space globally and avoid getting stuck in local minimums.
    -   Algorithm: `ParticleSwarm(n_particles=1000)`
    -   Bounds: Defined for masses and initial velocities to keep physical realism.

2.  **Local Refinement (Nelder-Mead)**:
    -   The best solution found by ParticleSwarm is used as the starting point for the Nelder-Mead algorithm.
    -   This step fine-tunes the parameters to converge to the precise minimum.

