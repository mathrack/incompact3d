# Structure of the repository

The present repository contains:

- The source code used to run the simulation (Makefile, .f90 and .prm files)
- A *raw* folder with raw statistics (1D text files). Our code outputs statistics averaged in the streamwise direction in a binary format. The statistics in the raw folder are averaged spanwise and preprocessed. As a result, the file *vphim1d.dat* contains $`\overline{v \phi} - \overline{v} \overline{\phi}`$.
- A scilab script (.sce) that reads the raw statistics and output quantities in wall-units.
- A *csv* folder with statistics in wall-units.
- A *xls* file with statistics in wall-units.
- The present README.md

# Configuration of the turbulent channel flow, dynamic part.

Here, $`[x,y,z]`$ and $`[1,2,3]`$ will be used for the streamwise, wall-normal and spanwise directions, respectively.

---

The domain is a parallelepiped: $`[0,0,0] \leq [x,y,z] \leq [25.6, 2, 8.52]`$. The mesh is streched in the wall-normal direction (`istret = 2` and `beta = 0.225`). Periodic boundary conditions are used in the directions $`x`$ and $`z`$. At $`y=0`$ and $`y=2`$, the velocity is null and the pressure satisfies an homogeneous Neumann boundary condition.

---

The momentum equation solved reads:
```math
\partial_t u_i = - \frac{\partial_j \left( u_i u_j \right) + u_j \partial_j u_i}{2} - \partial_i p + \nu \partial_{jj} u_i + f_i
```
The kinematic viscosity $`\nu`$ is the inverse of the bulk Reynolds number, which is equal to $`2280`$ here. The source term is present only in the streamwise direction. Its amplitude is exactly $`0.0042661405`$, which leads to a unit bulk velocity.

---

The time step is $`0.002`$. After the flow reached a statistically steady state, statistics were gathered for $`1,500,000`$ time steps.

# Configuration of the turbulent channel flow, thermal part.

The scalar conservation equation reads:
```math
\partial_t \phi = - \partial_j \left( \phi u_j \right) + \frac{\nu}{Pr} \partial_{jj} \phi + \frac{\nu}{Pr} u_x
```
The value of the Prandtl number is $`0.71`$. At $`y=0`$, $`\partial_y \phi = 1`$. At $`y=2`$, $`\partial_y \phi = -1`$.
