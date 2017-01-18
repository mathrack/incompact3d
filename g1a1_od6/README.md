# Structure of the repository

The present repository contains:

- The source code used to run the simulation (Makefile, .f90 and .prm files)
- A *raw* folder with raw statistics (1D text files). Our code outputs statistics averaged in the streamwise direction in a binary format. The statistics in the raw folder are averaged spanwise and preprocessed. As a result, the file *vphim1d.dat* contains $`\overline{v \phi} - \overline{v} \overline{\phi}`$.
- A scilab script (.sce) that reads the raw statistics and output quantities in wall-units.
- A *csv* folder with statistics in wall-units.
- A [xls file](/../raw/master/g1a1_od6/g1a1_od6.xls) with statistics in wall-units.
- The present README.md

# Configuration of the turbulent channel flow, dynamic part.

Here, $`[x,y,z]`$ and $`[1,2,3]`$ will be used for the streamwise, wall-normal and spanwise directions, respectively.

---

The fluid domain is a parallelepiped: $`[0,0,0] \leq [x,y,z] \leq [25.6, 2, 8.52]`$. The mesh is streched in the wall-normal direction (`istret = 2` and `beta = 0.225`). Periodic boundary conditions are used in the directions $`x`$ and $`z`$. At $`y=0`$ and $`y=2`$, the velocity is null and the pressure satisfies an homogeneous Neumann boundary condition. The number of nodes in the $`[x,y,z]`$ directions is $`[256, 193, 256]`$.

---

The solid domains are parallelepipeds located on top and on bottom of the fluid domain : $`[0,2,0] \leq [x,y,z] \leq [25.6, 3, 8.52]`$ and $`[0,-1,0] \leq [x,y,z] \leq [25.6, 0, 8.52]`$. In the streamwise and spanwise directions, the grid in the solid domain is identical to the fluid one. In the wall-normal direction, a Chebyshev grid with $`129`$ interior nodes is used.

---

The momentum equation solved reads:
```math
\partial_t u_i = - \frac{\partial_j \left( u_i u_j \right) + u_j \partial_j u_i}{2} - \partial_i p + \nu \partial_{jj} u_i + f_i
```
The kinematic viscosity $`\nu`$ is the inverse of the bulk Reynolds number, which is equal to $`2280`$ here. The source term is present only in the streamwise direction. Its amplitude is exactly $`0.0042661405`$, which leads to a unit bulk velocity.

---

The time step is $`0.002`$. After the flow reached a statistically steady state, statistics were gathered for $`1,500,000`$ time steps.

# Configuration of the turbulent channel flow, thermal part.

The scalar conservation equation in the fluid domain reads:
```math
\partial_t \phi = - \partial_j \left( \phi u_j \right) + \frac{\nu}{Pr} \partial_{jj} \phi + \frac{\nu}{Pr} u_x
```
The value of the Prandtl number is $`0.71`$. At the fluid-solid interfaces, the scalar satisfies:
```math
\phi = \phi_s \mbox{ and } \partial_y \phi = G_2 \partial_y \phi_s
```
Where $`G_2`$ is the ratio of solid-to-fluid thermal conductivities and $`\phi_s`$ the scalar in the solid domain. There, the scalar conservation equation reads:
```math
\partial_t \phi_s = \frac{\nu}{G Pr} \partial_{jj} \phi_s
```
Where $`G`$ is the ratio of fluid-to-solid thermal diffusivity. At $`y=-1`$, $`\partial_y \phi_s = G_2`$. At $`y=3`$, $`\partial_y \phi_s = -G_2`$.

Here, $`G = 1`$ and $`G_2 = 1`$.
