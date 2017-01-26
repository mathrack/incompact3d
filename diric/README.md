# Structure of the repository

The present repository contains:

- The source code used to run the simulation (Makefile, .f90 and .prm files)
- A *raw* folder with raw statistics (1D text files). Our code outputs statistics averaged in the streamwise direction in a binary format. The statistics in the raw folder are averaged spanwise and preprocessed. As a result, the file *vphim1d.dat* contains $`\overline{v \phi} - \overline{v} \overline{\phi}`$.
- A scilab script (.sce) that reads the raw statistics and output quantities in wall-units.
- A *csv* folder with statistics in wall-units.
- A [xls file](/../raw/master/diric/diric.xls) with statistics in wall-units.
- The present README.md

# Configuration of the turbulent channel flow, dynamic part.

Here, $`[x,y,z]`$ and $`[1,2,3]`$ will be used for the streamwise, wall-normal and spanwise directions, respectively.

---

The fluid domain is a parallelepiped: $`[0,0,0] \leq [x,y,z] \leq [25.6, 2, 8.52]`$. The mesh is streched in the wall-normal direction (`istret = 2` and `beta = 0.225`). Periodic boundary conditions are used in the directions $`x`$ and $`z`$. At $`y=0`$ and $`y=2`$, the velocity is null and the pressure satisfies an homogeneous Neumann boundary condition. The number of nodes in the $`[x,y,z]`$ directions is $`[256, 193, 256]`$.

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
The value of the Prandtl number is $`0.71`$. At $`y=0`$ and $`y=2`$, $`\phi = 0`$.

# Wall-units

Statistics in the *xls* file and in the *csv* folder are in wall-units. The conversion from computational units to wall-units is performed in the scilab script (.sce). This conversion is briefly described here. For further details, please consult a good book on Turbulence and (Computational) Fluid Mechanics. For instance *Turbulent flows* by S. B. Pope, *The theory of homogeneous turbulence* by G. K. Batchelor or *A first course in turbulence* by H. Tennekes and J. L. Lumley.

At the wall $`y=0`$, the friction velocity $`u_\tau`$ verifies:
```math
u_\tau = \sqrt{ \nu \partial_y \overline{U_x} \left( y=0 \right) }
```

And the friction temperature $`T_\tau`$ verifies:
```math
T_\tau = \frac{\overline{q_w}}{\rho \; C_p \; u_\tau} = \nu \frac{\partial_y \overline{\phi} \left( y=0 \right)}{Pr \; u_\tau}
```

The velocity is converted to wall-units when divided by $`u_\tau`$. The temperature is converted to wall-units when divided by $`T_\tau`$. Distances are converted to wall-units when multiplied by $`\frac{u_\tau}{\nu}`$. Application of dimensional analysis should easily allow one to convert time or pressure to wall-units.

For the budgets of the Reynolds stresses, please see equation (1) in [Mansour, Kim and Moin](https://doi.org/10.1017/S0022112088002885)

For the budgets of the turbulent heat fluxes, please see equation (12) in [Kozuka, Seki and Kawamura](http://dx.doi.org/10.1016/j.ijheatfluidflow.2009.02.023)

For the budget of the temperature variance, some look at the budget of $`\overline{\phi^2}`$ and some look at the budget of $`\frac{\overline{\phi'^2}}{2}`$, by analogy with $`k`$, the turbulent kinetic energy, which also contains a factor 2. Below is the budget equation of the latter:
```math
\partial_t \frac{\overline{\phi'^2}}{2} + \partial_k \left( \overline{u_k} \frac{\overline{\phi'^2}}{2} \right) = - \overline{u'_k \phi'} \partial_k \overline{\phi} -\partial_k \left( \overline{u'_k \frac{\phi'^2}{2}}\right) + \frac{1}{Pr} \partial_{kk} \frac{\overline{\phi'^2}}{2} - \frac{1}{Pr} \overline{\partial_k \phi' \partial_k \phi'}
```
