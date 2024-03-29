The present data and the associated source code are freely available under the GNU GPL v3 licence, EXCEPT for the *tex* files in the folders [IJHFF](/IJHFF/) and [IJHMT](/IJHMT/). Below are direct links to the different cases, followed by links to websites and publications related to the present work.

# Direct access to data and associated source code

## Imposed temperature case

The statistics are directly available in a [xls file](/../raw/master/diric/diric.xls) but it is recommended to visit the [directory](/diric/) to obtain a more detailed description of the case.

## Imposed heat flux case

The statistics are directly available in a [xls file](/../raw/master/neuma/neuma.xls) but it is recommended to visit the [directory](/neuma/) to obtain a more detailed description of the case.

## Conjugate heat-transfer cases

The statistics are directly available in xls files but it is recommended to visit the associated directories to obtain a more detailed description of the cases:

- Cases with the same thermal properties in the fluid and solid domains:
    - Regular grid: [xls file](/../raw/master/g1a1_od6/g1a1_od6.xls) and [repository](/g1a1_od6/)
    - Regular grid with spectral vanishing viscosity on the scalar: [xls file](/../raw/master/g1a1/g1a1.xls) and [repository](/g1a1/)
    - Finer grid: [xls file](/../raw/master/g1a1_refined/g1a1_refined.xls) and [repository](/g1a1_refined/)
- Cases with different thermal properties in the fluid and solid domains:
    - This [repository](/gxay/) contains the cases and their description.
    - [xls](/../raw/master/gxay/g05a2/g05a2.xls) file for $`G=\frac{1}{2}`$ and $`G_2=\frac{1}{2}`$
    - [xls](/../raw/master/gxay/g05a1/g05a1.xls) file for $`G=\frac{1}{2}`$ and $`G_2=1`$
    - [xls](/../raw/master/gxay/g05a05/g05a05.xls) file for $`G=\frac{1}{2}`$ and $`G_2=2`$
    - [xls](/../raw/master/gxay/g1a2/g1a2.xls) file for $`G=1`$ and $`G_2=\frac{1}{2}`$
    - [xls](/../raw/master/gxay/g1a05/g1a05.xls) file for $`G=1`$ and $`G_2=2`$
    - [xls](/../raw/master/gxay/g2a2/g2a2.xls) file for $`G=2`$ and $`G_2=\frac{1}{2}`$
    - [xls](/../raw/master/gxay/g2a1/g2a1.xls) file for $`G=2`$ and $`G_2=1`$
    - [xls](/../raw/master/gxay/g2a05/g2a05.xls) file for $`G=2`$ and $`G_2=2`$

## Cases with a Robin boundary condition

The statistics are directly available in xls files but it is recommended to visit the associated directories to obtain a more detailed description of the cases:

- Robin case reproducing the case $`K=\frac{1}{\sqrt{2}}`$: [xls file](/../raw/master/robin_27/robin_27.xls) and [repository](/robin_27/)
- Robin case reproducing the case $`G=1`$ and $`G_2=1`$: [xls file](/../raw/master/robin_19/robin.xls) and [repository](/robin_19/)
- Robin case reproducing the case $`K=\sqrt{2}`$: [xls file](/../raw/master/robin_14/robin_14.xls) and [repository](/robin_14/)

## Data associated with our publication in [International Journal of Heat and Fluid Flow](http://dx.doi.org/10.1016/j.ijheatfluidflow.2015.07.009)

We have identified 2 minor typos in our publication in IJHFF:

- $`G_2`$ is said to be the ratio of fluid-to-solid thermal conductivities while it is the ratio of solid-to-fluid thermal conductivities. Anyway, this ratio is equal to one in the publication.
- After equation (14), the values given for $`R`$ are erroneous as they were obtained using the dissipation rate of the temperature variance, $`\varepsilon_\theta`$, which includes the Prandtl number. Both number given should be multiplied by $`\sqrt{Pr} = \sqrt{0.71}`$. The resulting values are close to $`0.13`$, not $`0.16`$.

The [repository](/IJHFF/) contains:

- All the python scripts used to generate the figures in the publication.
- The *tex* and *bib* files used to obtain the preprint of the publication.
- A *jpg* picture sketching the computational domain.
- Several *xls* and *xlsx* files containing our results. Most of them are only symbolic links.

## Data associated with our publication in [International Journal of Heat and Mass Transfer](http://dx.doi.org/10.1016/j.ijheatmasstransfer.2017.04.005)

The [repository](/IJHMT/) contains:

- All the python scripts used to generate the figures in the publication.
- The *tex* and *bib* files used to obtain the preprint of the publication.
- A *jpg* picture sketching the computational domain. It is a symbolic link to the one in [IJHFF](/IJHFF/)
- Several *xls* and *xlsx* files containing our results. Most of them are only symbolic links.

# Bibliography

## Our roots

The present source code was first downloaded from the website of Incompact3d [the 2016-05-26 at 10:12 CEST = 08:12 UTC](http://www.incompact3d.com/uploads/5/8/7/2/58724623/channel.tar).

If using the software, you are kindly asked to cite:
- [High-order compact schemes for incompressible flows: A simple and efficient method with quasi-spectral accuracy](http://dx.doi.org/10.1016/j.jcp.2009.05.010)
- [Incompact3d: A powerful tool to tackle turbulence problems with up to O(10<sup>5</sup>) computational cores](http://dx.doi.org/10.1002/fld.2480)

## Preliminary work

Semi-implicit version of the code available here was originally implemented by T. Dairay:

- [LES of a turbulent jet impinging on a heated wall using high-order numerical schemes](http://dx.doi.org/10.1016/j.ijheatfluidflow.2014.08.001)
- [Direct numerical simulation of a turbulent jet impinging on a heated wall](http://dx.doi.org/10.1017/jfm.2014.715)
- [Simulation haute fidélité de l’aérothermique d’un jet en impact](https://tel.archives-ouvertes.fr/tel-01101235/)

## Our work

I have applied some modifications to the code during my PhD. All the publications associated with our work are available [online](https://hal.archives-ouvertes.fr/search/index/q/*/authFullName_s/C%C3%A9dric+Flageul):

- [On the discontinuity of the dissipation rate associated with the temperature variance at the fluid-solid interface for cases with conjugate heat transfer.](http://dx.doi.org/10.1016/j.ijheatmasstransfer.2017.04.005)
- [On the Discontinuity of $`\varepsilon_\theta`$ -the Dissipation Rate Associated with the Temperature Variance- at the Fluid-Solid Interface for Cases with Conjugate Heat Transfer.](https://hal.archives-ouvertes.fr/hal-01456930)
- [DNS of turbulent channel flow with conjugate heat transfer: Effect of thermal boundary conditions on the second moments and budgets.](http://dx.doi.org/10.1016/j.ijheatfluidflow.2015.07.009)
- [DNS of turbulent channel flow: can we imitate conjugate heat-transfer with a Robin boundary condition?](https://hal.archives-ouvertes.fr/hal-01323794v1)
- [Création de bases de données fines par simulation directe pour les effets de la turbulence sur les transferts thermiques pariétaux.](https://hal.archives-ouvertes.fr/tel-01321596v1)

# Remarks

There seems to be an issue with symbolic links. If you can not compile some cases, please check carefully each file and symbolic link.

# TODO

The module for the statistics would probably look much better if a derived type was used.

# Acknowledgements

The author and coworkers thank the French National Research Agency and EDF R&D for funding the study (CIFRE 2012/0047) and providing computational time on Zumbrota supercomputer (IBM - Blue-geneQ). We also thank Framasoft and the Institut Jozef Stefan for providing the gitlab services that host the present project at [https://repo.ijs.si/CFLAG/incompact3d](https://repo.ijs.si/CFLAG/incompact3d) and [https://framagit.org/CFLAG/incompact3d](https://framagit.org/CFLAG/incompact3d).
