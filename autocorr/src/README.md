The autocorrelations of the temperature and wall-normal heat flux were computed on the machine IVANOE at EDF R&D using the library FFTW *AND* on the machine ZUMBROTA at EDF R&D using the library ESSL. Some *f90* files are adapted to one library, some to the other.

---

Important: The autocorrelations are *NOT* normalized.

---

The autocorrelations in the fluid were taken at $`y^+=[0,5,15,150]`$. In the solid, they were taken at $`y^+=[0,5,15,75]`$:

- jj=[1,3,5,7] stands for the autocorrelations of the temperature at increasing distance from the wall.

- jj=[2,4,6,8] stands for the autocorrelations of the wall-normal heat-flux at increasing distance from the wall.

Some cross-correlations were also computed:

- j=1 : cross correlation between T(y=0) and dTdy(y=0)
- j=2 : T(y=0) and T(y=5)
- j=3 : T(y=0) and dTdy(y=5)
- j=4 : T(y=0) and T(y=15)
- j=5 : T(y=0) and dTdy(y=15)
- j=6 : dTdy(y=0) and T(y=5)
- j=7 : dTdy(y=0) and dTdy(y=5)
- j=8 : dTdy(y=0) and T(y=15)
- j=9 : dTdy(y=0) and dTdy(y=15)
