!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
!
!********************************************************************
!
subroutine parameter()
!
!********************************************************************
  
USE param
USE IBM 
USE variables
USE decomp_2d

implicit none

real(mytype) :: re, theta, cfl,cf2 
integer :: longueur ,impi,j
character :: a*80

#ifdef DOUBLE_PREC 
pi=dacos(-1.d0) 
#else
pi=acos(-1.)
#endif

twopi=2.*pi

1000 format(a,80x) 
1003 format(a,80x)
open(10,file='incompact3d.prm',status='unknown',form='formatted') 
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) xlx
read (10,*) yly 
read (10,*) zlz 
read (10,*) re 
read (10,*) sc
read (10,*) u1 
read (10,*) u2 
read (10,*) noise 
read (10,*) noise1
read (10,*) dt
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) nclx 
read (10,*) ncly 
read (10,*) nclz 
read (10,*) itype 
read (10,*) iin
read (10,*) ifirst
read (10,*) ilast
read (10,*) nscheme
read (10,*) istret
read (10,*) beta
read (10,*) iskew
read (10,*) iscalar
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) ilit 
read (10,*) isave
read (10,*) imodulo
read (10,1000) a 
read (10,1000) a 
read (10,1000) a 
read (10,*) ivirt
read (10,*) cex 
read (10,*) cey 
read (10,*) cez 
read (10,*) ra 
read (10,1000) a 
close(10) 
if (nrank==0) then
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,'======================Incompact3d=========================='
print *,'===Copyright (c) 2012 Eric Lamballais and Sylvain Laizet==='
print *,'eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com'
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,''
print *,''
print *,''
if (itype.eq.1) print *,'Constant flow field'
if (itype.eq.2) print *,'Channel flow'
if (itype.eq.3) print *,'Wake flow'
if (itype.eq.4) print *,'Mixing layer with splitter plate'
if (itype.eq.5) print *,'Channel flow'
if (itype.eq.6) print *,'Taylor Green vortices'
if (itype.eq.7) print *,'Cavity flow'
if (itype.eq.8) print *,'Flat plate Boundary layer'
if (itype.eq.9) print *,'Water tank'
write(*,1101) nx,ny,nz
write(*,1103) xlx,yly,zlz 
write(*,1102) nclx,ncly,nclz 
write(*,1104) u1,u2 
write(*,1105) re
write(*,1106) dt
if (nscheme.eq.1) print *,'Temporal scheme   : Adams-bashforth 4'
if (nscheme.eq.2) print *,'Temporal scheme   : Runge-Kutta 3'
if (nscheme.eq.3) print *,'Temporal scheme   : Runge-Kutta 4'
if (iscalar.eq.0) print *,'Passive scalar    : off'
if (iscalar.eq.1) then
   print *,'Passive scalar : on'
   write (*,1113) sc
endif
if (ivirt.eq.0) print *,'Immersed boundary : off'
if (ivirt.eq.1) then
   print *,'Immersed boundary : on old school'
   write(*,1107) cex,cey,cez
   write(*,1110) ra
endif
if (ivirt.eq.2) then
   print *,'Immersed boundary : on with Lagrangian Poly'
endif


 1101 format(' Spatial Resolution: (nx,ny,nz)=(',I4,',',I4,',',I4,')')
 1102 format(' Boundary condition: (nclx,ncly,nclz)=(',I1,',',I1,',',I1,')')
 1103 format(' Domain dimension  : (lx,ly,lz)=(',F6.1,',',F6.1,',',F6.1,')')
 1104 format(' High and low speed: u1=',F6.2,' and u2=',F6.2)
 1105 format(' Reynolds number Re: ',F15.8)
 1106 format(' Time step dt      : ',F15.8)
 1107 format(' Object centred at : (',F6.2,',',F6.2,',',F6.2,')')
 1110 format(' Object length     : ',F6.2)
 1113 format(' Schmidt number    : ',F6.2)
endif
xnu=1./re 
   
if (nclx==0) dx=xlx/nx 
if (nclx==1 .or. nclx==2) dx=xlx/(nx-1.) 
if (ncly==0) dy=yly/ny 
if (ncly==1.or.ncly==2) dy   =yly/(ny-1.) 
dx2=dx*dx
dy2=dy*dy
#ifndef TWOD
   if (nclz==0) dz=zlz/nz 
   if (nclz==1.or.nclz==2) dz=zlz/(nz-1.) 
   dz2=dz*dz
#endif

if (istret.eq.0) then
   do j=1,ny
      yp(j)=(j-1.)*dy
      ypi(j)=(j-0.5)*dy
   enddo
else
   call stretching()
endif


!******************************************************************
!
!**TIME ADVANCE***1=AB2***2=RK3***3=RK4C&K****************** 
!
!******************************************************************

adt(:)=0. ; bdt(:)=0. ; cdt(:)=0. ; gdt(:)=0.
if (nscheme==1) then!AB2
   iadvance_time=1 
   adt(1)=1.5*dt
   bdt(1)=-0.5*dt
   gdt(1)=adt(1)+bdt(1)
   gdt(3)=gdt(1)
endif
if (nscheme==2) then !RK3
   iadvance_time=3 
   adt(1)=(8./15.)*dt
   bdt(1)=0.
   gdt(1)=adt(1)
   adt(2)=(5./12.)*dt
   bdt(2)=(-17./60.)*dt
   gdt(2)=adt(2)+bdt(2)
   adt(3)=(3./4.)*dt
   bdt(3)=(-5./12.)*dt
   gdt(3)=adt(3)+bdt(3)
endif
if (nscheme==3) then !RK4 Carpenter and Kennedy  
   iadvance_time=5 
   adt(1)=0.
   adt(2)=-0.4178904745
   adt(3)=-1.192151694643
   adt(4)=-1.697784692471
   adt(5)=-1.514183444257
   bdt(1)=0.1496590219993
   bdt(2)=0.3792103129999
   bdt(3)=0.8229550293869
   bdt(4)=0.6994504559488
   bdt(5)=0.1530572479681
   gdt(1)=0.1496590219993*dt
   gdt(2)=0.220741935365*dt
   gdt(3)=0.25185480577*dt
   gdt(4)=0.33602636754*dt
   gdt(5)=0.041717869325*dt
endif

if (nscheme==4) then!AB3
   iadvance_time=1
   adt(1)= (23./12.)*dt
   bdt(1)=-(16./12.)*dt
   cdt(1)= ( 5./12.)*dt
   gdt(1)=adt(1)+bdt(1)+cdt(1)
   gdt(3)=gdt(1)
endif

return  
end subroutine parameter



