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

module variables

  use decomp_2d, only : mytype

#include "module_param_diff.f90"
!end module variables

!module filter
real(mytype), dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
real(mytype), dimension(nx,2) ::filax,filaxp
real(mytype), dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
real(mytype), dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
real(mytype), dimension(ny,2) ::filay,filayp
real(mytype), dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
real(mytype), dimension(nz) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
real(mytype), dimension(nz,2) ::filaz,filazp
real(mytype), dimension(nz) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
integer, dimension(200) :: idata

!module derivative
real(mytype), dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
real(mytype), dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
real(mytype), dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
real(mytype), dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
real(mytype), dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
real(mytype), dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
real(mytype), save, allocatable, dimension(:,:) :: sx,vx
real(mytype), save, allocatable, dimension(:,:) :: sy,vy
real(mytype), save, allocatable, dimension(:,:) :: sz,vz

!module scalar
real(mytype), dimension(nx) :: sfxt,scxt,sbxt,ssxt,swxt
real(mytype), dimension(nx) :: sfxpt,ssxpt,swxpt
real(mytype), dimension(ny) :: sfyt,scyt,sbyt,ssyt,swyt
real(mytype), dimension(ny) :: sfypt,ssypt,swypt
real(mytype), dimension(nz) :: sfzt,sczt,sbzt,sszt,swzt
real(mytype), dimension(nz) :: sfzpt,sszpt,swzpt

!module implicit
real(mytype), dimension(ny) :: aam,bbm,ccm,ddm,eem,ggm,hhm,wwm,zzm !!TIME IMPLICIT, ncl=2
real(mytype), dimension(ny) :: rrm,qqm,vvm,ssm !!TIME IMPLICIT (with HPL), ncl=2
real(mytype), dimension(ny) :: aam10,bbm10,ccm10,ddm10,eem10,ggm10,hhm10,wwm10,zzm10 !!TIME IMPLICIT, ncl=1, npaire=0
real(mytype), dimension(ny) :: rrm10,qqm10,vvm10,ssm10 !!TIME IMPLICIT (with HPL), ncl=1, npaire=0
real(mytype), dimension(ny) :: aam11,bbm11,ccm11,ddm11,eem11,ggm11,hhm11,wwm11,zzm11 !!TIME IMPLICIT, ncl=1, npaire=1
real(mytype), dimension(ny) :: rrm11,qqm11,vvm11,ssm11 !!TIME IMPLICIT (with HPL), ncl=1, npaire=1
real(mytype), dimension(ny) :: aam0,bbm0,ccm0,ddm0,eem0,ggm0,hhm0,wwm0,zzm0 !!TIME IMPLICIT, ncl=0
real(mytype), dimension(ny) :: rrm0,qqm0,vvm0,ssm0,l1m,l2m,l3m,u1m,u2m,u3m !!TIME IMPLICIT (with HPL), ncl=0
real(mytype), dimension(ny) :: aamt,bbmt,ccmt,ddmt,eemt,ggmt,hhmt,wwmt,zzmt !!TIME IMPLICIT SCALAR, ncl=2
real(mytype), dimension(ny) :: rrmt,qqmt,vvmt,ssmt !!TIME IMPLICIT SCALAR (with HPL), ncl=2
real(mytype), dimension(ny) :: aamt1,bbmt1,ccmt1,ddmt1,eemt1,ggmt1,hhmt1,wwmt1,zzmt1 !!TIME IMPLICIT SCALAR, ncl=1
real(mytype), dimension(ny) :: rrmt1,qqmt1,vvmt1,ssmt1 !!TIME IMPLICIT SCALAR (with HPL), ncl=1
real(mytype), dimension(ny) :: aamt0,bbmt0,ccmt0,ddmt0,eemt0,ggmt0,hhmt0,wwmt0,zzmt0 !!TIME IMPLICIT SCALAR, ncl=0
real(mytype), dimension(ny) :: rrmt0,qqmt0,vvmt0,ssmt0,l1mt,l2mt,l3mt,u1mt,u2mt,u3mt !!TIME IMPLICIT SCALAR (with HPL), ncl=0
!module implicit in x
real(mytype), dimension(nx) :: xaam0,xbbm0,xccm0,xddm0,xeem0,xggm0,xhhm0,xwwm0,xzzm0 !!TIME IMPLICIT, ncl=0
real(mytype), dimension(nx) :: xrrm0,xqqm0,xvvm0,xssm0,xl1m,xl2m,xl3m,xu1m,xu2m,xu3m !!TIME IMPLICIT (with HPL), ncl=0
!module implicit in z
real(mytype), dimension(nz) :: zaam0,zbbm0,zccm0,zddm0,zeem0,zggm0,zhhm0,zwwm0,zzzm0 !!TIME IMPLICIT, ncl=0
real(mytype), dimension(nz) :: zrrm0,zqqm0,zvvm0,zssm0,zl1m,zl2m,zl3m,zu1m,zu2m,zu3m !!TIME IMPLICIT (with HPL), ncl=0

!module pressure
real(mytype), save, allocatable, dimension(:,:) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn
real(mytype), save, allocatable, dimension(:,:) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn
real(mytype), save, allocatable, dimension(:,:) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn

!module solid_body
integer,parameter           :: nxfin=(nx-1)*10+1,nyfin=ny*10,nzfin=(nz-1)*10+1
integer,dimension(ny,nz)    :: nobjx
integer,dimension(nx,nz)    :: nobjy
integer,dimension(nx,ny)    :: nobjz
real(mytype),dimension(20,ny,nz) :: xi,xf
real(mytype),dimension(20,nx,nz) :: yi,yf
real(mytype),dimension(20,nx,ny) :: zi,zf


!module inflow
real(mytype), save, allocatable, dimension(:,:) :: bxx1,bxy1,bxz1,bxxn,bxyn,bxzn,bxo,byo,bzo
real(mytype), save, allocatable, dimension(:,:) :: byx1,byy1,byz1,byxn,byyn,byzn
real(mytype), save, allocatable, dimension(:,:) :: bzx1,bzy1,bzz1,bzxn,bzyn,bzzn

!module derpres
real(mytype),dimension(nxm) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
     cwxp6,csx6,cwx6,cifx6,cicx6,cisx6   
real(mytype),dimension(nxm) :: cibx6,cifxp6,cisxp6,ciwx6
real(mytype),dimension(nx) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
    cwi6,cifi6,cici6,cibi6,cifip6  
real(mytype),dimension(nx) :: cisip6,ciwip6,cisi6,ciwi6 
real(mytype),dimension(nym) :: cfy6,ccy6,cby6,cfyp6,csyp6,cwyp6,csy6 
real(mytype),dimension(nym) :: cwy6,cify6,cicy6,ciby6,cifyp6,cisyp6,&
     ciwyp6,cisy6,ciwy6 
real(mytype),dimension(ny) :: cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,&
     csi6y,cwi6y,cifi6y,cici6y  
real(mytype),dimension(ny) :: cibi6y,cifip6y,cisip6y,ciwip6y,cisi6y,ciwi6y  
real(mytype),dimension(nzm) :: cfz6,ccz6,cbz6,cfzp6,cszp6,cwzp6,csz6 
real(mytype),dimension(nzm) :: cwz6,cifz6,cicz6,cibz6,cifzp6,ciszp6,&
     ciwzp6,cisz6,ciwz6 
real(mytype),dimension(nz) :: cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,&
     csi6z,cwi6z,cifi6z,cici6z  
real(mytype),dimension(nz) :: cibi6z,cifip6z,cisip6z,ciwip6z,cisi6z,ciwi6z 

!module waves
complex(mytype), dimension(nz/2+1) :: zkz,zk2,ezs
complex(mytype), dimension(ny) :: yky,yk2,eys	
complex(mytype), dimension(nx) :: xkx,xk2,exs

!module mesh
real(mytype), dimension(ny) :: ppy,pp2y,pp4y
real(mytype), dimension(ny) :: ppyi,pp2yi,pp4yi
real(mytype), dimension(ny) :: yp,ypi
real(mytype), dimension(ny) :: yeta,yetai
real(mytype) :: alpha,beta
end module variables

module param

use decomp_2d, only : mytype

  integer, save :: nclx,ncly,nclz
  integer, save :: ifft, ivirt,istret,iforc_entree,iturb
  integer, save :: itype, iskew, iin, nscheme, ifirst, ilast, iles
  integer, save :: isave,ilit,idebmod, imodulo, idemarre, icommence, irecord
  integer, save :: iscalar
  integer, save :: nxboite, istat,iread,iadvance_time 
  real(mytype), save :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
  real(mytype), save :: dt,xnu,noise,noise1,pi,twopi,u1,u2,sc
  real(mytype), save :: t,xxk1,xxk2
  integer, save :: itr,itime
  character, save :: filesauve*80, filenoise*80, &
       nchamp*80,filepath*80, fileturb*80, filevisu*80 
  real(mytype), dimension(5), save :: adt,bdt,cdt,gdt
  integer, save :: iimplicit !!TIME IMPLICIT
  real(mytype), save :: xcst, xcst_pr !!TIME IMPLICIT
  !!
  !! Robin boundary condition on temperature
  !! alpha * T + beta * dT/dn = g
  !! alpha=1, beta=0 is dirichlet
  !! alpha=0, beta=1 is neumann
  !! 
  !! WARNING ATTENTION ACHTUNG WARNING ATTENTION ACHTUNG
  !!
  !! beta is the coefficient for NORMAL derivative :
  !!
  !! alpha_0*T(0) - beta_0*dTdy(0)=g_0
  !! alpha_n*T(L) + beta_n*dTdy(L)=g_n
  !!
  !! WARNING ATTENTION ACHTUNG WARNING ATTENTION ACHTUNG
  !!
  real(mytype), save :: alpha_0, beta_0, g_0, alpha_n, beta_n, g_n
end module param

module IBM

use decomp_2d, only : mytype

  real(mytype) :: cex,cey,cez,ra
end module IBM

module derivX

use decomp_2d, only : mytype

  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6,dicix6
  real(mytype) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(mytype) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(mytype) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(mytype) :: asmx,alsaix,asix,bsix,csix,alsa3x,as3x,bs3x,alsatx,astx,bstx 
  real(mytype) :: alsaixt,asixt,bsixt,csixt
end module derivX

module derivY

use decomp_2d, only : mytype

  real(mytype) :: alcaiy6,aciy6,bciy6
  real(mytype) :: ailcaiy6,aiciy6,biciy6,ciciy6,diciy6
  real(mytype) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(mytype) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(mytype) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(mytype) :: asmy,alsajy,asjy,bsjy,csjy,alsa3y,as3y,bs3y,alsaty,asty,bsty 
  real(mytype) :: alsajyt,asjyt,bsjyt,csjyt
end module derivY

module derivZ

use decomp_2d, only : mytype

  real(mytype) :: alcaiz6,aciz6,bciz6
  real(mytype) :: ailcaiz6,aiciz6,biciz6,ciciz6,diciz6
  real(mytype) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(mytype) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(mytype) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(mytype) :: asmz,alsakz,askz,bskz,cskz,alsa3z,as3z,bs3z,alsatz,astz,bstz
  real(mytype) :: alsakzt,askzt,bskzt,cskzt
end module derivZ


module parfiX

use decomp_2d, only : mytype

  real(mytype) :: fia1x, fib1x, fic1x, fid1x, fie1x, fia2x, fib2x, fic2x, fid2x
  real(mytype) :: fie2x, fia3x, fib3x, fic3x, fid3x, fie3x, fianx, fibnx, ficnx, fidnx
  real(mytype) :: fienx, fiamx, fibmx, ficmx, fidmx, fiemx, fiapx, fibpx, ficpx, fidpx
  real(mytype) :: fiepx, fiaix, fibix, ficix, fidix, fialx, fibex, fih1x, fih2x, fih3x,fih4x 
end module parfiX
!
module parfiY

use decomp_2d, only : mytype

  real(mytype) :: fia1y, fib1y, fic1y, fid1y, fie1y, fia2y, fib2y, fic2y, fid2y
  real(mytype) :: fie2y, fia3y, fib3y, fic3y, fid3y, fie3y, fiany, fibny, ficny, fidny
  real(mytype) :: fieny, fiamy, fibmy, ficmy, fidmy, fiemy, fiapy, fibpy, ficpy, fidpy
  real(mytype) :: fiepy, fiaiy, fibiy, ficiy, fidiy, fialy, fibey, fih1y, fih2y, fih3y,fih4y 
end module parfiY

module parfiZ

use decomp_2d, only : mytype

  real(mytype) :: fia1z, fib1z, fic1z, fid1z, fie1z, fia2z, fib2z, fic2z, fid2z
  real(mytype) :: fie2z, fia3z, fib3z, fic3z, fid3z, fie3z, fianz, fibnz, ficnz, fidnz
  real(mytype) :: fienz, fiamz, fibmz, ficmz, fidmz, fiemz, fiapz, fibpz, ficpz, fidpz
  real(mytype) :: fiepz, fiaiz, fibiz, ficiz, fidiz, fialz, fibez, fih1z, fih2z, fih3z,fih4z 
end module parfiZ



