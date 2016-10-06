
module user_stats

use decomp_2d, only : mytype, DECOMP_INFO

implicit none

interface update_user_stats
   module procedure update_user_stats_non_phi
   module procedure update_user_stats_oui_phi
end interface update_user_stats

logical, save :: bool_user_stat
integer, save :: beg_stat
TYPE(DECOMP_INFO) :: decomp_user_stats
integer, save :: xst1, xst2, xst3
integer, save :: xen1, xen2, xen3

! Ordre 1
real(mytype), save, allocatable, dimension(:,:,:) :: um,vm,wm
real(mytype), save, allocatable, dimension(:,:,:) :: dpdxm, dpdym, dpdzm
real(mytype), save, allocatable, dimension(:,:,:) :: dudym, dvdym, dwdym
real(mytype), save, allocatable, dimension(:,:,:) :: dudyym, dvdyym, dwdyym

! Ordre 2
real(mytype), save, allocatable, dimension(:,:,:) :: uum,vvm,wwm
real(mytype), save, allocatable, dimension(:,:,:) :: uvm,uwm,vwm
real(mytype), save, allocatable, dimension(:,:,:) :: udpdxm, udpdym, udpdzm
real(mytype), save, allocatable, dimension(:,:,:) :: vdpdxm, vdpdym, vdpdzm
real(mytype), save, allocatable, dimension(:,:,:) :: wdpdxm, wdpdym, wdpdzm

! Ordre 3
real(mytype), save, allocatable, dimension(:,:,:) :: uuum, uvvm, uwwm
real(mytype), save, allocatable, dimension(:,:,:) :: vuum, vvvm, vwwm
real(mytype), save, allocatable, dimension(:,:,:) :: wuum, wvvm, wwwm
real(mytype), save, allocatable, dimension(:,:,:) :: uvwm

! Epsilon pour tous les Rij
real(mytype), save, allocatable, dimension(:,:,:) :: dudx2m, dudy2m, dudz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dvdx2m, dvdy2m, dvdz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dwdx2m, dwdy2m, dwdz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dudvdxm, dudwdxm, dvdwdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudvdym, dudwdym, dvdwdym
real(mytype), save, allocatable, dimension(:,:,:) :: dudvdzm, dudwdzm, dvdwdzm

! Budget epsilon
! P0 : Production par convection avec vitesse moyenne, ok, cf ci-dessus
! P1 : Production par gradient vitesse moyenne : dUidXj*dUkdXj, ok, cf ci-dessus
! P2 : Production mixte, on a déjà les dudx², dudy², dudz², ajout termes croisés :
  ! uu, vv, ww
real(mytype), save, allocatable, dimension(:,:,:) :: dudxdudym, dudxdudzm, dudydudzm
real(mytype), save, allocatable, dimension(:,:,:) :: dvdxdvdym, dvdxdvdzm, dvdydvdzm
real(mytype), save, allocatable, dimension(:,:,:) :: dwdxdwdym, dwdxdwdzm, dwdydwdzm
  ! uv
real(mytype), save, allocatable, dimension(:,:,:) :: dudxdvdym, dudydvdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudxdvdzm, dudzdvdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudydvdzm, dudzdvdym
  ! uw
real(mytype), save, allocatable, dimension(:,:,:) :: dudxdwdym, dudydwdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudxdwdzm, dudzdwdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudydwdzm, dudzdwdym
  ! vw
real(mytype), save, allocatable, dimension(:,:,:) :: dvdxdwdym, dvdydwdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dvdxdwdzm, dvdzdwdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dvdydwdzm, dvdzdwdym

! P3 : Production du gradient (hessienne de la vitesse), ajout des termes manquants :
real(mytype), save, allocatable, dimension(:,:,:) :: vdudxm, vdudym, vdudzm
real(mytype), save, allocatable, dimension(:,:,:) :: vdvdxm, vdvdym, vdvdzm
real(mytype), save, allocatable, dimension(:,:,:) :: vdwdxm, vdwdym, vdwdzm
! P4 : Production turbulente, produits triples
  ! EpsXX
real(mytype), save, allocatable, dimension(:,:,:) :: dudx_dududxm, dudx_dududym, dudx_dududzm!dudx_eps11
real(mytype), save, allocatable, dimension(:,:,:) :: dudy_dudvdxm, dudy_dudvdym, dudy_dudvdzm!dudy_eps12
real(mytype), save, allocatable, dimension(:,:,:) :: dudz_dudwdxm, dudz_dudwdym, dudz_dudwdzm!dudz_eps13
  ! EpsYY
real(mytype), save, allocatable, dimension(:,:,:) :: dvdx_dvdudxm, dvdx_dvdudym, dvdx_dvdudzm!dvdx_eps12
real(mytype), save, allocatable, dimension(:,:,:) :: dvdy_dvdvdxm, dvdy_dvdvdym, dvdy_dvdvdzm!dvdy_eps22
real(mytype), save, allocatable, dimension(:,:,:) :: dvdz_dvdwdxm, dvdz_dvdwdym, dvdz_dvdwdzm!dvdz_eps23
  ! EpsZZ
real(mytype), save, allocatable, dimension(:,:,:) :: dwdx_dwdudxm, dwdx_dwdudym, dwdx_dwdudzm!dwdx_eps13
real(mytype), save, allocatable, dimension(:,:,:) :: dwdy_dwdvdxm, dwdy_dwdvdym, dwdy_dwdvdzm!dwdy_eps23
real(mytype), save, allocatable, dimension(:,:,:) :: dwdz_dwdwdxm, dwdz_dwdwdym, dwdz_dwdwdzm!dwdz_eps33
  ! EpsXY
real(mytype), save, allocatable, dimension(:,:,:) :: dudx_dudvdxm, dudx_dudvdym, dudx_dudvdzm!dudx_eps21
real(mytype), save, allocatable, dimension(:,:,:) :: dudy_dvdvdxm, dudy_dvdvdym, dudy_dvdvdzm!dudy_eps22
real(mytype), save, allocatable, dimension(:,:,:) :: dudz_dvdwdxm, dudz_dvdwdym, dudz_dvdwdzm!dudz_eps23
real(mytype), save, allocatable, dimension(:,:,:) :: dvdx_dududxm, dvdx_dududym, dvdx_dududzm!dvdx_eps11
real(mytype), save, allocatable, dimension(:,:,:) :: dvdy_dudvdxm, dvdy_dudvdym, dvdy_dudvdzm!dvdy_eps12
real(mytype), save, allocatable, dimension(:,:,:) :: dvdz_dudwdxm, dvdz_dudwdym, dvdz_dudwdzm!dvdz_eps13
! Turbulent transport, ajout des termes NON-NULS :
real(mytype), save, allocatable, dimension(:,:,:) :: vdudx2m, vdudy2m, vdudz2m ! EpsXX
real(mytype), save, allocatable, dimension(:,:,:) :: vdvdx2m, vdvdy2m, vdvdz2m ! EpsYY
real(mytype), save, allocatable, dimension(:,:,:) :: vdwdx2m, vdwdy2m, vdwdz2m ! EpsZZ
real(mytype), save, allocatable, dimension(:,:,:) :: vdudvdxm, vdudvdym, vdudvdzm ! EpsXY
! Pressure transport, ajout des termes NON-NULS
real(mytype), save, allocatable, dimension(:,:,:) :: dpdvdxm, dpdvdym, dpdvdzm ! EpsYY
real(mytype), save, allocatable, dimension(:,:,:) :: dpdudxm, dpdudym, dpdudzm ! EpsXY
! Pressure production sur EpsXY
real(mytype), save, allocatable, dimension(:,:,:) :: dpdxdudyxm, dpdydudyym, dpdzdudyzm
real(mytype), save, allocatable, dimension(:,:,:) :: dpdxdvdxxm, dpdydvdxym, dpdzdvdxzm
! Dissipation
  ! uu
real(mytype), save, allocatable, dimension(:,:,:) :: dudxx2m, dudyy2m, dudzz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dudxy2m, dudxz2m, dudyz2m
  ! vv
real(mytype), save, allocatable, dimension(:,:,:) :: dvdxx2m, dvdyy2m, dvdzz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dvdxy2m, dvdxz2m, dvdyz2m
  ! ww
real(mytype), save, allocatable, dimension(:,:,:) :: dwdxx2m, dwdyy2m, dwdzz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dwdxy2m, dwdxz2m, dwdyz2m
  ! uv
real(mytype), save, allocatable, dimension(:,:,:) :: dudvdxx2m, dudvdyy2m, dudvdzz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dudvdxy2m, dudvdxz2m, dudvdyz2m

! Scalaire
! Ordre 1
real(mytype), save, allocatable, dimension(:,:,:) :: phim
real(mytype), save, allocatable, dimension(:,:,:) :: dphidym
real(mytype), save, allocatable, dimension(:,:,:) :: dphidyym
real(mytype), save, allocatable, dimension(:,:,:) :: dphidyyym
! Ordre 2
real(mytype), save, allocatable, dimension(:,:,:) :: uphim, vphim, wphim, phiphim
real(mytype), save, allocatable, dimension(:,:,:) :: phidpdxm, phidpdym, phidpdzm
!
real(mytype), save, allocatable, dimension(:,:,:) :: dphidx2m, dphidy2m, dphidz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dudphidxm, dvdphidxm, dwdphidxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudphidym, dvdphidym, dwdphidym
real(mytype), save, allocatable, dimension(:,:,:) :: dudphidzm, dvdphidzm, dwdphidzm
real(mytype), save, allocatable, dimension(:,:,:) :: udtdxxm, udtdyym, udtdzzm
real(mytype), save, allocatable, dimension(:,:,:) :: vdtdxxm, vdtdyym, vdtdzzm
real(mytype), save, allocatable, dimension(:,:,:) :: wdtdxxm, wdtdyym, wdtdzzm
! Ordre 3
real(mytype), save, allocatable, dimension(:,:,:) :: uphi2m,vphi2m,wphi2m
real(mytype), save, allocatable, dimension(:,:,:) :: phiuum,phivvm,phiwwm
real(mytype), save, allocatable, dimension(:,:,:) :: phiuvm,phiuwm,phivwm
!
! Budget epsilon thermique
! P0 : Production par convection avec vitesse moyenne, OK
! P1 : Production par gradient vitesse moyenne : dUidXj*dTdXj, ok, cf ci-dessus
! P2 : Production mixte, gradU*gradT, ajout des termes croisés
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdtdym, dtdxdtdzm, dtdydtdzm ! TT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdudym, dtdydudxm ! uT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdudzm, dtdzdudxm
real(mytype), save, allocatable, dimension(:,:,:) :: dtdydudzm, dtdzdudym
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdvdym, dtdydvdxm ! vT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdvdzm, dtdzdvdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dtdydvdzm, dtdzdvdym
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdwdym, dtdydwdxm ! wT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxdwdzm, dtdzdwdxm
real(mytype), save, allocatable, dimension(:,:,:) :: dtdydwdzm, dtdzdwdym
! P3 : Production du gradient (hessienne vitesse et temp.), ajout termes manquants :
real(mytype), save, allocatable, dimension(:,:,:) :: vdtdxm, vdtdym, vdtdzm
! P4 : Production turbulente, produits triples
  ! EpsXT
real(mytype), save, allocatable, dimension(:,:,:) :: dudx_dudtdxm, dudx_dudtdym, dudx_dudtdzm!dudx_epsxT
real(mytype), save, allocatable, dimension(:,:,:) :: dudy_dvdtdxm, dudy_dvdtdym, dudy_dvdtdzm!dudy_epsyT
real(mytype), save, allocatable, dimension(:,:,:) :: dudz_dwdtdxm, dudz_dwdtdym, dudz_dwdtdzm!dudz_epszT
real(mytype), save, allocatable, dimension(:,:,:) ::               dtdx_dududym, dtdx_dududzm!dtdx_epsxx
real(mytype), save, allocatable, dimension(:,:,:) :: dtdy_dudvdxm,               dtdy_dudvdzm!dtdy_epsxy
real(mytype), save, allocatable, dimension(:,:,:) :: dtdz_dudwdxm, dtdz_dudwdym              !dtdz_epsxz
  ! EpxYT
real(mytype), save, allocatable, dimension(:,:,:) :: dvdx_dudtdxm, dvdx_dudtdym, dvdx_dudtdzm!dvdx_epsxT
real(mytype), save, allocatable, dimension(:,:,:) :: dvdy_dvdtdxm, dvdy_dvdtdym, dvdy_dvdtdzm!dvdy_epsyT
real(mytype), save, allocatable, dimension(:,:,:) :: dvdz_dwdtdxm, dvdz_dwdtdym, dvdz_dwdtdzm!dvdz_epszT
real(mytype), save, allocatable, dimension(:,:,:) ::               dtdx_dudvdym, dtdx_dudvdzm!dtdx_epsxy
real(mytype), save, allocatable, dimension(:,:,:) :: dtdy_dvdvdxm,               dtdy_dvdvdzm!dtdy_epsyy
real(mytype), save, allocatable, dimension(:,:,:) :: dtdz_dvdwdxm, dtdz_dvdwdym              !dtdz_epsyz
  ! EpsTT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdx_dudtdxm, dtdx_dudtdym, dtdx_dudtdzm!dtdx_epsxT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdy_dvdtdxm, dtdy_dvdtdym, dtdy_dvdtdzm!dtdy_epsyT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdz_dwdtdxm, dtdz_dwdtdym, dtdz_dwdtdzm!dtdz_epszT
! Turbulent transport
real(mytype), save, allocatable, dimension(:,:,:) :: vdudtdxm, vdudtdym, vdudtdzm!v_epsxT
real(mytype), save, allocatable, dimension(:,:,:) :: vdvdtdxm, vdvdtdym, vdvdtdzm!v_epsyT
real(mytype), save, allocatable, dimension(:,:,:) :: vdtdx2m, vdtdy2m, vdtdz2m!v_epsTT
! Pressure transport + production
real(mytype), save, allocatable, dimension(:,:,:) :: dpdtdxm, dpdtdym, dpdtdzm
real(mytype), save, allocatable, dimension(:,:,:) :: dpdtdxxm, dpdtdxym, dpdtdxzm ! EpsXT
real(mytype), save, allocatable, dimension(:,:,:) :: dpdtdyxm, dpdtdyym, dpdtdyzm ! EpsYT
! Dissipation
  ! xT
real(mytype), save, allocatable, dimension(:,:,:) :: dudtdxxm, dudtdyym, dudtdzzm
real(mytype), save, allocatable, dimension(:,:,:) :: dudtdxym, dudtdxzm, dudtdyzm
real(mytype), save, allocatable, dimension(:,:,:) :: dudx_dtdxxxm, dudx_dtdyxxm, dudx_dtdzxxm
real(mytype), save, allocatable, dimension(:,:,:) :: dudy_dtdxyym, dudy_dtdyyym, dudy_dtdzyym
real(mytype), save, allocatable, dimension(:,:,:) :: dudz_dtdxzzm, dudz_dtdyzzm, dudz_dtdzzzm
  ! yT
real(mytype), save, allocatable, dimension(:,:,:) :: dvdtdxxm, dvdtdyym, dvdtdzzm
real(mytype), save, allocatable, dimension(:,:,:) :: dvdtdxym, dvdtdxzm, dvdtdyzm
real(mytype), save, allocatable, dimension(:,:,:) :: dvdx_dtdxxxm, dvdx_dtdyxxm, dvdx_dtdzxxm
real(mytype), save, allocatable, dimension(:,:,:) :: dvdy_dtdxyym, dvdy_dtdyyym, dvdy_dtdzyym
real(mytype), save, allocatable, dimension(:,:,:) :: dvdz_dtdxzzm, dvdz_dtdyzzm, dvdz_dtdzzzm
  ! TT
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxx2m, dtdyy2m, dtdzz2m
real(mytype), save, allocatable, dimension(:,:,:) :: dtdxy2m, dtdxz2m, dtdyz2m

! Pression
real(mytype), save, allocatable, dimension(:,:,:) :: presm, pres2m
!real(mytype), save, allocatable, dimension(:,:,:) :: pp3m, pp32m

contains

subroutine allocate_user_stats(nx_user, ny_user, nz_user, phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : decomp_info_init, DECOMP_INFO, get_decomp_info
use param, only : iscalar

implicit none

integer, intent(in) :: nx_user, ny_user, nz_user
TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
TYPE(DECOMP_INFO) :: decomp_main

call decomp_info_init(nx_user, ny_user, nz_user, decomp_user_stats)

xst1=decomp_user_stats%xst(1)
xst2=decomp_user_stats%xst(2)
xst3=decomp_user_stats%xst(3)
xen1=decomp_user_stats%xst(1)
xen2=decomp_user_stats%xen(2)
xen3=decomp_user_stats%xen(3)

allocate(um     (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vm     (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wm     (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdxm  (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdym  (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdzm  (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudym  (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdym  (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdym  (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudyym (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdyym (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdyym (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uum    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vvm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wwm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uvm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uwm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vwm    (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(udpdxm (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(udpdym (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(udpdzm (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdpdxm (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdpdym (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdpdzm (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wdpdxm (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wdpdym (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wdpdzm (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uuum   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uvvm   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uwwm   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vuum   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vvvm   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vwwm   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wuum   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wvvm   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(wwwm   (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(uvwm   (xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudx2m (xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsXX
allocate(dudy2m (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudz2m (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdx2m (xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsYY
allocate(dvdy2m (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdz2m (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdx2m (xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsZZ
allocate(dwdy2m (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdz2m (xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsXY
allocate(dudvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsXZ
allocate(dudwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsYZ
allocate(dvdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudxdudym(xst1:xen1,xst2:xen2,xst3:xen3)) ! uu
allocate(dudxdudzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudydudzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdxdvdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! vv
allocate(dvdxdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdydvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdxdwdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! ww
allocate(dwdxdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdydwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudxdvdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! uv
allocate(dudydvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudxdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudzdvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudydvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudzdvdym(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudxdwdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! uw
allocate(dudydwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudxdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudzdwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudydwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudzdwdym(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dvdxdwdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! vw
allocate(dvdydwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdxdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdzdwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdydwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdzdwdym(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(vdudxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdudym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdudzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudx_dududxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dudx_eps11
allocate(dudx_dududym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudx_dududzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudy_dudvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dudy_eps12
allocate(dudy_dudvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudy_dudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudz_dudwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dudz_eps13
allocate(dudz_dudwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudz_dudwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dvdx_dvdudxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dvdx_eps12
allocate(dvdx_dvdudym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdx_dvdudzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdy_dvdvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dvdy_eps22
allocate(dvdy_dvdvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdy_dvdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdz_dvdwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dvdz_eps23
allocate(dvdz_dvdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdz_dvdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dwdx_dwdudxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dwdx_eps13
allocate(dwdx_dwdudym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdx_dwdudzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdy_dwdvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dwdy_eps23
allocate(dwdy_dwdvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdy_dwdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdz_dwdwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dwdz_eps33
allocate(dwdz_dwdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdz_dwdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudx_dudvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dudx_eps21
allocate(dudx_dudvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudx_dudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudy_dvdvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dudy_eps22
allocate(dudy_dvdvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudy_dvdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudz_dvdwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dudz_eps23
allocate(dudz_dvdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudz_dvdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdx_dududxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dvdx_eps11
allocate(dvdx_dududym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdx_dududzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdy_dudvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dvdy_eps12
allocate(dvdy_dudvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdy_dudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdz_dudwdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! dvdz_eps13
allocate(dvdz_dudwdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdz_dudwdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(vdudx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsXX
allocate(vdudy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdudz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdvdx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsYY
allocate(vdvdy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdvdz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdwdx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsZZ
allocate(vdwdy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdwdz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdudvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsXY
allocate(vdudvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(vdudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dpdvdxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsYY
allocate(dpdvdym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdudxm(xst1:xen1,xst2:xen2,xst3:xen3)) ! EpsXY
allocate(dpdudym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdudzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dpdxdudyxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdydudyym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdzdudyzm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdxdvdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdydvdxym(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dpdzdvdxzm(xst1:xen1,xst2:xen2,xst3:xen3))

allocate(dudxx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! uu
allocate(dudyy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudzz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudxy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudxz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudyz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdxx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! vv
allocate(dvdyy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdzz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdxy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdxz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dvdyz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdxx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! ww
allocate(dwdyy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdzz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdxy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdxz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dwdyz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdxx2m(xst1:xen1,xst2:xen2,xst3:xen3)) ! uv
allocate(dudvdyy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdzz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdxy2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdxz2m(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(dudvdyz2m(xst1:xen1,xst2:xen2,xst3:xen3))

um=0.; vm=0.; wm=0.;
dpdxm=0.; dpdym=0.; dpdzm=0.;
dudym=0.; dvdym=0.; dwdym=0.;
dudyym=0.; dvdyym=0.; dwdyym=0.;
uum=0.; vvm=0.; wwm=0.;
uvm=0.; uwm=0.; vwm=0.;
udpdxm=0.; udpdym=0.; udpdzm=0.;
vdpdxm=0.; vdpdym=0.; vdpdzm=0.;
wdpdxm=0.; wdpdym=0.; wdpdzm=0.;
uuum=0.; uvvm=0.; uwwm=0.;
vuum=0.; vvvm=0.; vwwm=0.;
wuum=0.; wvvm=0.; wwwm=0.
uvwm=0.;

dudx2m=0.; dudy2m=0.; dudz2m=0.;
dvdx2m=0.; dvdy2m=0.; dvdz2m=0.;
dwdx2m=0.; dwdy2m=0.; dwdz2m=0.;
dudvdxm=0.; dudwdxm=0.; dvdwdxm=0.;
dudvdym=0.; dudwdym=0.; dvdwdym=0.;
dudvdzm=0.; dudwdzm=0.; dvdwdzm=0.;

dudxdudym=0.; dudxdudzm=0.; dudydudzm=0.; ! uu
dvdxdvdym=0.; dvdxdvdzm=0.; dvdydvdzm=0.; ! vv
dwdxdwdym=0.; dwdxdwdzm=0.; dwdydwdzm=0.; ! ww

dudxdvdym=0.; dudydvdxm=0.; ! uv
dudxdvdzm=0.; dudzdvdxm=0.;
dudydvdzm=0.; dudzdvdym=0.;

dudxdwdym=0.; dudydwdxm=0.; ! uw
dudxdwdzm=0.; dudzdwdxm=0.;
dudydwdzm=0.; dudzdwdym=0.;

dvdxdwdym=0.; dvdydwdxm=0.; ! vw
dvdxdwdzm=0.; dvdzdwdxm=0.;
dvdydwdzm=0.; dvdzdwdym=0.;

vdudxm=0.; vdvdxm=0.; vdwdxm=0.;
vdudym=0.; vdvdym=0.; vdwdym=0.;
vdudzm=0.; vdvdzm=0.; vdwdzm=0.;

dudx_dududxm=0.; dudx_dududym=0.; dudx_dududzm=0.; ! dudx_eps11
dudy_dudvdxm=0.; dudy_dudvdym=0.; dudy_dudvdzm=0.; ! dudy_eps12
dudz_dudwdxm=0.; dudz_dudwdym=0.; dudz_dudwdzm=0.; ! dudz_eps13

dvdx_dvdudxm=0.; dvdx_dvdudym=0.; dvdx_dvdudzm=0.; ! dvdx_eps12
dvdy_dvdvdxm=0.; dvdy_dvdvdym=0.; dvdy_dvdvdzm=0.; ! dvdy_eps22
dvdz_dvdwdxm=0.; dvdz_dvdwdym=0.; dvdz_dvdwdzm=0.; ! dvdz_eps23

dwdx_dwdudxm=0.; dwdx_dwdudym=0.; dwdx_dwdudzm=0.; ! dwdx_eps13
dwdy_dwdvdxm=0.; dwdy_dwdvdym=0.; dwdy_dwdvdzm=0.; ! dwdy_eps23
dwdz_dwdwdxm=0.; dwdz_dwdwdym=0.; dwdz_dwdwdzm=0.; ! dwdz_eps33

dudx_dudvdxm=0.; dudx_dudvdym=0.; dudx_dudvdzm=0.; !dudx_eps21
dudy_dvdvdxm=0.; dudy_dvdvdym=0.; dudy_dvdvdzm=0.; !dudy_eps22
dudz_dvdwdxm=0.; dudz_dvdwdym=0.; dudz_dvdwdzm=0.; !dudz_eps23
dvdx_dududxm=0.; dvdx_dududym=0.; dvdx_dududzm=0.; !dvdx_eps11
dvdy_dudvdxm=0.; dvdy_dudvdym=0.; dvdy_dudvdzm=0.; !dvdy_eps12
dvdz_dudwdxm=0.; dvdz_dudwdym=0.; dvdz_dudwdzm=0.; !dvdz_eps13

vdudx2m=0.; vdudy2m=0.; vdudz2m=0.; ! EpsXX
vdvdx2m=0.; vdvdy2m=0.; vdvdz2m=0.; ! EpsYY
vdwdx2m=0.; vdwdy2m=0.; vdwdz2m=0.; ! EpsZZ
vdudvdxm=0.; vdudvdym=0.; vdudvdzm=0.; ! EpsXY

dpdvdxm=0.; dpdvdym=0.; dpdvdzm=0.; ! EpsYY
dpdudxm=0.; dpdudym=0.; dpdudzm=0.; ! EpsXY

dpdxdudyxm=0.; dpdydudyym=0.; dpdzdudyzm=0.;
dpdxdvdxxm=0.; dpdydvdxym=0.; dpdzdvdxzm=0.;

dudxx2m=0.; dudyy2m=0.; dudzz2m=0.; ! uu
dudxy2m=0.; dudxz2m=0.; dudyz2m=0.;
dvdxx2m=0.; dvdyy2m=0.; dvdzz2m=0.; ! vv
dvdxy2m=0.; dvdxz2m=0.; dvdyz2m=0.;
dwdxx2m=0.; dwdyy2m=0.; dwdzz2m=0.; ! ww
dwdxy2m=0.; dwdxz2m=0.; dwdyz2m=0.;
dudvdxx2m=0.; dudvdyy2m=0.; dudvdzz2m=0.; ! uv
dudvdxy2m=0.; dudvdxz2m=0.; dudvdyz2m=0.;

! Scalaire
if (iscalar.eq.1) then
  ! Ordre 1
  allocate(phim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dphidym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dphidyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dphidyyym(xst1:xen1,xst2:xen2,xst3:xen3))
  ! Ordre 2
  allocate(uphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(wphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phiphim(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phidpdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phidpdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phidpdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dphidx2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dphidy2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dphidz2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudphidxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdphidxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dwdphidxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudphidym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdphidym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dwdphidym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudphidzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdphidzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dwdphidzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(udtdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(udtdyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(udtdzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(wdtdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(wdtdyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(wdtdzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  ! Ordre 3
  allocate(uphi2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vphi2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(wphi2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phiuum(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phivvm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phiwwm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phiuvm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phiuwm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(phivwm(xst1:xen1,xst2:xen2,xst3:xen3))
  ! P2
  allocate(dtdxdtdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! TT
  allocate(dtdxdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdydtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxdudym(xst1:xen1,xst2:xen2,xst3:xen3)) ! uT
  allocate(dtdydudxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxdudzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzdudxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdydudzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzdudym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxdvdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! vT
  allocate(dtdydvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzdvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdydvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzdvdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxdwdym(xst1:xen1,xst2:xen2,xst3:xen3)) ! wT
  allocate(dtdydwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxdwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzdwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdydwdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(vdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dudx_dudtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudx_dudtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudx_dudtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudy_dvdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudy_dvdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudy_dvdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudz_dwdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudz_dwdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudz_dwdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dtdx_dududym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdx_dududzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dudvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dudwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dudwdym(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dvdx_dudtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdx_dudtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdx_dudtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdy_dvdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdy_dvdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdy_dvdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdz_dwdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdz_dwdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdz_dwdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dtdx_dudvdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdx_dudvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dvdvdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dvdvdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dvdwdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dvdwdym(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dtdx_dudtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdx_dudtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdx_dudtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dvdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dvdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdy_dvdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dwdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dwdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdz_dwdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(vdudtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdudtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdudtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdvdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdvdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdvdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdx2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdy2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(vdtdz2m(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dpdtdxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdxym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdxzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdyxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dpdtdyzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dudtdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudtdyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudtdzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudtdxym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudtdxzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudtdyzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dudx_dtdxxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudx_dtdyxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudx_dtdzxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudy_dtdxyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudy_dtdyyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudy_dtdzyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudz_dtdxzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudz_dtdyzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dudz_dtdzzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dvdtdxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdtdyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdtdzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdtdxym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdtdxzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdtdyzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dvdx_dtdxxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdx_dtdyxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdx_dtdzxxm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdy_dtdxyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdy_dtdyyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdy_dtdzyym(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdz_dtdxzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdz_dtdyzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dvdz_dtdzzzm(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  allocate(dtdxx2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdyy2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdzz2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxy2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdxz2m(xst1:xen1,xst2:xen2,xst3:xen3))
  allocate(dtdyz2m(xst1:xen1,xst2:xen2,xst3:xen3))
  !
  phim=0.;
  dphidym=0.;
  dphidyym=0.;
  dphidyyym=0.;
  uphim=0.; vphim=0.; wphim=0.; phiphim=0.;
  phidpdxm=0.; phidpdym=0.; phidpdzm=0.;
  dphidx2m=0.; dphidy2m=0.; dphidz2m=0.;
  dudphidxm=0.; dvdphidxm=0.; dwdphidxm=0.;
  dudphidym=0.; dvdphidym=0.; dwdphidym=0.;
  dudphidzm=0.; dvdphidzm=0.; dwdphidzm=0.;
  udtdxxm=0.; udtdyym=0.; udtdzzm=0.;
  vdtdxxm=0.; vdtdyym=0.; vdtdzzm=0.;
  wdtdxxm=0.; wdtdyym=0.; wdtdzzm=0.;
  uphi2m=0.; vphi2m=0.; wphi2m=0.;
  phiuum=0.; phivvm=0.; phiwwm=0.;
  phiuvm=0.; phiuwm=0.; phivwm=0.;
  !
  dtdxdtdym=0.; dtdxdtdzm=0.; dtdydtdzm=0.; ! TT
  dtdxdudym=0.; dtdydudxm=0.; ! uT
  dtdxdudzm=0.; dtdzdudxm=0.;
  dtdydudzm=0.; dtdzdudym=0.;
  dtdxdvdym=0.; dtdydvdxm=0.; ! vT
  dtdxdvdzm=0.; dtdzdvdxm=0.;
  dtdydvdzm=0.; dtdzdvdym=0.;
  dtdxdwdym=0.; dtdydwdxm=0.; ! wT
  dtdxdwdzm=0.; dtdzdwdxm=0.;
  dtdydwdzm=0.; dtdzdwdym=0.;
  !
  vdtdxm=0.; vdtdym=0.; vdtdzm=0.;
  !
  dudx_dudtdxm=0.; dudx_dudtdym=0.; dudx_dudtdzm=0.;
  dudy_dvdtdxm=0.; dudy_dvdtdym=0.; dudy_dvdtdzm=0.;
  dudz_dwdtdxm=0.; dudz_dwdtdym=0.; dudz_dwdtdzm=0.;
  dtdx_dududym=0.; dtdx_dududzm=0.;
  dtdy_dudvdxm=0.; dtdy_dudvdzm=0.;
  dtdz_dudwdxm=0.; dtdz_dudwdym=0.;
  !
  dvdx_dudtdxm=0.; dvdx_dudtdym=0.; dvdx_dudtdzm=0.;
  dvdy_dvdtdxm=0.; dvdy_dvdtdym=0.; dvdy_dvdtdzm=0.;
  dvdz_dwdtdxm=0.; dvdz_dwdtdym=0.; dvdz_dwdtdzm=0.;
  dtdx_dudvdym=0.; dtdx_dudvdzm=0.;
  dtdy_dvdvdxm=0.; dtdy_dvdvdzm=0.;
  dtdz_dvdwdxm=0.; dtdz_dvdwdym=0.;
  !
  dtdx_dudtdxm=0.; dtdx_dudtdym=0.; dtdx_dudtdzm=0.;
  dtdy_dvdtdxm=0.; dtdy_dvdtdym=0.; dtdy_dvdtdzm=0.;
  dtdz_dwdtdxm=0.; dtdz_dwdtdym=0.; dtdz_dwdtdzm=0.;
  !
  vdudtdxm=0.; vdudtdym=0.; vdudtdzm=0.;
  vdvdtdxm=0.; vdvdtdym=0.; vdvdtdzm=0.;
  vdtdx2m=0.; vdtdy2m=0.; vdtdz2m=0.;
  !
  dpdtdxm=0.; dpdtdym=0.; dpdtdzm=0.;
  dpdtdxxm=0.; dpdtdxym=0.; dpdtdxzm=0.;
  dpdtdyxm=0.; dpdtdyym=0.; dpdtdyzm=0.;
  !
  dudtdxxm=0.; dudtdyym=0.; dudtdzzm=0.;
  dudtdxym=0.; dudtdxzm=0.; dudtdyzm=0.;
  !
  dudx_dtdxxxm=0.; dudx_dtdyxxm=0.; dudx_dtdzxxm=0.;
  dudy_dtdxyym=0.; dudy_dtdyyym=0.; dudy_dtdzyym=0.;
  dudz_dtdxzzm=0.; dudz_dtdyzzm=0.; dudz_dtdzzzm=0.;
  !
  dvdtdxxm=0.; dvdtdyym=0.; dvdtdzzm=0.;
  dvdtdxym=0.; dvdtdxzm=0.; dvdtdyzm=0.;
  !
  dvdx_dtdxxxm=0.; dvdx_dtdyxxm=0.; dvdx_dtdzxxm=0.;
  dvdy_dtdxyym=0.; dvdy_dtdyyym=0.; dvdy_dtdzyym=0.;
  dvdz_dtdxzzm=0.; dvdz_dtdyzzm=0.; dvdz_dtdzzzm=0.;
  !
  dtdxx2m=0.; dtdyy2m=0.; dtdzz2m=0.;
  dtdxy2m=0.; dtdxz2m=0.; dtdyz2m=0.;
  !
endif

! Pression
allocate(presm(xst1:xen1,xst2:xen2,xst3:xen3))
allocate(pres2m(xst1:xen1,xst2:xen2,xst3:xen3))
!
presm=0.; pres2m=0.;
!

end subroutine allocate_user_stats

subroutine update_user_stats_non_phi(ux,uy,uz,dpdx,dpdy,dpdz,&
						dudx,dudy,dudz,&
						dvdx,dvdy,dvdz,&
						dwdx,dwdy,dwdz,&
						dudxx, dudyy, dudzz, dudxy, dudxz, dudyz,&
						dvdxx, dvdyy, dvdzz, dvdxy, dvdxz, dvdyz,&
						dwdxx, dwdyy, dwdzz, dwdxy, dwdxz, dwdyz,&
						     pres,pp3,&
                                                     phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : xsize, nrank, DECOMP_INFO
use param, only : itime

implicit none

TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux,uy,uz,&
						dpdx, dpdy, dpdz,&
						dudx, dudy, dudz,&
						dvdx, dvdy, dvdz,&
						dwdx, dwdy, dwdz,&
						dudxx, dudyy, dudzz, dudxy, dudxz, dudyz,&
						dvdxx, dvdyy, dvdzz, dvdxy, dvdxz, dvdyz,&
						dwdxx, dwdyy, dwdzz, dwdxy, dwdxz, dwdyz,&
						pres
real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),ph3%zst(3):ph3%zen(3)), intent(in) :: pp3

! Local variable
integer :: coeff

coeff=itime-beg_stat
if (coeff.le.0) then
  if (nrank.eq.0) print*,'ERROR : cannot start to compute statistics at ',itime,'. Wait till ',beg_stat
  return
endif

um = (real(coeff-1, kind=mytype)*um + fine_interpol(ux)) / real(coeff, kind=mytype)
vm = (real(coeff-1, kind=mytype)*vm + fine_interpol(uy)) / real(coeff, kind=mytype)
wm = (real(coeff-1, kind=mytype)*wm + fine_interpol(uz)) / real(coeff, kind=mytype)

dpdxm = (real(coeff-1, kind=mytype)*dpdxm + fine_interpol(dpdx)) / real(coeff, kind=mytype)
dpdym = (real(coeff-1, kind=mytype)*dpdym + fine_interpol(dpdy)) / real(coeff, kind=mytype)
dpdzm = (real(coeff-1, kind=mytype)*dpdzm + fine_interpol(dpdz)) / real(coeff, kind=mytype)

dudym = (real(coeff-1, kind=mytype)*dudym + fine_interpol(dudy)) / real(coeff, kind=mytype)
dvdym = (real(coeff-1, kind=mytype)*dvdym + fine_interpol(dvdy)) / real(coeff, kind=mytype)
dwdym = (real(coeff-1, kind=mytype)*dwdym + fine_interpol(dwdy)) / real(coeff, kind=mytype)

dudyym = (real(coeff-1, kind=mytype)*dudyym + fine_interpol(dudyy)) / real(coeff, kind=mytype)
dvdyym = (real(coeff-1, kind=mytype)*dvdyym + fine_interpol(dvdyy)) / real(coeff, kind=mytype)
dwdyym = (real(coeff-1, kind=mytype)*dwdyym + fine_interpol(dwdyy)) / real(coeff, kind=mytype)

uum = (real(coeff-1, kind=mytype)*uum + fine_interpol(ux*ux)) / real(coeff, kind=mytype)
vvm = (real(coeff-1, kind=mytype)*vvm + fine_interpol(uy*uy)) / real(coeff, kind=mytype)
wwm = (real(coeff-1, kind=mytype)*wwm + fine_interpol(uz*uz)) / real(coeff, kind=mytype)

uvm = (real(coeff-1, kind=mytype)*uvm + fine_interpol(ux*uy)) / real(coeff, kind=mytype)
uwm = (real(coeff-1, kind=mytype)*uwm + fine_interpol(ux*uz)) / real(coeff, kind=mytype)
vwm = (real(coeff-1, kind=mytype)*vwm + fine_interpol(uy*uz)) / real(coeff, kind=mytype)

udpdxm = (real(coeff-1, kind=mytype)*udpdxm + fine_interpol(ux*dpdx)) / real(coeff, kind=mytype)
udpdym = (real(coeff-1, kind=mytype)*udpdym + fine_interpol(ux*dpdy)) / real(coeff, kind=mytype)
udpdzm = (real(coeff-1, kind=mytype)*udpdzm + fine_interpol(ux*dpdz)) / real(coeff, kind=mytype)

vdpdxm = (real(coeff-1, kind=mytype)*vdpdxm + fine_interpol(uy*dpdx)) / real(coeff, kind=mytype)
vdpdym = (real(coeff-1, kind=mytype)*vdpdym + fine_interpol(uy*dpdy)) / real(coeff, kind=mytype)
vdpdzm = (real(coeff-1, kind=mytype)*vdpdzm + fine_interpol(uy*dpdz)) / real(coeff, kind=mytype)

wdpdxm = (real(coeff-1, kind=mytype)*wdpdxm + fine_interpol(uz*dpdx)) / real(coeff, kind=mytype)
wdpdym = (real(coeff-1, kind=mytype)*wdpdym + fine_interpol(uz*dpdy)) / real(coeff, kind=mytype)
wdpdzm = (real(coeff-1, kind=mytype)*wdpdzm + fine_interpol(uz*dpdz)) / real(coeff, kind=mytype)

uuum = (real(coeff-1, kind=mytype)*uuum + fine_interpol(ux*ux*ux)) / real(coeff, kind=mytype)
uvvm = (real(coeff-1, kind=mytype)*uvvm + fine_interpol(ux*uy*uy)) / real(coeff, kind=mytype)
uwwm = (real(coeff-1, kind=mytype)*uwwm + fine_interpol(ux*uz*uz)) / real(coeff, kind=mytype)

vuum = (real(coeff-1, kind=mytype)*vuum + fine_interpol(uy*ux*ux)) / real(coeff, kind=mytype)
vvvm = (real(coeff-1, kind=mytype)*vvvm + fine_interpol(uy*uy*uy)) / real(coeff, kind=mytype)
vwwm = (real(coeff-1, kind=mytype)*vwwm + fine_interpol(uy*uz*uz)) / real(coeff, kind=mytype)

wuum = (real(coeff-1, kind=mytype)*wuum + fine_interpol(uz*ux*ux)) / real(coeff, kind=mytype)
wvvm = (real(coeff-1, kind=mytype)*wvvm + fine_interpol(uz*uy*uy)) / real(coeff, kind=mytype)
wwwm = (real(coeff-1, kind=mytype)*wwwm + fine_interpol(uz*uz*uz)) / real(coeff, kind=mytype)

uvwm = (real(coeff-1, kind=mytype)*uvwm + fine_interpol(ux*uy*uz)) / real(coeff, kind=mytype)

dudx2m = (real(coeff-1, kind=mytype)*dudx2m + fine_interpol(dudx*dudx)) / real(coeff, kind=mytype)
dudy2m = (real(coeff-1, kind=mytype)*dudy2m + fine_interpol(dudy*dudy)) / real(coeff, kind=mytype)
dudz2m = (real(coeff-1, kind=mytype)*dudz2m + fine_interpol(dudz*dudz)) / real(coeff, kind=mytype)

dvdx2m = (real(coeff-1, kind=mytype)*dvdx2m + fine_interpol(dvdx*dvdx)) / real(coeff, kind=mytype)
dvdy2m = (real(coeff-1, kind=mytype)*dvdy2m + fine_interpol(dvdy*dvdy)) / real(coeff, kind=mytype)
dvdz2m = (real(coeff-1, kind=mytype)*dvdz2m + fine_interpol(dvdz*dvdz)) / real(coeff, kind=mytype)

dwdx2m = (real(coeff-1, kind=mytype)*dwdx2m + fine_interpol(dwdx*dwdx)) / real(coeff, kind=mytype)
dwdy2m = (real(coeff-1, kind=mytype)*dwdy2m + fine_interpol(dwdy*dwdy)) / real(coeff, kind=mytype)
dwdz2m = (real(coeff-1, kind=mytype)*dwdz2m + fine_interpol(dwdz*dwdz)) / real(coeff, kind=mytype)

dudvdxm= (real(coeff-1, kind=mytype)*dudvdxm+ real(1, kind=mytype)*fine_interpol(dudx*dvdx)) / real(coeff, kind=mytype)
dudwdxm= (real(coeff-1, kind=mytype)*dudwdxm+ real(1, kind=mytype)*fine_interpol(dudx*dwdx)) / real(coeff, kind=mytype)
dvdwdxm= (real(coeff-1, kind=mytype)*dvdwdxm+ real(1, kind=mytype)*fine_interpol(dvdx*dwdx)) / real(coeff, kind=mytype)

dudvdym= (real(coeff-1, kind=mytype)*dudvdym+ real(1, kind=mytype)*fine_interpol(dudy*dvdy)) / real(coeff, kind=mytype)
dudwdym= (real(coeff-1, kind=mytype)*dudwdym+ real(1, kind=mytype)*fine_interpol(dudy*dwdy)) / real(coeff, kind=mytype)
dvdwdym= (real(coeff-1, kind=mytype)*dvdwdym+ real(1, kind=mytype)*fine_interpol(dvdy*dwdy)) / real(coeff, kind=mytype)

dudvdzm= (real(coeff-1, kind=mytype)*dudvdzm+ real(1, kind=mytype)*fine_interpol(dudz*dvdz)) / real(coeff, kind=mytype)
dudwdzm= (real(coeff-1, kind=mytype)*dudwdzm+ real(1, kind=mytype)*fine_interpol(dudz*dwdz)) / real(coeff, kind=mytype)
dvdwdzm= (real(coeff-1, kind=mytype)*dvdwdzm+ real(1, kind=mytype)*fine_interpol(dvdz*dwdz)) / real(coeff, kind=mytype)

dudxdudym = (real(coeff-1, kind=mytype)*dudxdudym + fine_interpol(dudx*dudy)) / real(coeff, kind=mytype) ! uu
dudxdudzm = (real(coeff-1, kind=mytype)*dudxdudzm + fine_interpol(dudx*dudz)) / real(coeff, kind=mytype)
dudydudzm = (real(coeff-1, kind=mytype)*dudydudzm + fine_interpol(dudy*dudz)) / real(coeff, kind=mytype)
dvdxdvdym = (real(coeff-1, kind=mytype)*dvdxdvdym + fine_interpol(dvdx*dvdy)) / real(coeff, kind=mytype) ! vv
dvdxdvdzm = (real(coeff-1, kind=mytype)*dvdxdvdzm + fine_interpol(dvdx*dvdz)) / real(coeff, kind=mytype)
dvdydvdzm = (real(coeff-1, kind=mytype)*dvdydvdzm + fine_interpol(dvdy*dvdz)) / real(coeff, kind=mytype)
dwdxdwdym = (real(coeff-1, kind=mytype)*dwdxdwdym + fine_interpol(dwdx*dwdy)) / real(coeff, kind=mytype) ! ww
dwdxdwdzm = (real(coeff-1, kind=mytype)*dwdxdwdzm + fine_interpol(dwdx*dwdz)) / real(coeff, kind=mytype)
dwdydwdzm = (real(coeff-1, kind=mytype)*dwdydwdzm + fine_interpol(dwdy*dwdz)) / real(coeff, kind=mytype)

dudxdvdym = (real(coeff-1, kind=mytype)*dudxdvdym + fine_interpol(dudx*dvdy)) / real(coeff, kind=mytype) ! uv
dudydvdxm = (real(coeff-1, kind=mytype)*dudydvdxm + fine_interpol(dudy*dvdx)) / real(coeff, kind=mytype)
dudxdvdzm = (real(coeff-1, kind=mytype)*dudxdvdzm + fine_interpol(dudx*dvdz)) / real(coeff, kind=mytype)
dudzdvdxm = (real(coeff-1, kind=mytype)*dudzdvdxm + fine_interpol(dudz*dvdx)) / real(coeff, kind=mytype)
dudydvdzm = (real(coeff-1, kind=mytype)*dudydvdzm + fine_interpol(dudy*dvdz)) / real(coeff, kind=mytype)
dudzdvdym = (real(coeff-1, kind=mytype)*dudzdvdym + fine_interpol(dudz*dvdy)) / real(coeff, kind=mytype)

dudxdwdym = (real(coeff-1, kind=mytype)*dudxdwdym + fine_interpol(dudx*dwdy)) / real(coeff, kind=mytype) ! uw
dudydwdxm = (real(coeff-1, kind=mytype)*dudydwdxm + fine_interpol(dudy*dwdx)) / real(coeff, kind=mytype)
dudxdwdzm = (real(coeff-1, kind=mytype)*dudxdwdzm + fine_interpol(dudx*dwdz)) / real(coeff, kind=mytype)
dudzdwdxm = (real(coeff-1, kind=mytype)*dudzdwdxm + fine_interpol(dudz*dwdx)) / real(coeff, kind=mytype)
dudydwdzm = (real(coeff-1, kind=mytype)*dudydwdzm + fine_interpol(dudy*dwdz)) / real(coeff, kind=mytype)
dudzdwdym = (real(coeff-1, kind=mytype)*dudzdwdym + fine_interpol(dudz*dwdy)) / real(coeff, kind=mytype)

dvdxdwdym = (real(coeff-1, kind=mytype)*dvdxdwdym + fine_interpol(dvdx*dwdy)) / real(coeff, kind=mytype) ! vw
dvdydwdxm = (real(coeff-1, kind=mytype)*dvdydwdxm + fine_interpol(dvdy*dwdx)) / real(coeff, kind=mytype)
dvdxdwdzm = (real(coeff-1, kind=mytype)*dvdxdwdzm + fine_interpol(dvdx*dwdz)) / real(coeff, kind=mytype)
dvdzdwdxm = (real(coeff-1, kind=mytype)*dvdzdwdxm + fine_interpol(dvdz*dwdx)) / real(coeff, kind=mytype)
dvdydwdzm = (real(coeff-1, kind=mytype)*dvdydwdzm + fine_interpol(dvdy*dwdz)) / real(coeff, kind=mytype)
dvdzdwdym = (real(coeff-1, kind=mytype)*dvdzdwdym + fine_interpol(dvdz*dwdy)) / real(coeff, kind=mytype)

vdudxm = (real(coeff-1, kind=mytype)*vdudxm + fine_interpol(uy*dudx)) / real(coeff, kind=mytype)
vdudym = (real(coeff-1, kind=mytype)*vdudym + fine_interpol(uy*dudy)) / real(coeff, kind=mytype)
vdudzm = (real(coeff-1, kind=mytype)*vdudzm + fine_interpol(uy*dudz)) / real(coeff, kind=mytype)
vdvdxm = (real(coeff-1, kind=mytype)*vdvdxm + fine_interpol(uy*dvdx)) / real(coeff, kind=mytype)
vdvdym = (real(coeff-1, kind=mytype)*vdvdym + fine_interpol(uy*dvdy)) / real(coeff, kind=mytype)
vdvdzm = (real(coeff-1, kind=mytype)*vdvdzm + fine_interpol(uy*dvdz)) / real(coeff, kind=mytype)
vdwdxm = (real(coeff-1, kind=mytype)*vdwdxm + fine_interpol(uy*dwdx)) / real(coeff, kind=mytype)
vdwdym = (real(coeff-1, kind=mytype)*vdwdym + fine_interpol(uy*dwdy)) / real(coeff, kind=mytype)
vdwdzm = (real(coeff-1, kind=mytype)*vdwdzm + fine_interpol(uy*dwdz)) / real(coeff, kind=mytype)

dudx_dududxm = (real(coeff-1, kind=mytype)*dudx_dududxm + fine_interpol(dudx*dudx*dudx)) / real(coeff, kind=mytype) ! dudx_eps11
dudx_dududym = (real(coeff-1, kind=mytype)*dudx_dududym + fine_interpol(dudx*dudy*dudy)) / real(coeff, kind=mytype)
dudx_dududzm = (real(coeff-1, kind=mytype)*dudx_dududzm + fine_interpol(dudx*dudz*dudz)) / real(coeff, kind=mytype)
dudy_dudvdxm = (real(coeff-1, kind=mytype)*dudy_dudvdxm + fine_interpol(dudy*dudx*dvdx)) / real(coeff, kind=mytype) ! dudy_eps12
dudy_dudvdym = (real(coeff-1, kind=mytype)*dudy_dudvdym + fine_interpol(dudy*dudy*dvdy)) / real(coeff, kind=mytype)
dudy_dudvdzm = (real(coeff-1, kind=mytype)*dudy_dudvdzm + fine_interpol(dudy*dudz*dvdz)) / real(coeff, kind=mytype)
dudz_dudwdxm = (real(coeff-1, kind=mytype)*dudz_dudwdxm + fine_interpol(dudz*dudx*dwdx)) / real(coeff, kind=mytype) ! dudz_eps13
dudz_dudwdym = (real(coeff-1, kind=mytype)*dudz_dudwdym + fine_interpol(dudz*dudy*dwdy)) / real(coeff, kind=mytype)
dudz_dudwdzm = (real(coeff-1, kind=mytype)*dudz_dudwdzm + fine_interpol(dudz*dudz*dwdz)) / real(coeff, kind=mytype)

dvdx_dvdudxm = (real(coeff-1, kind=mytype)*dvdx_dvdudxm + fine_interpol(dvdx*dvdx*dudx)) / real(coeff, kind=mytype) ! dvdx_eps12
dvdx_dvdudym = (real(coeff-1, kind=mytype)*dvdx_dvdudym + fine_interpol(dvdx*dvdy*dudy)) / real(coeff, kind=mytype)
dvdx_dvdudzm = (real(coeff-1, kind=mytype)*dvdx_dvdudzm + fine_interpol(dvdx*dvdz*dudz)) / real(coeff, kind=mytype)
dvdy_dvdvdxm = (real(coeff-1, kind=mytype)*dvdy_dvdvdxm + fine_interpol(dvdy*dvdx*dvdx)) / real(coeff, kind=mytype) ! dvdy_eps22
dvdy_dvdvdym = (real(coeff-1, kind=mytype)*dvdy_dvdvdym + fine_interpol(dvdy*dvdy*dvdy)) / real(coeff, kind=mytype)
dvdy_dvdvdzm = (real(coeff-1, kind=mytype)*dvdy_dvdvdzm + fine_interpol(dvdy*dvdz*dvdz)) / real(coeff, kind=mytype)
dvdz_dvdwdxm = (real(coeff-1, kind=mytype)*dvdz_dvdwdxm + fine_interpol(dvdz*dvdx*dwdx)) / real(coeff, kind=mytype) ! dvdz_eps23
dvdz_dvdwdym = (real(coeff-1, kind=mytype)*dvdz_dvdwdym + fine_interpol(dvdz*dvdy*dwdy)) / real(coeff, kind=mytype)
dvdz_dvdwdzm = (real(coeff-1, kind=mytype)*dvdz_dvdwdzm + fine_interpol(dvdz*dvdz*dwdz)) / real(coeff, kind=mytype)

dwdx_dwdudxm = (real(coeff-1, kind=mytype)*dwdx_dwdudxm + fine_interpol(dwdx*dwdx*dudx)) / real(coeff, kind=mytype) ! dwdx_eps13
dwdx_dwdudym = (real(coeff-1, kind=mytype)*dwdx_dwdudym + fine_interpol(dwdx*dwdy*dudy)) / real(coeff, kind=mytype)
dwdx_dwdudzm = (real(coeff-1, kind=mytype)*dwdx_dwdudzm + fine_interpol(dwdx*dwdz*dudz)) / real(coeff, kind=mytype)
dwdy_dwdvdxm = (real(coeff-1, kind=mytype)*dwdy_dwdvdxm + fine_interpol(dwdy*dwdx*dvdx)) / real(coeff, kind=mytype) ! dwdy_eps23
dwdy_dwdvdym = (real(coeff-1, kind=mytype)*dwdy_dwdvdym + fine_interpol(dwdy*dwdy*dvdy)) / real(coeff, kind=mytype)
dwdy_dwdvdzm = (real(coeff-1, kind=mytype)*dwdy_dwdvdzm + fine_interpol(dwdy*dwdz*dvdz)) / real(coeff, kind=mytype)
dwdz_dwdwdxm = (real(coeff-1, kind=mytype)*dwdz_dwdwdxm + fine_interpol(dwdz*dwdx*dwdx)) / real(coeff, kind=mytype) ! dwdz_eps33
dwdz_dwdwdym = (real(coeff-1, kind=mytype)*dwdz_dwdwdym + fine_interpol(dwdz*dwdy*dwdy)) / real(coeff, kind=mytype)
dwdz_dwdwdzm = (real(coeff-1, kind=mytype)*dwdz_dwdwdzm + fine_interpol(dwdz*dwdz*dwdz)) / real(coeff, kind=mytype)

dudx_dudvdxm = (real(coeff-1, kind=mytype)*dudx_dudvdxm + fine_interpol(dudx*dudx*dvdx)) / real(coeff, kind=mytype) ! dudx_eps21
dudx_dudvdym = (real(coeff-1, kind=mytype)*dudx_dudvdym + fine_interpol(dudx*dudy*dvdy)) / real(coeff, kind=mytype)
dudx_dudvdzm = (real(coeff-1, kind=mytype)*dudx_dudvdzm + fine_interpol(dudx*dudz*dvdz)) / real(coeff, kind=mytype)
dudy_dvdvdxm = (real(coeff-1, kind=mytype)*dudy_dvdvdxm + fine_interpol(dudy*dvdx*dvdx)) / real(coeff, kind=mytype) ! dudy_eps22
dudy_dvdvdym = (real(coeff-1, kind=mytype)*dudy_dvdvdym + fine_interpol(dudy*dvdy*dvdy)) / real(coeff, kind=mytype)
dudy_dvdvdzm = (real(coeff-1, kind=mytype)*dudy_dvdvdzm + fine_interpol(dudy*dvdz*dvdz)) / real(coeff, kind=mytype)
dudz_dvdwdxm = (real(coeff-1, kind=mytype)*dudz_dvdwdxm + fine_interpol(dudz*dvdx*dwdx)) / real(coeff, kind=mytype) ! dudz_eps23
dudz_dvdwdym = (real(coeff-1, kind=mytype)*dudz_dvdwdym + fine_interpol(dudz*dvdy*dwdy)) / real(coeff, kind=mytype)
dudz_dvdwdzm = (real(coeff-1, kind=mytype)*dudz_dvdwdzm + fine_interpol(dudz*dvdz*dwdz)) / real(coeff, kind=mytype)
dvdx_dududxm = (real(coeff-1, kind=mytype)*dvdx_dududxm + fine_interpol(dvdx*dudx*dudx)) / real(coeff, kind=mytype) ! dvdx_eps11
dvdx_dududym = (real(coeff-1, kind=mytype)*dvdx_dududym + fine_interpol(dvdx*dudy*dudy)) / real(coeff, kind=mytype)
dvdx_dududzm = (real(coeff-1, kind=mytype)*dvdx_dududzm + fine_interpol(dvdx*dudz*dudz)) / real(coeff, kind=mytype)
dvdy_dudvdxm = (real(coeff-1, kind=mytype)*dvdy_dudvdxm + fine_interpol(dvdy*dudx*dvdx)) / real(coeff, kind=mytype) ! dvdy_eps12
dvdy_dudvdym = (real(coeff-1, kind=mytype)*dvdy_dudvdym + fine_interpol(dvdy*dudy*dvdy)) / real(coeff, kind=mytype)
dvdy_dudvdzm = (real(coeff-1, kind=mytype)*dvdy_dudvdzm + fine_interpol(dvdy*dudz*dvdz)) / real(coeff, kind=mytype)
dvdz_dudwdxm = (real(coeff-1, kind=mytype)*dvdz_dudwdxm + fine_interpol(dvdz*dudx*dwdx)) / real(coeff, kind=mytype) ! dvdz_eps13
dvdz_dudwdym = (real(coeff-1, kind=mytype)*dvdz_dudwdym + fine_interpol(dvdz*dudy*dwdy)) / real(coeff, kind=mytype)
dvdz_dudwdzm = (real(coeff-1, kind=mytype)*dvdz_dudwdzm + fine_interpol(dvdz*dudz*dwdz)) / real(coeff, kind=mytype)

vdudx2m = (real(coeff-1, kind=mytype)*vdudx2m + fine_interpol(uy*dudx*dudx)) / real(coeff, kind=mytype) ! EpsXX
vdudy2m = (real(coeff-1, kind=mytype)*vdudy2m + fine_interpol(uy*dudy*dudy)) / real(coeff, kind=mytype)
vdudz2m = (real(coeff-1, kind=mytype)*vdudz2m + fine_interpol(uy*dudz*dudz)) / real(coeff, kind=mytype)
vdvdx2m = (real(coeff-1, kind=mytype)*vdvdx2m + fine_interpol(uy*dvdx*dvdx)) / real(coeff, kind=mytype) ! EpsYY
vdvdy2m = (real(coeff-1, kind=mytype)*vdvdy2m + fine_interpol(uy*dvdy*dvdy)) / real(coeff, kind=mytype)
vdvdz2m = (real(coeff-1, kind=mytype)*vdvdz2m + fine_interpol(uy*dvdz*dvdz)) / real(coeff, kind=mytype)
vdwdx2m = (real(coeff-1, kind=mytype)*vdwdx2m + fine_interpol(uy*dwdx*dwdx)) / real(coeff, kind=mytype) ! EpsZZ
vdwdy2m = (real(coeff-1, kind=mytype)*vdwdy2m + fine_interpol(uy*dwdy*dwdy)) / real(coeff, kind=mytype)
vdwdz2m = (real(coeff-1, kind=mytype)*vdwdz2m + fine_interpol(uy*dwdz*dwdz)) / real(coeff, kind=mytype)
vdudvdxm = (real(coeff-1, kind=mytype)*vdudvdxm + fine_interpol(uy*dudx*dvdx)) / real(coeff, kind=mytype) ! EpsXY
vdudvdym = (real(coeff-1, kind=mytype)*vdudvdym + fine_interpol(uy*dudy*dvdy)) / real(coeff, kind=mytype)
vdudvdzm = (real(coeff-1, kind=mytype)*vdudvdzm + fine_interpol(uy*dudz*dvdz)) / real(coeff, kind=mytype)

dpdvdxm = (real(coeff-1, kind=mytype)*dpdvdxm + fine_interpol(dpdx*dvdx)) / real(coeff, kind=mytype) ! EpsYY
dpdvdym = (real(coeff-1, kind=mytype)*dpdvdym + fine_interpol(dpdy*dvdy)) / real(coeff, kind=mytype)
dpdvdzm = (real(coeff-1, kind=mytype)*dpdvdzm + fine_interpol(dpdz*dvdz)) / real(coeff, kind=mytype)
dpdudxm = (real(coeff-1, kind=mytype)*dpdudxm + fine_interpol(dpdx*dudx)) / real(coeff, kind=mytype) ! EpsXY
dpdudym = (real(coeff-1, kind=mytype)*dpdudym + fine_interpol(dpdy*dudy)) / real(coeff, kind=mytype)
dpdudzm = (real(coeff-1, kind=mytype)*dpdudzm + fine_interpol(dpdz*dudz)) / real(coeff, kind=mytype)

dpdxdudyxm = (real(coeff-1, kind=mytype)*dpdxdudyxm + fine_interpol(dpdx*dudxy)) / real(coeff, kind=mytype)
dpdydudyym = (real(coeff-1, kind=mytype)*dpdydudyym + fine_interpol(dpdy*dudyy)) / real(coeff, kind=mytype)
dpdzdudyzm = (real(coeff-1, kind=mytype)*dpdzdudyzm + fine_interpol(dpdz*dudyz)) / real(coeff, kind=mytype)
dpdxdvdxxm = (real(coeff-1, kind=mytype)*dpdxdvdxxm + fine_interpol(dpdx*dvdxx)) / real(coeff, kind=mytype)
dpdydvdxym = (real(coeff-1, kind=mytype)*dpdydvdxym + fine_interpol(dpdy*dvdxy)) / real(coeff, kind=mytype)
dpdzdvdxzm = (real(coeff-1, kind=mytype)*dpdzdvdxzm + fine_interpol(dpdz*dvdxz)) / real(coeff, kind=mytype)

dudxx2m = (real(coeff-1, kind=mytype)*dudxx2m + fine_interpol(dudxx*dudxx)) / real(coeff, kind=mytype) ! uu
dudyy2m = (real(coeff-1, kind=mytype)*dudyy2m + fine_interpol(dudyy*dudyy)) / real(coeff, kind=mytype)
dudzz2m = (real(coeff-1, kind=mytype)*dudzz2m + fine_interpol(dudzz*dudzz)) / real(coeff, kind=mytype)
dudxy2m = (real(coeff-1, kind=mytype)*dudxy2m + fine_interpol(dudxy*dudxy)) / real(coeff, kind=mytype)
dudxz2m = (real(coeff-1, kind=mytype)*dudxz2m + fine_interpol(dudxz*dudxz)) / real(coeff, kind=mytype)
dudyz2m = (real(coeff-1, kind=mytype)*dudyz2m + fine_interpol(dudyz*dudyz)) / real(coeff, kind=mytype)

dvdxx2m = (real(coeff-1, kind=mytype)*dvdxx2m + fine_interpol(dvdxx*dvdxx)) / real(coeff, kind=mytype) ! vv
dvdyy2m = (real(coeff-1, kind=mytype)*dvdyy2m + fine_interpol(dvdyy*dvdyy)) / real(coeff, kind=mytype)
dvdzz2m = (real(coeff-1, kind=mytype)*dvdzz2m + fine_interpol(dvdzz*dvdzz)) / real(coeff, kind=mytype)
dvdxy2m = (real(coeff-1, kind=mytype)*dvdxy2m + fine_interpol(dvdxy*dvdxy)) / real(coeff, kind=mytype)
dvdxz2m = (real(coeff-1, kind=mytype)*dvdxz2m + fine_interpol(dvdxz*dvdxz)) / real(coeff, kind=mytype)
dvdyz2m = (real(coeff-1, kind=mytype)*dvdyz2m + fine_interpol(dvdyz*dvdyz)) / real(coeff, kind=mytype)

dwdxx2m = (real(coeff-1, kind=mytype)*dwdxx2m + fine_interpol(dwdxx*dwdxx)) / real(coeff, kind=mytype) ! ww
dwdyy2m = (real(coeff-1, kind=mytype)*dwdyy2m + fine_interpol(dwdyy*dwdyy)) / real(coeff, kind=mytype)
dwdzz2m = (real(coeff-1, kind=mytype)*dwdzz2m + fine_interpol(dwdzz*dwdzz)) / real(coeff, kind=mytype)
dwdxy2m = (real(coeff-1, kind=mytype)*dwdxy2m + fine_interpol(dwdxy*dwdxy)) / real(coeff, kind=mytype)
dwdxz2m = (real(coeff-1, kind=mytype)*dwdxz2m + fine_interpol(dwdxz*dwdxz)) / real(coeff, kind=mytype)
dwdyz2m = (real(coeff-1, kind=mytype)*dwdyz2m + fine_interpol(dwdyz*dwdyz)) / real(coeff, kind=mytype)

dudvdxx2m = (real(coeff-1, kind=mytype)*dudvdxx2m + fine_interpol(dudxx*dvdxx)) / real(coeff, kind=mytype) ! uv
dudvdyy2m = (real(coeff-1, kind=mytype)*dudvdyy2m + fine_interpol(dudyy*dvdyy)) / real(coeff, kind=mytype)
dudvdzz2m = (real(coeff-1, kind=mytype)*dudvdzz2m + fine_interpol(dudzz*dvdzz)) / real(coeff, kind=mytype)
dudvdxy2m = (real(coeff-1, kind=mytype)*dudvdxy2m + fine_interpol(dudxy*dvdxy)) / real(coeff, kind=mytype)
dudvdxz2m = (real(coeff-1, kind=mytype)*dudvdxz2m + fine_interpol(dudxz*dvdxz)) / real(coeff, kind=mytype)
dudvdyz2m = (real(coeff-1, kind=mytype)*dudvdyz2m + fine_interpol(dudyz*dvdyz)) / real(coeff, kind=mytype)

presm  = (real(coeff-1, kind=mytype)*presm  + real(1, kind=mytype)*fine_interpol(pres)) / real(coeff, kind=mytype)
pres2m = (real(coeff-1, kind=mytype)*pres2m + fine_interpol(pres*pres)) / real(coeff, kind=mytype)

end subroutine update_user_stats_non_phi

subroutine update_user_stats_oui_phi(ux,uy,uz,dpdx,dpdy,dpdz,&
						dudx,dudy,dudz,&
						dvdx,dvdy,dvdz,&
						dwdx,dwdy,dwdz,&
						dudxx, dudyy, dudzz, dudxy, dudxz, dudyz,&
						dvdxx, dvdyy, dvdzz, dvdxy, dvdxz, dvdyz,&
						dwdxx, dwdyy, dwdzz, dwdxy, dwdxz, dwdyz,&
						pres,pp3,&
						phi,&
						dtdx,dtdy,dtdz,&
						dtdxx,dtdyy,dtdzz,dtdxy,dtdxz,dtdyz,&
						dtdxxx,dtdxyy,dtdxzz,&
						dtdyxx,dtdyyy,dtdyzz,&
						dtdzxx,dtdzyy,dtdzzz,&
						phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : xsize, nrank, DECOMP_INFO
use param, only : itime

implicit none

TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux,uy,uz,&
						dpdx, dpdy, dpdz,&
						dudx, dudy, dudz,&
						dvdx, dvdy, dvdz,&
						dwdx, dwdy, dwdz,&
						dudxx, dudyy, dudzz, dudxy, dudxz, dudyz,&
						dvdxx, dvdyy, dvdzz, dvdxy, dvdxz, dvdyz,&
						dwdxx, dwdyy, dwdzz, dwdxy, dwdxz, dwdyz,&
						pres
real(mytype), dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),ph3%zst(3):ph3%zen(3)), intent(in) :: pp3
real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: phi,&
						dtdx,dtdy,dtdz,&
						dtdxx,dtdyy,dtdzz,dtdxy,dtdxz,dtdyz,&
						dtdxxx,dtdxyy,dtdxzz,&
						dtdyxx,dtdyyy,dtdyzz,&
						dtdzxx,dtdzyy,dtdzzz

! Local variable
integer :: coeff

coeff=itime-beg_stat
if (coeff.le.0) then
  if (nrank.eq.0) print*,'ERROR : cannot start to compute statistics at ',itime,'. Wait till ',beg_stat
  return
endif

call update_user_stats(ux,uy,uz,dpdx,dpdy,dpdz,&
				dudx,dudy,dudz,&
				dvdx,dvdy,dvdz,&
				dwdx,dwdy,dwdz,&
				dudxx,dudyy,dudzz,dudxy,dudxz,dudyz,&
				dvdxx,dvdyy,dvdzz,dvdxy,dvdxz,dvdyz,&
				dwdxx,dwdyy,dwdzz,dwdxy,dwdxz,dwdyz,&
				pres,pp3,&
				phG,ph1,ph2,ph3,ph4)

phim = (real(coeff-1, kind=mytype)*phim + fine_interpol(phi)) / real(coeff, kind=mytype)
dphidym = (real(coeff-1, kind=mytype)*dphidym + fine_interpol(dtdy)) / real(coeff, kind=mytype)
dphidyym = (real(coeff-1, kind=mytype)*dphidyym + fine_interpol(dtdyy)) / real(coeff, kind=mytype)
dphidyyym = (real(coeff-1, kind=mytype)*dphidyyym + fine_interpol(dtdyyy)) / real(coeff, kind=mytype)

uphim = (real(coeff-1, kind=mytype)*uphim + fine_interpol(ux*phi)) / real(coeff, kind=mytype)
vphim = (real(coeff-1, kind=mytype)*vphim + fine_interpol(uy*phi)) / real(coeff, kind=mytype)
wphim = (real(coeff-1, kind=mytype)*wphim + fine_interpol(uz*phi)) / real(coeff, kind=mytype)
phiphim = (real(coeff-1, kind=mytype)*phiphim + fine_interpol(phi*phi)) / real(coeff, kind=mytype)

phidpdxm = (real(coeff-1, kind=mytype)*phidpdxm + fine_interpol(phi*dpdx)) / real(coeff, kind=mytype)
phidpdym = (real(coeff-1, kind=mytype)*phidpdym + fine_interpol(phi*dpdy)) / real(coeff, kind=mytype)
phidpdzm = (real(coeff-1, kind=mytype)*phidpdzm + fine_interpol(phi*dpdz)) / real(coeff, kind=mytype)

dphidx2m = (real(coeff-1, kind=mytype)*dphidx2m + fine_interpol(dtdx*dtdx)) / real(coeff, kind=mytype)
dphidy2m = (real(coeff-1, kind=mytype)*dphidy2m + fine_interpol(dtdy*dtdy)) / real(coeff, kind=mytype)
dphidz2m = (real(coeff-1, kind=mytype)*dphidz2m + fine_interpol(dtdz*dtdz)) / real(coeff, kind=mytype)
dudphidxm = (real(coeff-1, kind=mytype)*dudphidxm + fine_interpol(dudx*dtdx)) / real(coeff, kind=mytype)
dvdphidxm = (real(coeff-1, kind=mytype)*dvdphidxm + fine_interpol(dvdx*dtdx)) / real(coeff, kind=mytype)
dwdphidxm = (real(coeff-1, kind=mytype)*dwdphidxm + fine_interpol(dwdx*dtdx)) / real(coeff, kind=mytype)
dudphidym = (real(coeff-1, kind=mytype)*dudphidym + fine_interpol(dudy*dtdy)) / real(coeff, kind=mytype)
dvdphidym = (real(coeff-1, kind=mytype)*dvdphidym + fine_interpol(dvdy*dtdy)) / real(coeff, kind=mytype)
dwdphidym = (real(coeff-1, kind=mytype)*dwdphidym + fine_interpol(dwdy*dtdy)) / real(coeff, kind=mytype)
dudphidzm = (real(coeff-1, kind=mytype)*dudphidzm + fine_interpol(dudz*dtdz)) / real(coeff, kind=mytype)
dvdphidzm = (real(coeff-1, kind=mytype)*dvdphidzm + fine_interpol(dvdz*dtdz)) / real(coeff, kind=mytype)
dwdphidzm = (real(coeff-1, kind=mytype)*dwdphidzm + fine_interpol(dwdz*dtdz)) / real(coeff, kind=mytype)

udtdxxm = (real(coeff-1, kind=mytype)*udtdxxm + fine_interpol(ux*dtdxx)) / real(coeff, kind=mytype)
udtdyym = (real(coeff-1, kind=mytype)*udtdyym + fine_interpol(ux*dtdyy)) / real(coeff, kind=mytype)
udtdzzm = (real(coeff-1, kind=mytype)*udtdzzm + fine_interpol(ux*dtdzz)) / real(coeff, kind=mytype)

vdtdxxm = (real(coeff-1, kind=mytype)*vdtdxxm + fine_interpol(uy*dtdxx)) / real(coeff, kind=mytype)
vdtdyym = (real(coeff-1, kind=mytype)*vdtdyym + fine_interpol(uy*dtdyy)) / real(coeff, kind=mytype)
vdtdzzm = (real(coeff-1, kind=mytype)*vdtdzzm + fine_interpol(uy*dtdzz)) / real(coeff, kind=mytype)

wdtdxxm = (real(coeff-1, kind=mytype)*wdtdxxm + fine_interpol(uz*dtdxx)) / real(coeff, kind=mytype)
wdtdyym = (real(coeff-1, kind=mytype)*wdtdyym + fine_interpol(uz*dtdyy)) / real(coeff, kind=mytype)
wdtdzzm = (real(coeff-1, kind=mytype)*wdtdzzm + fine_interpol(uz*dtdzz)) / real(coeff, kind=mytype)

uphi2m = (real(coeff-1, kind=mytype)*uphi2m + fine_interpol(ux*phi*phi)) / real(coeff, kind=mytype)
vphi2m = (real(coeff-1, kind=mytype)*vphi2m + fine_interpol(uy*phi*phi)) / real(coeff, kind=mytype)
wphi2m = (real(coeff-1, kind=mytype)*wphi2m + fine_interpol(uz*phi*phi)) / real(coeff, kind=mytype)
phiuum = (real(coeff-1, kind=mytype)*phiuum + fine_interpol(ux*ux*phi)) / real(coeff, kind=mytype)
phivvm = (real(coeff-1, kind=mytype)*phivvm + fine_interpol(uy*uy*phi)) / real(coeff, kind=mytype)
phiwwm = (real(coeff-1, kind=mytype)*phiwwm + fine_interpol(uz*uz*phi)) / real(coeff, kind=mytype)
phiuvm = (real(coeff-1, kind=mytype)*phiuvm + fine_interpol(ux*uy*phi)) / real(coeff, kind=mytype)
phiuwm = (real(coeff-1, kind=mytype)*phiuwm + fine_interpol(ux*uz*phi)) / real(coeff, kind=mytype)
phivwm = (real(coeff-1, kind=mytype)*phivwm + fine_interpol(uy*uz*phi)) / real(coeff, kind=mytype)

dtdxdtdym = (real(coeff-1, kind=mytype)*dtdxdtdym + fine_interpol(dtdx*dtdy)) / real(coeff, kind=mytype)
dtdxdtdzm = (real(coeff-1, kind=mytype)*dtdxdtdzm + fine_interpol(dtdx*dtdz)) / real(coeff, kind=mytype)
dtdydtdzm = (real(coeff-1, kind=mytype)*dtdydtdzm + fine_interpol(dtdy*dtdz)) / real(coeff, kind=mytype)

dtdxdudym = (real(coeff-1, kind=mytype)*dtdxdudym + fine_interpol(dtdx*dudy)) / real(coeff, kind=mytype)
dtdydudxm = (real(coeff-1, kind=mytype)*dtdydudxm + fine_interpol(dtdy*dudx)) / real(coeff, kind=mytype)
dtdxdudzm = (real(coeff-1, kind=mytype)*dtdxdudzm + fine_interpol(dtdx*dudz)) / real(coeff, kind=mytype)
dtdzdudxm = (real(coeff-1, kind=mytype)*dtdzdudxm + fine_interpol(dtdz*dudx)) / real(coeff, kind=mytype)
dtdydudzm = (real(coeff-1, kind=mytype)*dtdydudzm + fine_interpol(dtdy*dudz)) / real(coeff, kind=mytype)
dtdzdudym = (real(coeff-1, kind=mytype)*dtdzdudym + fine_interpol(dtdz*dudy)) / real(coeff, kind=mytype)

dtdxdvdym = (real(coeff-1, kind=mytype)*dtdxdvdym + fine_interpol(dtdx*dvdy)) / real(coeff, kind=mytype)
dtdydvdxm = (real(coeff-1, kind=mytype)*dtdydvdxm + fine_interpol(dtdy*dvdx)) / real(coeff, kind=mytype)
dtdxdvdzm = (real(coeff-1, kind=mytype)*dtdxdvdzm + fine_interpol(dtdx*dvdz)) / real(coeff, kind=mytype)
dtdzdvdxm = (real(coeff-1, kind=mytype)*dtdzdvdxm + fine_interpol(dtdz*dvdx)) / real(coeff, kind=mytype)
dtdydvdzm = (real(coeff-1, kind=mytype)*dtdydvdzm + fine_interpol(dtdy*dvdz)) / real(coeff, kind=mytype)
dtdzdvdym = (real(coeff-1, kind=mytype)*dtdzdvdym + fine_interpol(dtdz*dvdy)) / real(coeff, kind=mytype)

dtdxdwdym = (real(coeff-1, kind=mytype)*dtdxdwdym + fine_interpol(dtdx*dwdy)) / real(coeff, kind=mytype)
dtdydwdxm = (real(coeff-1, kind=mytype)*dtdydwdxm + fine_interpol(dtdy*dwdx)) / real(coeff, kind=mytype)
dtdxdwdzm = (real(coeff-1, kind=mytype)*dtdxdwdzm + fine_interpol(dtdx*dwdz)) / real(coeff, kind=mytype)
dtdzdwdxm = (real(coeff-1, kind=mytype)*dtdzdwdxm + fine_interpol(dtdz*dwdx)) / real(coeff, kind=mytype)
dtdydwdzm = (real(coeff-1, kind=mytype)*dtdydwdzm + fine_interpol(dtdy*dwdz)) / real(coeff, kind=mytype)
dtdzdwdym = (real(coeff-1, kind=mytype)*dtdzdwdym + fine_interpol(dtdz*dwdy)) / real(coeff, kind=mytype)

vdtdxm = (real(coeff-1, kind=mytype)*vdtdxm + fine_interpol(uy*dtdx)) / real(coeff, kind=mytype)
vdtdym = (real(coeff-1, kind=mytype)*vdtdym + fine_interpol(uy*dtdy)) / real(coeff, kind=mytype)
vdtdzm = (real(coeff-1, kind=mytype)*vdtdzm + fine_interpol(uy*dtdz)) / real(coeff, kind=mytype)

dudx_dudtdxm = (real(coeff-1, kind=mytype)*dudx_dudtdxm + fine_interpol(dudx*dudx*dtdx)) / real(coeff, kind=mytype)
dudx_dudtdym = (real(coeff-1, kind=mytype)*dudx_dudtdym + fine_interpol(dudx*dudy*dtdy)) / real(coeff, kind=mytype)
dudx_dudtdzm = (real(coeff-1, kind=mytype)*dudx_dudtdzm + fine_interpol(dudx*dudz*dtdz)) / real(coeff, kind=mytype)
dudy_dvdtdxm = (real(coeff-1, kind=mytype)*dudy_dvdtdxm + fine_interpol(dudy*dvdx*dtdx)) / real(coeff, kind=mytype)
dudy_dvdtdym = (real(coeff-1, kind=mytype)*dudy_dvdtdym + fine_interpol(dudy*dvdy*dtdy)) / real(coeff, kind=mytype)
dudy_dvdtdzm = (real(coeff-1, kind=mytype)*dudy_dvdtdzm + fine_interpol(dudy*dvdz*dtdz)) / real(coeff, kind=mytype)
dudz_dwdtdxm = (real(coeff-1, kind=mytype)*dudz_dwdtdxm + fine_interpol(dudz*dwdx*dtdx)) / real(coeff, kind=mytype)
dudz_dwdtdym = (real(coeff-1, kind=mytype)*dudz_dwdtdym + fine_interpol(dudz*dwdy*dtdy)) / real(coeff, kind=mytype)
dudz_dwdtdzm = (real(coeff-1, kind=mytype)*dudz_dwdtdzm + fine_interpol(dudz*dwdz*dtdz)) / real(coeff, kind=mytype)

dtdx_dududym = (real(coeff-1, kind=mytype)*dtdx_dududym + fine_interpol(dtdx*dudy*dudy)) / real(coeff, kind=mytype)
dtdx_dududzm = (real(coeff-1, kind=mytype)*dtdx_dududzm + fine_interpol(dtdx*dudz*dudz)) / real(coeff, kind=mytype)
dtdy_dudvdxm = (real(coeff-1, kind=mytype)*dtdy_dudvdxm + fine_interpol(dtdy*dudx*dvdx)) / real(coeff, kind=mytype)
dtdy_dudvdzm = (real(coeff-1, kind=mytype)*dtdy_dudvdzm + fine_interpol(dtdy*dudz*dvdz)) / real(coeff, kind=mytype)
dtdz_dudwdxm = (real(coeff-1, kind=mytype)*dtdz_dudwdxm + fine_interpol(dtdz*dudx*dwdx)) / real(coeff, kind=mytype)
dtdz_dudwdym = (real(coeff-1, kind=mytype)*dtdz_dudwdym + fine_interpol(dtdz*dudy*dwdy)) / real(coeff, kind=mytype)

dvdx_dudtdxm = (real(coeff-1, kind=mytype)*dvdx_dudtdxm + fine_interpol(dvdx*dudx*dtdx)) / real(coeff, kind=mytype)
dvdx_dudtdym = (real(coeff-1, kind=mytype)*dvdx_dudtdym + fine_interpol(dvdx*dudy*dtdy)) / real(coeff, kind=mytype)
dvdx_dudtdzm = (real(coeff-1, kind=mytype)*dvdx_dudtdzm + fine_interpol(dvdx*dudz*dtdz)) / real(coeff, kind=mytype)
dvdy_dvdtdxm = (real(coeff-1, kind=mytype)*dvdy_dvdtdxm + fine_interpol(dvdy*dvdx*dtdx)) / real(coeff, kind=mytype)
dvdy_dvdtdym = (real(coeff-1, kind=mytype)*dvdy_dvdtdym + fine_interpol(dvdy*dvdy*dtdy)) / real(coeff, kind=mytype)
dvdy_dvdtdzm = (real(coeff-1, kind=mytype)*dvdy_dvdtdzm + fine_interpol(dvdy*dvdz*dtdz)) / real(coeff, kind=mytype)
dvdz_dwdtdxm = (real(coeff-1, kind=mytype)*dvdz_dwdtdxm + fine_interpol(dvdz*dwdx*dtdx)) / real(coeff, kind=mytype)
dvdz_dwdtdym = (real(coeff-1, kind=mytype)*dvdz_dwdtdym + fine_interpol(dvdz*dwdy*dtdy)) / real(coeff, kind=mytype)
dvdz_dwdtdzm = (real(coeff-1, kind=mytype)*dvdz_dwdtdzm + fine_interpol(dvdz*dwdz*dtdz)) / real(coeff, kind=mytype)

dtdx_dudvdym = (real(coeff-1, kind=mytype)*dtdx_dudvdym + fine_interpol(dtdx*dudy*dvdy)) / real(coeff, kind=mytype)
dtdx_dudvdzm = (real(coeff-1, kind=mytype)*dtdx_dudvdzm + fine_interpol(dtdx*dudz*dvdz)) / real(coeff, kind=mytype)
dtdy_dvdvdxm = (real(coeff-1, kind=mytype)*dtdy_dvdvdxm + fine_interpol(dtdy*dvdx*dvdx)) / real(coeff, kind=mytype)
dtdy_dvdvdzm = (real(coeff-1, kind=mytype)*dtdy_dvdvdzm + fine_interpol(dtdy*dvdz*dvdz)) / real(coeff, kind=mytype)
dtdz_dvdwdxm = (real(coeff-1, kind=mytype)*dtdz_dvdwdxm + fine_interpol(dtdz*dvdx*dwdx)) / real(coeff, kind=mytype)
dtdz_dvdwdym = (real(coeff-1, kind=mytype)*dtdz_dvdwdym + fine_interpol(dtdz*dvdy*dwdy)) / real(coeff, kind=mytype)

dtdx_dudtdxm = (real(coeff-1, kind=mytype)*dtdx_dudtdxm + fine_interpol(dtdx*dudx*dtdx)) / real(coeff, kind=mytype)
dtdx_dudtdym = (real(coeff-1, kind=mytype)*dtdx_dudtdym + fine_interpol(dtdx*dudy*dtdy)) / real(coeff, kind=mytype)
dtdx_dudtdzm = (real(coeff-1, kind=mytype)*dtdx_dudtdzm + fine_interpol(dtdx*dudz*dtdz)) / real(coeff, kind=mytype)
dtdy_dvdtdxm = (real(coeff-1, kind=mytype)*dtdy_dvdtdxm + fine_interpol(dtdy*dvdx*dtdx)) / real(coeff, kind=mytype)
dtdy_dvdtdym = (real(coeff-1, kind=mytype)*dtdy_dvdtdym + fine_interpol(dtdy*dvdy*dtdy)) / real(coeff, kind=mytype)
dtdy_dvdtdzm = (real(coeff-1, kind=mytype)*dtdy_dvdtdzm + fine_interpol(dtdy*dvdz*dtdz)) / real(coeff, kind=mytype)
dtdz_dwdtdxm = (real(coeff-1, kind=mytype)*dtdz_dwdtdxm + fine_interpol(dtdz*dwdx*dtdx)) / real(coeff, kind=mytype)
dtdz_dwdtdym = (real(coeff-1, kind=mytype)*dtdz_dwdtdym + fine_interpol(dtdz*dwdy*dtdy)) / real(coeff, kind=mytype)
dtdz_dwdtdzm = (real(coeff-1, kind=mytype)*dtdz_dwdtdzm + fine_interpol(dtdz*dwdz*dtdz)) / real(coeff, kind=mytype)

vdudtdxm = (real(coeff-1, kind=mytype)*vdudtdxm + fine_interpol(uy*dudx*dtdx)) / real(coeff, kind=mytype)
vdudtdym = (real(coeff-1, kind=mytype)*vdudtdym + fine_interpol(uy*dudy*dtdy)) / real(coeff, kind=mytype)
vdudtdzm = (real(coeff-1, kind=mytype)*vdudtdzm + fine_interpol(uy*dudz*dtdz)) / real(coeff, kind=mytype)
vdvdtdxm = (real(coeff-1, kind=mytype)*vdvdtdxm + fine_interpol(uy*dvdx*dtdx)) / real(coeff, kind=mytype)
vdvdtdym = (real(coeff-1, kind=mytype)*vdvdtdym + fine_interpol(uy*dvdy*dtdy)) / real(coeff, kind=mytype)
vdvdtdzm = (real(coeff-1, kind=mytype)*vdvdtdzm + fine_interpol(uy*dvdz*dtdz)) / real(coeff, kind=mytype)
vdtdx2m = (real(coeff-1, kind=mytype)*vdtdx2m + fine_interpol(uy*dtdx*dtdx)) / real(coeff, kind=mytype)
vdtdy2m = (real(coeff-1, kind=mytype)*vdtdy2m + fine_interpol(uy*dtdy*dtdy)) / real(coeff, kind=mytype)
vdtdz2m = (real(coeff-1, kind=mytype)*vdtdz2m + fine_interpol(uy*dtdz*dtdz)) / real(coeff, kind=mytype)

dpdtdxm = (real(coeff-1, kind=mytype)*dpdtdxm + fine_interpol(dpdx*dtdx)) / real(coeff, kind=mytype)
dpdtdym = (real(coeff-1, kind=mytype)*dpdtdym + fine_interpol(dpdy*dtdy)) / real(coeff, kind=mytype)
dpdtdzm = (real(coeff-1, kind=mytype)*dpdtdzm + fine_interpol(dpdz*dtdz)) / real(coeff, kind=mytype)
dpdtdxxm = (real(coeff-1, kind=mytype)*dpdtdxxm + fine_interpol(dpdx*dtdxx)) / real(coeff, kind=mytype)
dpdtdxym = (real(coeff-1, kind=mytype)*dpdtdxym + fine_interpol(dpdx*dtdxy)) / real(coeff, kind=mytype)
dpdtdxzm = (real(coeff-1, kind=mytype)*dpdtdxzm + fine_interpol(dpdx*dtdxz)) / real(coeff, kind=mytype)
dpdtdyxm = (real(coeff-1, kind=mytype)*dpdtdyxm + fine_interpol(dpdy*dtdxy)) / real(coeff, kind=mytype)
dpdtdyym = (real(coeff-1, kind=mytype)*dpdtdyym + fine_interpol(dpdy*dtdyy)) / real(coeff, kind=mytype)
dpdtdyzm = (real(coeff-1, kind=mytype)*dpdtdyzm + fine_interpol(dpdy*dtdyz)) / real(coeff, kind=mytype)

dudtdxxm = (real(coeff-1, kind=mytype)*dudtdxxm + fine_interpol(dudxx*dtdxx)) / real(coeff, kind=mytype)
dudtdyym = (real(coeff-1, kind=mytype)*dudtdyym + fine_interpol(dudyy*dtdyy)) / real(coeff, kind=mytype)
dudtdzzm = (real(coeff-1, kind=mytype)*dudtdzzm + fine_interpol(dudzz*dtdzz)) / real(coeff, kind=mytype)
dudtdxym = (real(coeff-1, kind=mytype)*dudtdxym + fine_interpol(dudxy*dtdxy)) / real(coeff, kind=mytype)
dudtdxzm = (real(coeff-1, kind=mytype)*dudtdxzm + fine_interpol(dudxz*dtdxz)) / real(coeff, kind=mytype)
dudtdyzm = (real(coeff-1, kind=mytype)*dudtdyzm + fine_interpol(dudyz*dtdyz)) / real(coeff, kind=mytype)

dudx_dtdxxxm = (real(coeff-1, kind=mytype)*dudx_dtdxxxm + fine_interpol(dudx*dtdxxx)) / real(coeff, kind=mytype)
dudx_dtdyxxm = (real(coeff-1, kind=mytype)*dudx_dtdyxxm + fine_interpol(dudx*dtdyxx)) / real(coeff, kind=mytype)
dudx_dtdzxxm = (real(coeff-1, kind=mytype)*dudx_dtdzxxm + fine_interpol(dudx*dtdzxx)) / real(coeff, kind=mytype)
dudy_dtdxyym = (real(coeff-1, kind=mytype)*dudy_dtdxyym + fine_interpol(dudy*dtdxyy)) / real(coeff, kind=mytype)
dudy_dtdyyym = (real(coeff-1, kind=mytype)*dudy_dtdyyym + fine_interpol(dudy*dtdyyy)) / real(coeff, kind=mytype)
dudy_dtdzyym = (real(coeff-1, kind=mytype)*dudy_dtdzyym + fine_interpol(dudy*dtdzyy)) / real(coeff, kind=mytype)
dudz_dtdxzzm = (real(coeff-1, kind=mytype)*dudz_dtdxzzm + fine_interpol(dudz*dtdxzz)) / real(coeff, kind=mytype)
dudz_dtdyzzm = (real(coeff-1, kind=mytype)*dudz_dtdyzzm + fine_interpol(dudz*dtdyzz)) / real(coeff, kind=mytype)
dudz_dtdzzzm = (real(coeff-1, kind=mytype)*dudz_dtdzzzm + fine_interpol(dudz*dtdzzz)) / real(coeff, kind=mytype)

dvdtdxxm = (real(coeff-1, kind=mytype)*dvdtdxxm + fine_interpol(dvdxx*dtdxx)) / real(coeff, kind=mytype)
dvdtdyym = (real(coeff-1, kind=mytype)*dvdtdyym + fine_interpol(dvdyy*dtdyy)) / real(coeff, kind=mytype)
dvdtdzzm = (real(coeff-1, kind=mytype)*dvdtdzzm + fine_interpol(dvdzz*dtdzz)) / real(coeff, kind=mytype)
dvdtdxym = (real(coeff-1, kind=mytype)*dvdtdxym + fine_interpol(dvdxy*dtdxy)) / real(coeff, kind=mytype)
dvdtdxzm = (real(coeff-1, kind=mytype)*dvdtdxzm + fine_interpol(dvdxz*dtdxz)) / real(coeff, kind=mytype)
dvdtdyzm = (real(coeff-1, kind=mytype)*dvdtdyzm + fine_interpol(dvdyz*dtdyz)) / real(coeff, kind=mytype)

dvdx_dtdxxxm = (real(coeff-1, kind=mytype)*dvdx_dtdxxxm + fine_interpol(dvdx*dtdxxx)) / real(coeff, kind=mytype)
dvdx_dtdyxxm = (real(coeff-1, kind=mytype)*dvdx_dtdyxxm + fine_interpol(dvdx*dtdyxx)) / real(coeff, kind=mytype)
dvdx_dtdzxxm = (real(coeff-1, kind=mytype)*dvdx_dtdzxxm + fine_interpol(dvdx*dtdzxx)) / real(coeff, kind=mytype)
dvdy_dtdxyym = (real(coeff-1, kind=mytype)*dvdy_dtdxyym + fine_interpol(dvdy*dtdxyy)) / real(coeff, kind=mytype)
dvdy_dtdyyym = (real(coeff-1, kind=mytype)*dvdy_dtdyyym + fine_interpol(dvdy*dtdyyy)) / real(coeff, kind=mytype)
dvdy_dtdzyym = (real(coeff-1, kind=mytype)*dvdy_dtdzyym + fine_interpol(dvdy*dtdzyy)) / real(coeff, kind=mytype)
dvdz_dtdxzzm = (real(coeff-1, kind=mytype)*dvdz_dtdxzzm + fine_interpol(dvdz*dtdxzz)) / real(coeff, kind=mytype)
dvdz_dtdyzzm = (real(coeff-1, kind=mytype)*dvdz_dtdyzzm + fine_interpol(dvdz*dtdyzz)) / real(coeff, kind=mytype)
dvdz_dtdzzzm = (real(coeff-1, kind=mytype)*dvdz_dtdzzzm + fine_interpol(dvdz*dtdzzz)) / real(coeff, kind=mytype)

dtdxx2m = (real(coeff-1, kind=mytype)*dtdxx2m + fine_interpol(dtdxx*dtdxx)) / real(coeff, kind=mytype)
dtdyy2m = (real(coeff-1, kind=mytype)*dtdyy2m + fine_interpol(dtdyy*dtdyy)) / real(coeff, kind=mytype)
dtdzz2m = (real(coeff-1, kind=mytype)*dtdzz2m + fine_interpol(dtdzz*dtdzz)) / real(coeff, kind=mytype)
dtdxy2m = (real(coeff-1, kind=mytype)*dtdxy2m + fine_interpol(dtdxy*dtdxy)) / real(coeff, kind=mytype)
dtdxz2m = (real(coeff-1, kind=mytype)*dtdxz2m + fine_interpol(dtdxz*dtdxz)) / real(coeff, kind=mytype)
dtdyz2m = (real(coeff-1, kind=mytype)*dtdyz2m + fine_interpol(dtdyz*dtdyz)) / real(coeff, kind=mytype)

end subroutine update_user_stats_oui_phi

subroutine pre_update_user_stats(phG,ph1,ph2,ph3,ph4)

  use decomp_2d_io, only : decomp_2d_write_plane
  use decomp_2d, only : xsize, ysize, zsize, DECOMP_INFO, get_decomp_info
  use param, only : istret, iscalar, iimplicit
  use variables
  use var

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
  logical :: myflag
  integer :: i,j,k

  ! Pressure interpolation
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: tmp_pres
  real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: myta3,mydi3
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: myta2
  real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: mytb2,mydi2
  real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: myta1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: mytb1
  TYPE(DECOMP_INFO) :: decomp_main
  character(len=10) :: timer

  real(mytype),dimension(ysize(1),8,ysize(3)) :: output_tdyt

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dudx, dudy, dudz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dvdx, dvdy, dvdz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dwdx, dwdy, dwdz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dudxx,dudyy,dudzz,dudxy,dudxz,dudyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dvdxx,dvdyy,dvdzz,dvdxy,dvdxz,dvdyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dwdxx,dwdyy,dwdzz,dwdxy,dwdxz,dwdyz

  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dtdx, dtdy, dtdz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dtdxx,dtdyy,dtdzz,dtdxy,dtdxz,dtdyz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dtdxxx,dtdxyy,dtdxzz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dtdyxx,dtdyyy,dtdyzz
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dtdzxx,dtdzyy,dtdzzz

TYPE(DECOMP_INFO) :: phglob

call get_decomp_info(phglob)

  myflag=.false.
  if (iimplicit.eq.1) myflag=.true.
  if (myflag) iimplicit=0

  ! Compute laplacian temperature (budget uPhi, vPhi, wPhi) and grad(phi)
  if (iscalar.eq.1) then

    ! Compute laplacian phi
    !
    call transpose_x_to_y(phi1,phi2,phglob)
    call transpose_y_to_z(phi2,phi3)
    !
    call derx(dtdx,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derxxt(dtdxx,phi1,di1,sx,sfxt,ssxt,swxt,xsize(1),xsize(2),xsize(3),0)
    call derx(dtdxxx,dtdxx,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    !
    call dery (tb2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call transpose_y_to_x(tb2,dtdy) ! dTdy
    call derx(dtdxy,dtdy,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! dTdxy
    call derxxt(dtdyxx,dtdy,di1,sx,sfxt,ssxt,swxt,xsize(1),xsize(2),xsize(3),0) ! dTdyxx
    call deryyt(ta2,phi2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)
    if (istret.ne.0) then 
      do j=1,ysize(2)
        ta2(:,j,:)=ta2(:,j,:)*pp2y(j)-pp4y(j)*tb2(:,j,:) ! dyy dans ta2
      enddo
    endif
    call transpose_y_to_x(ta2,dtdyy) ! dTdyy
    call derx(dtdxyy,dtdyy,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! dTdxyy
    call dery (tc2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call transpose_y_to_x(tc2,dtdyyy) ! dTdyyy
    call transpose_y_to_z(ta2,ta3)
    call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(tb3,ta2)
    call transpose_y_to_x(ta2,dtdzyy)
    !
    call derz (ta3,phi3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(ta3,ta2) ! dTdz dans tb2
    call transpose_y_to_x(ta2,dtdz) ! dTdz
    call derxxt(dtdzxx,dtdz,di1,sx,sfxt,ssxt,swxt,xsize(1),xsize(2),xsize(3),0)
    call dery (tc2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call transpose_y_to_x(tc2,dtdyz) ! dTdyz
    call derx(dtdxz,dtdz,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derzzt (ta3,phi3,di3,sz,sfzt,sszt,swzt,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(ta3,ta2) ! dzz dans ta2
    call transpose_y_to_x(ta2,dtdzz) ! dTdzz
    call derx(dtdxzz,dtdzz,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call dery (tc2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call transpose_y_to_x(tc2,dtdyzz)
    call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call transpose_z_to_y(tb3,ta2)
    call transpose_y_to_x(ta2,dtdzzz)

if (mod(itime,20).eq.0) then
154 format('t',I9.9)
  write(timer,154) itime
  output_tdyt(:,1,:)=phi2(:,1,:)
  output_tdyt(:,2,:)= tb2(:,1,:)
  output_tdyt(:,3,:)=phi2(:,11,:)
  output_tdyt(:,4,:)= tb2(:,11,:)
  output_tdyt(:,5,:)=phi2(:,30,:)
  output_tdyt(:,6,:)= tb2(:,30,:)
  output_tdyt(:,7,:)=phi2(:,96,:)
  output_tdyt(:,8,:)= tb2(:,96,:)
  call decomp_2d_write_plane(2,output_tdyt,5,8,'slices/tdyt_flu_'//timer//'.dat')
endif

  endif

  ! Compute velocity gradient
  call derx (dudx,ux1,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
  call derx (dvdx,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (dwdx,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  ! ta1=dudx / tb1=dvdx / tc1=dwdx
  call transpose_x_to_y(ux1,ux2)
  call transpose_x_to_y(uy1,uy2)
  call transpose_x_to_y(uz1,uz2)
  call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call transpose_y_to_x(ta2,dudy)!dudy
  call transpose_y_to_x(tb2,dvdy)!dvdy
  call transpose_y_to_x(tc2,dwdy)!dwdy
! ta2=dudy / tb2=dvdy / tc2=dwdy
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(uz2,uz3)
  call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tc3,uz3,di3,sz,ffz ,fsz ,fwz ,zsize(1),zsize(2),zsize(3),0)
  ! ta3=dudz / tb3=dvdz / tc3=dwdz
  ! Back in y
  call transpose_z_to_y(ta3,td2)!dudz
  call transpose_z_to_y(tb3,te2)!dvdz
  call transpose_z_to_y(tc3,tf2)!dwdz
  ! Back to x
  call transpose_y_to_x(td2,dudz)!dudz
  call transpose_y_to_x(te2,dvdz)!dvdz
  call transpose_y_to_x(tf2,dwdz)!dwdz

  ! Compute velocity second derivatives
  call derxx (dudxx,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
  call derx (dudxy,dudy,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
  call derx (dudxz,dudz,di1,sx,ffx ,fsx ,fwx ,xsize(1),xsize(2),xsize(3),0)
  call derxx (dvdxx,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
  call derx (dvdxy,dvdy,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (dvdxz,dvdz,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derxx (dwdxx,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
  call derx (dwdxy,dwdy,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (dwdxz,dwdz,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !
  call dery (tg2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call transpose_y_to_x(tg2,dudyz)!dudyz
  call dery (tg2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call transpose_y_to_x(tg2,dvdyz)!dvdyz
  call dery (tg2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call transpose_y_to_x(tg2,dwdyz)!dwdyz
  !
  call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
  do j=1,ysize(2)
    td2(:,j,:)=td2(:,j,:)*pp2y(j)-pp4y(j)*ta2(:,j,:)
  enddo
  call transpose_y_to_x(td2,dudyy)!dudyy
  call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
  do j=1,ysize(2)
    te2(:,j,:)=te2(:,j,:)*pp2y(j)-pp4y(j)*tb2(:,j,:)
  enddo
  call transpose_y_to_x(te2,dvdyy)!dvdyy
  call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
  do j=1,ysize(2)
    tf2(:,j,:)=tf2(:,j,:)*pp2y(j)-pp4y(j)*tc2(:,j,:)
  enddo
  call transpose_y_to_x(tf2,dwdyy)!dwdyy
  !
  call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
  call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
  call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
  call transpose_z_to_y(ta3,td2)!dudzz
  call transpose_z_to_y(tb3,te2)!dvdzz
  call transpose_z_to_y(tc3,tf2)!dwdzz
  call transpose_y_to_x(td2,dudzz)!dudzz
  call transpose_y_to_x(te2,dvdzz)!dvdzz
  call transpose_y_to_x(tf2,dwdzz)!dwdzz

  if (.true.) then
    !WORK Z-PENCILS
    call interiz6(myta3,pp3,mydi3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    !WORK Y-PENCILS
    call transpose_z_to_y(myta3,myta2,ph3) !nxm nym nz
    call interiy6(mytb2,myta2,mydi2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    !WORK X-PENCILS
    call transpose_y_to_x(mytb2,myta1,ph2) !nxm ny nz
    call interi6(tmp_pres,myta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    !The pressure field on the main mesh is in tmp_pres
  endif

  if (myflag) iimplicit=1

  if (iscalar.eq.1) then
    call update_user_stats(ux1,uy1,uz1,px1,py1,pz1,&
				dudx,dudy,dudz,& ! Grad ux
				dvdx,dvdy,dvdz,& ! Grad uy
				dwdx,dwdy,dwdz,& ! Grad uz
				dudxx,dudyy,dudzz,dudxy,dudxz,dudyz,&
				dvdxx,dvdyy,dvdzz,dvdxy,dvdxz,dvdyz,&
				dwdxx,dwdyy,dwdzz,dwdxy,dwdxz,dwdyz,&
				tmp_pres,pp3,&
				phi1,&
				dtdx,dtdy,dtdz,& ! Grad phi
				dtdxx,dtdyy,dtdzz,dtdxy,dtdxz,dtdyz,& ! Laplacian phi
				dtdxxx,dtdxyy,dtdxzz,&
				dtdyxx,dtdyyy,dtdyzz,&
				dtdzxx,dtdzyy,dtdzzz,&
				phG,ph1,ph2,ph3,ph4)
  else
    call update_user_stats(ux1,uy1,uz1,px1,py1,pz1,&
				dudx,dudy,dudz,& ! Grad ux
				dvdx,dvdy,dvdz,& ! Grad uy
				dwdx,dwdy,dwdz,& ! Grad uz
				dudxx,dudyy,dudzz,dudxy,dudxz,dudyz,&
				dvdxx,dvdyy,dvdzz,dvdxy,dvdxz,dvdyz,&
				dwdxx,dwdyy,dwdzz,dwdxy,dwdxz,dwdyz,&
				tmp_pres,pp3,&
				phG,ph1,ph2,ph3,ph4)
  endif

end subroutine pre_update_user_stats

subroutine write_user_stats(phG,ph1,ph2,ph3,ph4)

use decomp_2d, only : DECOMP_INFO, get_decomp_info
use decomp_2d_io, only : decomp_2d_write_plane, decomp_2d_write_one
use param, only : iscalar

implicit none

TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4
TYPE(DECOMP_INFO) :: decomp_main

call decomp_2d_write_plane(1,um,4,1,'um.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vm,4,1,'vm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,wm,4,1,'wm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dpdxm,4,1,'dpdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdym,4,1,'dpdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdzm,4,1,'dpdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudym,4,1,'dudym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdym,4,1,'dvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdym,4,1,'dwdym.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudyym,4,1,'dudyym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdyym,4,1,'dvdyym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdyym,4,1,'dwdyym.dat',decomp_user_stats)

call decomp_2d_write_plane(1,uum,4,1,'uum.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vvm,4,1,'vvm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,wwm,4,1,'wwm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,uvm,4,1,'uvm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,uwm,4,1,'uwm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vwm,4,1,'vwm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,udpdxm,4,1,'udpdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,udpdym,4,1,'udpdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,udpdzm,4,1,'udpdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdpdxm,4,1,'vdpdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdpdym,4,1,'vdpdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdpdzm,4,1,'vdpdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,wdpdxm,4,1,'wdpdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,wdpdym,4,1,'wdpdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,wdpdzm,4,1,'wdpdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,uuum,4,1,'uuum.dat',decomp_user_stats)
call decomp_2d_write_plane(1,uvvm,4,1,'uvvm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,uwwm,4,1,'uwwm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vuum,4,1,'vuum.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vvvm,4,1,'vvvm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vwwm,4,1,'vwwm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,wuum,4,1,'wuum.dat',decomp_user_stats)
call decomp_2d_write_plane(1,wvvm,4,1,'wvvm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,wwwm,4,1,'wwwm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,uvwm,4,1,'uvwm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudx2m,4,1,'dudx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudy2m,4,1,'dudy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudz2m,4,1,'dudz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdx2m,4,1,'dvdx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdy2m,4,1,'dvdy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdz2m,4,1,'dvdz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dwdx2m,4,1,'dwdx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdy2m,4,1,'dwdy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdz2m,4,1,'dwdz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudvdxm,4,1,'dudvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudwdxm,4,1,'dudwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdwdxm,4,1,'dvdwdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudvdym,4,1,'dudvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudwdym,4,1,'dudwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdwdym,4,1,'dvdwdym.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudvdzm,4,1,'dudvdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudwdzm,4,1,'dudwdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdwdzm,4,1,'dvdwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudxdudym,4,1,'dudxdudym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudxdudzm,4,1,'dudxdudzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudydudzm,4,1,'dudydudzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdxdvdym,4,1,'dvdxdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdxdvdzm,4,1,'dvdxdvdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdydvdzm,4,1,'dvdydvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dwdxdwdym,4,1,'dwdxdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdxdwdzm,4,1,'dwdxdwdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdydwdzm,4,1,'dwdydwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudxdvdym,4,1,'dudxdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudydvdxm,4,1,'dudydvdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudxdvdzm,4,1,'dudxdvdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudzdvdxm,4,1,'dudzdvdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudydvdzm,4,1,'dudydvdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudzdvdym,4,1,'dudzdvdym.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudxdwdym,4,1,'dudxdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudydwdxm,4,1,'dudydwdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudxdwdzm,4,1,'dudxdwdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudzdwdxm,4,1,'dudzdwdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudydwdzm,4,1,'dudydwdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudzdwdym,4,1,'dudzdwdym.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdxdwdym,4,1,'dvdxdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdydwdxm,4,1,'dvdydwdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdxdwdzm,4,1,'dvdxdwdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdzdwdxm,4,1,'dvdzdwdxm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdydwdzm,4,1,'dvdydwdzm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdzdwdym,4,1,'dvdzdwdym.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdudxm,4,1,'vdudxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdudym,4,1,'vdudym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdudzm,4,1,'vdudzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdvdxm,4,1,'vdvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdvdym,4,1,'vdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdvdzm,4,1,'vdvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdwdxm,4,1,'vdwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdwdym,4,1,'vdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdwdzm,4,1,'vdwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudx_dududxm,4,1,'dudx_dududxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudx_dududym,4,1,'dudx_dududym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudx_dududzm,4,1,'dudx_dududzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudy_dudvdxm,4,1,'dudy_dudvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudy_dudvdym,4,1,'dudy_dudvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudy_dudvdzm,4,1,'dudy_dudvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudz_dudwdxm,4,1,'dudz_dudwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudz_dudwdym,4,1,'dudz_dudwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudz_dudwdzm,4,1,'dudz_dudwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdx_dvdudxm,4,1,'dvdx_dvdudxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdx_dvdudym,4,1,'dvdx_dvdudym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdx_dvdudzm,4,1,'dvdx_dvdudzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdy_dvdvdxm,4,1,'dvdy_dvdvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdy_dvdvdym,4,1,'dvdy_dvdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdy_dvdvdzm,4,1,'dvdy_dvdvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdz_dvdwdxm,4,1,'dvdz_dvdwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdz_dvdwdym,4,1,'dvdz_dvdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdz_dvdwdzm,4,1,'dvdz_dvdwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dwdx_dwdudxm,4,1,'dwdx_dwdudxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdx_dwdudym,4,1,'dwdx_dwdudym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdx_dwdudzm,4,1,'dwdx_dwdudzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dwdy_dwdvdxm,4,1,'dwdy_dwdvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdy_dwdvdym,4,1,'dwdy_dwdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdy_dwdvdzm,4,1,'dwdy_dwdvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dwdz_dwdwdxm,4,1,'dwdz_dwdwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdz_dwdwdym,4,1,'dwdz_dwdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdz_dwdwdzm,4,1,'dwdz_dwdwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudx_dudvdxm,4,1,'dudx_dudvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudx_dudvdym,4,1,'dudx_dudvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudx_dudvdzm,4,1,'dudx_dudvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudy_dvdvdxm,4,1,'dudy_dvdvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudy_dvdvdym,4,1,'dudy_dvdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudy_dvdvdzm,4,1,'dudy_dvdvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudz_dvdwdxm,4,1,'dudz_dvdwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudz_dvdwdym,4,1,'dudz_dvdwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudz_dvdwdzm,4,1,'dudz_dvdwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdx_dududxm,4,1,'dvdx_dududxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdx_dududym,4,1,'dvdx_dududym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdx_dududzm,4,1,'dvdx_dududzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdy_dudvdxm,4,1,'dvdy_dudvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdy_dudvdym,4,1,'dvdy_dudvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdy_dudvdzm,4,1,'dvdy_dudvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdz_dudwdxm,4,1,'dvdz_dudwdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdz_dudwdym,4,1,'dvdz_dudwdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdz_dudwdzm,4,1,'dvdz_dudwdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdudx2m,4,1,'vdudx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdudy2m,4,1,'vdudy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdudz2m,4,1,'vdudz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdvdx2m,4,1,'vdvdx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdvdy2m,4,1,'vdvdy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdvdz2m,4,1,'vdvdz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdwdx2m,4,1,'vdwdx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdwdy2m,4,1,'vdwdy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdwdz2m,4,1,'vdwdz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,vdudvdxm,4,1,'vdudvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdudvdym,4,1,'vdudvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,vdudvdzm,4,1,'vdudvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dpdvdxm,4,1,'dpdvdxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdvdym,4,1,'dpdvdym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdvdzm,4,1,'dpdvdzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dpdudxm,4,1,'dpdudxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdudym,4,1,'dpdudym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdudzm,4,1,'dpdudzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dpdxdudyxm,4,1,'dpdxdudyxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdydudyym,4,1,'dpdydudyym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdzdudyzm,4,1,'dpdzdudyzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dpdxdvdxxm,4,1,'dpdxdvdxxm.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdydvdxym,4,1,'dpdydvdxym.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dpdzdvdxzm,4,1,'dpdzdvdxzm.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudxx2m,4,1,'dudxx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudyy2m,4,1,'dudyy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudzz2m,4,1,'dudzz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudxy2m,4,1,'dudxy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudxz2m,4,1,'dudxz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudyz2m,4,1,'dudyz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dvdxx2m,4,1,'dvdxx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdyy2m,4,1,'dvdyy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdzz2m,4,1,'dvdzz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdxy2m,4,1,'dvdxy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdxz2m,4,1,'dvdxz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dvdyz2m,4,1,'dvdyz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dwdxx2m,4,1,'dwdxx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdyy2m,4,1,'dwdyy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdzz2m,4,1,'dwdzz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdxy2m,4,1,'dwdxy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdxz2m,4,1,'dwdxz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dwdyz2m,4,1,'dwdyz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudvdxx2m,4,1,'dudvdxx2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudvdyy2m,4,1,'dudvdyy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudvdzz2m,4,1,'dudvdzz2m.dat',decomp_user_stats)

call decomp_2d_write_plane(1,dudvdxy2m,4,1,'dudvdxy2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudvdxz2m,4,1,'dudvdxz2m.dat',decomp_user_stats)
call decomp_2d_write_plane(1,dudvdyz2m,4,1,'dudvdyz2m.dat',decomp_user_stats)


if (iscalar.eq.1) then
  call decomp_2d_write_plane(1,phim,4,1,'phim.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dphidym,4,1,'dphidym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dphidyym,4,1,'dphidyym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dphidyyym,4,1,'dphidyyym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,uphim,4,1,'uphim.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vphim,4,1,'vphim.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,wphim,4,1,'wphim.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phiphim,4,1,'phiphim.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,phidpdxm,4,1,'phidpdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phidpdym,4,1,'phidpdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phidpdzm,4,1,'phidpdzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dphidx2m,4,1,'dphidx2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dphidy2m,4,1,'dphidy2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dphidz2m,4,1,'dphidz2m.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dudphidxm,4,1,'dudphidxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdphidxm,4,1,'dvdphidxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dwdphidxm,4,1,'dwdphidxm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dudphidym,4,1,'dudphidym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdphidym,4,1,'dvdphidym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dwdphidym,4,1,'dwdphidym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dudphidzm,4,1,'dudphidzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdphidzm,4,1,'dvdphidzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dwdphidzm,4,1,'dwdphidzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,udtdxxm,4,1,'udtdxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,udtdyym,4,1,'udtdyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,udtdzzm,4,1,'udtdzzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,vdtdxxm,4,1,'vdtdxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdyym,4,1,'vdtdyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdzzm,4,1,'vdtdzzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,wdtdxxm,4,1,'wdtdxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,wdtdyym,4,1,'wdtdyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,wdtdzzm,4,1,'wdtdzzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,uphi2m,4,1,'uphi2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vphi2m,4,1,'vphi2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,wphi2m,4,1,'wphi2m.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,phiuum,4,1,'phiuum.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phivvm,4,1,'phivvm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phiwwm,4,1,'phiwwm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,phiuvm,4,1,'phiuvm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phiuwm,4,1,'phiuwm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,phivwm,4,1,'phivwm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdxdtdym,4,1,'dtdxdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdxdtdzm,4,1,'dtdxdtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydtdzm,4,1,'dtdydtdzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdxdudym,4,1,'dtdxdudym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydudxm,4,1,'dtdydudxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdxdudzm,4,1,'dtdxdudzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzdudxm,4,1,'dtdzdudxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydudzm,4,1,'dtdydudzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzdudym,4,1,'dtdzdudym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdxdvdym,4,1,'dtdxdvdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydvdxm,4,1,'dtdydvdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdxdvdzm,4,1,'dtdxdvdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzdvdxm,4,1,'dtdzdvdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydvdzm,4,1,'dtdydvdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzdvdym,4,1,'dtdzdvdym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdxdwdym,4,1,'dtdxdwdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydwdxm,4,1,'dtdydwdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdxdwdzm,4,1,'dtdxdwdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzdwdxm,4,1,'dtdzdwdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdydwdzm,4,1,'dtdydwdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzdwdym,4,1,'dtdzdwdym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,vdtdxm,4,1,'vdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdym,4,1,'vdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdzm,4,1,'vdtdzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dudx_dudtdxm,4,1,'dudx_dudtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudx_dudtdym,4,1,'dudx_dudtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudx_dudtdzm,4,1,'dudx_dudtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudy_dvdtdxm,4,1,'dudy_dvdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudy_dvdtdym,4,1,'dudy_dvdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudy_dvdtdzm,4,1,'dudy_dvdtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudz_dwdtdxm,4,1,'dudz_dwdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudz_dwdtdym,4,1,'dudz_dwdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudz_dwdtdzm,4,1,'dudz_dwdtdzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdx_dududym,4,1,'dtdx_dududym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdx_dududzm,4,1,'dtdx_dududzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dudvdxm,4,1,'dtdy_dudvdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dudvdzm,4,1,'dtdy_dudvdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dudwdxm,4,1,'dtdz_dudwdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dudwdym,4,1,'dtdz_dudwdym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dvdx_dudtdxm,4,1,'dvdx_dudtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdx_dudtdym,4,1,'dvdx_dudtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdx_dudtdzm,4,1,'dvdx_dudtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdy_dvdtdxm,4,1,'dvdy_dvdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdy_dvdtdym,4,1,'dvdy_dvdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdy_dvdtdzm,4,1,'dvdy_dvdtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdz_dwdtdxm,4,1,'dvdz_dwdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdz_dwdtdym,4,1,'dvdz_dwdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdz_dwdtdzm,4,1,'dvdz_dwdtdzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdx_dudvdym,4,1,'dtdx_dudvdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdx_dudvdzm,4,1,'dtdx_dudvdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dvdvdxm,4,1,'dtdy_dvdvdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dvdvdzm,4,1,'dtdy_dvdvdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dvdwdxm,4,1,'dtdz_dvdwdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dvdwdym,4,1,'dtdz_dvdwdym.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdx_dudtdxm,4,1,'dtdx_dudtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdx_dudtdym,4,1,'dtdx_dudtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdx_dudtdzm,4,1,'dtdx_dudtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dvdtdxm,4,1,'dtdy_dvdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dvdtdym,4,1,'dtdy_dvdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdy_dvdtdzm,4,1,'dtdy_dvdtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dwdtdxm,4,1,'dtdz_dwdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dwdtdym,4,1,'dtdz_dwdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdz_dwdtdzm,4,1,'dtdz_dwdtdzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,vdudtdxm,4,1,'vdudtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdudtdym,4,1,'vdudtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdudtdzm,4,1,'vdudtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdvdtdxm,4,1,'vdvdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdvdtdym,4,1,'vdvdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdvdtdzm,4,1,'vdvdtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdx2m,4,1,'vdtdx2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdy2m,4,1,'vdtdy2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,vdtdz2m,4,1,'vdtdz2m.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dpdtdxm,4,1,'dpdtdxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdym,4,1,'dpdtdym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdzm,4,1,'dpdtdzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdxxm,4,1,'dpdtdxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdxym,4,1,'dpdtdxym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdxzm,4,1,'dpdtdxzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdyxm,4,1,'dpdtdyxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdyym,4,1,'dpdtdyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dpdtdyzm,4,1,'dpdtdyzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dudtdxxm,4,1,'dudtdxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudtdyym,4,1,'dudtdyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudtdzzm,4,1,'dudtdzzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudtdxym,4,1,'dudtdxym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudtdxzm,4,1,'dudtdxzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudtdyzm,4,1,'dudtdyzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dudx_dtdxxxm,4,1,'dudx_dtdxxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudx_dtdyxxm,4,1,'dudx_dtdyxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudx_dtdzxxm,4,1,'dudx_dtdzxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudy_dtdxyym,4,1,'dudy_dtdxyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudy_dtdyyym,4,1,'dudy_dtdyyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudy_dtdzyym,4,1,'dudy_dtdzyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudz_dtdxzzm,4,1,'dudz_dtdxzzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudz_dtdyzzm,4,1,'dudz_dtdyzzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dudz_dtdzzzm,4,1,'dudz_dtdzzzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dvdtdxxm,4,1,'dvdtdxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdtdyym,4,1,'dvdtdyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdtdzzm,4,1,'dvdtdzzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdtdxym,4,1,'dvdtdxym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdtdxzm,4,1,'dvdtdxzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdtdyzm,4,1,'dvdtdyzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dvdx_dtdxxxm,4,1,'dvdx_dtdxxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdx_dtdyxxm,4,1,'dvdx_dtdyxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdx_dtdzxxm,4,1,'dvdx_dtdzxxm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdy_dtdxyym,4,1,'dvdy_dtdxyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdy_dtdyyym,4,1,'dvdy_dtdyyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdy_dtdzyym,4,1,'dvdy_dtdzyym.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdz_dtdxzzm,4,1,'dvdz_dtdxzzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdz_dtdyzzm,4,1,'dvdz_dtdyzzm.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dvdz_dtdzzzm,4,1,'dvdz_dtdzzzm.dat',decomp_user_stats)
  !
  call decomp_2d_write_plane(1,dtdxx2m,4,1,'dtdxx2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdyy2m,4,1,'dtdyy2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdzz2m,4,1,'dtdzz2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdxy2m,4,1,'dtdxy2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdxz2m,4,1,'dtdxz2m.dat',decomp_user_stats)
  call decomp_2d_write_plane(1,dtdyz2m,4,1,'dtdyz2m.dat',decomp_user_stats)
  !
endif

call decomp_2d_write_plane(1,presm ,4,1,'presm.dat' ,decomp_user_stats)
call decomp_2d_write_plane(1,pres2m,4,1,'pres2m.dat',decomp_user_stats)

end subroutine write_user_stats

function fine_interpol(var)

  USE decomp_2d, only : mytype, xsize

  implicit none

  real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: var
  real(mytype), dimension(xst1:xen1,xst2:xen2,xst3:xen3) :: fine_interpol
  integer :: i,j,k

  fine_interpol=0.

  do k=xst3,xen3
  do j=xst2,xen2
  do i=1,xsize(1)
    fine_interpol(xst1,j,k)=fine_interpol(xst1,j,k) + &
			 var(i,j-xst2+1,k-xst3+1)/xsize(1)
  enddo
  enddo
  enddo

end function fine_interpol

subroutine read_user_stats(phG,ph1,ph2,ph3,ph4)
  !
  ! Read average when using a restart
  !
  use decomp_2d_io, only : decomp_2d_read_one
  use decomp_2d, only : DECOMP_INFO, get_decomp_info
  use decomp_2d_io, only : decomp_2d_read_one
  use param, only : iscalar

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

  logical :: file_exist
  character(len=30) :: filename
  TYPE(DECOMP_INFO) :: decomp_main

  filename='um.dat'
  call read_one_user_stat(filename,um)
  filename='vm.dat'
  call read_one_user_stat(filename,vm)
  filename='wm.dat'
  call read_one_user_stat(filename,wm)

  filename='dpdxm.dat'
  call read_one_user_stat(filename,dpdxm)
  filename='dpdym.dat'
  call read_one_user_stat(filename,dpdym)
  filename='dpdzm.dat'
  call read_one_user_stat(filename,dpdzm)

  filename='dudym.dat'
  call read_one_user_stat(filename,dudym)
  filename='dvdym.dat'
  call read_one_user_stat(filename,dvdym)
  filename='dwdym.dat'
  call read_one_user_stat(filename,dwdym)

  filename='dudyym.dat'
  call read_one_user_stat(filename,dudyym)
  filename='dvdyym.dat'
  call read_one_user_stat(filename,dvdyym)
  filename='dwdyym.dat'
  call read_one_user_stat(filename,dwdyym)

  filename='uum.dat'
  call read_one_user_stat(filename,uum)
  filename='vvm.dat'
  call read_one_user_stat(filename,vvm)
  filename='wwm.dat'
  call read_one_user_stat(filename,wwm)

  filename='uvm.dat'
  call read_one_user_stat(filename,uvm)
  filename='uwm.dat'
  call read_one_user_stat(filename,uwm)
  filename='vwm.dat'
  call read_one_user_stat(filename,vwm)

  filename='udpdxm.dat'
  call read_one_user_stat(filename,udpdxm)
  filename='udpdym.dat'
  call read_one_user_stat(filename,udpdym)
  filename='udpdzm.dat'
  call read_one_user_stat(filename,udpdzm)

  filename='vdpdxm.dat'
  call read_one_user_stat(filename,vdpdxm)
  filename='vdpdym.dat'
  call read_one_user_stat(filename,vdpdym)
  filename='vdpdzm.dat'
  call read_one_user_stat(filename,vdpdzm)

  filename='wdpdxm.dat'
  call read_one_user_stat(filename,wdpdxm)
  filename='wdpdym.dat'
  call read_one_user_stat(filename,wdpdym)
  filename='wdpdzm.dat'
  call read_one_user_stat(filename,wdpdzm)

  filename='uuum.dat'
  call read_one_user_stat(filename,uuum)
  filename='uvvm.dat'
  call read_one_user_stat(filename,uvvm)
  filename='uwwm.dat'
  call read_one_user_stat(filename,uwwm)

  filename='vuum.dat'
  call read_one_user_stat(filename,vuum)
  filename='vvvm.dat'
  call read_one_user_stat(filename,vvvm)
  filename='vwwm.dat'
  call read_one_user_stat(filename,vwwm)

  filename='wuum.dat'
  call read_one_user_stat(filename,wuum)
  filename='wvvm.dat'
  call read_one_user_stat(filename,wvvm)
  filename='wwwm.dat'
  call read_one_user_stat(filename,wwwm)

  filename='uvwm.dat'
  call read_one_user_stat(filename,uvwm)

  filename='dudx2m.dat'
  call read_one_user_stat(filename,dudx2m)
  filename='dudy2m.dat'
  call read_one_user_stat(filename,dudy2m)
  filename='dudz2m.dat'
  call read_one_user_stat(filename,dudz2m)

  filename='dvdx2m.dat'
  call read_one_user_stat(filename,dvdx2m)
  filename='dvdy2m.dat'
  call read_one_user_stat(filename,dvdy2m)
  filename='dvdz2m.dat'
  call read_one_user_stat(filename,dvdz2m)

  filename='dwdx2m.dat'
  call read_one_user_stat(filename,dwdx2m)
  filename='dwdy2m.dat'
  call read_one_user_stat(filename,dwdy2m)
  filename='dwdz2m.dat'
  call read_one_user_stat(filename,dwdz2m)

  filename='dudvdxm.dat'
  call read_one_user_stat(filename,dudvdxm)
  filename='dudwdxm.dat'
  call read_one_user_stat(filename,dudwdxm)
  filename='dvdwdxm.dat'
  call read_one_user_stat(filename,dvdwdxm)

  filename='dudvdym.dat'
  call read_one_user_stat(filename,dudvdym)
  filename='dudwdym.dat'
  call read_one_user_stat(filename,dudwdym)
  filename='dvdwdym.dat'
  call read_one_user_stat(filename,dvdwdym)

  filename='dudvdzm.dat'
  call read_one_user_stat(filename,dudvdzm)
  filename='dudwdzm.dat'
  call read_one_user_stat(filename,dudwdzm)
  filename='dvdwdzm.dat'
  call read_one_user_stat(filename,dvdwdzm)

  filename='dudxdudym.dat'
  call read_one_user_stat(filename,dudxdudym)
  filename='dudxdudzm.dat'
  call read_one_user_stat(filename,dudxdudzm)
  filename='dudydudzm.dat'
  call read_one_user_stat(filename,dudydudzm)

  filename='dvdxdvdym.dat'
  call read_one_user_stat(filename,dvdxdvdym)
  filename='dvdxdvdzm.dat'
  call read_one_user_stat(filename,dvdxdvdzm)
  filename='dvdydvdzm.dat'
  call read_one_user_stat(filename,dvdydvdzm)

  filename='dwdxdwdym.dat'
  call read_one_user_stat(filename,dwdxdwdym)
  filename='dwdxdwdzm.dat'
  call read_one_user_stat(filename,dwdxdwdzm)
  filename='dwdydwdzm.dat'
  call read_one_user_stat(filename,dwdydwdzm)

  filename='dudxdvdym.dat'
  call read_one_user_stat(filename,dudxdvdym)
  filename='dudydvdxm.dat'
  call read_one_user_stat(filename,dudydvdxm)

  filename='dudxdvdzm.dat'
  call read_one_user_stat(filename,dudxdvdzm)
  filename='dudzdvdxm.dat'
  call read_one_user_stat(filename,dudzdvdxm)

  filename='dudydvdzm.dat'
  call read_one_user_stat(filename,dudydvdzm)
  filename='dudzdvdym.dat'
  call read_one_user_stat(filename,dudzdvdym)

  filename='dudxdwdym.dat'
  call read_one_user_stat(filename,dudxdwdym)
  filename='dudydwdxm.dat'
  call read_one_user_stat(filename,dudydwdxm)

  filename='dudxdwdzm.dat'
  call read_one_user_stat(filename,dudxdwdzm)
  filename='dudzdwdxm.dat'
  call read_one_user_stat(filename,dudzdwdxm)

  filename='dudydwdzm.dat'
  call read_one_user_stat(filename,dudydwdzm)
  filename='dudzdwdym.dat'
  call read_one_user_stat(filename,dudzdwdym)

  filename='dvdxdwdym.dat'
  call read_one_user_stat(filename,dvdxdwdym)
  filename='dvdydwdxm.dat'
  call read_one_user_stat(filename,dvdydwdxm)

  filename='dvdxdwdzm.dat'
  call read_one_user_stat(filename,dvdxdwdzm)
  filename='dvdzdwdxm.dat'
  call read_one_user_stat(filename,dvdzdwdxm)

  filename='dvdydwdzm.dat'
  call read_one_user_stat(filename,dvdydwdzm)
  filename='dvdzdwdym.dat'
  call read_one_user_stat(filename,dvdzdwdym)

  filename='vdudxm.dat'
  call read_one_user_stat(filename,vdudxm)
  filename='vdvdxm.dat'
  call read_one_user_stat(filename,vdvdxm)
  filename='vdwdxm.dat'
  call read_one_user_stat(filename,vdwdxm)
  filename='vdudym.dat'
  call read_one_user_stat(filename,vdudym)
  filename='vdvdym.dat'
  call read_one_user_stat(filename,vdvdym)
  filename='vdwdym.dat'
  call read_one_user_stat(filename,vdwdym)
  filename='vdudzm.dat'
  call read_one_user_stat(filename,vdudzm)
  filename='vdvdzm.dat'
  call read_one_user_stat(filename,vdvdzm)
  filename='vdwdzm.dat'
  call read_one_user_stat(filename,vdwdzm)

  filename='dudx_dududxm.dat'
  call read_one_user_stat(filename,dudx_dududxm)
  filename='dudx_dududym.dat'
  call read_one_user_stat(filename,dudx_dududym)
  filename='dudx_dududzm.dat'
  call read_one_user_stat(filename,dudx_dududzm)

  filename='dudy_dudvdxm.dat'
  call read_one_user_stat(filename,dudy_dudvdxm)
  filename='dudy_dudvdym.dat'
  call read_one_user_stat(filename,dudy_dudvdym)
  filename='dudy_dudvdzm.dat'
  call read_one_user_stat(filename,dudy_dudvdzm)

  filename='dudz_dudwdxm.dat'
  call read_one_user_stat(filename,dudz_dudwdxm)
  filename='dudz_dudwdym.dat'
  call read_one_user_stat(filename,dudz_dudwdym)
  filename='dudz_dudwdzm.dat'
  call read_one_user_stat(filename,dudz_dudwdzm)

  filename='dvdx_dvdudxm.dat'
  call read_one_user_stat(filename,dvdx_dvdudxm)
  filename='dvdx_dvdudym.dat'
  call read_one_user_stat(filename,dvdx_dvdudym)
  filename='dvdx_dvdudzm.dat'
  call read_one_user_stat(filename,dvdx_dvdudzm)

  filename='dvdy_dvdvdxm.dat'
  call read_one_user_stat(filename,dvdy_dvdvdxm)
  filename='dvdy_dvdvdym.dat'
  call read_one_user_stat(filename,dvdy_dvdvdym)
  filename='dvdy_dvdvdzm.dat'
  call read_one_user_stat(filename,dvdy_dvdvdzm)

  filename='dvdz_dvdwdxm.dat'
  call read_one_user_stat(filename,dvdz_dvdwdxm)
  filename='dvdz_dvdwdym.dat'
  call read_one_user_stat(filename,dvdz_dvdwdym)
  filename='dvdz_dvdwdzm.dat'
  call read_one_user_stat(filename,dvdz_dvdwdzm)

  filename='dwdx_dwdudxm.dat'
  call read_one_user_stat(filename,dwdx_dwdudxm)
  filename='dwdx_dwdudym.dat'
  call read_one_user_stat(filename,dwdx_dwdudym)
  filename='dwdx_dwdudzm.dat'
  call read_one_user_stat(filename,dwdx_dwdudzm)

  filename='dwdy_dwdvdxm.dat'
  call read_one_user_stat(filename,dwdy_dwdvdxm)
  filename='dwdy_dwdvdym.dat'
  call read_one_user_stat(filename,dwdy_dwdvdym)
  filename='dwdy_dwdvdzm.dat'
  call read_one_user_stat(filename,dwdy_dwdvdzm)

  filename='dwdz_dwdwdxm.dat'
  call read_one_user_stat(filename,dwdz_dwdwdxm)
  filename='dwdz_dwdwdym.dat'
  call read_one_user_stat(filename,dwdz_dwdwdym)
  filename='dwdz_dwdwdzm.dat'
  call read_one_user_stat(filename,dwdz_dwdwdzm)

  filename='dudx_dudvdxm.dat'
  call read_one_user_stat(filename,dudx_dudvdxm)
  filename='dudx_dudvdym.dat'
  call read_one_user_stat(filename,dudx_dudvdym)
  filename='dudx_dudvdzm.dat'
  call read_one_user_stat(filename,dudx_dudvdzm)

  filename='dudy_dvdvdxm.dat'
  call read_one_user_stat(filename,dudy_dvdvdxm)
  filename='dudy_dvdvdym.dat'
  call read_one_user_stat(filename,dudy_dvdvdym)
  filename='dudy_dvdvdzm.dat'
  call read_one_user_stat(filename,dudy_dvdvdzm)

  filename='dudz_dvdwdxm.dat'
  call read_one_user_stat(filename,dudz_dvdwdxm)
  filename='dudz_dvdwdym.dat'
  call read_one_user_stat(filename,dudz_dvdwdym)
  filename='dudz_dvdwdzm.dat'
  call read_one_user_stat(filename,dudz_dvdwdzm)

  filename='dvdx_dududxm.dat'
  call read_one_user_stat(filename,dvdx_dududxm)
  filename='dvdx_dududym.dat'
  call read_one_user_stat(filename,dvdx_dududym)
  filename='dvdx_dududzm.dat'
  call read_one_user_stat(filename,dvdx_dududzm)

  filename='dvdy_dudvdxm.dat'
  call read_one_user_stat(filename,dvdy_dudvdxm)
  filename='dvdy_dudvdym.dat'
  call read_one_user_stat(filename,dvdy_dudvdym)
  filename='dvdy_dudvdzm.dat'
  call read_one_user_stat(filename,dvdy_dudvdzm)

  filename='dvdz_dudwdxm.dat'
  call read_one_user_stat(filename,dvdz_dudwdxm)
  filename='dvdz_dudwdym.dat'
  call read_one_user_stat(filename,dvdz_dudwdym)
  filename='dvdz_dudwdzm.dat'
  call read_one_user_stat(filename,dvdz_dudwdzm)

  filename='vdudx2m.dat'
  call read_one_user_stat(filename,vdudx2m)
  filename='vdudy2m.dat'
  call read_one_user_stat(filename,vdudy2m)
  filename='vdudz2m.dat'
  call read_one_user_stat(filename,vdudz2m)

  filename='vdvdx2m.dat'
  call read_one_user_stat(filename,vdvdx2m)
  filename='vdvdy2m.dat'
  call read_one_user_stat(filename,vdvdy2m)
  filename='vdvdz2m.dat'
  call read_one_user_stat(filename,vdvdz2m)

  filename='vdwdx2m.dat'
  call read_one_user_stat(filename,vdwdx2m)
  filename='vdwdy2m.dat'
  call read_one_user_stat(filename,vdwdy2m)
  filename='vdwdz2m.dat'
  call read_one_user_stat(filename,vdwdz2m)

  filename='vdudvdxm.dat'
  call read_one_user_stat(filename,vdudvdxm)
  filename='vdudvdym.dat'
  call read_one_user_stat(filename,vdudvdym)
  filename='vdudvdzm.dat'
  call read_one_user_stat(filename,vdudvdzm)

  filename='dpdvdxm.dat'
  call read_one_user_stat(filename,dpdvdxm)
  filename='dpdvdym.dat'
  call read_one_user_stat(filename,dpdvdym)
  filename='dpdvdzm.dat'
  call read_one_user_stat(filename,dpdvdzm)

  filename='dpdudxm.dat'
  call read_one_user_stat(filename,dpdudxm)
  filename='dpdudym.dat'
  call read_one_user_stat(filename,dpdudym)
  filename='dpdudzm.dat'
  call read_one_user_stat(filename,dpdudzm)

  filename='dpdxdudyxm.dat'
  call read_one_user_stat(filename,dpdxdudyxm)
  filename='dpdydudyym.dat'
  call read_one_user_stat(filename,dpdydudyym)
  filename='dpdzdudyzm.dat'
  call read_one_user_stat(filename,dpdzdudyzm)

  filename='dpdxdvdxxm.dat'
  call read_one_user_stat(filename,dpdxdvdxxm)
  filename='dpdydvdxym.dat'
  call read_one_user_stat(filename,dpdydvdxym)
  filename='dpdzdvdxzm.dat'
  call read_one_user_stat(filename,dpdzdvdxzm)

  filename='dudxx2m.dat'
  call read_one_user_stat(filename,dudxx2m)
  filename='dudyy2m.dat'
  call read_one_user_stat(filename,dudyy2m)
  filename='dudzz2m.dat'
  call read_one_user_stat(filename,dudzz2m)
  filename='dudxy2m.dat'
  call read_one_user_stat(filename,dudxy2m)
  filename='dudxz2m.dat'
  call read_one_user_stat(filename,dudxz2m)
  filename='dudyz2m.dat'
  call read_one_user_stat(filename,dudyz2m)

  filename='dvdxx2m.dat'
  call read_one_user_stat(filename,dvdxx2m)
  filename='dvdyy2m.dat'
  call read_one_user_stat(filename,dvdyy2m)
  filename='dvdzz2m.dat'
  call read_one_user_stat(filename,dvdzz2m)
  filename='dvdxy2m.dat'
  call read_one_user_stat(filename,dvdxy2m)
  filename='dvdxz2m.dat'
  call read_one_user_stat(filename,dvdxz2m)
  filename='dvdyz2m.dat'
  call read_one_user_stat(filename,dvdyz2m)

  filename='dwdxx2m.dat'
  call read_one_user_stat(filename,dwdxx2m)
  filename='dwdyy2m.dat'
  call read_one_user_stat(filename,dwdyy2m)
  filename='dwdzz2m.dat'
  call read_one_user_stat(filename,dwdzz2m)
  filename='dwdxy2m.dat'
  call read_one_user_stat(filename,dwdxy2m)
  filename='dwdxz2m.dat'
  call read_one_user_stat(filename,dwdxz2m)
  filename='dwdyz2m.dat'
  call read_one_user_stat(filename,dwdyz2m)

  filename='dudvdxx2m.dat'
  call read_one_user_stat(filename,dudvdxx2m)
  filename='dudvdyy2m.dat'
  call read_one_user_stat(filename,dudvdyy2m)
  filename='dudvdzz2m.dat'
  call read_one_user_stat(filename,dudvdzz2m)
  filename='dudvdxy2m.dat'
  call read_one_user_stat(filename,dudvdxy2m)
  filename='dudvdxz2m.dat'
  call read_one_user_stat(filename,dudvdxz2m)
  filename='dudvdyz2m.dat'
  call read_one_user_stat(filename,dudvdyz2m)

  if (iscalar.eq.1) then
    filename='phim.dat'
    call read_one_user_stat(filename,phim)
    !
    filename='dphidym.dat'
    call read_one_user_stat(filename,dphidym)
    !
    filename='dphidyym.dat'
    call read_one_user_stat(filename,dphidyym)
    !
    filename='dphidyyym.dat'
    call read_one_user_stat(filename,dphidyyym)
    !
    filename='uphim.dat'
    call read_one_user_stat(filename,uphim)
    filename='vphim.dat'
    call read_one_user_stat(filename,vphim)
    filename='wphim.dat'
    call read_one_user_stat(filename,wphim)
    filename='phiphim.dat'
    call read_one_user_stat(filename,phiphim)
    !
    filename='phidpdxm.dat'
    call read_one_user_stat(filename,phidpdxm)
    filename='phidpdym.dat'
    call read_one_user_stat(filename,phidpdym)
    filename='phidpdzm.dat'
    call read_one_user_stat(filename,phidpdzm)
    !
    filename='dphidx2m.dat'
    call read_one_user_stat(filename,dphidx2m)
    filename='dphidy2m.dat'
    call read_one_user_stat(filename,dphidy2m)
    filename='dphidz2m.dat'
    call read_one_user_stat(filename,dphidz2m)
    !
    filename='dudphidxm.dat'
    call read_one_user_stat(filename,dudphidxm)
    filename='dvdphidxm.dat'
    call read_one_user_stat(filename,dvdphidxm)
    filename='dwdphidxm.dat'
    call read_one_user_stat(filename,dwdphidxm)
    !
    filename='dudphidym.dat'
    call read_one_user_stat(filename,dudphidym)
    filename='dvdphidym.dat'
    call read_one_user_stat(filename,dvdphidym)
    filename='dwdphidym.dat'
    call read_one_user_stat(filename,dwdphidym)
    !
    filename='dudphidzm.dat'
    call read_one_user_stat(filename,dudphidzm)
    filename='dvdphidzm.dat'
    call read_one_user_stat(filename,dvdphidzm)
    filename='dwdphidzm.dat'
    call read_one_user_stat(filename,dwdphidzm)
    !
    filename='udtdxxm.dat'
    call read_one_user_stat(filename,udtdxxm)
    filename='udtdyym.dat'
    call read_one_user_stat(filename,udtdyym)
    filename='udtdzzm.dat'
    call read_one_user_stat(filename,udtdzzm)
    !
    filename='vdtdxxm.dat'
    call read_one_user_stat(filename,vdtdxxm)
    filename='vdtdyym.dat'
    call read_one_user_stat(filename,vdtdyym)
    filename='vdtdzzm.dat'
    call read_one_user_stat(filename,vdtdzzm)
    !
    filename='wdtdxxm.dat'
    call read_one_user_stat(filename,wdtdxxm)
    filename='wdtdyym.dat'
    call read_one_user_stat(filename,wdtdyym)
    filename='wdtdzzm.dat'
    call read_one_user_stat(filename,wdtdzzm)
    !
    filename='uphi2m.dat'
    call read_one_user_stat(filename,uphi2m)
    filename='vphi2m.dat'
    call read_one_user_stat(filename,vphi2m)
    filename='wphi2m.dat'
    call read_one_user_stat(filename,wphi2m)
    !
    filename='phiuum.dat'
    call read_one_user_stat(filename,phiuum)
    filename='phivvm.dat'
    call read_one_user_stat(filename,phivvm)
    filename='phiwwm.dat'
    call read_one_user_stat(filename,phiwwm)
    !
    filename='phiuvm.dat'
    call read_one_user_stat(filename,phiuvm)
    filename='phiuwm.dat'
    call read_one_user_stat(filename,phiuwm)
    filename='phivwm.dat'
    call read_one_user_stat(filename,phivwm)
    !
    filename='dtdxdtdym.dat'
    call read_one_user_stat(filename,dtdxdtdym)
    filename='dtdxdtdzm.dat'
    call read_one_user_stat(filename,dtdxdtdzm)
    filename='dtdydtdzm.dat'
    call read_one_user_stat(filename,dtdydtdzm)
    !
    filename='dtdxdudym.dat'
    call read_one_user_stat(filename,dtdxdudym)
    filename='dtdydudxm.dat'
    call read_one_user_stat(filename,dtdydudxm)
    filename='dtdxdudzm.dat'
    call read_one_user_stat(filename,dtdxdudzm)
    filename='dtdzdudxm.dat'
    call read_one_user_stat(filename,dtdzdudxm)
    filename='dtdydudzm.dat'
    call read_one_user_stat(filename,dtdydudzm)
    filename='dtdzdudym.dat'
    call read_one_user_stat(filename,dtdzdudym)
    !
    filename='dtdxdvdym.dat'
    call read_one_user_stat(filename,dtdxdvdym)
    filename='dtdydvdxm.dat'
    call read_one_user_stat(filename,dtdydvdxm)
    filename='dtdxdvdzm.dat'
    call read_one_user_stat(filename,dtdxdvdzm)
    filename='dtdzdvdxm.dat'
    call read_one_user_stat(filename,dtdzdvdxm)
    filename='dtdydvdzm.dat'
    call read_one_user_stat(filename,dtdydvdzm)
    filename='dtdzdvdym.dat'
    call read_one_user_stat(filename,dtdzdvdym)
    !
    filename='dtdxdwdym.dat'
    call read_one_user_stat(filename,dtdxdwdym)
    filename='dtdydwdxm.dat'
    call read_one_user_stat(filename,dtdydwdxm)
    filename='dtdxdwdzm.dat'
    call read_one_user_stat(filename,dtdxdwdzm)
    filename='dtdzdwdxm.dat'
    call read_one_user_stat(filename,dtdzdwdxm)
    filename='dtdydwdzm.dat'
    call read_one_user_stat(filename,dtdydwdzm)
    filename='dtdzdwdym.dat'
    call read_one_user_stat(filename,dtdzdwdym)
    !
    filename='vdtdxm.dat'
    call read_one_user_stat(filename,vdtdxm)
    filename='vdtdym.dat'
    call read_one_user_stat(filename,vdtdym)
    filename='vdtdzm.dat'
    call read_one_user_stat(filename,vdtdzm)
    !
    filename='dudx_dudtdxm.dat'
    call read_one_user_stat(filename,dudx_dudtdxm)
    filename='dudx_dudtdym.dat'
    call read_one_user_stat(filename,dudx_dudtdym)
    filename='dudx_dudtdzm.dat'
    call read_one_user_stat(filename,dudx_dudtdzm)
    filename='dudy_dvdtdxm.dat'
    call read_one_user_stat(filename,dudy_dvdtdxm)
    filename='dudy_dvdtdym.dat'
    call read_one_user_stat(filename,dudy_dvdtdym)
    filename='dudy_dvdtdzm.dat'
    call read_one_user_stat(filename,dudy_dvdtdzm)
    filename='dudz_dwdtdxm.dat'
    call read_one_user_stat(filename,dudz_dwdtdxm)
    filename='dudz_dwdtdym.dat'
    call read_one_user_stat(filename,dudz_dwdtdym)
    filename='dudz_dwdtdzm.dat'
    call read_one_user_stat(filename,dudz_dwdtdzm)
    !
    filename='dtdx_dududym.dat'
    call read_one_user_stat(filename,dtdx_dududym)
    filename='dtdx_dududzm.dat'
    call read_one_user_stat(filename,dtdx_dududzm)
    filename='dtdy_dudvdxm.dat'
    call read_one_user_stat(filename,dtdy_dudvdxm)
    filename='dtdy_dudvdzm.dat'
    call read_one_user_stat(filename,dtdy_dudvdzm)
    filename='dtdz_dudwdxm.dat'
    call read_one_user_stat(filename,dtdz_dudwdxm)
    filename='dtdz_dudwdym.dat'
    call read_one_user_stat(filename,dtdz_dudwdym)
    !
    filename='dvdx_dudtdxm.dat'
    call read_one_user_stat(filename,dvdx_dudtdxm)
    filename='dvdx_dudtdym.dat'
    call read_one_user_stat(filename,dvdx_dudtdym)
    filename='dvdx_dudtdzm.dat'
    call read_one_user_stat(filename,dvdx_dudtdzm)
    filename='dvdy_dvdtdxm.dat'
    call read_one_user_stat(filename,dvdy_dvdtdxm)
    filename='dvdy_dvdtdym.dat'
    call read_one_user_stat(filename,dvdy_dvdtdym)
    filename='dvdy_dvdtdzm.dat'
    call read_one_user_stat(filename,dvdy_dvdtdzm)
    filename='dvdz_dwdtdxm.dat'
    call read_one_user_stat(filename,dvdz_dwdtdxm)
    filename='dvdz_dwdtdym.dat'
    call read_one_user_stat(filename,dvdz_dwdtdym)
    filename='dvdz_dwdtdzm.dat'
    call read_one_user_stat(filename,dvdz_dwdtdzm)
    !
    filename='dtdx_dudvdym.dat'
    call read_one_user_stat(filename,dtdx_dudvdym)
    filename='dtdx_dudvdzm.dat'
    call read_one_user_stat(filename,dtdx_dudvdzm)
    filename='dtdy_dvdvdxm.dat'
    call read_one_user_stat(filename,dtdy_dvdvdxm)
    filename='dtdy_dvdvdzm.dat'
    call read_one_user_stat(filename,dtdy_dvdvdzm)
    filename='dtdz_dvdwdxm.dat'
    call read_one_user_stat(filename,dtdz_dvdwdxm)
    filename='dtdz_dvdwdym.dat'
    call read_one_user_stat(filename,dtdz_dvdwdym)
    !
    filename='dtdx_dudtdxm.dat'
    call read_one_user_stat(filename,dtdx_dudtdxm)
    filename='dtdx_dudtdym.dat'
    call read_one_user_stat(filename,dtdx_dudtdym)
    filename='dtdx_dudtdzm.dat'
    call read_one_user_stat(filename,dtdx_dudtdzm)
    filename='dtdy_dvdtdxm.dat'
    call read_one_user_stat(filename,dtdy_dvdtdxm)
    filename='dtdy_dvdtdym.dat'
    call read_one_user_stat(filename,dtdy_dvdtdym)
    filename='dtdy_dvdtdzm.dat'
    call read_one_user_stat(filename,dtdy_dvdtdzm)
    filename='dtdz_dwdtdxm.dat'
    call read_one_user_stat(filename,dtdz_dwdtdxm)
    filename='dtdz_dwdtdym.dat'
    call read_one_user_stat(filename,dtdz_dwdtdym)
    filename='dtdz_dwdtdzm.dat'
    call read_one_user_stat(filename,dtdz_dwdtdzm)
    !
    filename='vdudtdxm.dat'
    call read_one_user_stat(filename,vdudtdxm)
    filename='vdudtdym.dat'
    call read_one_user_stat(filename,vdudtdym)
    filename='vdudtdzm.dat'
    call read_one_user_stat(filename,vdudtdzm)
    filename='vdvdtdxm.dat'
    call read_one_user_stat(filename,vdvdtdxm)
    filename='vdvdtdym.dat'
    call read_one_user_stat(filename,vdvdtdym)
    filename='vdvdtdzm.dat'
    call read_one_user_stat(filename,vdvdtdzm)
    filename='vdtdx2m.dat'
    call read_one_user_stat(filename,vdtdx2m)
    filename='vdtdy2m.dat'
    call read_one_user_stat(filename,vdtdy2m)
    filename='vdtdz2m.dat'
    call read_one_user_stat(filename,vdtdz2m)
    !
    filename='dpdtdxm.dat'
    call read_one_user_stat(filename,dpdtdxm)
    filename='dpdtdym.dat'
    call read_one_user_stat(filename,dpdtdym)
    filename='dpdtdzm.dat'
    call read_one_user_stat(filename,dpdtdzm)
    filename='dpdtdxxm.dat'
    call read_one_user_stat(filename,dpdtdxxm)
    filename='dpdtdxym.dat'
    call read_one_user_stat(filename,dpdtdxym)
    filename='dpdtdxzm.dat'
    call read_one_user_stat(filename,dpdtdxzm)
    filename='dpdtdyxm.dat'
    call read_one_user_stat(filename,dpdtdyxm)
    filename='dpdtdyym.dat'
    call read_one_user_stat(filename,dpdtdyym)
    filename='dpdtdyzm.dat'
    call read_one_user_stat(filename,dpdtdyzm)
    !
    filename='dudtdxxm.dat'
    call read_one_user_stat(filename,dudtdxxm)
    filename='dudtdyym.dat'
    call read_one_user_stat(filename,dudtdyym)
    filename='dudtdzzm.dat'
    call read_one_user_stat(filename,dudtdzzm)
    filename='dudtdxym.dat'
    call read_one_user_stat(filename,dudtdxym)
    filename='dudtdxzm.dat'
    call read_one_user_stat(filename,dudtdxzm)
    filename='dudtdyzm.dat'
    call read_one_user_stat(filename,dudtdyzm)
    !
    filename='dudx_dtdxxxm.dat'
    call read_one_user_stat(filename,dudx_dtdxxxm)
    filename='dudx_dtdyxxm.dat'
    call read_one_user_stat(filename,dudx_dtdyxxm)
    filename='dudx_dtdzxxm.dat'
    call read_one_user_stat(filename,dudx_dtdzxxm)
    filename='dudy_dtdxyym.dat'
    call read_one_user_stat(filename,dudy_dtdxyym)
    filename='dudy_dtdyyym.dat'
    call read_one_user_stat(filename,dudy_dtdyyym)
    filename='dudy_dtdzyym.dat'
    call read_one_user_stat(filename,dudy_dtdzyym)
    filename='dudz_dtdxzzm.dat'
    call read_one_user_stat(filename,dudz_dtdxzzm)
    filename='dudz_dtdyzzm.dat'
    call read_one_user_stat(filename,dudz_dtdyzzm)
    filename='dudz_dtdzzzm.dat'
    call read_one_user_stat(filename,dudz_dtdzzzm)
    !
    filename='dvdtdxxm.dat'
    call read_one_user_stat(filename,dvdtdxxm)
    filename='dvdtdyym.dat'
    call read_one_user_stat(filename,dvdtdyym)
    filename='dvdtdzzm.dat'
    call read_one_user_stat(filename,dvdtdzzm)
    filename='dvdtdxym.dat'
    call read_one_user_stat(filename,dvdtdxym)
    filename='dvdtdxzm.dat'
    call read_one_user_stat(filename,dvdtdxzm)
    filename='dvdtdyzm.dat'
    call read_one_user_stat(filename,dvdtdyzm)
    !
    filename='dvdx_dtdxxxm.dat'
    call read_one_user_stat(filename,dvdx_dtdxxxm)
    filename='dvdx_dtdyxxm.dat'
    call read_one_user_stat(filename,dvdx_dtdyxxm)
    filename='dvdx_dtdzxxm.dat'
    call read_one_user_stat(filename,dvdx_dtdzxxm)
    filename='dvdy_dtdxyym.dat'
    call read_one_user_stat(filename,dvdy_dtdxyym)
    filename='dvdy_dtdyyym.dat'
    call read_one_user_stat(filename,dvdy_dtdyyym)
    filename='dvdy_dtdzyym.dat'
    call read_one_user_stat(filename,dvdy_dtdzyym)
    filename='dvdz_dtdxzzm.dat'
    call read_one_user_stat(filename,dvdz_dtdxzzm)
    filename='dvdz_dtdyzzm.dat'
    call read_one_user_stat(filename,dvdz_dtdyzzm)
    filename='dvdz_dtdzzzm.dat'
    call read_one_user_stat(filename,dvdz_dtdzzzm)
    !
    filename='dtdxx2m.dat'
    call read_one_user_stat(filename,dtdxx2m)
    filename='dtdyy2m.dat'
    call read_one_user_stat(filename,dtdyy2m)
    filename='dtdzz2m.dat'
    call read_one_user_stat(filename,dtdzz2m)
    filename='dtdxy2m.dat'
    call read_one_user_stat(filename,dtdxy2m)
    filename='dtdxz2m.dat'
    call read_one_user_stat(filename,dtdxz2m)
    filename='dtdyz2m.dat'
    call read_one_user_stat(filename,dtdyz2m)
    !
  endif

  filename='presm.dat'
  call read_one_user_stat(filename,presm)
  filename='pres2m.dat'
  call read_one_user_stat(filename,pres2m)

end subroutine read_user_stats

subroutine read_one_user_stat(filename,var)

  USE decomp_2d, only : mytype, xsize, nrank
  use decomp_2d_io, only : decomp_2d_read_plane

  use mpi

  implicit none

  real(mytype), dimension(xst1:xen1,xst2:xen2,xst3:xen3), intent(out) :: var
  character(len=30) :: filename

  logical :: file_exist
  real(mytype), dimension(xst1:xst1,xst2:xen2,xst3:xen3) :: tmp
  integer :: i,j,k

  INQUIRE(FILE=filename, EXIST=file_exist)

  if (file_exist) then
    !
    call decomp_2d_read_plane(1,tmp,filename,decomp_user_stats)
    !
    do k=xst3,xen3
      do j=xst2,xen2
        do i=xst1,xen1
          var(i,j,k)=tmp(xst1,j,k)
        end do
      end do
    end do
    !
    if (nrank.eq.0) print*,'user_stats : reading average for file ',trim(filename),' done'
    !
  else
    !
    var = 0.
    !
    if (nrank.eq.0) print*,'user_stats : average file ',trim(filename),' is not existing, statistics may be incorrect'
    !
  endif

end subroutine read_one_user_stat

end module user_stats
