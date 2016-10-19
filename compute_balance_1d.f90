subroutine compute_balance_1d()

  use decomp_2d, only : mytype
  use param, only : iscalar
  use user_stats
  ! Pression sur grille décalée
!  use decomp_2d, only : transpose_z_to_y, transpose_y_to_x
!  use variables, only : ph3

  implicit none

  integer :: i,j,k,n
  ! Ordre 1
  real(mytype), dimension(xst2:xen2) :: u,v,w
  real(mytype), dimension(xst2:xen2) :: dpdx, dpdy, dpdz
  ! Ordre 2
  real(mytype), dimension(xst2:xen2) :: uu,vv,ww
  real(mytype), dimension(xst2:xen2) :: uv,uw,vw
  real(mytype), dimension(xst2:xen2) :: udpdx, udpdy, udpdz
  real(mytype), dimension(xst2:xen2) :: vdpdx, vdpdy, vdpdz
  real(mytype), dimension(xst2:xen2) :: wdpdx, wdpdy, wdpdz
  real(mytype), dimension(xst2:xen2) :: dudx2, dudy2, dudz2
  real(mytype), dimension(xst2:xen2) :: dvdx2, dvdy2, dvdz2
  real(mytype), dimension(xst2:xen2) :: dwdx2, dwdy2, dwdz2
  real(mytype), dimension(xst2:xen2) :: dudvdx, dudwdx, dvdwdx
  real(mytype), dimension(xst2:xen2) :: dudvdy, dudwdy, dvdwdy
  real(mytype), dimension(xst2:xen2) :: dudvdz, dudwdz, dvdwdz
  ! Ordre 3
  real(mytype), dimension(xst2:xen2) :: uuu, uvv, uww
  real(mytype), dimension(xst2:xen2) :: vuu, vvv, vww
  real(mytype), dimension(xst2:xen2) :: wuu, wvv, www
  real(mytype), dimension(xst2:xen2) :: uvw

  ! Scalaire
  ! Ordre 1
  real(mytype), dimension(xst2:xen2) :: phi
  real(mytype), dimension(xst2:xen2) :: dphidx, dphidy, dphidz
  real(mytype), dimension(xst2:xen2) :: deltaphi
  ! Ordre 2
  real(mytype), dimension(xst2:xen2) :: uphi, vphi, wphi, phiphi
  real(mytype), dimension(xst2:xen2) :: dphidx2, dphidy2, dphidz2
  real(mytype), dimension(xst2:xen2) :: phidpdx, phidpdy, phidpdz
  real(mytype), dimension(xst2:xen2) :: dudphidx, dvdphidx, dwdphidx
  real(mytype), dimension(xst2:xen2) :: dudphidy, dvdphidy, dwdphidy
  real(mytype), dimension(xst2:xen2) :: dudphidz, dvdphidz, dwdphidz
  real(mytype), dimension(xst2:xen2) :: udeltaphi, vdeltaphi, wdeltaphi
  ! Ordre 3
  real(mytype), dimension(xst2:xen2) :: uphi2,vphi2,wphi2
  real(mytype), dimension(xst2:xen2) :: phiuu,phivv,phiww
  real(mytype), dimension(xst2:xen2) :: phiuv,phiuw,phivw

  ! Pression
  real(mytype), dimension(xst2:xen2) :: pres, pres2

  ! Pression sur grille décalée
!  real(mytype), dimension(ph3%xst(2):ph3%xen(2)) :: presi, pres2i
!  real(mytype), dimension(ph3%xst(1):ph3%xen(1),ph3%xst(2):ph3%xen(2),ph3%xst(3):ph3%xen(3)) :: pp3mx, pp32mx
!  real(mytype), dimension(ph3%yst(1):ph3%yen(1),ph3%yst(2):ph3%yen(2),ph3%yst(3):ph3%yen(3)) :: pp3my

  ! Gradient vitesse
  real(mytype), dimension(xst2:xen2) :: dudx, dvdx, dwdx
  real(mytype), dimension(xst2:xen2) :: dudy, dvdy, dwdy
  real(mytype), dimension(xst2:xen2) :: dudz, dvdz, dwdz
  dudx=0.; dvdx=0.; dwdx=0.;
  dudz=0.; dvdz=0.; dwdz=0.;

  ! Ordre 1
  if (.true.) then
!
  ! Init Ordre 1
  u=0.;		v=0.;		w=0.
  dpdx=0.;	dpdy=0.;	dpdz=0.;
  dudy=0.;      dvdy=0.;        dwdy=0.;
!
  ! Compute Average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    u(j)	=u(j)	+um(i,j,k)
    v(j)	=v(j)	+vm(i,j,k)
    w(j)	=w(j)	+wm(i,j,k)
    dpdx(j)	=dpdx(j)+dpdxm(i,j,k)
    dpdy(j)	=dpdy(j)+dpdym(i,j,k)
    dpdz(j)	=dpdz(j)+dpdzm(i,j,k)
    dudy(j)     =dudy(j)+dudym(i,j,k)
    dvdy(j)     =dvdy(j)+dvdym(i,j,k)
    dwdy(j)     =dwdy(j)+dwdym(i,j,k)
  end do
  end do
  end do
!
  ! Average ordre 1
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  u=u/n;	v=v/n;		w=w/n;
  dpdx=dpdx/n;	dpdy=dpdy/n;	dpdz=dpdz/n;
  dudy=dudy/n;  dvdy=dvdy/n;    dwdy=dwdy/n;
!
  ! Write data
  call write_myaverage_real_1d(u   ,xen2-xst2+1,n,xst2,xen2,   'um1d')
  call write_myaverage_real_1d(v   ,xen2-xst2+1,n,xst2,xen2,   'vm1d')
  call write_myaverage_real_1d(w   ,xen2-xst2+1,n,xst2,xen2,   'wm1d')
  call write_myaverage_real_1d(dpdx,xen2-xst2+1,n,xst2,xen2,'dpdxm1d')
  call write_myaverage_real_1d(dpdy,xen2-xst2+1,n,xst2,xen2,'dpdym1d')
  call write_myaverage_real_1d(dpdz,xen2-xst2+1,n,xst2,xen2,'dpdzm1d')
  call write_myaverage_real_1d(dudy,xen2-xst2+1,n,xst2,xen2,'dudym1d')
  call write_myaverage_real_1d(dvdy,xen2-xst2+1,n,xst2,xen2,'dvdym1d')
  call write_myaverage_real_1d(dwdy,xen2-xst2+1,n,xst2,xen2,'dwdym1d')
!
  endif ! Ordre 1

  ! Ordre 2
  if (.true.) then
!
  ! Init ordre 2
  uu=0.; 	 vv=0.; 	ww=0.; 
  uv=0.; 	 uw=0.; 	vw=0.; 
  udpdx=0.; 	 udpdy=0.; 	udpdz=0.; 
  vdpdx=0.; 	 vdpdy=0.; 	vdpdz=0.; 
  wdpdx=0.; 	 wdpdy=0.; 	wdpdz=0.; 
  dudx2=0.; 	 dudy2=0.; 	dudz2=0.; 
  dvdx2=0.; 	 dvdy2=0.; 	dvdz2=0.; 
  dwdx2=0.; 	 dwdy2=0.; 	dwdz2=0.; 
  dudvdx=0.; 	 dudwdx=0.; 	dvdwdx=0.; 
  dudvdy=0.; 	 dudwdy=0.; 	dvdwdy=0.; 
  dudvdz=0.; 	 dudwdz=0.; 	dvdwdz=0.;
!
  ! Compute average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    uu(j)	=uu(j)		+uum(i,j,k)
    vv(j)	=vv(j)		+vvm(i,j,k)
    ww(j)	=ww(j)		+wwm(i,j,k)
    uv(j)	=uv(j)		+uvm(i,j,k)
    uw(j)	=uw(j)		+uwm(i,j,k)
    vw(j)	=vw(j)		+vwm(i,j,k)
!
    udpdx(j)	=udpdx(j)	+udpdxm(i,j,k)
    udpdy(j)	=udpdy(j)	+udpdym(i,j,k)
    udpdz(j)	=udpdz(j)	+udpdzm(i,j,k)
!
    vdpdx(j)	=vdpdx(j)	+vdpdxm(i,j,k)
    vdpdy(j)	=vdpdy(j)	+vdpdym(i,j,k)
    vdpdz(j)	=vdpdz(j)	+vdpdzm(i,j,k)
!
    wdpdx(j)	=wdpdx(j)	+wdpdxm(i,j,k)
    wdpdy(j)	=wdpdy(j)	+wdpdym(i,j,k)
    wdpdz(j)	=wdpdz(j)	+wdpdzm(i,j,k)
!
    dudx2(j)	=dudx2(j)	+dudx2m(i,j,k)
    dudy2(j)	=dudy2(j)	+dudy2m(i,j,k)
    dudz2(j)	=dudz2(j)	+dudz2m(i,j,k)
!
    dvdx2(j)	=dvdx2(j)	+dvdx2m(i,j,k)
    dvdy2(j)	=dvdy2(j)	+dvdy2m(i,j,k)
    dvdz2(j)	=dvdz2(j)	+dvdz2m(i,j,k)
!
    dwdx2(j)	=dwdx2(j)	+dwdx2m(i,j,k)
    dwdy2(j)	=dwdy2(j)	+dwdy2m(i,j,k)
    dwdz2(j)	=dwdz2(j)	+dwdz2m(i,j,k)
!
    dudvdx(j)	 =dudvdx(j)	+dudvdxm(i,j,k)
    dudwdx(j)	 =dudwdx(j)	+dudwdxm(i,j,k)
    dvdwdx(j)	 =dvdwdx(j)	+dvdwdxm(i,j,k)
!
    dudvdy(j)	 =dudvdy(j)	+dudvdym(i,j,k)
    dudwdy(j)	 =dudwdy(j)	+dudwdym(i,j,k)
    dvdwdy(j)	 =dvdwdy(j)	+dvdwdym(i,j,k)
!
    dudvdz(j)	 =dudvdz(j)	+dudvdzm(i,j,k)
    dudwdz(j)	 =dudwdz(j)	+dudwdzm(i,j,k)
    dvdwdz(j)	 =dvdwdz(j)	+dvdwdzm(i,j,k)
  end do
  end do
  end do
!
  ! Average ordre 2
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  uu=uu/n;		vv=vv/n;		ww=ww/n;
  uv=uv/n;		uw=uw/n;		vw=vw/n;
  udpdx=udpdx/n;	udpdy=udpdy/n;		udpdz=udpdz/n;
  vdpdx=vdpdx/n;	vdpdy=vdpdy/n;		vdpdz=vdpdz/n;
  wdpdx=wdpdx/n;	wdpdy=wdpdy/n;		wdpdz=wdpdz/n;
  dudx2=dudx2/n;	dudy2=dudy2/n;		dudz2=dudz2/n;
  dvdx2=dvdx2/n;	dvdy2=dvdy2/n;		dvdz2=dvdz2/n;
  dwdx2=dwdx2/n;	dwdy2=dwdy2/n;		dwdz2=dwdz2/n;
  dudvdx=dudvdx/n;	dudwdx=dudwdx/n;	dvdwdx=dvdwdx/n;
  dudvdy=dudvdy/n;	dudwdy=dudwdy/n;	dvdwdy=dvdwdy/n;
  dudvdz=dudvdz/n;	dudwdz=dudwdz/n;	dvdwdz=dvdwdz/n;
!
  ! Write data
  call write_myaverage_real_1d(uu-u*u,xen2-xst2+1,n,xst2,xen2,'uum1d')
  call write_myaverage_real_1d(vv-v*v,xen2-xst2+1,n,xst2,xen2,'vvm1d')
  call write_myaverage_real_1d(ww-w*w,xen2-xst2+1,n,xst2,xen2,'wwm1d')
  call write_myaverage_real_1d(uv-u*v,xen2-xst2+1,n,xst2,xen2,'uvm1d')
  call write_myaverage_real_1d(uw-u*w,xen2-xst2+1,n,xst2,xen2,'uwm1d')
  call write_myaverage_real_1d(vw-v*w,xen2-xst2+1,n,xst2,xen2,'vwm1d')
!
  call write_myaverage_real_1d(udpdx-u*dpdx,xen2-xst2+1,n,xst2,xen2,'udpdxm1d')
  call write_myaverage_real_1d(udpdy-u*dpdy,xen2-xst2+1,n,xst2,xen2,'udpdym1d')
  call write_myaverage_real_1d(udpdz-u*dpdz,xen2-xst2+1,n,xst2,xen2,'udpdzm1d')
  call write_myaverage_real_1d(vdpdx-v*dpdx,xen2-xst2+1,n,xst2,xen2,'vdpdxm1d')
  call write_myaverage_real_1d(vdpdy-v*dpdy,xen2-xst2+1,n,xst2,xen2,'vdpdym1d')
  call write_myaverage_real_1d(vdpdz-v*dpdz,xen2-xst2+1,n,xst2,xen2,'vdpdzm1d')
  call write_myaverage_real_1d(wdpdx-w*dpdx,xen2-xst2+1,n,xst2,xen2,'wdpdxm1d')
  call write_myaverage_real_1d(wdpdy-w*dpdy,xen2-xst2+1,n,xst2,xen2,'wdpdym1d')
  call write_myaverage_real_1d(wdpdz-w*dpdz,xen2-xst2+1,n,xst2,xen2,'wdpdzm1d')
!
  call write_myaverage_real_1d(dudx2-dudx*dudx,xen2-xst2+1,n,xst2,xen2,'dudx2m1d')
  call write_myaverage_real_1d(dudy2-dudy*dudy,xen2-xst2+1,n,xst2,xen2,'dudy2m1d')
  call write_myaverage_real_1d(dudz2-dudz*dudz,xen2-xst2+1,n,xst2,xen2,'dudz2m1d')
  call write_myaverage_real_1d(dvdx2-dvdx*dvdx,xen2-xst2+1,n,xst2,xen2,'dvdx2m1d')
  call write_myaverage_real_1d(dvdy2-dvdy*dvdy,xen2-xst2+1,n,xst2,xen2,'dvdy2m1d')
  call write_myaverage_real_1d(dvdz2-dvdz*dvdz,xen2-xst2+1,n,xst2,xen2,'dvdz2m1d')
  call write_myaverage_real_1d(dwdx2-dwdx*dwdx,xen2-xst2+1,n,xst2,xen2,'dwdx2m1d')
  call write_myaverage_real_1d(dwdy2-dwdy*dwdy,xen2-xst2+1,n,xst2,xen2,'dwdy2m1d')
  call write_myaverage_real_1d(dwdz2-dwdz*dwdz,xen2-xst2+1,n,xst2,xen2,'dwdz2m1d')
!
  call write_myaverage_real_1d(dudvdx-dudx*dvdx,xen2-xst2+1,n,xst2,xen2,'dudvdxm1d')
  call write_myaverage_real_1d(dudvdy-dudy*dvdy,xen2-xst2+1,n,xst2,xen2,'dudvdym1d')
  call write_myaverage_real_1d(dudvdz-dudz*dvdz,xen2-xst2+1,n,xst2,xen2,'dudvdzm1d')
  call write_myaverage_real_1d(dudwdx-dudx*dwdx,xen2-xst2+1,n,xst2,xen2,'dudwdxm1d')
  call write_myaverage_real_1d(dudwdy-dudy*dwdy,xen2-xst2+1,n,xst2,xen2,'dudwdym1d')
  call write_myaverage_real_1d(dudwdz-dudz*dwdz,xen2-xst2+1,n,xst2,xen2,'dudwdzm1d')
  call write_myaverage_real_1d(dvdwdx-dvdx*dwdx,xen2-xst2+1,n,xst2,xen2,'dvdwdxm1d')
  call write_myaverage_real_1d(dvdwdy-dvdy*dwdy,xen2-xst2+1,n,xst2,xen2,'dvdwdym1d')
  call write_myaverage_real_1d(dvdwdz-dvdz*dwdz,xen2-xst2+1,n,xst2,xen2,'dvdwdzm1d')!
!
  endif ! Ordre 2

  ! Ordre 3
  if (.true.) then
!
  ! Init ordre 3
  uuu=0; uvv=0; uww=0;
  vuu=0; vvv=0; vww=0;
  wuu=0; wvv=0; www=0;
  uvw=0;
!
  ! Compute average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    uuu(j)	=uuu(j)		+uuum(i,j,k)
    uvv(j)	=uvv(j)		+uvvm(i,j,k)
    uww(j)	=uww(j)		+uwwm(i,j,k)
!
    vuu(j)	=vuu(j)		+vuum(i,j,k)
    vvv(j)	=vvv(j)		+vvvm(i,j,k)
    vww(j)	=vww(j)		+vwwm(i,j,k)
!
    wuu(j)	=wuu(j)		+wuum(i,j,k)
    wvv(j)	=wvv(j)		+wvvm(i,j,k)
    www(j)	=www(j)		+wwwm(i,j,k)
!
    uvw(j)	=uvw(j)		+uvwm(i,j,k)
  end do
  end do
  end do
!
  ! Average ordre 3
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  uuu=uuu/n; uvv=uvv/n; uww=uww/n;
  vuu=vuu/n; vvv=vvv/n; vww=vww/n;
  wuu=wuu/n; wvv=wvv/n; www=www/n;
  uvw=uvw/n;
!
  ! Write data
  call write_myaverage_real_1d(uuu-u*uu-2*u*uu,xen2-xst2+1,n,xst2,xen2,'uuum1d')
  call write_myaverage_real_1d(uvv-u*vv-2*v*uv,xen2-xst2+1,n,xst2,xen2,'uvvm1d')
  call write_myaverage_real_1d(uww-u*ww-2*w*uw,xen2-xst2+1,n,xst2,xen2,'uwwm1d')
!
  call write_myaverage_real_1d(vuu-v*uu-2*u*uv,xen2-xst2+1,n,xst2,xen2,'vuum1d')
  call write_myaverage_real_1d(vvv-v*vv-2*v*vv,xen2-xst2+1,n,xst2,xen2,'vvvm1d')
  call write_myaverage_real_1d(vww-v*ww-2*w*vw,xen2-xst2+1,n,xst2,xen2,'vwwm1d')
!
  call write_myaverage_real_1d(wuu-w*uu-2*u*uw,xen2-xst2+1,n,xst2,xen2,'wuum1d')
  call write_myaverage_real_1d(wvv-w*vv-2*v*vw,xen2-xst2+1,n,xst2,xen2,'wvvm1d')
  call write_myaverage_real_1d(www-w*ww-2*w*ww,xen2-xst2+1,n,xst2,xen2,'wwwm1d')
!
  call write_myaverage_real_1d(uvw-u*vw-v*uw-w*uv,xen2-xst2+1,n,xst2,xen2,'uvwm1d')
!
  endif ! Ordre 3

  ! Scalaire
  if (iscalar.eq.1) then

  ! Ordre 1
  if (.true.) then
!
  ! Init ordre 1
  phi=0.;
  dphidy=0.;
  deltaphi=0.;
!
  ! Compute Average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    phi(j)	=phi(j)		+phim(i,j,k)
    dphidy(j)	=dphidy(j)	+dphidym(i,j,k)
    deltaphi(j)	=deltaphi(j)	+dphidyym(i,j,k)
  end do
  end do
  end do
!
  ! Average ordre 1
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  phi=phi/n;
  dphidy=dphidy/n;
  deltaphi=deltaphi/n;
!
  ! Write data
  call write_myaverage_real_1d(phi	,xen2-xst2+1,n,xst2,xen2,     'phim1d')
  call write_myaverage_real_1d(dphidy	,xen2-xst2+1,n,xst2,xen2,  'dphidym1d')
  call write_myaverage_real_1d(deltaphi	,xen2-xst2+1,n,xst2,xen2,'deltaphim1d')
!
  dphidx=0.; dphidz=0.;
!
  endif ! Ordre 1
!
  ! Ordre 2
  if (.true.) then
!
  ! Init ordre 2
  uphi=0.; vphi=0.; wphi=0.; phiphi=0.;
  dphidx2=0.; dphidy2=0.; dphidz2=0.;
  phidpdx=0.; phidpdy=0.; phidpdz=0.;
  dudphidx=0.; dvdphidx=0.; dwdphidx=0.;
  dudphidy=0.; dvdphidy=0.; dwdphidy=0.;
  dudphidz=0.; dvdphidz=0.; dwdphidz=0.;
  udeltaphi=0.; vdeltaphi=0.; wdeltaphi=0.;
!
  ! Compute Average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    uphi(j)		=uphi(j)	+uphim(i,j,k)
    vphi(j)		=vphi(j)	+vphim(i,j,k)
    wphi(j)		=wphi(j)	+wphim(i,j,k)
    phiphi(j)		=phiphi(j)	+phiphim(i,j,k)
    dphidx2(j)		=dphidx2(j)	+dphidx2m(i,j,k)
    dphidy2(j)		=dphidy2(j)	+dphidy2m(i,j,k)
    dphidz2(j)		=dphidz2(j)	+dphidz2m(i,j,k)
    phidpdx(j)		=phidpdx(j)	+phidpdxm(i,j,k)
    phidpdy(j)		=phidpdy(j)	+phidpdym(i,j,k)
    phidpdz(j)		=phidpdz(j)	+phidpdzm(i,j,k)
    dudphidx(j)		=dudphidx(j)	+dudphidxm(i,j,k)
    dvdphidx(j)		=dvdphidx(j)	+dvdphidxm(i,j,k)
    dwdphidx(j)		=dwdphidx(j)	+dwdphidxm(i,j,k)
    dudphidy(j)		=dudphidy(j)	+dudphidym(i,j,k)
    dvdphidy(j)		=dvdphidy(j)	+dvdphidym(i,j,k)
    dwdphidy(j)		=dwdphidy(j)	+dwdphidym(i,j,k)
    dudphidz(j)		=dudphidz(j)	+dudphidzm(i,j,k)
    dvdphidz(j)		=dvdphidz(j)	+dvdphidzm(i,j,k)
    dwdphidz(j)		=dwdphidz(j)	+dwdphidzm(i,j,k)
    udeltaphi(j)	=udeltaphi(j)&
 +udtdxxm(i,j,k)+udtdyym(i,j,k)+udtdzzm(i,j,k)
    vdeltaphi(j)	=vdeltaphi(j)&
 +vdtdxxm(i,j,k)+vdtdyym(i,j,k)+vdtdzzm(i,j,k)
    wdeltaphi(j)	=wdeltaphi(j)&
 +wdtdxxm(i,j,k)+wdtdyym(i,j,k)+wdtdzzm(i,j,k)
  end do
  end do
  end do
!
  ! Average ordre 2
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  uphi=uphi/n; vphi=vphi/n; wphi=wphi/n; phiphi=phiphi/n;
  dphidx2=dphidx2/n; dphidy2=dphidy2/n; dphidz2=dphidz2/n;
  phidpdx=phidpdx/n; phidpdy=phidpdy/n; phidpdz=phidpdz/n;
  dudphidx=dudphidx/n; dvdphidx=dvdphidx/n; dwdphidx=dwdphidx/n;
  dudphidy=dudphidy/n; dvdphidy=dvdphidy/n; dwdphidy=dwdphidy/n;
  dudphidz=dudphidz/n; dvdphidz=dvdphidz/n; dwdphidz=dwdphidz/n;
  udeltaphi=udeltaphi/n; vdeltaphi=vdeltaphi/n; wdeltaphi=wdeltaphi/n;
!
  ! Write data
  call write_myaverage_real_1d(uphi - u*phi		,xen2-xst2+1,n,xst2,xen2,      'uphim1d')
  call write_myaverage_real_1d(vphi - v*phi		,xen2-xst2+1,n,xst2,xen2,      'vphim1d')
  call write_myaverage_real_1d(wphi - w*phi		,xen2-xst2+1,n,xst2,xen2,      'wphim1d')
  call write_myaverage_real_1d(phiphi - phi*phi		,xen2-xst2+1,n,xst2,xen2,    'phiphim1d')
  call write_myaverage_real_1d(dphidx2 - dphidx*dphidx	,xen2-xst2+1,n,xst2,xen2,  'dphidx2m1d')
  call write_myaverage_real_1d(dphidy2 - dphidy*dphidy	,xen2-xst2+1,n,xst2,xen2,  'dphidy2m1d')
  call write_myaverage_real_1d(dphidz2 - dphidz*dphidz	,xen2-xst2+1,n,xst2,xen2,  'dphidz2m1d')
  call write_myaverage_real_1d(phidpdx - phi*dpdx	,xen2-xst2+1,n,xst2,xen2,  'phidpdxm1d')
  call write_myaverage_real_1d(phidpdy - phi*dpdy	,xen2-xst2+1,n,xst2,xen2,  'phidpdym1d')
  call write_myaverage_real_1d(phidpdz - phi*dpdz	,xen2-xst2+1,n,xst2,xen2,  'phidpdzm1d')
  call write_myaverage_real_1d(dudphidx - dudx*dphidx	,xen2-xst2+1,n,xst2,xen2, 'dudphidxm1d')
  call write_myaverage_real_1d(dvdphidx - dvdx*dphidx	,xen2-xst2+1,n,xst2,xen2, 'dvdphidxm1d')
  call write_myaverage_real_1d(dwdphidx - dwdx*dphidx	,xen2-xst2+1,n,xst2,xen2, 'dwdphidxm1d')
  call write_myaverage_real_1d(dudphidy - dudy*dphidy	,xen2-xst2+1,n,xst2,xen2, 'dudphidym1d')
  call write_myaverage_real_1d(dvdphidy - dvdy*dphidy	,xen2-xst2+1,n,xst2,xen2, 'dvdphidym1d')
  call write_myaverage_real_1d(dwdphidy - dwdy*dphidy	,xen2-xst2+1,n,xst2,xen2, 'dwdphidym1d')
  call write_myaverage_real_1d(dudphidz - dudz*dphidz	,xen2-xst2+1,n,xst2,xen2, 'dudphidzm1d')
  call write_myaverage_real_1d(dvdphidz - dvdz*dphidz	,xen2-xst2+1,n,xst2,xen2, 'dvdphidzm1d')
  call write_myaverage_real_1d(dwdphidz - dwdz*dphidz	,xen2-xst2+1,n,xst2,xen2, 'dwdphidzm1d')
  call write_myaverage_real_1d(udeltaphi - u*deltaphi	,xen2-xst2+1,n,xst2,xen2,'udeltaphim1d')
  call write_myaverage_real_1d(vdeltaphi - v*deltaphi	,xen2-xst2+1,n,xst2,xen2,'vdeltaphim1d')
  call write_myaverage_real_1d(wdeltaphi - w*deltaphi	,xen2-xst2+1,n,xst2,xen2,'wdeltaphim1d')
!
  endif ! Ordre 2

  ! Ordre 3
  if (.true.) then
!
  ! Init ordre 3
  uphi2=0.; vphi2=0.; wphi2=0.; 
  phiuu=0.; phivv=0.; phiww=0.; 
  phiuv=0.; phiuw=0.; phivw=0.; 
!
  ! Compute Average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    uphi2(j) = uphi2(j) + uphi2m(i,j,k)
    vphi2(j) = vphi2(j) + vphi2m(i,j,k)
    wphi2(j) = wphi2(j) + wphi2m(i,j,k)
    phiuu(j) = phiuu(j) + phiuum(i,j,k)
    phivv(j) = phivv(j) + phivvm(i,j,k)
    phiww(j) = phiww(j) + phiwwm(i,j,k)
    phiuv(j) = phiuv(j) + phiuvm(i,j,k)
    phiuw(j) = phiuw(j) + phiuwm(i,j,k)
    phivw(j) = phivw(j) + phivwm(i,j,k)
  end do
  end do
  end do
!
  ! Average ordre 3
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  uphi2=uphi2/n; vphi2=vphi2/n; wphi2=wphi2/n; 
  phiuu=phiuu/n; phivv=phivv/n; phiww=phiww/n; 
  phiuv=phiuv/n; phiuw=phiuw/n; phivw=phivw/n; 
!
  ! Write data
  call write_myaverage_real_1d(uphi2 - u*phiphi -2*phi*uphi  ,xen2-xst2+1,n,xst2,xen2,'uphi2m1d')
  call write_myaverage_real_1d(vphi2 - v*phiphi -2*phi*vphi  ,xen2-xst2+1,n,xst2,xen2,'vphi2m1d')
  call write_myaverage_real_1d(wphi2 - w*phiphi -2*phi*wphi  ,xen2-xst2+1,n,xst2,xen2,'wphi2m1d')
  call write_myaverage_real_1d(phiuu - phi*uu -2*u*uphi      ,xen2-xst2+1,n,xst2,xen2,'phiuum1d')
  call write_myaverage_real_1d(phivv - phi*vv -2*v*vphi      ,xen2-xst2+1,n,xst2,xen2,'phivvm1d')
  call write_myaverage_real_1d(phiww - phi*ww -2*w*wphi      ,xen2-xst2+1,n,xst2,xen2,'phiwwm1d')
  call write_myaverage_real_1d(phiuv - phi*uv -u*vphi -v*uphi,xen2-xst2+1,n,xst2,xen2,'phiuvm1d')
  call write_myaverage_real_1d(phiuw - phi*uw -u*wphi -w*uphi,xen2-xst2+1,n,xst2,xen2,'phiuwm1d')
  call write_myaverage_real_1d(phivw - phi*vw -v*wphi -w*vphi,xen2-xst2+1,n,xst2,xen2,'phivwm1d')
!
  endif ! Ordre 3

  endif ! Scalaire

  ! Pression
  if (.true.) then !(.true.) then
!
  ! Init
  pres=0.; pres2=0.;
!
  ! Compute average
  do k=xst3,xen3
  do j=xst2,xen2
  do i=xst1,xen1
    pres(j)  = pres(j)  + presm(i,j,k)
    pres2(j) = pres2(j) + pres2m(i,j,k)
  end do
  end do
  end do
!
  ! Average
  n= (xen3-xst3+1)*(xen1-xst1+1)
  !
  pres=pres/n; pres2=pres2/n;
!
  ! Write data
  call write_myaverage_real_1d(pres             ,xen2-xst2+1,n,xst2,xen2, 'presm1d')
  call write_myaverage_real_1d(pres2 - pres*pres,xen2-xst2+1,n,xst2,xen2,'pres2m1d')
!
  endif

! Pression sur maillage décalé
!
!if (.true.) then
!
  ! Init
!  presi=0.; pres2i=0.;
!
  ! From z-pencil to x-pencil
!  call transpose_z_to_y(pp3m, pp3my,ph3)
!  call transpose_y_to_x(pp3my,pp3mx,ph3)
!  call transpose_z_to_y(pp32m, pp3my,ph3)
!  call transpose_y_to_x(pp3my,pp32mx,ph3)
!
  ! Compute average
!  do k=ph3%xst(3),ph3%xen(3)
!  do j=ph3%xst(2),ph3%xen(2)
!  do i=ph3%xst(1),ph3%xen(1)
!    presi(j)  = presi(j)  + pp3mx(i,j,k)
!    pres2i(j) = pres2i(j) + pp32mx(i,j,k)
!  end do
!  end do
!  end do
!
  ! Average
!  n=(ph3%xen(3)-ph3%xst(3)+1)*(ph3%xen(1)-ph3%xst(1)+1)
  !
!  presi=presi/n; pres2i=pres2i/n;
!
  ! Write data
!  call write_myaverage_real_1di(presi               ,ph3%xen(2)-ph3%xst(2)+1,n,ph3%xst(2),ph3%xen(2), 'presim1d')
!  call write_myaverage_real_1di(pres2i - presi*presi,ph3%xen(2)-ph3%xst(2)+1,n,ph3%xst(2),ph3%xen(2),'pres2im1d')
!
!endif
!

end subroutine compute_balance_1d

subroutine write_myaverage_real_1d(my_var,var_size,average_size,jst,jen,filename)

use decomp_2d, only : mytype, nrank, real_type
use variables, only : yp
use MPI

implicit none

integer, intent(in) :: var_size,average_size,jst,jen
real(mytype), intent(in) :: my_var(jst:jen)
character(len=*), intent(IN) :: filename

integer, parameter :: num_unit_average=1111
integer :: total_size, jmax, jmin
integer :: ierror, i
integer, allocatable, dimension(:)      :: compteur, tmp_int
real(mytype), allocatable, dimension(:) :: glob_var, tmp_dbe

! Compute total size
call MPI_AllReduce(jst,jmin,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierror)
call MPI_AllReduce(jen,jmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierror) 

! Init
total_size = jmax - jmin + 1
allocate(glob_var(jmin:jmax),tmp_dbe(jmin:jmax))
glob_var = 0.
tmp_dbe = 0.
allocate(compteur(jmin:jmax),tmp_int(jmin:jmax))
compteur = 0
tmp_int = 0

! Affect
glob_var(jst:jen) = average_size * my_var(jst:jen)
compteur(jst:jen) = average_size

! Regroup data on rank 0
call MPI_Reduce(compteur,tmp_int,total_size,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
call MPI_Reduce(glob_var,tmp_dbe,total_size,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
if (nrank.eq.0) then
  compteur = tmp_int
  glob_var = tmp_dbe
endif
deallocate(tmp_int,tmp_dbe)

! Rank 0 compute the average then write data in text format
if (nrank.eq.0) then
  ! Compute average
  if (any(compteur.eq.0)) print *,'ERROR : Division by zero when averaging data, value replaced by 0'
  where (compteur.ne.0)
    glob_var = glob_var / (1.*compteur)
  elsewhere
    glob_var = 0.
  end where
  ! Write average
  open(num_unit_average,file=filename//'.dat',position='APPEND')
  do i=jmin,jmax
    write(num_unit_average,*)yp(i), glob_var(i)
  enddo
  close(num_unit_average)
endif

deallocate(glob_var,compteur)

end subroutine write_myaverage_real_1d
