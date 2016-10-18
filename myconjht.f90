module conjugate_ht
!
! Module for conjugate heat transfer
!
use decomp_2d, only : mytype, DECOMP_INFO

type(decomp_info), save :: mydecomp_bot, mydecomp_top
logical, save :: bool_conjugate_ht, bool_sol_stats
integer, save :: beg_stat_sol
integer, save :: ny_sol_bot, ny_sol_top
real(mytype), parameter :: sol_imp_var=0. ! 1=explicit, 0=implicit, 0.5=Crank-Nicolson. 0 recommended
real(mytype), save :: ly_sol_bot, ly_sol_top
real(mytype), save :: repr_sol_bot, repr_sol_top
real(mytype), save :: fluxratio_bot, fluxratio_top
real(mytype), save, allocatable, dimension(:) :: yp_bot, yp_top
real(mytype), save, allocatable, dimension(:,:) :: val_bot, val_top
real(mytype), save, allocatable, dimension(:,:) :: flux_bot, flux_top
real(mytype), save, allocatable, dimension(:,:) :: d2pn_bot, d2pn_top
real(mytype), save, allocatable, dimension(:,:) :: invdiffy_bot, invdiffy_top
real(mytype), save, allocatable, dimension(:,:,:) :: cl_bot, cl_top
real(mytype), save, allocatable, dimension(:,:,:) :: temp_bot, temp_top
real(mytype), save, allocatable, dimension(:,:,:) :: bot_dm1, top_dm1
real(mytype), save, allocatable, dimension(:,:,:) :: bot_dm2, top_dm2
real(mytype), save, allocatable, dimension(:,:,:) :: bot_dm3, top_dm3
!
! stats
!
real(mytype), save, allocatable, dimension(:,:,:) :: tempm_1_bot, tempm_2_bot, tdeltat_bot
real(mytype), save, allocatable, dimension(:,:,:) :: tempm_1_top, tempm_2_top, tdeltat_top
!
! second derivative in x & z
!
real(mytype), save :: fpi2ts
real(mytype), save :: alsaixts, asixts, bsixts, csixts
real(mytype), save :: alsakzts, askzts, bskzts, cskzts
real(mytype), save, allocatable, dimension(:) :: sfxts, scxts, sbxts, ssxts, swxts
real(mytype), save, allocatable, dimension(:) :: sfxpts, ssxpts, swxpts
real(mytype), save, allocatable, dimension(:) :: sfzts, sczts, sbzts, sszts, swzts
real(mytype), save, allocatable, dimension(:) :: sfzpts, sszpts, swzpts
!
! first derivative in y
!
real(mytype), save, allocatable, dimension(:,:) :: mat_dery

contains

subroutine conjugate_ht_init()
!
! Subroutine to initialize the module
!
  use decomp_2d, only : ysize, nrank, ystart
  use param, only : yly, ilit, twopi, dz, dx, dt

  implicit none

  ! Local variables
  integer :: i,j,k
  real(mytype) :: ytmp
  real(mytype) :: x, z, mysx, mysy, mysz
  real(mytype), dimension(ysize(1),ny_sol_bot+1+2,ysize(3)) :: ydiff_bot
  real(mytype), dimension(ysize(1),ny_sol_top+1+2,ysize(3)) :: ydiff_top

  if (nrank.eq.0) print *,'  User module for conjugate heat transfer :'
  if (nrank.eq.0) print *,'    - ',ny_sol_bot,' & ',ny_sol_top,' Chebyshev nodes top & bottom'

  ! ny_bot et ny_top : ordre de l'interpolation
  ! ny_bot+1 et ny_top+1 : nombres de points de chebyshev
  ! +2 points pour les bords
  !
  allocate(yp_bot(ny_sol_bot+1+2), yp_top(ny_sol_top+1+2))
  call cheb_nodes(yp_bot,ny_sol_bot,-ly_sol_bot,real(0,kind=mytype))
  call cheb_nodes(yp_top,ny_sol_top,yly,yly+ly_sol_top)
  if (nrank.eq.0) print *,'    - Chebyshev nodes computed'
  if (nrank.eq.0) then
   open(21,file='yp_bot.dat',form='formatted')
   do j=1,ny_sol_bot+1+2
    write(21,*) yp_bot(j)
   enddo
   close(21)
   open(21,file='yp_top.dat',form='formatted')
   do j=1,ny_sol_top+1+2
    write(21,*) yp_top(j)
   enddo
   close(21)
  endif
  !
  ! Compute coefficients for value at boundary
  !
  allocate(val_bot(2,ny_sol_bot+1),val_top(2,ny_sol_top+1))
  call cheb_edge(val_bot,yp_bot,ny_sol_bot)
  call cheb_edge(val_top,yp_top,ny_sol_top)
  if (nrank.eq.0) print *,'    - Boundary-value matrix computed'
  !
  ! Compute coefficients for flux at boundary
  !
  allocate(flux_bot(2,ny_sol_bot+1),flux_top(2,ny_sol_top+1))
  call cheb_flux(flux_bot,yp_bot,ny_sol_bot)
  call cheb_flux(flux_top,yp_top,ny_sol_top)
  if (nrank.eq.0) print *,'    - Boundary-flux matrix computed'
  !
  ! Compute coefficients of second derivative
  ! on Chebyshev nodes (y-direction)
  !
  allocate(d2pn_bot(ny_sol_bot+1,ny_sol_bot+1), d2pn_top(ny_sol_top+1,ny_sol_top+1))
  call cheb_der2(d2pn_bot,yp_bot,ny_sol_bot)
  call cheb_der2(d2pn_top,yp_top,ny_sol_top)
  if (nrank.eq.0) print *,'    - Y-diffusion matrix computed'
  if (nrank.eq.0) print *,'    - Computing inverse of matrix'
  allocate(invdiffy_bot(ny_sol_bot-1,ny_sol_bot-1))
  i = my_ydiff_mat_inv(d2pn_bot,ny_sol_bot,invdiffy_bot,flux_bot(1,:),flux_bot(2,:),repr_sol_bot)
  if (nrank.eq.0) print *,'      o Bottom matrix inverted : ', i
  allocate(invdiffy_top(ny_sol_top-1,ny_sol_top-1))
  i = my_ydiff_mat_inv(d2pn_top,ny_sol_top,invdiffy_top,flux_top(1,:),flux_top(2,:),repr_sol_top)
  if (nrank.eq.0) print *,'      o Top matrix inverted : ', i
  !
  allocate(cl_bot(ysize(1),1,ysize(3)))
  allocate(cl_top(ysize(1),1,ysize(3)))
  cl_bot=0.
  cl_top=0.
  !
  allocate(temp_bot(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(temp_top(ysize(1),ny_sol_top+1+2,ysize(3)))
  temp_bot=0.
  temp_top=0.
  allocate(bot_dm1(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(top_dm1(ysize(1),ny_sol_top+1+2,ysize(3)))
  bot_dm1=0.
  top_dm1=0.
  allocate(bot_dm2(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(top_dm2(ysize(1),ny_sol_top+1+2,ysize(3)))
  bot_dm2=0.
  top_dm2=0.
  allocate(bot_dm3(ysize(1),ny_sol_bot+1+2,ysize(3)))
  allocate(top_dm3(ysize(1),ny_sol_top+1+2,ysize(3)))
  bot_dm3=0.
  top_dm3=0.
  !
  ! 2DECOMP for bottom and top solid domain
  !
  call my2decomp_solide()
  !
  ! Initial value for temperature
  !
  if (ilit.eq.0) then
    !
    do k=1,ysize(3)
      z=(k-1+ystart(3)-1)*dz
      mysz=sin(twopi*z)
    do j=1,ny_sol_bot+1+2
      ytmp=yp_bot(j)
      mysy=sin(twopi*ytmp)
    do i=1,ysize(1)
      x=(i-1+ystart(1)-1)*dx
      mysx=sin(twopi*x)
      temp_bot(i,j,k)=ytmp*fluxratio_bot
    enddo
    enddo
    enddo
    !
    do k=1,ysize(3)
      z=(k-1+ystart(3)-1)*dz
      mysz=sin(twopi*z)
    do j=1,ny_sol_top+1+2
      ytmp=yp_top(j)
      mysy=sin(twopi*ytmp)
    do i=1,ysize(1)
      x=(i-1+ystart(1)-1)*dx
      mysx=sin(twopi*x)
      temp_top(i,j,k)=-(ytmp-2.)*fluxratio_top
    enddo
    enddo
    enddo
    !
    ! Compute current xz diffusion
    call xzdiff_temp_solide(temp_bot,bot_dm1,ny_sol_bot,mydecomp_bot)
    call xzdiff_temp_solide(temp_top,top_dm1,ny_sol_top,mydecomp_top)
    ! Compute current y diffusion
    call ydiff_temp_solide(ydiff_bot,temp_bot,d2pn_bot,ny_sol_bot)
    call ydiff_temp_solide(ydiff_top,temp_top,d2pn_top,ny_sol_top)
    ! Diffusion non explicite
    bot_dm1=(bot_dm1+sol_imp_var*ydiff_bot)/repr_sol_bot
    top_dm1=(top_dm1+sol_imp_var*ydiff_top)/repr_sol_top
    if (nrank.eq.0) print *,'    - Solid temperature initialized'
    ! Compute min/max & check convergence
    if (nrank.eq.0) print *,'Solide bas : '
    call test_sol_min_max(temp_bot,yp_bot,ny_sol_bot,fluxratio_bot,repr_sol_bot)
    if (nrank.eq.0) print *,'Solide haut : '
    call test_sol_min_max(temp_top,yp_top,ny_sol_top,fluxratio_top,repr_sol_top)
    !
  else
    !
    call solide_restart(.false.)
    !
    if (nrank.eq.0) print *,'    - Solid temperature restart ok'
    !
  endif
  !
  if (nrank.eq.0) print *,'  Conjugate heat transfer initialized'
  !
end subroutine conjugate_ht_init

subroutine cheb_nodes(y,n,a,b)
!
! Store chebyshev nodes in y
!
  use param, only : pi

  implicit none

  integer, intent(in) :: n
  real(mytype), intent(in) :: a, b
  real(mytype), dimension(n+1+2), intent(out) :: y

  ! Local variables
  integer :: i

  y(1)=a
  do i=2,n+2
    y(i)=(a+b)/2.d0+(b-a)*cos(((2*(n-(i-2))+1)*1.d0/(2*n+2))*pi)/2.d0
  enddo
  y(n+1+2)=b

end subroutine cheb_nodes

subroutine cheb_der2(der2,y,n)
!
! Compute second derivative of (Chebyshev) polynoms
! At y(i) nodes
! Line for a node
! Col. for a polynom
!
  use decomp_2d, only : nproc, nrank, real_type
  use MPI

  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(n+1,n+1), intent(out) :: der2

  ! Local variables
  integer :: i, j, k, l, code
  real(mytype), dimension(n+1) :: d2li, mysum, myprod

  ! Loop on polynoms
  if (nproc.ge.n+1) then
    !
    ! Procs 0 to n compute a column
    !
    if (nrank.le.n) then
      i=nrank+1
      d2li=0.
      do l=1,n+1
      if (l.ne.i) then
        mysum=0.
        do k=1,n+1
        if ((k.ne.i).and.(k.ne.l)) then
          myprod=1.
          do j=1,n+1
          if ((j.ne.i).and.(j.ne.k).and.(j.ne.l)) then
            myprod=myprod*(y(2:n+2)-y(j+1))/(y(i+1)-y(j+1))
          endif
          enddo
          mysum=mysum+myprod/(y(i+1)-y(k+1))
        endif
        enddo
        d2li=d2li+mysum/(y(i+1)-y(l+1))
      endif
      enddo
      der2(:,i)=d2li
    endif
    !
    ! then information is broadcasted to all procs
    !
    do i=1,n+1
      call MPI_BCAST(der2(:,i),n+1,real_type,i-1,MPI_COMM_WORLD,code)
    enddo
    !
  else
    !
    do i=1,n+1
      d2li=0.
      do l=1,n+1
      if (l.ne.i) then
        mysum=0.
        do k=1,n+1
        if ((k.ne.i).and.(k.ne.l)) then
          myprod=1.
          do j=1,n+1
          if ((j.ne.i).and.(j.ne.k).and.(j.ne.l)) then
            myprod=myprod*(y(2:n+2)-y(j+1))/(y(i+1)-y(j+1))
          endif
          enddo
          mysum=mysum+myprod/(y(i+1)-y(k+1))
        endif
        enddo
        d2li=d2li+mysum/(y(i+1)-y(l+1))
      endif
      enddo
      der2(:,i)=d2li
    enddo
    !
  endif

end subroutine cheb_der2

subroutine cheb_flux(flux,y,n)
!
! Compute First derivative of chebyshev polynoms
! At boundary nodes
! First line for derivative at first node
! Second line for derivative at last node
!
  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(2,n+1), intent(out) :: flux

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(2) :: d1li, myprod, tmpy

  tmpy(1)=y(1)
  tmpy(2)=y(n+1+2)
  
  do i=1,n+1
    d1li=0.
    do k=1,n+1
    if (k.ne.i) then
      myprod=1.
      do j=1,n+1
      if ((j.ne.i).and.(j.ne.k)) then
        myprod=myprod*(tmpy-y(j+1))/(y(i+1)-y(j+1))
      endif
      enddo
      d1li=d1li+myprod/(y(i+1)-y(k+1))
    endif
    enddo
    flux(:,i)=d1li
  enddo

end subroutine cheb_flux

subroutine cheb_edge(val,y,n)
!
! Compute polynom value at boundary nodes
! First line for value at first node
! Second line for value at last node
!
  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(2,n+1), intent(out) :: val

  ! Local variables
  integer :: i, j
  real(mytype), dimension(2) :: li, tmpy

  tmpy(1)=y(1)
  tmpy(2)=y(n+1+2)

  do i=1,n+1
    li=1.
    do j=1,n+1
    if (j.ne.i) then
      li=li*(tmpy-y(j+1))/(y(i+1)-y(j+1))
    endif
    enddo
    val(:,i)=li
  enddo

end subroutine cheb_edge

subroutine update_temp_solide()
!
! Compute temperature in the solid
! with the boundary flux imposed from the fluid
!
  use decomp_2d, only : ysize, nrank
  use param, only : dt, twopi, itime, ilit, t

  implicit none

  ! Local variables
  integer :: i, j, k
  real(mytype) :: a, b, c
  real(mytype), dimension(ysize(1),ny_sol_bot+1+2,ysize(3)) :: ydiff_bot
  real(mytype), dimension(ysize(1),ny_sol_top+1+2,ysize(3)) :: ydiff_top

  ! Switch from Euler1 to AB3
  if ((ilit.eq.0).and.(itime.eq.1)) then
    a=   dt
    b=   0.
    c=   0.
  elseif ((ilit.eq.0).and.(itime.eq.2)) then
    a=   2.*dt
    b= - dt
    c=   0.
  else
    a=   3.*dt
    b= - 3.*dt
    c=   dt
  endif

  ! Full RHS in temporary array ydiff_***
  ydiff_bot = temp_bot &
            + a*bot_dm1 &
            + b*bot_dm2 &
            + c*bot_dm3
  ydiff_top = temp_top &
            + a*top_dm1 &
            + b*top_dm2 &
            + c*top_dm3

  ! Boundary coundition : flux imposed with fluid and temperature imposed with exterior
  ! Derivate in solid is derivate in fluid * fluxratio
  ydiff_bot(:,1,:)=fluxratio_bot
  ydiff_bot(:,ny_sol_bot+1+2,:)=cl_bot(:,1,:)*fluxratio_bot
  ydiff_top(:,1,:)=cl_top(:,1,:)*fluxratio_top
  ydiff_top(:,ny_sol_top+1+2,:)=-fluxratio_top

  ! Use matrix inversion to update temperature
  ! New temperature in temp_bot and temp_top
  call my_new_solid_temp(temp_bot,ydiff_bot,invdiffy_bot,&
           d2pn_bot,val_bot,ny_sol_bot,flux_bot(1,:),flux_bot(2,:),repr_sol_bot)
  call my_new_solid_temp(temp_top,ydiff_top,invdiffy_top,&
           d2pn_top,val_top,ny_sol_top,flux_top(1,:),flux_top(2,:),repr_sol_top)

  ! Compute min/max & check convergence if required
  if (nrank.eq.0) print *,'Solide bas : '
  call test_sol_min_max(temp_bot,yp_bot,ny_sol_bot,fluxratio_bot,repr_sol_bot)
  if (nrank.eq.0) print *,'Solide haut : '
  call test_sol_min_max(temp_top,yp_top,ny_sol_top,fluxratio_top,repr_sol_top)

  ! Update history
  bot_dm3=bot_dm2
  top_dm3=top_dm2
  bot_dm2=bot_dm1
  top_dm2=top_dm1

  ! Compute current xz diffusion
  ! (bot_dm1 and top_dm1 are updated)
  call xzdiff_temp_solide(temp_bot,bot_dm1,ny_sol_bot,mydecomp_bot)
  call xzdiff_temp_solide(temp_top,top_dm1,ny_sol_top,mydecomp_top)
  ! Compute current y diffusion
  ! (stored in temporary array ydiff_***)
  call ydiff_temp_solide(ydiff_bot,temp_bot,d2pn_bot,ny_sol_bot)
  call ydiff_temp_solide(ydiff_top,temp_top,d2pn_top,ny_sol_top)
  ! Diffusion non explicite
  bot_dm1=(bot_dm1+sol_imp_var*ydiff_bot)/repr_sol_bot
  top_dm1=(top_dm1+sol_imp_var*ydiff_top)/repr_sol_top

end subroutine update_temp_solide

subroutine ydiff_temp_solide(diff,temp,d2pn,n)
!
! Compute diffusion in y-direction
! Only at chebyshev nodes
!
  use decomp_2d, only : ysize

  implicit none

  integer, intent(in) :: n
  real(mytype), dimension(n+1,n+1), intent(in) :: d2pn
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: temp
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(out) :: diff

  ! Local variables
  integer :: j,l

  diff=0.
  do j=2,n+2
    do l=1,n+1
      diff(:,j,:)=diff(:,j,:)+d2pn(j-1,l)*temp(:,l+1,:)
    enddo
  enddo

end subroutine ydiff_temp_solide

integer function my_ydiff_mat_inv(matin,n,matout,cl_b,cl_t,repr)
!
! Compute the inverse of y diffusion
! With cl_* boundary condition
! Put it in matout
!
  use param, only : dt

  implicit none

  integer, intent(in) :: n
  real(mytype), intent(in) :: repr
  real(mytype), dimension(1,n+1), intent(in) :: cl_b, cl_t
  real(mytype), dimension(n+1,n+1), intent(in) :: matin
  real(mytype), dimension(n-1,n-1), intent(out) :: matout

  ! Local variables
  integer :: info
  integer :: i,j,k,l
  integer, dimension(n-1) :: ipiv
  real(mytype), dimension(2,2) :: math
  real(mytype), dimension(2,n-1) :: matg
  real(mytype), dimension(2,n+1) :: mycl

  real(mytype), dimension(n+1,n+1) :: mattmp
  real(mytype), dimension(:), allocatable :: work

  ! Lapack routines to inverse matrix
  external DGETRF ! LU DECOMP, double prec
  external SGETRF ! LU DECOMP, simple prec
  external DGETRI ! Matrix inversion, double prec
  external SGETRI ! Matrix inversion, simple prec

  my_ydiff_mat_inv=0
  matout=0.
  mattmp=0.

  mycl(1,:)=cl_b(1,:)
  mycl(2,:)=cl_t(1,:)

  ! Compute math (inverse of a 2x2 matrix)
  math(1,1)=mycl(2,n+1)
  math(2,2)=mycl(1,1)
  math(1,2)=-mycl(1,n+1)
  math(2,1)=-mycl(2,1)
  math=math/(mycl(1,1)*mycl(2,n+1)-mycl(1,n+1)*mycl(2,1))

  ! Compute matg
  matg=0.
  do i=1,2
  do j=1,n-1
  do l=1,2
    matg(i,j) = matg(i,j) - math(i,l)*mycl(l,j+1)
  enddo
  enddo
  enddo

  ! Compute matout (matrix to inverse)
  mattmp=0.
  do i=1,n+1
    mattmp(i,i)=1.
  enddo
  ! Diffusion is not explicit
  do i=1,n+1
  do j=1,n+1
    mattmp(i,j)=mattmp(i,j)-(1.-sol_imp_var)*dt*matin(i,j)/repr
  enddo
  enddo
  do i=1,n-1
  do j=1,n-1
    matout(i,j) = mattmp(i+1,j+1) &
                + mattmp(i+1,1)*matg(1,j) &
                + mattmp(i+1,n+1)*matg(2,j)
  enddo
  enddo

  ! Call lapack to inverse the matrix
  if (allocated(work)) deallocate(work)
  allocate(work(1))
#ifdef DOUBLE_PREC
  call DGETRF( n-1, n-1, matout, n-1, ipiv, info )
  CALL DGETRI( n-1, matout, n-1, ipiv, work, -1, info)
#else
  call SGETRF( n-1, n-1, matout, n-1, ipiv, info )
  CALL SGETRI( n-1, matout, n-1, ipiv, work, -1, info)
#endif
  i=work(1)
  info=1
  do while (info.ne.0)
    j=i
    if (allocated(work)) deallocate(work)
    allocate(work(j),stat=info)
    i=j/2
  end do
#ifdef DOUBLE_PREC
  CALL DGETRI( n-1, matout, n-1, ipiv, work, j, info)
#else
  CALL SGETRI( n-1, matout, n-1, ipiv, work, j, info)
#endif
  if (allocated(work)) deallocate(work)

  my_ydiff_mat_inv = info

  ! Terminus
  return

end function my_ydiff_mat_inv

subroutine my_new_solid_temp(temp,rhs,invdiffy,d2pn,val,n,cl_b,cl_t,repr)
!
! Update temperature temp :
! temp = invdiffy * rhs
!
  use decomp_2d, only : ysize
  use param, only : dt

  implicit none

  integer, intent(in) :: n
  real(mytype), intent(in) :: repr
  real(mytype), dimension(2,n+1), intent(in) :: val
  real(mytype), dimension(1,n+1), intent(in) :: cl_b, cl_t
  real(mytype), dimension(n+1,n+1), intent(in) :: d2pn
  real(mytype), dimension(n-1,n-1), intent(in) :: invdiffy
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: rhs
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(out) :: temp

  ! Local variables
  integer :: i,j,k,l
  real(mytype), dimension(2,2) :: math
  real(mytype), dimension(2,n-1) :: matg
  real(mytype), dimension(2,n+1) :: mycl

  temp=0.

  mycl(1,:)=cl_b(1,:)
  mycl(2,:)=cl_t(1,:)

  ! Compute math (inverse of a 2x2 matrix)
  math(1,1)=mycl(2,n+1)
  math(2,2)=mycl(1,1)
  math(1,2)=-mycl(1,n+1)
  math(2,1)=-mycl(2,1)
  math=math/(mycl(1,1)*mycl(2,n+1)-mycl(1,n+1)*mycl(2,1))
  ! Compute matg
  matg=0.
  do i=1,2
  do j=1,n-1
  do l=1,2
    matg(i,j) = matg(i,j) - math(i,l)*mycl(l,j+1)
  enddo
  enddo
  enddo

  ! Compute temp inside (j=3:n+1)
  do k=1,ysize(3)
  do j=3,n+1
  do i=1,ysize(1)
    do l=1,n-1
      temp(i,j,k) = temp(i,j,k) + invdiffy(j-2,l)*( rhs(i,l+2,k) &
            + ((1.-sol_imp_var)/repr)*dt*d2pn(l+1,  1)* &
                  ( math(1,1)*rhs(i,1,k)+math(1,2)*rhs(i,n+3,k) ) &
            + ((1.-sol_imp_var)/repr)*dt*d2pn(l+1,n+1)* &
                  ( math(2,1)*rhs(i,1,k)+math(2,2)*rhs(i,n+3,k) ) )
    enddo
  enddo
  enddo
  enddo

  ! Compute temp at boundary chebyshev nodes
  ! Using boundary condition
  ! j=2 & j=n+2
  !
  ! Non-homogeneous boundary-condition from fluid
  do k=1,ysize(3)
  do i=1,ysize(1)
    temp(i,2,k)   = math(1,1)*rhs(i,1,k) + math(1,2)*rhs(i,n+3,k)
    temp(i,n+2,k) = math(2,1)*rhs(i,1,k) + math(2,2)*rhs(i,n+3,k)
  enddo
  enddo
  !
  ! Homogeneous boundary condition from solid
  do k=1,ysize(3)
  do i=1,ysize(1)
    do j=1,n-1
      temp(i,2,k)   = temp(i,2,k) &
                    + matg(1,j)*temp(i,j+2,k)
      temp(i,n+2,k) = temp(i,n+2,k) & 
                    + matg(2,j)*temp(i,j+2,k)
    enddo
  enddo
  enddo

  ! Compute temp at boundary nodes
  ! Using chebyshev coefficients
  do k=1,ysize(3)
  do j=2,n+2
  do i=1,ysize(1)
    temp(i,1,k)   = temp(i,1,k) &
                  + temp(i,j,k)*val(1,j-1)
    temp(i,n+3,k) = temp(i,n+3,k) &
                  + temp(i,j,k)*val(2,j-1)
  enddo
  enddo
  enddo

end subroutine my_new_solid_temp

subroutine xzdiff_temp_solide(temp,dm1,n,mydeco)
!
! Compute diffusion of temp
! In x and z direction
! and put it in dm1
!
  use decomp_2d, only : ysize, nproc, nrank, ystart, yend, real_type, &
                        alloc_x, alloc_y, alloc_z, &
                        transpose_x_to_y, transpose_y_to_x, &
                        transpose_z_to_y, transpose_y_to_z
  use param, only : dx2, dz2, itime, ifirst, pi
  use variables, only : nx, nz

  implicit none

  integer, intent(in) :: n
  type(decomp_info), intent(in) :: mydeco
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: temp
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(out) :: dm1

  ! Local variables
  integer :: i,j,k
  real(mytype), dimension(mydeco%xsz(1),mydeco%xsz(2),mydeco%xsz(3)) :: dtdxx
  real(mytype), dimension(mydeco%zsz(1),mydeco%zsz(2),mydeco%zsz(3)) :: dtdzz
  real(mytype), dimension(mydeco%ysz(1),mydeco%ysz(2),mydeco%ysz(3)) :: dtdxyz
  real(mytype), dimension(mydeco%xsz(1),mydeco%xsz(2),mydeco%xsz(3)) :: temp1
  real(mytype), dimension(mydeco%zsz(1),mydeco%zsz(2),mydeco%zsz(3)) :: temp3
  real(mytype), dimension(mydeco%xsz(1),mydeco%xsz(2),mydeco%xsz(3)) :: di1
  real(mytype), dimension(mydeco%zsz(1),mydeco%zsz(2),mydeco%zsz(3)) :: di3
  real(mytype), dimension(mydeco%xsz(2),mydeco%xsz(3)) :: sx
  real(mytype), dimension(mydeco%zsz(1),mydeco%zsz(2)) :: sz

  ! Init diffusion (dm1) and boundary values
  dm1=0.

  if (.not.allocated(sfxts)) then ! init x & z second derivative

    ! allocate
    allocate(sfxts(nx), scxts(nx), sbxts(nx), ssxts(nx), swxts(nx))
    allocate(sfxpts(nx), ssxpts(nx), swxpts(nx))
    allocate(sfzts(nz), sczts(nz), sbzts(nz), sszts(nz), swzts(nz))
    allocate(sfzpts(nz), sszpts(nz), swzpts(nz))
    ! init
    sfxts=0.; scxts=0.; sbxts=0.; ssxts=0.; swxts=0.
    sfxpts=0.; ssxpts=0.; swxpts=0.
    sfzts=0.; sczts=0.; sbzts=0.; sszts=0.; swzts=0.
    sfzpts=0.; sszpts=0.; swzpts=0.

    ! wavenumber at cut-off
    fpi2ts=4.
    ! compact finite difference coefficients in x
    alsaixts=(45.*fpi2ts*pi*pi-272.)/(2.*(45.*fpi2ts*pi*pi-208.))
    asixts  =((6.-9.*alsaixts)/4.)/dx2
    bsixts  =((-3.+24.*alsaixts)/5.)/(4.*dx2)
    csixts  =((2.-11.*alsaixts)/20.)/(9.*dx2)
    ! compact finite difference coefficients in z
    alsakzts=(45.*fpi2ts*pi*pi-272.)/(2.*(45.*fpi2ts*pi*pi-208.))
    askzts  =((6.-9.*alsakzts)/4.)/dz2
    bskzts  =((-3.+24.*alsakzts)/5.)/(4.*dz2)
    cskzts  =((2.-11.*alsakzts)/20.)/(9.*dz2)

    ! nclx=0 only implemented
    sfxts(1)   =alsaixts
    sfxts(2)   =alsaixts
    sfxts(nx-2)=alsaixts
    sfxts(nx-1)=alsaixts
    sfxts(nx)  =0.
    scxts(1)   =2.
    scxts(2)   =1.
    scxts(nx-2)=1.
    scxts(nx-1)=1.
    scxts(nx  )=1.+alsaixts*alsaixts
    sbxts(1)   =alsaixts
    sbxts(2)   =alsaixts
    sbxts(nx-2)=alsaixts
    sbxts(nx-1)=alsaixts
    sbxts(nx  )=0.
    do i=3,nx-3
       sfxts(i)=alsaixts
       scxts(i)=1.
       sbxts(i)=alsaixts
    enddo

    ! nclz=0 only implemented
    sfzts(1)   =alsakzts
    sfzts(2)   =alsakzts
    sfzts(nz-2)=alsakzts
    sfzts(nz-1)=alsakzts
    sfzts(nz)  =0.
    sczts(1)   =2.
    sczts(2)   =1.
    sczts(nz-2)=1.
    sczts(nz-1)=1.
    sczts(nz  )=1.+alsakzts*alsakzts
    sbzts(1)   =alsakzts
    sbzts(2)   =alsakzts
    sbzts(nz-2)=alsakzts
    sbzts(nz-1)=alsakzts
    sbzts(nz  )=0.
    do k=3,nz-3
       sfzts(k)=alsakzts
       sczts(k)=1.
       sbzts(k)=alsakzts
    enddo

    call prepare (sbxts,scxts,sfxts,ssxts,swxts,nx)
    call prepare (sbzts,sczts,sfzts,sszts,swzts,nz)

  endif ! init x & z second derivative

    ! transpose y to x & y to z
    call transpose_y_to_x(temp,temp1,mydeco)
    call transpose_y_to_z(temp,temp3,mydeco)

    ! x & z second derivative
    call derxxts(dtdxx,temp1,di1,sx,sfxts,ssxts,swxts,mydeco%xsz(1),mydeco%xsz(2),mydeco%xsz(3),0)
    call derzzts(dtdzz,temp3,di3,sz,sfzts,sszts,swzts,mydeco%zsz(1),mydeco%zsz(2),mydeco%zsz(3),0)

    ! gather data in dm1 (y-pencil)
    call transpose_x_to_y(dtdxx,dm1,mydeco)
    call transpose_z_to_y(dtdzz,dtdxyz,mydeco)
    dm1=dm1+dtdxyz

end subroutine xzdiff_temp_solide

subroutine test_sol_min_max(temp,y,n,ratio,repr)
!
! Compute minimum & maximum value on both sides of the solid
! Can be used to check convergence with objective function
!
  use decomp_2d, only : ysize, nrank, real_type, ystart
  use param, only : t, twopi, dt, dx, dz
  use variables, only : nx, nz
  use MPI
!
  implicit none
!
  integer, intent(in) :: n
  real(mytype), intent(in) :: ratio, repr
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(ysize(1),n+1+2,ysize(3)), intent(in) :: temp

  ! Local variables
  integer :: i,j,k,l
  integer, dimension(3) :: mymax
  real(mytype) :: x, z, mysx, mysy, mysz
  real(mytype) :: phimax, phimin, phimax1, phimin1
!  real(mytype), dimension(ysize(1),n+1+2,ysize(3)) :: temp_obj

  ! Print min/max
  phimax=maxval(temp)
  phimin=minval(temp)
  call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,i)
  call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,i)
  if (nrank.eq.0) print *,'Temperature min-max : ',phimin1, phimax1

end subroutine test_sol_min_max

subroutine solide_restart(will_write)
!
! Subroutine called to read or write restart file
! Use MPI-IO and 2DECOMP to read / write
!
  use decomp_2d_io, only : decomp_2d_write_var, decomp_2d_read_var, decomp_2d_write_one
  use param, only : itime, isave
  use MPI

  implicit none

  logical, intent(in) :: will_write

  ! Local variable
  integer :: fh, ierror
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  character(len=3) :: compteur

  if (will_write) then ! write restart file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'solide.dat', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_var(fh,disp,2,temp_bot,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,bot_dm1 ,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,bot_dm2 ,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,bot_dm3 ,mydecomp_bot)
    call decomp_2d_write_var(fh,disp,2,temp_top,mydecomp_top)
    call decomp_2d_write_var(fh,disp,2,top_dm1 ,mydecomp_top)
    call decomp_2d_write_var(fh,disp,2,top_dm2 ,mydecomp_top)
    call decomp_2d_write_var(fh,disp,2,top_dm3 ,mydecomp_top)
    call MPI_FILE_CLOSE(fh,ierror)
    write(compteur,'(I3.3)') itime/isave
    call decomp_2d_write_one(2,temp_bot,'temp_bot'//compteur,mydecomp_bot)
    call decomp_2d_write_one(2,temp_top,'temp_top'//compteur,mydecomp_top)
  else ! read restart file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'solide.dat', &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_var(fh,disp,2,temp_bot,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,bot_dm1 ,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,bot_dm2 ,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,bot_dm3 ,mydecomp_bot)
    call decomp_2d_read_var(fh,disp,2,temp_top,mydecomp_top)
    call decomp_2d_read_var(fh,disp,2,top_dm1 ,mydecomp_top)
    call decomp_2d_read_var(fh,disp,2,top_dm2 ,mydecomp_top)
    call decomp_2d_read_var(fh,disp,2,top_dm3 ,mydecomp_top)
    call MPI_FILE_CLOSE(fh,ierror)
  endif

end subroutine solide_restart

subroutine update_solide_stats()
!
!  Subroutine updating statistics
!
  use decomp_2d, only : nrank, ysize
  use param, only : itime

  implicit none

  ! Local variables
  integer :: coeff
  real(mytype), dimension(ysize(1),ny_sol_bot+1+2,ysize(3)) :: deltab
  real(mytype), dimension(ysize(1),ny_sol_top+1+2,ysize(3)) :: deltat

  coeff=itime-beg_stat_sol
  if (coeff.le.0) then
    if (nrank.eq.0) print*,'ERROR : cannot start to compute solid statistics now, skip.'
    return
  endif

  tempm_1_bot=((coeff-1)*1./coeff)*tempm_1_bot + (1./coeff)*temp_bot
  tempm_2_bot=((coeff-1)*1./coeff)*tempm_2_bot + (1./coeff)*temp_bot**2
  if (.true.) then
    call ydiff_temp_solide(deltab,temp_bot,d2pn_bot,ny_sol_bot)
    deltab=bot_dm1*repr_sol_bot+(1.-sol_imp_var)*deltab
    tdeltat_bot=((coeff-1)*1./coeff)*tdeltat_bot + (1./coeff)*deltab
  endif
  if (mod(itime,20).eq.0) then
    call ydery_temp_sol(temp_bot,yp_bot,ny_sol_bot,mydecomp_bot)
  endif

  tempm_1_top=((coeff-1)*1./coeff)*tempm_1_top + (1./coeff)*temp_top
  tempm_2_top=((coeff-1)*1./coeff)*tempm_2_top + (1./coeff)*temp_top**2
  if (.true.) then
    call ydiff_temp_solide(deltat,temp_top,d2pn_top,ny_sol_top)
    deltat=top_dm1*repr_sol_top+(1.-sol_imp_var)*deltat
    tdeltat_top=((coeff-1)*1./coeff)*tdeltat_top + (1./coeff)*deltat
  endif

end subroutine update_solide_stats

subroutine solide_stats_restart(will_write)
!
!  Subroutine to read/write solid statistics
!
  use decomp_2d_io, only : decomp_2d_write_one, decomp_2d_read_one
  use MPI

  implicit none

  logical, intent(in) :: will_write

  ! Local variable
  logical :: file_exist
  integer :: fh, ierror
  integer (kind=MPI_OFFSET_KIND) :: filesize, disp

  if (will_write) then ! write restart file
    call decomp_2d_write_one(2,tempm_1_bot,'sol_bot_tempm.dat' ,mydecomp_bot)
    call decomp_2d_write_one(2,tempm_2_bot,'sol_bot_temp2m.dat',mydecomp_bot)
    call decomp_2d_write_one(2,tdeltat_bot,'sol_bot_td2tm.dat',mydecomp_bot)
    call decomp_2d_write_one(2,tempm_1_top,'sol_top_tempm.dat' ,mydecomp_top)
    call decomp_2d_write_one(2,tempm_2_top,'sol_top_temp2m.dat',mydecomp_top)
    call decomp_2d_write_one(2,tdeltat_top,'sol_top_td2tm.dat',mydecomp_top)
  else ! read restart file
    INQUIRE(FILE='sol_bot_tempm.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_1_bot,'sol_bot_tempm.dat' ,mydecomp_bot)
    else
      tempm_1_bot=0.
    endif
    INQUIRE(FILE='sol_bot_temp2m.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_2_bot,'sol_bot_temp2m.dat',mydecomp_bot)
    else
      tempm_2_bot=0.
    endif
    INQUIRE(FILE='sol_bot_td2tm.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tdeltat_bot,'sol_bot_td2tm.dat',mydecomp_bot)
    else
      tdeltat_bot=0.
    endif
    INQUIRE(FILE='sol_top_tempm.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_1_top,'sol_top_tempm.dat',mydecomp_top)
    else
      tempm_1_top=0.
    endif
    INQUIRE(FILE='sol_top_temp2m.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tempm_2_top,'sol_top_temp2m.dat',mydecomp_top)
    else
      tempm_2_top=0.
    endif
    INQUIRE(FILE='sol_top_td2tm.dat', EXIST=file_exist)
    if (file_exist) then
      call decomp_2d_read_one(2,tdeltat_top,'sol_top_td2tm.dat',mydecomp_top)
    else
      tdeltat_top=0.
    endif
  endif

end subroutine solide_stats_restart

subroutine allocate_solide_stats()
!
!  Subroutine to allocate solid statistics
!
  use decomp_2d, only : ysize

  implicit none

  ! bottom
  allocate( tempm_1_bot(ysize(1),ny_sol_bot+1+2,ysize(3)) )
  allocate( tempm_2_bot(ysize(1),ny_sol_bot+1+2,ysize(3)) )
  allocate( tdeltat_bot(ysize(1),ny_sol_bot+1+2,ysize(3)) )
  tempm_1_bot=0.
  tempm_2_bot=0.
  tdeltat_bot=0.

  ! top
  allocate( tempm_1_top(ysize(1),ny_sol_top+1+2,ysize(3)) )
  allocate( tempm_2_top(ysize(1),ny_sol_top+1+2,ysize(3)) )
  allocate( tdeltat_top(ysize(1),ny_sol_top+1+2,ysize(3)) )
  tempm_1_top=0.
  tempm_2_top=0.
  tdeltat_top=0.

end subroutine allocate_solide_stats

subroutine my2decomp_solide()
!
! 2DECOMP information for MPI IO & X-diff & Z-diff
!
  use decomp_2d, only : get_decomp_info, decomp_info_init

  implicit none

  type(decomp_info) :: decomp_main

  call get_decomp_info(decomp_main)

  call decomp_info_init(decomp_main%xsz(1),ny_sol_bot+1+2, decomp_main%zsz(3),mydecomp_bot)
  call decomp_info_init(decomp_main%xsz(1),ny_sol_top+1+2, decomp_main%zsz(3),mydecomp_top)

if (decomp_main%ysz(1).ne.mydecomp_bot%ysz(1)) then
  print *,'Error in mydecomp_bot : ',mydecomp_bot%ysz
endif
if (decomp_main%ysz(3).ne.mydecomp_bot%ysz(3)) then
  print *,'Error in mydecomp_bot : ',mydecomp_bot%ysz
endif
if (decomp_main%ysz(1).ne.mydecomp_top%ysz(1)) then
  print *,'Error in mydecomp_top : ',mydecomp_top%ysz
endif
if (decomp_main%ysz(3).ne.mydecomp_top%ysz(3)) then
  print *,'Error in mydecomp_top : ',mydecomp_top%ysz
endif

end subroutine my2decomp_solide

subroutine derxxts(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire) 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx
real(mytype), dimension(ny,nz) :: sx
real(mytype),  dimension(nx):: sfx,ssx,swx 

   do k=1,nz
   do j=1,ny
      tx(1,j,k)=asixts*(ux(2,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx  ,j,k))&
           +bsixts*(ux(3,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx-1,j,k))&
           +csixts*(ux(4,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx-2,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=asixts*(ux(3,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(1   ,j,k))&
           +bsixts*(ux(4,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(nx  ,j,k))&
           +csixts*(ux(5,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(nx-1,j,k))
      rx(2,j,k)=0.
      tx(3,j,k)=asixts*(ux(4,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(2 ,j,k))&
           +bsixts*(ux(5,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(1 ,j,k))&
           +csixts*(ux(6,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(nx,j,k))
      rx(3,j,k)=0.
      do i=4,nx-3
         tx(i,j,k)=asixts*(ux(i+1,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-1,j,k))&
              +bsixts*(ux(i+2,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-2,j,k))&
              +csixts*(ux(i+3,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-3,j,k))
         rx(i,j,k)=0.
      enddo
      tx(nx-2,j,k)=asixts*(ux(nx-1,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-3,j,k))&
           +bsixts*(ux(nx  ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-4,j,k))&
           +csixts*(ux(1   ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-5,j,k))
      rx(nx-2,j,k)=0.
      tx(nx-1,j,k)=asixts*(ux(nx  ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-2,j,k))&
           +bsixts*(ux(1   ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-3,j,k))&
           +csixts*(ux(2   ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-4,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=asixts*(ux(1 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-1,j,k))&
           +bsixts*(ux(2 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-2,j,k))&
           +csixts*(ux(3 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-3,j,k))
      rx(nx  ,j,k)=alsaixts
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
      enddo
      sx(j,k)=(   tx(1,j,k)-alsaixts*tx(nx,j,k))/&
           (1.+rx(1,j,k)-alsaixts*rx(nx,j,k))
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
   enddo
   enddo

return  

end subroutine derxxts

subroutine derzzts(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire) 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: sz 
real(mytype), dimension(nz) :: sfz,ssz,swz

   do j=1,ny
   do i=1,nx
      tz(i,j,1)=askzts*(uz(i,j,2)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz  ))&
           +bskzts*(uz(i,j,3)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz-1))&
           +cskzts*(uz(i,j,4)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz-2))
      rz(i,j,1)=-1.
      tz(i,j,2)=askzts*(uz(i,j,3)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,1 ))&
           +bskzts*(uz(i,j,4)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,nz))&
           +cskzts*(uz(i,j,5)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,nz-1))
      rz(i,j,2)=0.
      tz(i,j,3)=askzts*(uz(i,j,4)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,2 ))&
           +bskzts*(uz(i,j,5)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,1 ))&
           +cskzts*(uz(i,j,6)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,nz))
      rz(i,j,3)=0.
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=askzts*(uz(i,j,k+1)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-1))&
           +bskzts*(uz(i,j,k+2)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-2))&
           +cskzts*(uz(i,j,k+3)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-3))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-2)=askzts*(uz(i,j,nz-1)-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-3))&
           +bskzts*(uz(i,j,nz  )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-4))&
           +cskzts*(uz(i,j,1   )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-5))
      rz(i,j,nz-2)=0.
      tz(i,j,nz-1)=askzts*(uz(i,j,nz  )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-2))&
           +bskzts*(uz(i,j,1   )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-3))&
           +cskzts*(uz(i,j,2   )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-4))
      rz(i,j,nz-1)=0.
      tz(i,j,nz  )=askzts*(uz(i,j,1 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-1))&
           +bskzts*(uz(i,j,2 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-2))&
           +cskzts*(uz(i,j,3 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-3))
      rz(i,j,nz  )=alsakzts
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      rz(i,j,nz)=rz(i,j,nz)*swz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(   tz(i,j,1)-alsakzts*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alsakzts*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo

return

end subroutine derzzts

subroutine ydery_temp_sol(temp,y,n,decomp)

  use decomp_2d_io, only : decomp_2d_write_plane
  use param, only : itime

  implicit none

  integer, parameter :: nplan_sol=4
  ! IO variables
  type(decomp_info), intent(in) :: decomp
  integer, intent(in) :: n
  real(mytype), dimension(n+1+2), intent(in) :: y
  real(mytype), dimension(decomp%ysz(1),n+1+2,decomp%ysz(3)), intent(in) :: temp
  ! Local variables
  integer :: i,j,k
  integer, dimension(nplan_sol) :: jtarget
  real(mytype), dimension(nplan_sol) :: ytarget, d1li, myprod
  real(mytype), dimension(decomp%ysz(1),2*nplan_sol,decomp%ysz(3)) :: output
  character(len=10) :: timer

  output=0.

  jtarget(1) = 131 ! y=0 pour n=128
  jtarget(2) = 115 ! y+=-5 pour n=128
  jtarget(3) = 104 ! y+=-15 pour n=128
  jtarget(4) = 65 ! demi-hauteur pour n=128
  do j=1,nplan_sol
    ytarget(j)=y(jtarget(j))
  enddo

  if (.not.allocated(mat_dery)) then
    allocate(mat_dery(nplan_sol,n+1))
    do i=1,n+1
      d1li=0.
      do k=1,n+1
      if (k.ne.i) then
        myprod=1.
        do j=1,n+1
        if ((j.ne.i).and.(j.ne.k)) then
          myprod=myprod*(ytarget-y(j+1))/(y(i+1)-y(j+1))
        endif
        enddo
        d1li=d1li+myprod/(y(i+1)-y(k+1))
      endif
      enddo
      mat_dery(:,i)=d1li
    enddo
  endif

  ! temperature & normal derivative
  do j=1,nplan_sol
    output(:,2*j-1,:) = temp(:,jtarget(j),:)
    output(:,2*j,:) = 0.
    do i=1,n+1
      output(:,2*j,:) = output(:,2*j,:) + &
        mat_dery(j,i)*temp(:,i+1,:)
    enddo
  enddo

  ! Write to disk
154 format('t',I9.9)
  write(timer,154) itime
  call decomp_2d_write_plane(2,output,5,nplan_sol*2,'slices/tdyt_sol_'//timer//'.dat',decomp)
!  j=1
!  call decomp_2d_write_plane(2,output,5,j,'slices/temp_sol_j131_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/dy_t_sol_j131_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/temp_sol_j115_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/dy_t_sol_j115_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/temp_sol_j104_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/dy_t_sol_j104_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/temp_sol_j065_'//timer//'.dat',decomp); j=j+1
!  call decomp_2d_write_plane(2,output,5,j,'slices/dy_t_sol_j065_'//timer//'.dat',decomp)

end subroutine ydery_temp_sol

end module conjugate_ht
