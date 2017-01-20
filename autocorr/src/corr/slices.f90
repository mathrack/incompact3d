!***********************************************************
!
program inc_to_para
!BE CAREFUL WITH FORTRAN COMPILER
!USE OPTION ifort -assume byterecl
!***********************************************************
use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f03'

 logical :: file_found

 ! Variables pour lecture fichier incompact3d, real(8) pour double precision
 integer,parameter :: nx=256,ny=8,nz=256, ndt=75000!,nt=81920!ndt=49995, nt=52488
 real(8), parameter :: lx=25.6, ly=2., lz=8.52, lt=ndt*20.d0*0.002d0
 real(8),dimension(ny) :: yp
 real(8) :: dx, dz
 integer :: i,j,k,jj
 integer(8) :: count,count2

 ! Variables pour gestion liste fichiers
 integer :: nca, l
 character(70):: nf

 ! Variables pour FFTW
  integer ( kind = 8 ) :: plan, plan2! type(C_PTR) :: plan, plan2
 real(8), allocatable :: in_data(:,:,:,:), read_data(:,:,:)
 real(8), allocatable :: io_pad(:,:), corr2d(:,:)
 complex(8), allocatable :: out_data(:,:)

allocate(read_data(nx,ny,nz)); read_data=0.d0
allocate(in_data(nx,ny,nz,ndt)); in_data=0.d0
allocate(io_pad(nx,nz)); io_pad=0.d0
allocate(corr2d(nx,nz)); corr2d=0.d0
allocate(out_data(nx/2+1,nz)); out_data=cmplx(0.d0,0.d0)

 nf = ''

 ! dx et dz fonction des CL
 dx=lx/nx!(nx-1) ! nx pour périodique, nx-1 sinon
 dz=lz/nz!(nz-1) ! nz pour périodique, nz-1 sinon

  ! Init fftw
  call dfftw_plan_dft_r2c_2d_(plan,  nx,nz, io_pad, out_data, FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_2d_(plan2, nx,nz, out_data, io_pad, FFTW_ESTIMATE)

  ! Lecture données binaires
 do l=20,20*ndt,20

991 format('../raw/robin_14/tdyt_flu_t',I9.9,'.dat')
  write(nf,991) l

  ! Lecture fichier
  write(*,*) 'Lecture fichier ',trim(nf)
  OPEN(21,FILE=trim(nf),FORM='UNFORMATTED',&
     ACCESS='DIRECT', RECL=8, STATUS='OLD', CONVERT='big_endian')
  COUNT = 1
  DO k=1,nz
   DO j=1,ny
    DO i=1,nx
     READ(21,REC=COUNT) read_data(i,j,k)
     COUNT = COUNT + 1
    ENDDO
   ENDDO
  ENDDO
  CLOSE(21)
  write(*,*) 'Max & Min : ',maxval(read_data),minval(read_data)

!  ! First snapshot for diric has a bug...
!  if (l.eq.20) then
!    read_data(:,1,:)=0.d0
!  endif

  if (maxval(read_data)-minval(read_data).ge.123456789.d0) then
    write(*,*) 'Error in file ',trim(nf)
    write(*,*) 'Abort'
    return
  endif

  in_data(:,:,:,l/20)=read_data

 enddo

 do jj=1,ny
  corr2d=0.d0
  do k=1,ndt
   ! Padding
   io_pad=0.d0
   io_pad(1:nx,1:nz)=in_data(1:nx,jj,1:nz,k)
   ! Compute FFT
   call dfftw_execute_(plan)
   out_data=cmplx(real(out_data)**2+aimag(out_data)**2,0.d0)
   out_data(1,1)=0.d0
   call dfftw_execute_(plan2)
   io_pad=io_pad/( real(nx,kind=8)*real(nz,kind=8) )
   corr2d=corr2d+io_pad/ndt
  enddo

! Write 2D FFT to file
992 format('./robin_14_flu_jj_',I1.1,'_nx_',I3.3,'_nz_',I3.3,'_nt_',I9.9,'.dat')
  write(nf,992) jj,nx,nz,ndt
  OPEN(21,FILE=trim(nf),FORM='UNFORMATTED',&
     ACCESS='DIRECT', RECL=8, STATUS='REPLACE')
  COUNT = 1
  DO j=1,nz
   DO i=1,nx
    WRITE(21,REC=COUNT) corr2d(i,j)
    COUNT = COUNT + 1
   ENDDO
  ENDDO
  CLOSE(21)

993 format('./1dx_robin_14_flu_jj_',I1.1,'_nx_',I3.3,'_nz_',I3.3,'_nt_',I9.9,'.dat')
994 format('./1dz_robin_14_flu_jj_',I1.1,'_nx_',I3.3,'_nz_',I3.3,'_nt_',I9.9,'.dat')
  write(nf,993) jj,nx,nz,ndt
  OPEN(21,FILE=trim(nf),FORM='FORMATTED',STATUS='REPLACE')
  DO i=1,nx
    WRITE(21,*) corr2d(i,1)
  ENDDO
  CLOSE(21)

  write(nf,994) jj,nx,nz,ndt
  OPEN(21,FILE=trim(nf),FORM='FORMATTED',STATUS='REPLACE')
  DO j=1,nz
    WRITE(21,*) corr2d(1,j)
  ENDDO
  CLOSE(21)

 enddo

call dfftw_destroy_plan_(plan)
call dfftw_destroy_plan_(plan2)

deallocate(in_data)
deallocate(io_pad)
deallocate(out_data)
deallocate(corr2d)

  write(*,*) 'Terminus. FFT3d : [',nx,',',nz,',',2*ndt,']'

end program inc_to_para

