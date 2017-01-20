!***********************************************************
!
program inc_to_para
!BE CAREFUL WITH FORTRAN COMPILER
!USE OPTION ifort -assume byterecl
!***********************************************************
!use, intrinsic :: iso_c_binding

implicit none

 logical :: file_found

 ! Variables pour lecture fichier incompact3d, real(8) pour double precision
 integer,parameter :: nx=512,ny=8,nz=512, ndt=75000!,nt=81920!ndt=49995, nt=52488
 real(8), parameter :: lx=25.6, ly=2., lz=8.52, lt=ndt*20.d0*0.002d0
 real(8),dimension(ny) :: yp
 real(8) :: dx, dz
 integer :: i,j,k,jj
 integer(8) :: count,count2

 ! Variables pour gestion liste fichiers
 integer :: nca, l
 character(70):: nf

 ! Variables pour FFTW
 real(8), allocatable :: read_data(:,:,:)
 real(8), allocatable :: corr2d(:,:,:)
 real(8) :: io_pad(nx+2,nz)
 complex(8) :: out_data(nx/2+1,nz)
 equivalence(io_pad,out_data)
 complex(8) :: tmp_spec(nx/2+1,9,nz)

 ! Variables pour ESSL
 integer :: initrc, initcr
 ! x = input = io_pad OU out_data
 integer :: inc2xrc, inc2xcr
 ! y = output = out_data OU io_pad
 integer :: inc2yrc, inc2ycr
 integer :: n1rc, n1cr
 integer :: n2rc, n2cr
 integer :: isignrc, isigncr
 real(8) :: scalerc, scalecr
 real(8), allocatable :: aux1rc(:), aux1cr(:)
 integer :: naux1rc, naux1cr
 real(8), allocatable :: aux2rc(:), aux2cr(:)
 integer :: naux2rc, naux2cr

 external DRCFT2 ! real to complex
 external DCRFT2 ! complex to real

 ! Init tmp_spec
 tmp_spec=cmplx(0.d0,0.d0)

print *,'Debut programme autocorrelation'

allocate(read_data(nx,ny,nz)); read_data=0.d0
allocate(corr2d(nx,9,nz)); corr2d=0.d0

print *,'init essl R2C'

 ! Init ESSL, Real to Complex
 initrc = 1
 inc2xrc = nx+2
 inc2yrc = nx/2+1
 n1rc = nx
 n2rc = nz
 isignrc = 1
 scalerc = 1.d0
 naux1rc = 62000
 allocate(aux1rc(naux1rc))
 naux2rc = 82000
 allocate(aux2rc(naux2rc))
 call DRCFT2(initrc,io_pad,inc2xrc,out_data,inc2yrc,n1rc,n2rc,isignrc,scalerc,aux1rc,naux1rc,aux2rc,naux2rc)
 initrc=0

print *,'done'
print *,'init ESSL C2R'
 ! Init ESSL, Complex to Real
 initcr = 1
 inc2xcr = nx/2+1
 inc2ycr = nx+2
 n1cr = nx
 n2cr = nz
 isigncr = -1
 scalecr = 1.d0
 naux1cr = 62000
 allocate(aux1cr(naux1cr))
 naux2cr = 82000 
 allocate(aux2cr(naux2cr)) 
 call DCRFT2(initcr,out_data,inc2xcr,io_pad,inc2ycr,n1cr,n2cr,isigncr,scalecr,aux1cr,naux1cr,aux2cr,naux2cr)
 initcr=0
print *,'done'

 nf = ''

 ! dx et dz fonction des CL
 dx=lx/nx!(nx-1) ! nx pour périodique, nx-1 sinon
 dz=lz/nz!(nz-1) ! nz pour périodique, nz-1 sinon

  ! Lecture données binaires
 do l=20,20*ndt,20

991 format('./data/conj/tdyt_flu_t',I9.9,'.dat')
  write(nf,991) l

  ! Lecture fichier
  write(*,*) 'Lecture fichier ',trim(nf)
  OPEN(21,FILE=trim(nf),FORM='UNFORMATTED',&
     ACCESS='DIRECT', RECL=8, STATUS='OLD')!, CONVERT='big_endian')
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

  ! First snapshot for diric has a bug...
!  if (l.eq.20) then
!    read_data(:,1,:)=0.d0
!  endif

  if (maxval(read_data)-minval(read_data).ge.123456789.d0) then
    write(*,*) 'Error in file ',trim(nf)
    write(*,*) 'Abort'
    return
  endif

  do jj=1,6
   ! Padding
   io_pad=0.d0
   io_pad(1:nx,1:nz)=read_data(1:nx,jj,1:nz)
   if (jj.le.6) then
    ! Compute FFT
    call DRCFT2(initrc,io_pad,inc2xrc,out_data,inc2yrc,n1rc,n2rc,isignrc,scalerc,aux1rc,naux1rc,aux2rc,naux2rc)
    out_data(1,1)=0.d0
    if (jj.eq.1) then
     tmp_spec(1:nx/2+1,1,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,2,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,3,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,4,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,5,1:nz)=out_data(1:nx/2+1,1:nz)
    endif
    if (jj.eq.2) then
     tmp_spec(1:nx/2+1,1,1:nz)=tmp_spec(1:nx/2+1,1,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
     tmp_spec(1:nx/2+1,6,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,7,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,8,1:nz)=out_data(1:nx/2+1,1:nz)
     tmp_spec(1:nx/2+1,9,1:nz)=out_data(1:nx/2+1,1:nz)
    endif
    if (jj.eq.3) then
     tmp_spec(1:nx/2+1,2,1:nz)=tmp_spec(1:nx/2+1,2,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
     tmp_spec(1:nx/2+1,6,1:nz)=tmp_spec(1:nx/2+1,6,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
    endif
    if (jj.eq.4) then
     tmp_spec(1:nx/2+1,3,1:nz)=tmp_spec(1:nx/2+1,3,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
     tmp_spec(1:nx/2+1,7,1:nz)=tmp_spec(1:nx/2+1,7,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
    endif
    if (jj.eq.5) then
     tmp_spec(1:nx/2+1,4,1:nz)=tmp_spec(1:nx/2+1,4,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
     tmp_spec(1:nx/2+1,8,1:nz)=tmp_spec(1:nx/2+1,8,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
    endif
    if (jj.eq.6) then
     tmp_spec(1:nx/2+1,5,1:nz)=tmp_spec(1:nx/2+1,5,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
     tmp_spec(1:nx/2+1,9,1:nz)=tmp_spec(1:nx/2+1,9,1:nz)*conjg(out_data(1:nx/2+1,1:nz))
    endif
   endif
  enddo

  do jj=1,9
   out_data=tmp_spec(1:nx/2+1,jj,1:nz)
   call DCRFT2(initcr,out_data,inc2xcr,io_pad,inc2ycr,n1cr,n2cr,isigncr,scalecr,aux1cr,naux1cr,aux2cr,naux2cr)
   io_pad(1:nx,1:nz)=io_pad(1:nx,1:nz)/( real(nx,kind=8)*real(nz,kind=8) )
   corr2d(:,jj,:)=corr2d(:,jj,:)+io_pad(1:nx,1:nz)/ndt
  enddo

 enddo

 do jj=1,9

! Write 2D FFT to file
992 format('./corrconj_flu_jj_',I1.1,'_nx_',I3.3,'_nz_',I3.3,'_nt_',I9.9,'.dat')
  write(nf,992) jj,nx,nz,ndt
  OPEN(21,FILE=trim(nf),FORM='UNFORMATTED',&
     ACCESS='DIRECT', RECL=8, STATUS='REPLACE')
  COUNT = 1
  DO j=1,nz
   DO i=1,nx
    WRITE(21,REC=COUNT) corr2d(i,jj,j)
    COUNT = COUNT + 1
   ENDDO
  ENDDO
  CLOSE(21)

993 format('./1dx_corrconj_flu_jj_',I1.1,'_nx_',I3.3,'_nz_',I3.3,'_nt_',I9.9,'.dat')
994 format('./1dz_corrconj_flu_jj_',I1.1,'_nx_',I3.3,'_nz_',I3.3,'_nt_',I9.9,'.dat')
  write(nf,993) jj,nx,nz,ndt
  OPEN(21,FILE=trim(nf),FORM='FORMATTED',STATUS='REPLACE')
  DO i=1,nx
    WRITE(21,*) corr2d(i,jj,1)
  ENDDO
  CLOSE(21)

  write(nf,994) jj,nx,nz,ndt
  OPEN(21,FILE=trim(nf),FORM='FORMATTED',STATUS='REPLACE')
  DO j=1,nz
    WRITE(21,*) corr2d(1,jj,j)
  ENDDO
  CLOSE(21)

 enddo

  write(*,*) 'Terminus fluide. FFT3d : [',nx,',',nz,',',2*ndt,']'

end program inc_to_para

