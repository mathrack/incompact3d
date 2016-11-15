module checksum

   use MPI
   use decomp_2d, only : mytype, nrank, &
       xsize, ysize, zsize, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x

   implicit none

   private ! Make everything private unless declared public

   real(mytype), parameter :: xx=1.
   integer, parameter, public :: chksum_size = size(transfer(xx,(/0,0/)))
   logical, save :: chksum_is_working

   ! Temporary work variables / arrays
   integer :: code
   integer, dimension(chksum_size) :: tmprchk
   integer, dimension(chksum_size) :: chkr1, chkr2, chkr3

   public :: init_chksum, chksum, equal_chksum

   contains

   !
   ! Function to compute the checksum of a real 3D array var
   !
   function chksum(var,nx,ny,nz)
      integer, intent(in) :: nx, ny, nz
      real(mytype), dimension(nx,ny,nz), intent(in) :: var
      integer, dimension(chksum_size) :: chksum

      tmprchk = sum(transfer(var,(/0,0/)))
      call MPI_ALLREDUCE(tmprchk,chksum,chksum_size,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,code)

   end function chksum

   !
   ! Subroutine to make sure input arrays have the same checksum
   ! First / second / third array are in X / Y / Z pencil
   ! If switch is provided, reference array is var3.
   !    Otherwise, reference array is var1
   !
   subroutine equal_chksum(var1, var2, var3, switch)
      real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: var1
      real(mytype), dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: var2
      real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(inout) :: var3
      logical, optional, intent(in) :: switch

      if (chksum_is_working) then ! compute checksums
         chkr1 = chksum(var1,xsize(1),xsize(2),xsize(3))
         chkr2 = chksum(var2,ysize(1),ysize(2),ysize(3))
         chkr3 = chksum(var3,zsize(1),zsize(2),zsize(3))
      else ! generate checksums
         chkr1 = 1
         chkr2 = 2
         chkr3 = 3
      endif

      if (present(switch)) then
         if (any(chkr3.ne.chkr2)) call transpose_z_to_y(var3,var2)
         if (any(chkr3.ne.chkr1)) call transpose_y_to_x(var2,var1)
      else
         if (any(chkr1.ne.chkr2)) call transpose_x_to_y(var1,var2)
         if (any(chkr1.ne.chkr3)) call transpose_y_to_z(var2,var3)
      endif

   end subroutine equal_chksum

   !
   ! Subroutine used to check we have a working checksum
   !
   subroutine init_chksum(var1,var2,var3)
      real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(out) :: var1
      real(mytype), dimension(ysize(1),ysize(2),ysize(3)), intent(out) :: var2
      real(mytype), dimension(zsize(1),zsize(2),zsize(3)), intent(out) :: var3

      ! Same random data inside all arrays
      call random_number(var1)
      call transpose_x_to_y(var1,var2)
      call transpose_y_to_z(var2,var3)

      ! Compute checksums
      chkr1 = chksum(var1,xsize(1),xsize(2),xsize(3))
      chkr2 = chksum(var2,ysize(1),ysize(2),ysize(3))
      chkr3 = chksum(var3,zsize(1),zsize(2),zsize(3))

      ! Check checksums
      if (any(chkr1.ne.chkr2).or.any(chkr1.ne.chkr3)) then
         chksum_is_working = .false.
         if (nrank.eq.0) print *,'Checksums based on integer overflow do not work'
      else
         chksum_is_working = .true.
         if (nrank.eq.0) print *,'Checksums based on integer overflow seem to work'
      endif

   end subroutine init_chksum


end module checksum
