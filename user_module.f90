#define my_mod_stats

#ifdef my_mod_stats
#include "mystats.f90"
#endif

module user_specific
!
! Module reserved for user_specific functions
!

contains

subroutine module_user_init(phG,ph1,ph2,ph3,ph4)
  !
  ! Call before main time loop
  !
  use decomp_2d, only : nrank, nx_global, ny_global, nz_global, DECOMP_INFO, mytype
  use param, only : ifirst, ilast, xlx, yly, zlz, ilit
  use variables, only : yp

#ifdef my_mod_stats
  use user_stats
#endif

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

  integer :: j
  real(mytype) :: xx, yy, zz

#ifdef my_mod_stats
  bool_user_stat=.true.
#endif

#ifdef my_mod_stats
  if (bool_user_stat) then
    beg_stat=0 ! Statistics are gathered from itime > beg_stat
    call allocate_user_stats(nx_global, ny_global, nz_global, phG,ph1,ph2,ph3,ph4)
    call read_user_stats(phG,ph1,ph2,ph3,ph4)
  endif
#endif

end subroutine module_user_init

subroutine module_user_write(phG,ph1,ph2,ph3,ph4)
  !
  ! Call at the end of every time step
  !
  use decomp_2d, only : DECOMP_INFO
  use var, only : ux1,uy1,uz1,phi1 ! For probes
  use param, only : itime

#ifdef my_mod_stats
  use user_stats, only : bool_user_stat, beg_stat, &
                         pre_update_user_stats
#endif

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

#ifdef my_mod_stats
  if (bool_user_stat) then
    if (itime.gt.beg_stat) then
      call pre_update_user_stats(phG,ph1,ph2,ph3,ph4)
    endif
  endif
#endif

end subroutine module_user_write

subroutine module_user_post(phG,ph1,ph2,ph3,ph4)
  !
  ! Call at the end of the program
  !
  use decomp_2d, only : DECOMP_INFO

#ifdef my_mod_stats
  use user_stats, only : bool_user_stat, write_user_stats
#endif

  implicit none

  TYPE(DECOMP_INFO), intent(in) :: phG,ph1,ph2,ph3,ph4

#ifdef my_mod_stats
  if (bool_user_stat) call write_user_stats(phG,ph1,ph2,ph3,ph4)
#endif

end subroutine module_user_post

end module user_specific
