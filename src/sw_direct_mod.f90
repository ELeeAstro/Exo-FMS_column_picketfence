!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!
! sw: Direct beam component only
!     Pros: Very very fast
!     Cons: Zero scattering
!!!

module sw_direct_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private :: sw_direct_beam
  public :: sw_direct

contains

  subroutine sw_direct(nlev, tau_e, mu_z, Finc, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlev), intent(in) :: mu_z
    real(dp), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Shortwave flux calculation
    if ((mu_z(nlev) > 0.0_dp) .and. (Finc > 0.0_dp)) then
      call sw_direct_beam(nlev, Finc, mu_z(:), tau_e(:), sw_down(:), sw_up(:))
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_direct

  subroutine sw_direct_beam(nlev, F0_in, mu_in, tau_in, flx_down, flx_up)
    implicit none

    integer, intent(in) :: nlev
    real(dp), intent(in) :: F0_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in

    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    integer :: k
    real(dp), dimension(nlev) :: cum_trans

    !! Calculate direct beam component only

    ! Check mu_z for spherical corrections
    if (mu_in(nlev) == mu_in(1)) then
      ! No zenith correction, use regular method
      flx_down(:) = F0_in* mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
    else
      ! Zenith angle correction, use cumulative transmission function
      cum_trans(1) = tau_in(1)/mu_in(1)
      do k = 1, nlev-1
        cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
      end do
      do k = 1, nlev
        flx_down(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
      end do
    end if

    ! Gas giant lower boundary - nothing is reflected back up
    flx_up(:) = 0.0_dp


  end subroutine sw_direct_beam


end module sw_direct_mod
