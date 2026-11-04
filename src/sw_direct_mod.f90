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
    real(dp), intent(in) :: mu_z, Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Shortwave flux calculation
    if ((mu_z > 0.0_dp) .and. (Finc > 0.0_dp)) then
      call sw_direct_beam(nlev, Finc, mu_z, tau_e, sw_down, sw_up)
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
    real(dp), intent(in) :: F0_in, mu_in
    real(dp), dimension(nlev), intent(in) :: tau_in

    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Calculate direct beam component only
    flx_down(:) = F0_in*mu_in * exp(-(tau_in(:) - tau_in(1))/mu_in)

    ! Gas giant lower boundary - nothing is reflected back up
    flx_up(:) = 0.0_dp


  end subroutine sw_direct_beam


end module sw_direct_mod
