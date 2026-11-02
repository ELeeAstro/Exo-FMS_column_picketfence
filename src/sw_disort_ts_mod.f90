!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Pass a local copy of nlev to TWOSTR
!
! sw: Two-stream DISORT method with modifications from Xianyu Tan
!     Pros: Multiple scattering method, independent check of the sw_Toon89 scheme (same answer)
!     Cons: Two-stream only, ~3x slower than sw_Toon89
!!!

module sw_disort_ts_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private
  public :: sw_disort_ts

contains

  subroutine sw_disort_ts(nlay, nlev, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev), intent(in) :: mu_z
    real(dp), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: i, ntau

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc
    real(dp), dimension(maxcly) :: ggg, ssalb
    real(dp), dimension(maxulv) :: sw_net_d, sw_up_d, sw_down_d, utau
    real(dp) :: umu0, fbeam, Tint
    logical :: planck

    !! Shortwave flux calculation
    if ((mu_z(nlev) > 0.0_dp) .and. (Finc > 0.0_dp)) then

      planck = .False.
      umu0 =  mu_z(nlev)
      sw_net_d(:) = 0.0_dp
      sw_up_d(:) = 0.0_dp
      sw_down_d(:) = 0.0_dp

      !! Initialise arrays
      ggg(:) = 0.0_dp
      ssalb(:) = 0.0_dp
      utau(:) = 0.0_dp
      dtauc(:) = 0.0_dp
      Te_0(:) = 0.0_dp

      fbeam = Finc
      wvnmlo = 0.0_dp
      wvnmhi = 0.0_dp
      Tint = 0.0_dp

      ggg(1:nlay) = gg(:)
      ssalb(1:nlay) = ssa(:)
      utau(1:nlev) = tau_e(:)
      do i = 1, nlay
        dtauc(i) = (tau_e(i+1) - tau_e(i))
      end do
      ! Gas giant lower boundary (albedo = 0), TWOSTR overwrites ntau so we pass a copy of nlev
      ntau = nlev
      call call_twostr(nlay,0.0_dp,Te_0,ggg,ssalb,dtauc,ntau,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0, &
        & sw_net_d(:),sw_up_d(:),sw_down_d(:))
    else
      sw_net_d(:) = 0.0_dp
      sw_up_d(:) = 0.0_dp
      sw_down_d(:) = 0.0_dp
    end if

    sw_net(:) = sw_net_d(1:nlev)
    sw_up(:) = sw_up_d(1:nlev)
    sw_down(:) = sw_down_d(1:nlev)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_disort_ts

end module sw_disort_ts_mod
