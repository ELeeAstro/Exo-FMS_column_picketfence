! Solvers ported from Exo-FMS_column_semigrey (76be3de, November 2026).
! Each call solves one band. The driver supplies the stellar band flux and
! multiplies the entire longwave solution by Beta_IR after solving.
module picket_radiation_mod
  use, intrinsic :: iso_fortran_env, only: real64
  use sw_direct_mod, only : sw_direct
  use sw_adding_mod, only : sw_adding
  use sw_SDA_mod, only : sw_SDA
  use sw_Toon89_mod, only : sw_Toon89
  use sw_SH2_mod, only : sw_SH2
  use sw_SH4_mod, only : sw_SH4
  use sw_AD_mod, only : sw_AD
  use sw_Liou88_mod, only : sw_Liou88
  use sw_disort_ts_mod, only : sw_disort_ts
  use sw_disort_mod, only : sw_disort
  !use sw_Feautrier_mod, only : sw_Feautrier

  use lw_AA_E_mod, only : lw_AA_E
  use lw_AA_L_mod, only : lw_AA_L
  use lw_sc_linear_mod, only : lw_sc_linear
  !use lw_sc_scat_mod, only : lw_sc_scat   ! lw_sc_scat disabled (in development)
  use lw_sc_scat_ric_mod, only : lw_sc_scat_ric
  !use lw_sc_scat_iter_mod, only : lw_sc_scat_iter   ! lw_sc_scat_iter disabled (in development)
  use lw_VIM_E_mod, only : lw_VIM_E
  use lw_VIM_L_mod, only : lw_VIM_L
  use lw_Toon89_mod, only : lw_Toon89
  use lw_SH2_mod, only : lw_SH2
  use lw_SH4_mod, only : lw_SH4
  use lw_disort_ts_mod, only : lw_disort_ts
  use lw_disort_mod, only : lw_disort
  use lw_AD_mod, only : lw_AD
  use lw_DDA_mod, only : lw_DDA
  use lw_Fu97_mod, only : lw_Fu97
  use lw_TAS_mod, only : lw_TAS
  !use lw_Feautrier_mod, only : lw_Feautrier   ! lw_Feautrier disabled (in development)
  use lw_DFE_mod, only : lw_DFE

  implicit none
  private
  integer, parameter :: dp = real64
  public :: solve_sw_band, solve_lw_band
contains
  subroutine solve_sw_band(sw_scheme, nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, &
    & sw_up, sw_down, sw_net, asr)
    character(len=*), intent(in) :: sw_scheme
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: tau_Ve(nlev), mu_z, Finc, sw_a(nlay), sw_g(nlay)
    real(dp), intent(out) :: sw_up(nlev), sw_down(nlev), sw_net(nlev), asr
    select case(sw_scheme)
    case('sw_direct')
      ! Direct beam only method with no scattering
      call sw_direct(nlev, tau_Ve, mu_z, Finc, sw_up, sw_down, sw_net, asr)
    case('sw_adding')
      ! Approximate adding method with approximate scattering
      call sw_adding(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_SDA')
      ! Spherical harmonic doubling adding (SDA) four stream method with multiple scattering
      call sw_SDA(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_Toon89')
      ! Toon89 shortwave method with multiple scattering
      call sw_Toon89(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_disort_ts')
      ! Two stream disort method with multiple scattering
      call sw_disort_ts(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_disort')
      ! n-stream disort method with multiple scattering - benchmark method
      call sw_disort(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_SH2')
      ! Two stream spherical harmonic method with multiple scattering (PICASO)
      call sw_SH2(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_SH4')
      ! Four stream spherical harmonic method with multiple scattering (PICASO)
      call sw_SH4(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_AD')
      ! Adding-doubling method with multiple scattering (AddingDoublingRT)
      call sw_AD(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_Liou88')
      ! Delta-four-stream discrete ordinates method with multiple scattering (Liou et al. 1988)
      call sw_Liou88(nlay, nlev, tau_Ve, mu_z, Finc, sw_a, sw_g, sw_up, sw_down, sw_net, asr)
    case('sw_Feautier','sw_Feautrier')
      print*, 'Unsupported sw_scheme (not implemented): ', trim(sw_scheme)
      error stop
    case('none','None')
      sw_up(:) = 0.0_dp
      sw_down(:) = 0.0_dp
      sw_net(:) = 0.0_dp
      asr = 0.0_dp
    case default
      print*, 'Invalid sw_scheme: ', trim(sw_scheme)
      error stop
    end select

  end subroutine solve_sw_band

  subroutine solve_lw_band(lw_scheme, nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, &
    & lw_up, lw_down, lw_net, olr)
    character(len=*), intent(in) :: lw_scheme
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Tl(nlay), pl(nlay), pe(nlev), tau_IRe(nlev), lw_a(nlay), lw_g(nlay), Tint
    real(dp), intent(out) :: lw_up(nlev), lw_down(nlev), lw_net(nlev), olr
    select case(lw_scheme)
    case('lw_AA_E')
      ! Absorption approximation exponential in tau (AA_E) with approximate scattering
      call lw_AA_E(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr) 
    case('lw_AA_L')
      ! Absorption approximation linear in tau (AA_L) with approximate scattering
      call lw_AA_L(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)     
    case('lw_sc_linear')
      ! Short characteristics (sc) with linear interpolants with no scattering
      call lw_sc_linear(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)     
    case('lw_sc_scat')
      ! Short characteristics (sc) discrete ordinates with multiple scattering, direct block-tridiagonal solve
      ! - disabled (in development)
      print*, 'Unsupported lw_scheme (not implemented or disabled): ', trim(lw_scheme)
      error stop
    case('lw_sc_scat_ric')
      ! Same SC equations solved by Riccati invariant imbedding
      call lw_sc_scat_ric(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_sc_scat_iter')
      ! Iterative SC: ALI, Ng-accelerated ALI, symmetric GS, or preconditioned GMRES
      ! - disabled (in development)
      print*, 'Unsupported lw_scheme (not implemented or disabled): ', trim(lw_scheme)
      error stop
    case('lw_VIM_E')
      ! Variational Iteration Method (VIM) with approximate scattering
      call lw_VIM_E(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_VIM_L')
      ! Variational Iteration Method (VIM), linear in tau Planck function
      call lw_VIM_L(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Toon89')
      ! Toon89 longwave method with multiple scattering
      call lw_Toon89(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_SH2')
      ! Two-stream spherical harmonic method (PICASO) with multiple scattering
      call lw_SH2(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_SH4')
      ! Four-stream spherical harmonic method (PICASO) with multiple scattering
      call lw_SH4(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort_ts')
      ! Two stream disort method with multiple scattering
      call lw_disort_ts(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort')
      ! n-stream disort method with multiple scattering - benchmark method
      call lw_disort(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_AD')
      ! Adding-doubling method with multiple scattering
      call lw_AD(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_DDA')
      ! Four-stream discrete ordinates adding method (4DDA) with multiple scattering
      call lw_DDA(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Fu97')
      ! Delta-two- and four-stream combination method (Fu et al. 1997) with multiple scattering
      call lw_Fu97(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_TAS')
      ! Chou et al. (1999) scaling with the Tang et al. (2018) adjustment (TAS)
      call lw_TAS(nlay, nlev, Tl, pl, pe, tau_IRe, lw_a, lw_g, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Feautrier')
      ! Feautrier method with ALI and Ng acceleration (petitRADTRANS), isotropic scattering - disabled (in development)
      print*, 'Unsupported lw_scheme (not implemented or disabled): ', trim(lw_scheme)
      error stop
    case('lw_DFE')
      ! Absorption-only Discontinuous Finite Element (DFE) method
      if (any(abs(lw_a(:)) > 1.0e-12_dp)) then
        print*, 'lw_DFE is absorption-only; set lw_ac = 0.'
        error stop
      end if
      call lw_DFE(nlay, nlev, Tl, pl, pe, tau_IRe, Tint, lw_up, lw_down, lw_net, olr)
    case('none','None')
      lw_up(:) = 0.0_dp
      lw_down(:) = 0.0_dp
      lw_net(:) = 0.0_dp
      olr = 0.0_dp
    case default
      print*, 'Invalid lw_scheme: ', trim(lw_scheme)
      error stop
    end select

  end subroutine solve_lw_band
end module picket_radiation_mod
