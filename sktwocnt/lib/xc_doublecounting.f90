#:include 'common.fypp'

!> General XC double-counting for DFTB Slater-Koster generation.
module xc_doublecounting

  use, intrinsic :: iso_c_binding, only : c_size_t
  use common_accuracy, only : dp
  use xcfunctionals, only : xcFunctional
#:if LIBXC_VERSION_MAJOR == 6 or LIBXC_VERSION_MAJOR == 7
  use xc_f03_lib_m, only : xc_f03_func_t, xc_f03_lda_exc_vxc, xc_f03_gga_exc_vxc,&
      & xc_f03_mgga_exc_vxc
#:endif

  implicit none
  private

  public :: getXcDoubleCounting

contains

  !> Integrate the XC double-counting density int(v_xc rho) - E_xc over a grid.
  function getXcDoubleCounting(iXC, xcfunc_x, xcfunc_c, xcfunc_xc, camAlpha, camBeta,&
      & rho, sigma, tau, weight) result(edc)

    !> libxc functional identifier (xcFunctional enum, identical to twocnt.f90)
    integer, intent(in) :: iXC

    !> libxc functional handles: separate exchange/correlation and combined xc
    type(xc_f03_func_t), intent(in) :: xcfunc_x, xcfunc_c, xcfunc_xc

    !> global-hybrid / CAM mixing fractions (for the manual PBE0 / CAMY-PBEh assembly)
    real(dp), intent(in) :: camAlpha, camBeta

    !> superimposed density, |grad rho|^2 and kinetic energy density on the grid (4pi-norm.)
    real(dp), intent(in), contiguous :: rho(:), sigma(:), tau(:)

    !> grid quadrature weights (for \int ... dV)
    real(dp), intent(in) :: weight(:)

    !> resulting XC double-counting energy (Hartree)
    real(dp) :: edc

    integer(c_size_t) :: np
    integer :: n
    real(dp), allocatable :: ex(:), ec(:), exc(:)
    real(dp), allocatable :: vx(:), vc(:), vxc(:)
    real(dp), allocatable :: vxsigma(:), vcsigma(:), vxcsigma(:)
    real(dp), allocatable :: vxtau(:), vctau(:), vxctau(:)
    real(dp), allocatable :: lapl(:), vxlapl(:), vclapl(:), vxclapl(:)
    real(dp), allocatable :: e_dens(:), ivxc_dens(:)
    ! local contiguous copies (libxc element-passing needs guaranteed-contiguous arrays)
    real(dp), allocatable :: rl(:), sl(:), tl(:)

    n = size(rho)
    np = n
    allocate(e_dens(n), ivxc_dens(n))
    e_dens(:) = 0.0_dp
    ivxc_dens(:) = 0.0_dp
    rl = rho
    sl = sigma
    tl = tau

#:if LIBXC_VERSION_MAJOR == 6 or LIBXC_VERSION_MAJOR == 7
    allocate(ex(n), ec(n), exc(n)); ex = 0.0_dp; ec = 0.0_dp; exc = 0.0_dp
    allocate(vx(n), vc(n), vxc(n)); vx = 0.0_dp; vc = 0.0_dp; vxc = 0.0_dp

    select case (iXC)

    ! ---- LDA ------------------------------------------------------------------
    case (xcFunctional%LDA_PW91)
      call xc_f03_lda_exc_vxc(xcfunc_x, np, rl(1), ex(1), vx(1))
      call xc_f03_lda_exc_vxc(xcfunc_c, np, rl(1), ec(1), vc(1))
      e_dens(:) = (ex + ec) * rho
      ivxc_dens(:) = (vx + vc) * rho

    ! ---- GGA (separate exchange/correlation) ----------------------------------
    case (xcFunctional%GGA_PBE96, xcFunctional%GGA_BLYP, xcFunctional%LCY_PBE96,&
        & xcFunctional%GGA_revPBE, xcFunctional%GGA_RPBE, xcFunctional%LC_PBE, xcFunctional%GGA_OPBE)
      allocate(vxsigma(n), vcsigma(n)); vxsigma = 0.0_dp; vcsigma = 0.0_dp
      call xc_f03_gga_exc_vxc(xcfunc_x, np, rl(1), sl(1), ex(1), vx(1), vxsigma(1))
      call xc_f03_gga_exc_vxc(xcfunc_c, np, rl(1), sl(1), ec(1), vc(1), vcsigma(1))
      e_dens(:) = (ex + ec) * rho
      ivxc_dens(:) = (vx + vc) * rho + 2.0_dp * (vxsigma + vcsigma) * sigma

    ! ---- GGA (combined xc: B3LYP, CAMY-B3LYP, B97-D, B97-x, HSE, CAM, ...) -----
    case (xcFunctional%HYB_B3LYP, xcFunctional%CAMY_B3LYP, xcFunctional%GGA_B97D,&
        & xcFunctional%HYB_B97_2, xcFunctional%HYB_B97_3, xcFunctional%WB97X_V,&
        & xcFunctional%HSE06, xcFunctional%LC_WPBE, xcFunctional%CAM_B3LYP, xcFunctional%CAM_PBEH,&
        & xcFunctional%WHPBE0, xcFunctional%HSE12, xcFunctional%GGA_B97_3c, xcFunctional%HYB_B97,&
        & xcFunctional%HYB_B97_1, xcFunctional%HYB_B97_K, xcFunctional%WB97, xcFunctional%WB97X,&
        & xcFunctional%WB97X_D, xcFunctional%WB97X_D3, xcFunctional%HYB_O3LYP)
      allocate(vxcsigma(n)); vxcsigma = 0.0_dp
      call xc_f03_gga_exc_vxc(xcfunc_xc, np, rl(1), sl(1), exc(1), vxc(1), vxcsigma(1))
      e_dens(:) = exc * rho
      ivxc_dens(:) = vxc * rho + 2.0_dp * vxcsigma * sigma

    ! ---- PBE0 / revPBE0 (manual (1-alpha)*X + C) ------------------------------
    case (xcFunctional%HYB_PBE0, xcFunctional%HYB_revPBE0)
      allocate(vxsigma(n), vcsigma(n)); vxsigma = 0.0_dp; vcsigma = 0.0_dp
      call xc_f03_gga_exc_vxc(xcfunc_x, np, rl(1), sl(1), ex(1), vx(1), vxsigma(1))
      call xc_f03_gga_exc_vxc(xcfunc_c, np, rl(1), sl(1), ec(1), vc(1), vcsigma(1))
      e_dens(:) = ((1.0_dp - camAlpha) * ex + ec) * rho
      ivxc_dens(:) = ((1.0_dp - camAlpha) * vx + vc) * rho&
          & + 2.0_dp * ((1.0_dp - camAlpha) * vxsigma + vcsigma) * sigma

    ! ---- meta-GGA (separate x/c: r2SCAN, M06-L, TPSS, MN15L, ...) --------------
    case (xcFunctional%MGGA_r2SCAN, xcFunctional%MGGA_M06L, xcFunctional%HMGGA_MN15,&
        & xcFunctional%HMGGA_M06_2X, xcFunctional%MGGA_TPSS, xcFunctional%MGGA_MN15L,&
        & xcFunctional%HMGGA_M06, xcFunctional%HMGGA_CF22D)
      allocate(vxsigma(n), vcsigma(n), lapl(n), vxlapl(n), vclapl(n), vxtau(n), vctau(n))
      vxsigma = 0.0_dp; vcsigma = 0.0_dp; lapl = 0.0_dp; vxlapl = 0.0_dp; vclapl = 0.0_dp
      vxtau = 0.0_dp; vctau = 0.0_dp
      call xc_f03_mgga_exc_vxc(xcfunc_x, np, rl(1), sl(1), lapl(1), tl(1), ex(1), vx(1),&
          & vxsigma(1), vxlapl(1), vxtau(1))
      call xc_f03_mgga_exc_vxc(xcfunc_c, np, rl(1), sl(1), lapl(1), tl(1), ec(1), vc(1),&
          & vcsigma(1), vclapl(1), vctau(1))
      e_dens(:) = (ex + ec) * rho
      ivxc_dens(:) = (vx + vc) * rho + 2.0_dp * (vxsigma + vcsigma) * sigma&
          & + (vxtau + vctau) * tau

    ! ---- meta-GGA (combined xc: B97M, r2SCANh, r2SCAN0, PW6B95, TPSSh, ...) ----
    case (xcFunctional%MGGA_B97M, xcFunctional%HMGGA_r2SCANh, xcFunctional%HMGGA_r2SCAN0,&
        & xcFunctional%HMGGA_PW6B95, xcFunctional%WB97M_V, xcFunctional%MGGA_TASK,&
        & xcFunctional%HMGGA_TPSSh, xcFunctional%HMGGA_r2SCAN50)
      allocate(vxcsigma(n), lapl(n), vxclapl(n), vxctau(n))
      vxcsigma = 0.0_dp; lapl = 0.0_dp; vxclapl = 0.0_dp; vxctau = 0.0_dp
      call xc_f03_mgga_exc_vxc(xcfunc_xc, np, rl(1), sl(1), lapl(1), tl(1), exc(1), vxc(1),&
          & vxcsigma(1), vxclapl(1), vxctau(1))
      e_dens(:) = exc * rho
      ivxc_dens(:) = vxc * rho + 2.0_dp * vxcsigma * sigma + vxctau * tau

    case default
      ! LCY-BNL/LC-BNL and CAMY-PBEh: assemble like twocnt.f90 if needed;
      ! error stop "getXcDoubleCounting: functional not yet wired (see twocnt.f90 dispatch)"
    end select
#:else
    error stop "getXcDoubleCounting requires LIBXC"
#:endif

    edc = sum(weight * (ivxc_dens - e_dens))

  end function getXcDoubleCounting

end module xc_doublecounting
