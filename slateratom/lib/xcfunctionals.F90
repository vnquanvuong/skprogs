#:include 'common.fypp'

!> Module related to supported xc-functionals of the slateratom code.
module xcfunctionals

  use, intrinsic :: iso_c_binding, only : c_size_t
  use common_accuracy, only : dp
  use common_constants, only : rec4pi
  use utilities, only : zeroOutCpotOfEmptyDensitySpinChannels
#:if LIBXC_VERSION_MAJOR == 6
  use xc_f03_lib_m, only : xc_f03_func_t, xc_f03_func_init, xc_f03_func_end, xc_f03_lda_exc_vxc,&
      & xc_f03_gga_exc_vxc, xc_f03_mgga_exc_vxc, xc_f03_func_set_ext_params,&
      & xc_f03_func_set_dens_threshold, XC_LDA_X,&
      & XC_LDA_X_YUKAWA, XC_LDA_C_PW, XC_GGA_X_PBE, XC_GGA_X_B88, XC_GGA_X_SFAT_PBE,&
      & XC_HYB_GGA_XC_B3LYP, XC_HYB_GGA_XC_CAMY_B3LYP, XC_GGA_C_PBE, XC_GGA_C_LYP,&
      & XC_MGGA_X_R2SCAN, XC_MGGA_C_R2SCAN, XC_POLARIZED,&
      & XC_GGA_XC_B97_D, XC_HYB_GGA_XC_B97_2, XC_HYB_GGA_XC_B97_3, XC_MGGA_X_M06_L,&
      & XC_MGGA_C_M06_L, XC_MGGA_XC_B97M_V, XC_HYB_MGGA_XC_R2SCANH, XC_HYB_MGGA_XC_R2SCAN0,&
      & XC_HYB_MGGA_XC_PW6B95, XC_HYB_MGGA_X_MN15, XC_MGGA_C_MN15, XC_HYB_MGGA_X_M06_2X,&
      & XC_MGGA_C_M06_2X, XC_HYB_GGA_XC_WB97X_V, XC_HYB_MGGA_XC_WB97M_V,&
      & XC_GGA_X_PBE_R, XC_GGA_X_RPBE, XC_MGGA_X_TPSS, XC_MGGA_C_TPSS, XC_MGGA_X_TASK,&
      & XC_MGGA_X_MN15_L, XC_MGGA_C_MN15_L, XC_HYB_MGGA_XC_TPSSH, XC_HYB_MGGA_XC_R2SCAN50,&
      & XC_HYB_MGGA_X_M06, XC_MGGA_C_M06, XC_HYB_MGGA_X_CF22D, XC_MGGA_C_CF22D,&
      & XC_HYB_GGA_XC_HSE06, XC_HYB_GGA_XC_LC_WPBE, XC_HYB_GGA_XC_CAM_B3LYP,&
      & XC_HYB_GGA_XC_CAM_PBEH, XC_HYB_GGA_XC_WHPBE0, XC_GGA_X_WPBEH, XC_LDA_X_ERF,&
      & XC_HYB_GGA_XC_HSE12,&
      & XC_GGA_XC_B97_3C, XC_HYB_GGA_XC_B97, XC_HYB_GGA_XC_B97_1, XC_HYB_GGA_XC_B97_K,&
      & XC_HYB_GGA_XC_WB97, XC_HYB_GGA_XC_WB97X, XC_HYB_GGA_XC_WB97X_D, XC_HYB_GGA_XC_WB97X_D3,&
      & XC_HYB_GGA_XC_O3LYP, XC_GGA_X_OPTX
#:elif LIBXC_VERSION_MAJOR == 7
  use xc_f03_lib_m, only : xc_f03_func_t, xc_f03_func_init, xc_f03_func_end, xc_f03_lda_exc_vxc,&
      & xc_f03_gga_exc_vxc, xc_f03_mgga_exc_vxc, xc_f03_func_set_ext_params,&
      & xc_f03_func_set_dens_threshold, XC_POLARIZED
  use xc_f03_funcs_m, only : XC_LDA_X, XC_LDA_X_YUKAWA, XC_LDA_C_PW, XC_GGA_X_PBE, XC_GGA_X_B88,&
      & XC_GGA_X_SFAT_PBE, XC_HYB_GGA_XC_B3LYP, XC_HYB_GGA_XC_CAMY_B3LYP, XC_GGA_C_PBE,&
      & XC_GGA_C_LYP, XC_MGGA_X_R2SCAN, XC_MGGA_C_R2SCAN,&
      & XC_GGA_XC_B97_D, XC_HYB_GGA_XC_B97_2, XC_HYB_GGA_XC_B97_3, XC_MGGA_X_M06_L,&
      & XC_MGGA_C_M06_L, XC_MGGA_XC_B97M_V, XC_HYB_MGGA_XC_R2SCANH, XC_HYB_MGGA_XC_R2SCAN0,&
      & XC_HYB_MGGA_XC_PW6B95, XC_HYB_MGGA_X_MN15, XC_MGGA_C_MN15, XC_HYB_MGGA_X_M06_2X,&
      & XC_MGGA_C_M06_2X, XC_HYB_GGA_XC_WB97X_V, XC_HYB_MGGA_XC_WB97M_V,&
      & XC_GGA_X_PBE_R, XC_GGA_X_RPBE, XC_MGGA_X_TPSS, XC_MGGA_C_TPSS, XC_MGGA_X_TASK,&
      & XC_MGGA_X_MN15_L, XC_MGGA_C_MN15_L, XC_HYB_MGGA_XC_TPSSH, XC_HYB_MGGA_XC_R2SCAN50,&
      & XC_HYB_MGGA_X_M06, XC_MGGA_C_M06, XC_HYB_MGGA_X_CF22D, XC_MGGA_C_CF22D,&
      & XC_HYB_GGA_XC_HSE06, XC_HYB_GGA_XC_LC_WPBE, XC_HYB_GGA_XC_CAM_B3LYP,&
      & XC_HYB_GGA_XC_CAM_PBEH, XC_HYB_GGA_XC_WHPBE0, XC_GGA_X_WPBEH, XC_LDA_X_ERF,&
      & XC_HYB_GGA_XC_HSE12,&
      & XC_GGA_XC_B97_3C, XC_HYB_GGA_XC_B97, XC_HYB_GGA_XC_B97_1, XC_HYB_GGA_XC_B97_K,&
      & XC_HYB_GGA_XC_WB97, XC_HYB_GGA_XC_WB97X, XC_HYB_GGA_XC_WB97X_D, XC_HYB_GGA_XC_WB97X_D3,&
      & XC_HYB_GGA_XC_O3LYP, XC_GGA_X_OPTX
#:endif

  implicit none
  private

  public :: xcFunctional, radial_divergence
  public :: getExcVxc_LDA_PW91
  public :: getExcVxc_GGA_PBE96, getExcVxc_GGA_BLYP
  public :: getExcVxc_LCY_PBE96, getExcVxc_LCY_BNL
  public :: getExcVxc_HYB_B3LYP, getExcVxc_HYB_PBE0
  public :: getExcVxc_CAMY_B3LYP, getExcVxc_CAMY_PBEh
  public :: getExcVxc_MGGA, getExcVxc_GGA_combined, getExcVxc_MGGA_byNr, getExcVxc_GGA_byNr
  public :: getExcVxc_LC_byNr


  interface libxcVxcToInternalVxc
    module procedure :: libxcVxcToInternalVxc_joined
    module procedure :: libxcVxcToInternalVxc_separate
  end interface libxcVxcToInternalVxc


  !> Enumerates available xc-functionals.
  type :: TXcFunctionalsEnum

    !> Hartree-Fock exchange (for spherically symmetric problems)
    !! Strictly speaking not an xc-functional but nevertheless listed here.
    integer :: HF_Exchange = 0

    !> X-Alpha
    integer :: X_Alpha = 1

    !=== LDA ==================================================================
    !> LDA-PW91
    integer :: LDA_PW91 = 2

    !=== GGA (pure, no exact exchange) ====================
    !> PBE
    integer :: GGA_PBE96 = 3
    integer :: GGA_revPBE = 4
    integer :: GGA_RPBE = 5
    !> BLYP
    integer :: GGA_BLYP = 6
    !> B97 (Grimme B97-D; B97-3c composite — pair with D3(/gCP) in DFTB+)
    integer :: GGA_B97D = 7
    integer :: GGA_B97_3c = 8
    !> OPTX
    integer :: GGA_OPBE = 9

    !=== meta-GGA (pure) =================================
    !> SCAN
    integer :: MGGA_r2SCAN = 10
    !> B97 (semilocal part of B97M-V; pair with D3 in DFTB+)
    integer :: MGGA_B97M = 11
    !> TPSS
    integer :: MGGA_TPSS = 12
    !> TASK (exchange-only meta-GGA, designed for band gaps)
    integer :: MGGA_TASK = 13
    !> Minnesota
    integer :: MGGA_M06L = 14
    integer :: MGGA_MN15L = 15

    !=== Global hybrid + GGA (camAlpha * full exact exchange) =====
    !> PBE
    integer :: HYB_PBE0 = 16
    !> BLYP
    integer :: HYB_B3LYP = 17
    !> B97 (B97, B97-1, B97-2, B97-3, B97-K)
    integer :: HYB_B97 = 18
    integer :: HYB_B97_1 = 19
    integer :: HYB_B97_2 = 20
    integer :: HYB_B97_3 = 21
    integer :: HYB_B97_K = 22
    !> revPBE
    integer :: HYB_revPBE0 = 23
    !> OPTX
    integer :: HYB_O3LYP = 24

    !=== Global hybrid + meta-GGA ================================
    !> SCAN (r2SCANh 10%, r2SCAN0 25%, r2SCAN50 50%)
    integer :: HMGGA_r2SCANh = 25
    integer :: HMGGA_r2SCAN0 = 26
    integer :: HMGGA_r2SCAN50 = 27
    !> PW6B95 (28%)
    integer :: HMGGA_PW6B95 = 28
    !> TPSS (TPSSh 10%)
    integer :: HMGGA_TPSSh = 29
    !> Minnesota (M06 27%, M06-2X 54%, MN15 44%)
    integer :: HMGGA_M06 = 30
    integer :: HMGGA_M06_2X = 31
    integer :: HMGGA_MN15 = 32
    !> CF22D (46.28%; pair with D4 in DFTB+)
    integer :: HMGGA_CF22D = 33

    !=== Range-separated + GGA (pure LC: camAlpha=0, only long-range HF) ======
    !> Yukawa-screened LC (kernel exact over STOs)
    integer :: LCY_PBE96 = 34
    integer :: LCY_BNL = 35
    !> erf-screened LC: wPBE / PBE
    integer :: LC_WPBE = 36
    integer :: LC_PBE = 37
    !> erf-screened LC: BNL
    integer :: LC_BNL = 38
    !> erf-screened LC: B97
    integer :: WB97 = 39

    !=== Range-separated + meta-GGA (pure LC): none implemented yet ===========

    !=== Range-separated + global hybrid + GGA (camAlpha>0: full + LR HF) =====
    !> HSE (screened, camAlpha=-camBeta)
    integer :: HSE06 = 40
    integer :: HSE12 = 41
    !> Yukawa CAM
    integer :: CAMY_B3LYP = 42
    integer :: CAMY_PBEh = 43
    !> erf CAM
    integer :: CAM_B3LYP = 44
    integer :: CAM_PBEH = 45
    integer :: WHPBE0 = 46
    !> erf B97 (wB97X CHG-2008; wB97X-D +D2; wB97X-D3 +D3; wB97X-V VV10 dropped)
    integer :: WB97X = 47
    integer :: WB97X_D = 48
    integer :: WB97X_D3 = 49
    integer :: WB97X_V = 50

    !=== Range-separated + global hybrid + meta-GGA (camAlpha>0) ==============
    !> erf B97 (wB97M-V / wB97M-D3; VV10 dropped)
    integer :: WB97M_V = 51

  contains

    procedure :: isLDA => TXcFunctionalsEnum_isLDA
    procedure :: isGGA => TXcFunctionalsEnum_isGGA
    procedure :: isMGGA => TXcFunctionalsEnum_isMGGA
    procedure :: isGlobalHybrid => TXcFunctionalsEnum_isGlobalHybrid
    procedure :: isLongRangeCorrected => TXcFunctionalsEnum_isLongRangeCorrected
    procedure :: isCAMY => TXcFunctionalsEnum_isCAMY
    procedure :: isRangeSepErf => TXcFunctionalsEnum_isRangeSepErf
    procedure :: isNotImplemented => TXcFunctionalsEnum_isNotImplemented

  end type TXcFunctionalsEnum


  !> Container for enumerated xc-functional types.
  type(TXcFunctionalsEnum), parameter :: xcFunctional = TXcFunctionalsEnum()


contains

  pure function TXcFunctionalsEnum_isLDA(this, xcnr) result(isLDA)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to an LDA functional
    logical :: isLDA

    isLDA = .false.

    if (xcnr == this%LDA_PW91) isLDA = .true.

  end function TXcFunctionalsEnum_isLDA


  pure function TXcFunctionalsEnum_isGGA(this, xcnr) result(isGGA)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a GGA functional
    logical :: isGGA

    isGGA = .false.

    if (xcnr == this%GGA_PBE96 .or. xcnr == this%GGA_BLYP .or. xcnr == this%GGA_B97D&
        & .or. xcnr == this%GGA_revPBE .or. xcnr == this%GGA_RPBE&
        & .or. xcnr == this%GGA_B97_3c .or. xcnr == this%GGA_OPBE) isGGA = .true.

  end function TXcFunctionalsEnum_isGGA


  pure function TXcFunctionalsEnum_isMGGA(this, xcnr) result(isMGGA)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a meta-GGA functional
    logical :: isMGGA

    isMGGA = .false.

    if (xcnr == this%MGGA_r2SCAN .or. xcnr == this%MGGA_M06L .or. xcnr == this%MGGA_B97M&
        & .or. xcnr == this%HMGGA_r2SCANh .or. xcnr == this%HMGGA_r2SCAN0&
        & .or. xcnr == this%HMGGA_PW6B95 .or. xcnr == this%HMGGA_MN15&
        & .or. xcnr == this%HMGGA_M06_2X .or. xcnr == this%WB97M_V&
        & .or. xcnr == this%MGGA_TPSS .or. xcnr == this%MGGA_TASK .or. xcnr == this%MGGA_MN15L&
        & .or. xcnr == this%HMGGA_TPSSh .or. xcnr == this%HMGGA_r2SCAN50 .or. xcnr == this%HMGGA_M06&
        & .or. xcnr == this%HMGGA_CF22D) isMGGA = .true.

  end function TXcFunctionalsEnum_isMGGA


  pure function TXcFunctionalsEnum_isLongRangeCorrected(this, xcnr) result(isLongRangeCorrected)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a long-range corrected functional
    logical :: isLongRangeCorrected

    isLongRangeCorrected = .false.

    if (xcnr == this%LCY_PBE96 .or. xcnr == this%LCY_BNL) then
      isLongRangeCorrected = .true.
    end if

  end function TXcFunctionalsEnum_isLongRangeCorrected


  pure function TXcFunctionalsEnum_isGlobalHybrid(this, xcnr) result(isGlobalHybrid)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a global hybrid functional
    logical :: isGlobalHybrid

    isGlobalHybrid = .false.

    if (xcnr == this%HYB_PBE0 .or. xcnr == this%HYB_B3LYP .or. xcnr == this%HYB_B97_2&
        & .or. xcnr == this%HYB_B97_3 .or. xcnr == this%HMGGA_r2SCANh&
        & .or. xcnr == this%HMGGA_r2SCAN0 .or. xcnr == this%HMGGA_PW6B95&
        & .or. xcnr == this%HMGGA_MN15 .or. xcnr == this%HMGGA_M06_2X&
        & .or. xcnr == this%HYB_revPBE0 .or. xcnr == this%HMGGA_TPSSh&
        & .or. xcnr == this%HMGGA_r2SCAN50 .or. xcnr == this%HMGGA_M06&
        & .or. xcnr == this%HMGGA_CF22D&
        & .or. xcnr == this%HYB_B97 .or. xcnr == this%HYB_B97_1&
        & .or. xcnr == this%HYB_B97_K .or. xcnr == this%HYB_O3LYP) then
      isGlobalHybrid = .true.
    end if

  end function TXcFunctionalsEnum_isGlobalHybrid


  pure function TXcFunctionalsEnum_isCAMY(this, xcnr) result(isCamy)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a general CAM functional
    logical :: isCamy

    isCamy = .false.

    if (xcnr == this%CAMY_B3LYP .or. xcnr == this%CAMY_PBEh) then
      isCamy = .true.
    end if

  end function TXcFunctionalsEnum_isCAMY


  pure function TXcFunctionalsEnum_isRangeSepErf(this, xcnr) result(isRangeSepErf)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional uses erf-based range separation realised through a sum of Yukawa
    !! kernels (wB97X, wB97M). The exact-exchange assembly mirrors the CAMY path
    !! (camAlpha*K_full + camBeta*K_erfLR), with K_erfLR = sum_i c_i K_LRYukawa(beta_i*omega).
    logical :: isRangeSepErf

    isRangeSepErf = .false.

    if (xcnr == this%WB97X_V .or. xcnr == this%WB97M_V .or. xcnr == this%HSE06&
        & .or. xcnr == this%LC_WPBE .or. xcnr == this%LC_PBE .or. xcnr == this%LC_BNL&
        & .or. xcnr == this%CAM_B3LYP .or. xcnr == this%CAM_PBEH .or. xcnr == this%WHPBE0&
        & .or. xcnr == this%HSE12 .or. xcnr == this%WB97 .or. xcnr == this%WB97X&
        & .or. xcnr == this%WB97X_D .or. xcnr == this%WB97X_D3) then
      isRangeSepErf = .true.
    end if

  end function TXcFunctionalsEnum_isRangeSepErf


  pure function TXcFunctionalsEnum_isNotImplemented(this, xcnr) result(isNotImplemented)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index is not in the expected range
    logical :: isNotImplemented

    isNotImplemented = .false.

    if (xcnr < 0 .or. xcnr > 51) then
      isNotImplemented = .true.
    end if

  end function TXcFunctionalsEnum_isNotImplemented


  !> Calculates exc and vxc for the LDA-PW91 xc-functional.
  subroutine getExcVxc_LDA_PW91(rho, exc, vxc)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    rhor = transpose(rho) * rec4pi

    allocate(ex(nn))
    allocate(ec(nn))
    allocate(vx(2, nn))
    allocate(vc(2, nn))

    call xc_f03_func_init(xcfunc_x, XC_LDA_X, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_LDA_C_PW, XC_POLARIZED)

    ! exchange
    call xc_f03_lda_exc_vxc(xcfunc_x, nn, rhor(1, 1), ex(1), vx(1, 1))

    ! correlation
    call xc_f03_lda_exc_vxc(xcfunc_c, nn, rhor(1, 1), ec(1), vc(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LDA_PW91


  !> Calculates exc and vxc for the GGA-PBE96 xc-functional.
  subroutine getExcVxc_GGA_PBE96(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_GGA_PBE96


  !> Calculates exc, vxc for a GGA via libxc ids. If cId<=0, xId is a COMBINED xc functional
  !! (e.g. GGA_XC_B97_D); otherwise xId/cId are separate exchange/correlation functionals.
  !! exScale scales the (separate) semilocal EXCHANGE (=1-alpha for a global hybrid GGA without a libxc
  !! HYB id, e.g. revPBE0; =1 for pure functionals and for libxc HYB combined ids that scale internally).
  subroutine getExcVxc_GGA_combined(xId, cId, exScale, abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> libxc id of the GGA exchange (or combined xc) functional
    integer, intent(in) :: xId

    !> libxc id of the correlation functional (cId<=0 => xId is a combined xc functional)
    integer, intent(in) :: cId

    !> scaling of the semilocal exchange (1-alpha for a manual global hybrid GGA, 1 otherwise)
    real(dp), intent(in) :: exScale

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma
    real(dp), allocatable :: vxsigma(:,:), vcsigma(:,:)

    nn = size(rho, dim=1)
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn)); ex(:) = 0.0_dp
    allocate(ec(nn)); ec(:) = 0.0_dp
    allocate(vx(2, nn)); vx(:,:) = 0.0_dp
    allocate(vc(2, nn)); vc(:,:) = 0.0_dp
    allocate(vxsigma(3, nn)); vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn)); vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, xId, XC_POLARIZED)
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    if (cId > 0) then
      call xc_f03_func_init(xcfunc_c, cId, XC_POLARIZED)
      call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
      call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)
      call xc_f03_func_end(xcfunc_c)
    end if

    exc(:) = exScale * ex + ec
    vxc(:,:) = transpose(exScale * vx + vc)
    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, exScale * vxsigma, vcsigma, vxc)

    call xc_f03_func_end(xcfunc_x)

  end subroutine getExcVxc_GGA_combined


  !> Dispatches an erf long-range-corrected xc-functional id (LC-PBE, LC-BNL) to the right short-range
  !! exchange and calls getExcVxc_LC_erf with the range-separation parameter omega.
  subroutine getExcVxc_LC_byNr(xcnr, abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    real(dp), intent(in) :: abcissa(:), dz, dzdr(:), rho(:,:), drho(:,:)
    real(dp), intent(in), allocatable :: sigma(:,:)
    real(dp), intent(in) :: omega
    real(dp), intent(out) :: exc(:), vxc(:,:)

    if (xcnr == xcFunctional%LC_PBE) then
      ! LC-PBE: short-range omega-PBE exchange (WPBEH, single omega ext-param) + PBE correlation
      call getExcVxc_LC_erf(XC_GGA_X_WPBEH, .false., abcissa, dz, dzdr, rho, drho, sigma, omega,&
          & exc, vxc)
    elseif (xcnr == xcFunctional%LC_BNL) then
      call getExcVxc_LC_erf(XC_LDA_X_ERF, .true., abcissa, dz, dzdr, rho, drho, sigma, omega, exc,&
          & vxc)
    end if

  end subroutine getExcVxc_LC_byNr


  !> Calculates exc, vxc for the SEMILOCAL part of an erf long-range-corrected functional: a short-range
  !! exchange functional xId_sr (its range-separation parameter set to omega) plus PBE correlation. With
  !! tLDAx the short-range exchange is LDA-level (XC_LDA_X_ERF -> LC-BNL); otherwise GGA-level
  !! (XC_GGA_X_WPBEH -> LC-PBE). The 100% long-range HF exchange is added in the Hamiltonian. This is
  !! the erf analogue of getExcVxc_LCY_PBE96 / getExcVxc_LCY_BNL (which use Yukawa screening).
  subroutine getExcVxc_LC_erf(xId_sr, tLDAx, abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)

    !> libxc id of the (omega-dependent) short-range exchange functional
    integer, intent(in) :: xId_sr

    !> true, if xId_sr is an LDA-level exchange (LC-BNL); false for a GGA-level one (LC-PBE)
    logical, intent(in) :: tLDAx

    real(dp), intent(in) :: abcissa(:), dz, dzdr(:), rho(:,:), drho(:,:)
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    real(dp), intent(out) :: exc(:), vxc(:,:)

    real(dp), allocatable :: rhor(:,:), ex(:), ec(:), vx(:,:), vc(:,:), vxsigma(:,:), vcsigma(:,:)
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c
    integer(c_size_t) :: nn

    nn = size(rho, dim=1)
    allocate(rhor(2, nn)); rhor(:,:) = transpose(rho) * rec4pi
    allocate(ex(nn)); ex(:) = 0.0_dp
    allocate(ec(nn)); ec(:) = 0.0_dp
    allocate(vx(2, nn)); vx(:,:) = 0.0_dp
    allocate(vc(2, nn)); vc(:,:) = 0.0_dp
    allocate(vxsigma(3, nn)); vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn)); vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, xId_sr, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_x, [omega])
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! short-range exchange (LDA- or GGA-level)
    if (tLDAx) then
      call xc_f03_lda_exc_vxc(xcfunc_x, nn, rhor(1, 1), ex(1), vx(1, 1))
    else
      call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))
    end if

    ! PBE correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)
    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LC_erf


  !> Calculates exc, vxc and the orbital-dependent tau-potential for a tau-dependent meta-GGA, given
  !! its libxc exchange/correlation ids (xId, cId) and exact-exchange scale -- the generic driver for
  !! every meta-GGA in the enum (r2SCAN, B97M, TPSS, TASK, M06-L, MN15-L and the hybrid meta-GGAs).
  !! The rho- and sigma-dependence yields a multiplicative potential (vxc, folded in via
  !! libxcVxcToInternalVxc as for a GGA), while the tau-dependence yields a generalized-Kohn-Sham
  !! operator -1/2 div(vtau grad) whose matrix elements are assembled in dft_exc_matrixelement.
  !! These meta-GGAs are tau-dependent but not Laplacian-dependent, so lapl/vlapl are dummies.
  subroutine getExcVxc_MGGA(xId, cId, exScale, abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc,&
      & vtau)

    !> libxc id of the (meta-)GGA exchange functional; if cId<=0, xId is a COMBINED xc functional
    integer, intent(in) :: xId

    !> libxc id of the correlation functional (cId<=0 => xId is a combined xc functional)
    integer, intent(in) :: cId

    !> scaling of the semilocal EXCHANGE (=1-alpha for a global hybrid meta-GGA, 1 for pure)
    real(dp), intent(in) :: exScale

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> kinetic energy density on grid
    real(dp), intent(in) :: tau(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !> orbital-dependent tau potential on grid
    real(dp), intent(out) :: vtau(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! kinetic energy density in libxc compatible format, i.e. tau/(4pi)
    real(dp), allocatable :: rtau(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange/correlation)
    real(dp), allocatable :: vxsigma(:,:), vcsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of tau (exchange/correlation)
    real(dp), allocatable :: vxtau(:,:), vctau(:,:)

    !! laplacian of the density and its potential derivatives (dummies; r2SCAN does not use them)
    real(dp), allocatable :: lapl(:,:), vxlapl(:,:), vclapl(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi
    allocate(rtau(2, nn))
    rtau(:,:) = transpose(tau) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    allocate(vxtau(2, nn))
    vxtau(:,:) = 0.0_dp
    allocate(vctau(2, nn))
    vctau(:,:) = 0.0_dp

    allocate(lapl(2, nn))
    lapl(:,:) = 0.0_dp
    allocate(vxlapl(2, nn))
    vxlapl(:,:) = 0.0_dp
    allocate(vclapl(2, nn))
    vclapl(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, xId, XC_POLARIZED)
    call xc_f03_func_set_dens_threshold(xcfunc_x, getDensThreshold())

    ! exchange (or, if cId<=0, the combined xc functional)
    call xc_f03_mgga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), lapl(1, 1), rtau(1, 1), ex(1),&
        & vx(1, 1), vxsigma(1, 1), vxlapl(1, 1), vxtau(1, 1))

    if (cId > 0) then
      ! separate correlation functional
      call xc_f03_func_init(xcfunc_c, cId, XC_POLARIZED)
      call xc_f03_func_set_dens_threshold(xcfunc_c, getDensThreshold())
      call xc_f03_mgga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), lapl(1, 1), rtau(1, 1), ec(1),&
          & vc(1, 1), vcsigma(1, 1), vclapl(1, 1), vctau(1, 1))
      call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)
      call xc_f03_func_end(xcfunc_c)
    end if

    ! scale the semilocal exchange by exScale (=1-alpha for a global hybrid meta-GGA, 1 for pure);
    ! for a combined functional (cId<=0) exScale=1 and libxc has already applied any internal scaling
    exc(:) = exScale * ex + ec
    vxc(:,:) = transpose(exScale * vx + vc)
    vtau(:,:) = transpose(exScale * vxtau + vctau)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, exScale * vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)

  end subroutine getExcVxc_MGGA


  !> libxc density floor: points with rho below it return zero xc, guarding the numerically fragile.
  function getDensThreshold() result(thr)

    !> resulting density threshold
    real(dp) :: thr

    character(len=64) :: buf
    integer :: stat, ln

    thr = 1.0e-9_dp
    call get_environment_variable("SLATERATOM_DENSTHR", buf, ln, stat)
    if (stat == 0 .and. ln > 0) then
      read(buf, *, iostat=stat) thr
      if (stat /= 0 .or. thr <= 0.0_dp) thr = 1.0e-9_dp
    end if

  end function getDensThreshold


  !> Dispatches a meta-GGA xc-functional id (xcnr) to libxc functional ids and calls getExcVxc_MGGA.
  !! Hybrid meta-GGAs use the libxc HYB functional for the semilocal part (libxc applies the internal
  !! exchange scaling); the alpha*exact-exchange is added separately in the Hamiltonian (isGlobalHybrid).
  subroutine getExcVxc_MGGA_byNr(xcnr, abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    real(dp), intent(in) :: abcissa(:), dz, dzdr(:), rho(:,:), drho(:,:)
    real(dp), intent(in), allocatable :: sigma(:,:)
    real(dp), intent(in) :: tau(:,:)
    real(dp), intent(out) :: exc(:), vxc(:,:), vtau(:,:)

    if (xcnr == xcFunctional%MGGA_r2SCAN) then
      call getExcVxc_MGGA(XC_MGGA_X_R2SCAN, XC_MGGA_C_R2SCAN, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%MGGA_M06L) then
      call getExcVxc_MGGA(XC_MGGA_X_M06_L, XC_MGGA_C_M06_L, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%MGGA_B97M) then
      ! combined semilocal functional (VV10 nonlocal part is NOT evaluated by libxc -> pair with D3)
      call getExcVxc_MGGA(XC_MGGA_XC_B97M_V, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma, tau,&
          & exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_r2SCANh) then
      call getExcVxc_MGGA(XC_HYB_MGGA_XC_R2SCANH, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_r2SCAN0) then
      call getExcVxc_MGGA(XC_HYB_MGGA_XC_R2SCAN0, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_PW6B95) then
      call getExcVxc_MGGA(XC_HYB_MGGA_XC_PW6B95, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_MN15) then
      call getExcVxc_MGGA(XC_HYB_MGGA_X_MN15, XC_MGGA_C_MN15, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_M06_2X) then
      call getExcVxc_MGGA(XC_HYB_MGGA_X_M06_2X, XC_MGGA_C_M06_2X, 1.0_dp, abcissa, dz, dzdr, rho,&
          & drho, sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%WB97M_V) then
      ! wB97M-V / wB97M-D3: combined erf range-separated semilocal meta-GGA (libxc bakes in the
      ! short-range exchange scaling at omega; the alpha*K_full + beta*K_erfLR exact exchange is added
      ! in the Hamiltonian). VV10 nonlocal part is NOT evaluated by libxc (drop it / pair with D3).
      call getExcVxc_MGGA(XC_HYB_MGGA_XC_WB97M_V, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%MGGA_TPSS) then
      call getExcVxc_MGGA(XC_MGGA_X_TPSS, XC_MGGA_C_TPSS, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%MGGA_TASK) then
      ! TASK is an exchange-only meta-GGA (no correlation), matching the bare libxc/pyscf name
      call getExcVxc_MGGA(XC_MGGA_X_TASK, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma, tau, exc,&
          & vxc, vtau)
    elseif (xcnr == xcFunctional%MGGA_MN15L) then
      call getExcVxc_MGGA(XC_MGGA_X_MN15_L, XC_MGGA_C_MN15_L, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_TPSSh) then
      call getExcVxc_MGGA(XC_HYB_MGGA_XC_TPSSH, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_r2SCAN50) then
      call getExcVxc_MGGA(XC_HYB_MGGA_XC_R2SCAN50, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_M06) then
      call getExcVxc_MGGA(XC_HYB_MGGA_X_M06, XC_MGGA_C_M06, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    elseif (xcnr == xcFunctional%HMGGA_CF22D) then
      ! CF22D global hybrid meta-GGA (46.28% HF added in the Hamiltonian); pair with D4 in DFTB+
      call getExcVxc_MGGA(XC_HYB_MGGA_X_CF22D, XC_MGGA_C_CF22D, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, tau, exc, vxc, vtau)
    end if

  end subroutine getExcVxc_MGGA_byNr


  !> Dispatches a GGA xc-functional id (xcnr) to its libxc x/c ids and calls getExcVxc_GGA_combined.
  !! Generic driver for GGA-based functionals: pure GGAs and the GGA-based global-hybrid and erf/
  !! Yukawa range-separated functionals. libxc returns the semilocal part (applying any internal
  !! exchange scaling); the exact-exchange admixture is added in the Hamiltonian.
  subroutine getExcVxc_GGA_byNr(xcnr, abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    real(dp), intent(in) :: abcissa(:), dz, dzdr(:), rho(:,:), drho(:,:)
    real(dp), intent(in), allocatable :: sigma(:,:)
    real(dp), intent(out) :: exc(:), vxc(:,:)

    if (xcnr == xcFunctional%GGA_B97D) then
      call getExcVxc_GGA_combined(XC_GGA_XC_B97_D, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho, sigma,&
          & exc, vxc)
    elseif (xcnr == xcFunctional%HYB_B97_2) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_B97_2, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HYB_B97_3) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_B97_3, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%WB97X_V) then
      ! wB97X-V: combined erf range-separated semilocal hybrid GGA (libxc bakes in the short-range
      ! exchange scaling at omega; the alpha*K_full + beta*K_erfLR exact exchange is added in the
      ! Hamiltonian). VV10 nonlocal part is NOT evaluated by libxc (drop it / pair with D3).
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_WB97X_V, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%GGA_revPBE) then
      call getExcVxc_GGA_combined(XC_GGA_X_PBE_R, XC_GGA_C_PBE, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%GGA_RPBE) then
      call getExcVxc_GGA_combined(XC_GGA_X_RPBE, XC_GGA_C_PBE, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HYB_revPBE0) then
      ! revPBE0 global hybrid: 25% HF (added in the Hamiltonian) + 0.75*revPBE exchange + revPBE corr.
      call getExcVxc_GGA_combined(XC_GGA_X_PBE_R, XC_GGA_C_PBE, 0.75_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    ! erf range-separated GGAs: libxc returns the range-separated SEMILOCAL xc (its built-in omega; the
    ! short-range PBE exchange already reduced by the hybrid fraction); the exact exchange
    ! camAlpha*K_full + camBeta*K_erfLR is added via the Hamiltonian (isRangeSepErf).
    elseif (xcnr == xcFunctional%HSE06) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_HSE06, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%LC_WPBE) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_LC_WPBE, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%CAM_B3LYP) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_CAM_B3LYP, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%CAM_PBEH) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_CAM_PBEH, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%WHPBE0) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_WHPBE0, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HSE12) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_HSE12, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    ! B97 (all combined libxc GGA semilocal; global hybrids/range-sep add the exact exchange
    ! via the Hamiltonian with their camAlpha/camBeta/omega)
    elseif (xcnr == xcFunctional%GGA_B97_3c) then
      call getExcVxc_GGA_combined(XC_GGA_XC_B97_3C, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HYB_B97) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_B97, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HYB_B97_1) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_B97_1, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HYB_B97_K) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_B97_K, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%WB97) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_WB97, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%WB97X) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_WB97X, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%WB97X_D) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_WB97X_D, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%WB97X_D3) then
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_WB97X_D3, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%HYB_O3LYP) then
      ! HYB_O3LYP global hybrid GGA: libxc returns the scaled semilocal part; 11.61% HF added in Hamiltonian
      call getExcVxc_GGA_combined(XC_HYB_GGA_XC_O3LYP, -1, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    elseif (xcnr == xcFunctional%GGA_OPBE) then
      ! GGA_OPBE pure GGA: OPTX exchange + PBE correlation (no HF), like revPBE
      call getExcVxc_GGA_combined(XC_GGA_X_OPTX, XC_GGA_C_PBE, 1.0_dp, abcissa, dz, dzdr, rho, drho,&
          & sigma, exc, vxc)
    end if

  end subroutine getExcVxc_GGA_byNr


  !> Calculates exc and vxc for the GGA-BLYP xc-functional.
  subroutine getExcVxc_GGA_BLYP(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_B88, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_LYP, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_GGA_BLYP


  !> Calculates exc and vxc for the LCY-PBE96 xc-functional.
  subroutine getExcVxc_LCY_PBE96(abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_SFAT_PBE, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_x, [omega])
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LCY_PBE96


  !> Calculates exc and vxc for the LCY-BNL xc-functional.
  subroutine getExcVxc_LCY_BNL(abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_LDA_X_YUKAWA, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_x, [omega])
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_lda_exc_vxc(xcfunc_x, nn, rhor(1, 1), ex(1), vx(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LCY_BNL


  !> Calculates exc and vxc for the HYB-PBE0 xc-functional.
  subroutine getExcVxc_HYB_PBE0(abcissa, dz, dzdr, rho, drho, sigma, camAlpha, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp
    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    ! build PBE0 functional
    vxcsigma(:,:) = (1.0_dp - camAlpha) * vxsigma + vcsigma
    vxc(:,:) = transpose((1.0_dp - camAlpha) * vx + vc)
    exc(:) = (1.0_dp - camAlpha) * ex + ec

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_HYB_PBE0


  !> Calculates exc and vxc for the HYB-B3LYP xc-functional.
  subroutine getExcVxc_HYB_B3LYP(abcissa, dz, dzdr, rho, drho, sigma, camAlpha, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_xc

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exc energy density on grid
    real(dp), allocatable :: exc_tmp(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vxc_tmp(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(exc_tmp(nn))
    exc_tmp(:) = 0.0_dp

    allocate(vxc_tmp(2, nn))
    vxc_tmp(:,:) = 0.0_dp

    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_B3LYP, XC_POLARIZED)
    ! Adjustable fraction of Fock-type exchange, otherwise standard parametrization taken from
    ! J. Phys. Chem. 1994, 98, 45, 11623-11627; DOI: 10.1021/j100096a001
    call xc_f03_func_set_ext_params(xcfunc_xc, [camAlpha, 0.72_dp, 0.81_dp])

    ! exchange + correlation
    call xc_f03_gga_exc_vxc(xcfunc_xc, nn, rhor(1, 1), sigma(1, 1), exc_tmp(1), vxc_tmp(1, 1),&
        & vxcsigma(1, 1))
    exc(:) = exc_tmp
    vxc(:,:) = transpose(vxc_tmp)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_xc)

  end subroutine getExcVxc_HYB_B3LYP


  !> Calculates exc and vxc for the CAMY-B3LYP xc-functional.
  subroutine getExcVxc_CAMY_B3LYP(abcissa, dz, dzdr, rho, drho, sigma, omega, camAlpha, camBeta,&
      & exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! exc energy density on grid
    real(dp), allocatable :: exc_tmp(:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_xc

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vxc_tmp(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(exc_tmp(nn))
    exc_tmp(:) = 0.0_dp

    allocate(vxc_tmp(2, nn))
    vxc_tmp(:,:) = 0.0_dp

    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_CAMY_B3LYP, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_xc, [0.81_dp, camAlpha + camBeta, -camBeta, omega])

    ! exchange + correlation
    call xc_f03_gga_exc_vxc(xcfunc_xc, nn, rhor(1, 1), sigma(1, 1), exc_tmp(1), vxc_tmp(1, 1),&
        & vxcsigma(1, 1))
    exc(:) = exc_tmp
    vxc(:,:) = transpose(vxc_tmp)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_xc)

  end subroutine getExcVxc_CAMY_B3LYP


  !> Calculates exc and vxc for the CAMY-PBEh xc-functional.
  subroutine getExcVxc_CAMY_PBEh(abcissa, dz, dzdr, rho, drho, sigma, omega, camAlpha, camBeta,&
      & exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_xsr, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ex_sr(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vx_sr(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:), vxsigma_sr(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ex_sr(nn))
    ex_sr(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp

    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vx_sr(2, nn))
    vx_sr(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vxsigma_sr(3, nn))
    vxsigma_sr(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp
    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    ! short-range exchange
    call xc_f03_func_init(xcfunc_xsr, XC_GGA_X_SFAT_PBE, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_xsr, [omega])

    ! full-range exchange
    call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)

    ! correlation
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! short-range exchange
    call xc_f03_gga_exc_vxc(xcfunc_xsr, nn, rhor(1, 1), sigma(1, 1), ex_sr(1), vx_sr(1, 1),&
        & vxsigma_sr(1, 1))
    ! full-range exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1),&
        & vxsigma(1, 1))
    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    ! build CAMY-PBEh functional
    vxcsigma(:,:) = camBeta * vxsigma_sr + (1.0_dp - (camAlpha + camBeta)) * vxsigma + vcsigma
    vxc(:,:) = transpose(camBeta * vx_sr + (1.0_dp - (camAlpha + camBeta)) * vx + vc)
    exc(:) = camBeta * ex_sr + (1.0_dp - (camAlpha + camBeta)) * ex + ec

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_xsr)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_CAMY_PBEh


  !> Converts libXC vxc to our internal representation.
  subroutine libxcVxcToInternalVxc_joined(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> first partial derivative of the energy per unit volume in terms of sigma (xc)
    real(dp), intent(in) :: vxcsigma(:,:)

    !> xc potential on grid
    real(dp), intent(inout) :: vxc(:,:)

    !! temporary storage
    real(dp), allocatable :: tmpv1(:), tmpv2(:)

    !! auxiliary variables
    integer :: iSpin, iSpin2, iSigma

    allocate(tmpv1(size(drho, dim=1)))
    allocate(tmpv2(size(drho, dim=1)))

    ! derivative of E vs. grad n
    do iSpin = 1, 2
      ! the other spin
      iSpin2 = 3 - iSpin
      ! 1 for spin up, 3 for spin down
      iSigma = 2 * iSpin - 1
      tmpv1(:) = vxcsigma(iSigma, :) * drho(:, iSpin) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - 2.0_dp * tmpv2
      tmpv1(:) = vxcsigma(2, :) * drho(:, iSpin2) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - tmpv2
    end do

  end subroutine libxcVxcToInternalVxc_joined


  !> Converts libXC vxc to our internal representation.
  subroutine libxcVxcToInternalVxc_separate(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), intent(in) :: vxsigma(:,:)

    !> first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), intent(in) :: vcsigma(:,:)

    !> xc potential on grid
    real(dp), intent(inout) :: vxc(:,:)

    !! temporary storage
    real(dp), allocatable :: tmpv1(:), tmpv2(:)

    !! auxiliary variables
    integer :: iSpin, iSpin2, iSigma

    allocate(tmpv1(size(drho, dim=1)))
    allocate(tmpv2(size(drho, dim=1)))

    ! derivative of E vs. grad n
    do iSpin = 1, 2
      ! the other spin
      iSpin2 = 3 - iSpin
      ! 1 for spin up, 3 for spin down
      iSigma = 2 * iSpin - 1
      tmpv1(:) = (vxsigma(iSigma, :) + vcsigma(iSigma, :)) * drho(:, iSpin) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - 2.0_dp * tmpv2
      tmpv1(:) = (vxsigma(2, :) +  vcsigma(2, :)) * drho(:, iSpin2) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - tmpv2
    end do

  end subroutine libxcVxcToInternalVxc_separate


  !>
  pure subroutine radial_divergence(ff, rr, dr, rdiv, jacobi)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: rr(:)
    real(dp), intent(in) :: dr
    real(dp), intent(out) :: rdiv(:)
    real(dp), intent(in), optional :: jacobi(:)

    call derive1_5(ff, dr, rdiv, jacobi)
    rdiv(:) = rdiv + 2.0_dp / rr * ff

  end subroutine radial_divergence


  !>
  pure subroutine derive(ff, dx, jacobi)
    real(dp), intent(inout) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(in), optional :: jacobi(:)

    real(dp), allocatable :: tmp1(:)
    integer :: nn

    nn = size(ff)
    allocate(tmp1(nn))
    tmp1(:) = ff
    ff(2:nn - 1) = (ff(3:nn) - ff(1:nn - 2)) / (2.0 * dx)
    ff(1) = (tmp1(2) - tmp1(1)) / dx
    ff(nn) = (tmp1(nn) - tmp1(nn - 1)) / dx
    if (present(jacobi)) then
      ff = ff * jacobi
    end if

  end subroutine derive


  !>
  pure subroutine derive1_5(ff, dx, dfdx, dudx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: dfdx(:)
    real(dp), intent(in), optional :: dudx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) =&
        reshape([&
        & -25.0_dp,  48.0_dp, -36.0_dp,  16.0_dp, -3.0_dp,&
        &  -3.0_dp, -10.0_dp,  18.0_dp,  -6.0_dp,  1.0_dp,&
        &   1.0_dp,  -8.0_dp,   0.0_dp,   8.0_dp, -1.0_dp,&
        &  -1.0_dp,   6.0_dp, -18.0_dp,  10.0_dp,  3.0_dp,&
        &   3.0_dp,  -16.0_dp, 36.0_dp, -48.0_dp, 25.0_dp], [np, np])

    integer :: ngrid
    integer :: ii

    ngrid = size(ff)
    do ii = 1, nleft
      dfdx(ii) = dot_product(coeffs(:, ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      dfdx(ii) = dot_product(coeffs(:, imiddle), ff(ii - nleft:ii + nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      dfdx(ii) = dot_product(coeffs(:, np - (ngrid - ii)), ff(ngrid - np + 1:ngrid))
    end do

    if (present(dudx)) then
      dfdx = dfdx * (dudx / (dxprefac * dx))
    else
      dfdx = dfdx / (dxprefac * dx)
    end if

  end subroutine derive1_5


  !>
  pure subroutine derive2_5(ff, dx, d2fdx2, dudx, d2udx2, dfdx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: d2fdx2(:)
    real(dp), intent(in), optional :: dudx(:), d2udx2(:)
    real(dp), intent(out), target, optional :: dfdx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) = &
        reshape([ &
        &  35.0_dp, -104.0_dp, 114.0_dp, -56.0_dp, 11.0_dp, &
        &  11.0_dp, -20.0_dp, 6.0_dp, 4.0_dp, -1.0_dp, &
        &  -1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp, &
        &  -1.0_dp, 4.0_dp, 6.0_dp, -20.0_dp, 11.0_dp, &
        &  11.0_dp, -56.0_dp, 114.0_dp, -104.0_dp, 35.0_dp], [np, np])

    integer :: ngrid
    integer :: ii
    real(dp), allocatable, target :: dfdxlocal(:)
    real(dp), pointer :: pdfdx(:)

    ngrid = size(ff)
    if (present(dfdx)) then
      pdfdx => dfdx
    elseif (present(d2udx2)) then
      allocate(dfdxlocal(ngrid))
      pdfdx => dfdxlocal
    end if

    do ii = 1, nleft
      d2fdx2(ii) = dot_product(coeffs(:, ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      d2fdx2(ii) = dot_product(coeffs(:, imiddle), ff(ii - nleft:ii + nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      d2fdx2(ii) = dot_product(coeffs(:, np - (ngrid - ii)), ff(ngrid - np + 1:ngrid))
    end do

    if (present(dudx)) then
      d2fdx2 = d2fdx2 * (dudx * dudx / (dxprefac * dx * dx))
    else
      d2fdx2 = d2fdx2 / (dxprefac * dx * dx)
    end if

    if (present(d2udx2) .or. present(dfdx)) then
      call derive1_5(ff, dx, pdfdx)
      if (present(d2udx2)) then
        d2fdx2 = d2fdx2 + pdfdx * d2udx2
      end if
      if (present(dfdx) .and. present(dudx)) then
        dfdx = dfdx * dudx
      end if
    end if

  end subroutine derive2_5

end module xcfunctionals
