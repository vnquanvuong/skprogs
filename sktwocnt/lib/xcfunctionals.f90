!> Module related to supported xc-functionals of the sktwocnt code.
module xcfunctionals

  implicit none
  private

  public :: xcFunctional


  !> Enumerates available xc-functionals.
  type :: TXcFunctionalsEnum

    !=== LDA =================================================================
    !> LDA-PW91
    integer :: LDA_PW91 = 1

    !=== GGA (pure) ======================================
    !> PBE
    integer :: GGA_PBE96 = 2
    integer :: GGA_revPBE = 3
    integer :: GGA_RPBE = 4
    !> BLYP
    integer :: GGA_BLYP = 5
    !> B97 (B97-D; B97-3c composite)
    integer :: GGA_B97D = 6
    integer :: GGA_B97_3c = 7
    !> OPTX
    integer :: GGA_OPBE = 8

    !=== meta-GGA (pure) =================================
    !> SCAN
    integer :: MGGA_r2SCAN = 9
    !> B97
    integer :: MGGA_B97M = 10
    !> TPSS
    integer :: MGGA_TPSS = 11
    !> TASK
    integer :: MGGA_TASK = 12
    !> Minnesota
    integer :: MGGA_M06L = 13
    integer :: MGGA_MN15L = 14

    !=== Global hybrid + GGA =====================================
    !> PBE
    integer :: HYB_PBE0 = 15
    !> BLYP
    integer :: HYB_B3LYP = 16
    !> B97 (B97, B97-1, B97-2, B97-3, B97-K)
    integer :: HYB_B97 = 17
    integer :: HYB_B97_1 = 18
    integer :: HYB_B97_2 = 19
    integer :: HYB_B97_3 = 20
    integer :: HYB_B97_K = 21
    !> revPBE
    integer :: HYB_revPBE0 = 22
    !> OPTX
    integer :: HYB_O3LYP = 23

    !=== Global hybrid + meta-GGA ================================
    !> SCAN
    integer :: HMGGA_r2SCANh = 24
    integer :: HMGGA_r2SCAN0 = 25
    integer :: HMGGA_r2SCAN50 = 26
    !> PW6B95
    integer :: HMGGA_PW6B95 = 27
    !> TPSS
    integer :: HMGGA_TPSSh = 28
    !> Minnesota (M06, M06-2X, MN15)
    integer :: HMGGA_M06 = 29
    integer :: HMGGA_M06_2X = 30
    integer :: HMGGA_MN15 = 31
    !> CF22D
    integer :: HMGGA_CF22D = 32

    !=== Range-separated + GGA (pure LC: camAlpha=0) =========================
    !> Yukawa-screened LC
    integer :: LCY_PBE96 = 33
    integer :: LCY_BNL = 34
    !> erf-screened LC: wPBE / PBE
    integer :: LC_WPBE = 35
    integer :: LC_PBE = 36
    !> erf-screened LC: BNL
    integer :: LC_BNL = 37
    !> erf-screened LC: B97
    integer :: WB97 = 38

    !=== Range-separated + global hybrid + GGA (camAlpha>0) ===================
    !> HSE (screened)
    integer :: HSE06 = 39
    integer :: HSE12 = 40
    !> Yukawa CAM
    integer :: CAMY_B3LYP = 41
    integer :: CAMY_PBEh = 42
    !> erf CAM
    integer :: CAM_B3LYP = 43
    integer :: CAM_PBEH = 44
    integer :: WHPBE0 = 45
    !> erf B97 (wB97X; wB97X-D; wB97X-D3; wB97X-V VV10 dropped)
    integer :: WB97X = 46
    integer :: WB97X_D = 47
    integer :: WB97X_D3 = 48
    integer :: WB97X_V = 49

    !=== Range-separated + global hybrid + meta-GGA (camAlpha>0) ==============
    !> erf B97 (wB97M-V / wB97M-D3; VV10 dropped)
    integer :: WB97M_V = 50

  contains

    procedure :: isLDA => TXcFunctionalsEnum_isLDA
    procedure :: isGGA => TXcFunctionalsEnum_isGGA
    procedure :: isMGGA => TXcFunctionalsEnum_isMGGA
    procedure :: isCombinedXc => TXcFunctionalsEnum_isCombinedXc
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


  pure function TXcFunctionalsEnum_isCombinedXc(this, xcnr) result(isCombinedXc)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    logical :: isCombinedXc

    isCombinedXc = .false.

    if (xcnr == this%HYB_B3LYP .or. xcnr == this%CAMY_B3LYP .or. xcnr == this%GGA_B97D&
        & .or. xcnr == this%HYB_B97_2 .or. xcnr == this%HYB_B97_3 .or. xcnr == this%MGGA_B97M&
        & .or. xcnr == this%HMGGA_r2SCANh .or. xcnr == this%HMGGA_r2SCAN0&
        & .or. xcnr == this%HMGGA_PW6B95 .or. xcnr == this%WB97X_V&
        & .or. xcnr == this%WB97M_V .or. xcnr == this%MGGA_TASK&
        & .or. xcnr == this%HMGGA_TPSSh .or. xcnr == this%HMGGA_r2SCAN50&
        & .or. xcnr == this%HSE06 .or. xcnr == this%LC_WPBE&
        & .or. xcnr == this%CAM_B3LYP .or. xcnr == this%CAM_PBEH&
        & .or. xcnr == this%WHPBE0 .or. xcnr == this%HSE12&
        & .or. xcnr == this%GGA_B97_3c&
        & .or. xcnr == this%HYB_B97 .or. xcnr == this%HYB_B97_1&
        & .or. xcnr == this%HYB_B97_K&
        & .or. xcnr == this%WB97 .or. xcnr == this%WB97X&
        & .or. xcnr == this%WB97X_D .or. xcnr == this%WB97X_D3&
        & .or. xcnr == this%HYB_O3LYP) isCombinedXc = .true.

  end function TXcFunctionalsEnum_isCombinedXc


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

    !> True, if xc-functional index corresponds to a general CAMY functional
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
    !! kernels (wB97X, wB97M): exact exchange = camAlpha*K_full + camBeta*sum_i c_i K_LRYukawa(beta_i*omega).
    logical :: isRangeSepErf

    isRangeSepErf = .false.

    if (xcnr == this%WB97X_V .or. xcnr == this%WB97M_V .or. xcnr == this%HSE06&
        & .or. xcnr == this%LC_WPBE .or. xcnr == this%LC_PBE .or. xcnr == this%LC_BNL&
        & .or. xcnr == this%CAM_B3LYP .or. xcnr == this%CAM_PBEH .or. xcnr == this%WHPBE0&
        & .or. xcnr == this%HSE12&
        & .or. xcnr == this%WB97 .or. xcnr == this%WB97X .or. xcnr == this%WB97X_D&
        & .or. xcnr == this%WB97X_D3) then
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

    if (xcnr < 0 .or. xcnr > 50) then
      isNotImplemented = .true.
    end if

  end function TXcFunctionalsEnum_isNotImplemented

end module xcfunctionals
