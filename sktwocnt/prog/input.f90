!> Module that handles input parsing of configuration and raw data.
module input

  use common_accuracy, only : dp
  use common_globalenv, only : stdOut, abortProgram
  use gridorbital, only : TGridorb2_init
  use twocnt, only : TTwocntInp, TAtomdata
  use xcfunctionals, only : xcFunctional

  implicit none
  private

  public :: readInput

  !> maximum line length of sktwocnt.in file
  integer, parameter :: maxlen = 1024

  !> expected line format when reading sktwocnt.in file
  character(len=*), parameter :: lineformat = "(A1024)"

  !> comment string
  character, parameter :: comment = "#"


contains

  !> Reads and extracts relevant information from 'sktwocnt.in' file.
  subroutine readInput(inp, fname)

    !> instance of parsed input for twocnt
    type(TTwocntInp), intent(out) :: inp

    !> filename
    character(len=*), intent(in) :: fname

    !! file identifier
    integer :: fp

    !! current line index
    integer :: iLine

    !! character buffer
    character(len=maxlen) :: line, buffer1, buffer2

    !! error status
    integer :: iErr

    !! xc-functional type (integer identifier; see the xcFunctional enum in the xcfunctionals module
    !! for the authoritative name <-> id mapping)
    integer :: iXC

    !! potential data columns, summed up in order to receive the total atomic potential
    integer, allocatable :: potcomps(:)

    !! true, if radial grid-orbital 1st/2nd derivative shall be read
    logical :: tReadRadDerivs

    inp%tLC = .false.
    inp%tCam = .false.
    inp%tGlobalHybrid = .false.
    inp%tRangeSepErf = .false.
    inp%mYukawa = 0

    open(newunit=fp, file=fname, form="formatted", action="read")
    ! general part
    iLine = 0

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) buffer1, buffer2, iXC
    call checkerror_(fname, line, iLine, iErr)
    if ((buffer1 /= "hetero") .and. (buffer1 /= "homo")) then
      call error_("Wrong interaction (must be 'hetero' or 'homo')!", fname, line, iLine)
    end if
    inp%tHetero = (buffer1 == "hetero")

    select case (buffer2)
    case("potential")
      inp%tDensitySuperpos = .false.
    case("density")
      inp%tDensitySuperpos = .true.
    case default
      call error_("Wrong superposition mode (must be 'potential' or 'density')!", fname, line,&
          & iline)
    end select

    select case (iXC)
    case(xcFunctional%LDA_PW91)
      ! LDA-PW91
    case(xcFunctional%GGA_PBE96)
      ! GGA-PBE96
    case(xcFunctional%GGA_BLYP)
      ! GGA-BLYP
    case(xcFunctional%MGGA_r2SCAN)
      ! MGGA-r2SCAN (pure meta-GGA, no range separation)
    case(xcFunctional%LCY_PBE96)
      ! LCY-PBE96 (purely long-range corrected)
      inp%tLC = .true.
    case(xcFunctional%LCY_BNL)
      ! LCY-BNL (purely long-range corrected)
      inp%tLC = .true.
    case(xcFunctional%HYB_PBE0)
      ! PBE0 (global hybrid)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HYB_B3LYP)
      ! B3LYP (global hybrid)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%CAMY_B3LYP)
      ! CAMY-B3LYP (general CAM form)
      inp%tCam = .true.
    case(xcFunctional%CAMY_PBEh)
      ! CAMY-PBEh (general CAM form)
      inp%tCam = .true.
    case(xcFunctional%GGA_B97D)
      ! GGA-B97-D (combined semilocal GGA)
    case(xcFunctional%MGGA_M06L)
      ! MGGA-M06-L (pure meta-GGA)
    case(xcFunctional%MGGA_B97M)
      ! MGGA-B97M (semilocal meta-GGA)
    case(xcFunctional%HYB_B97_2)
      ! B97-2 (global hybrid GGA)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HYB_B97_3)
      ! B97-3 (global hybrid GGA)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_r2SCANh)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_r2SCAN0)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_PW6B95)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_MN15)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_M06_2X)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%WB97X_V)
      ! wB97X-V (erf range-separated hybrid GGA, Yukawa-sum exchange)
      inp%tRangeSepErf = .true.
    case(xcFunctional%WB97M_V)
      ! wB97M-V / wB97M-D3 (erf range-separated hybrid meta-GGA, Yukawa-sum exchange)
      inp%tRangeSepErf = .true.
    case(xcFunctional%GGA_revPBE)
      ! revPBE (pure GGA)
    case(xcFunctional%GGA_RPBE)
      ! RPBE (pure GGA)
    case(xcFunctional%MGGA_TPSS)
      ! TPSS (pure meta-GGA)
    case(xcFunctional%MGGA_TASK)
      ! TASK (exchange-only meta-GGA)
    case(xcFunctional%MGGA_MN15L)
      ! MN15-L (pure meta-GGA)
    case(xcFunctional%HYB_revPBE0)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_TPSSh)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_r2SCAN50)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_M06)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HMGGA_CF22D)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%HSE06, xcFunctional%LC_WPBE, xcFunctional%LC_PBE, xcFunctional%LC_BNL,&
        & xcFunctional%CAM_B3LYP, xcFunctional%CAM_PBEH, xcFunctional%WHPBE0, xcFunctional%HSE12)
      ! erf range-separated GGAs (screened/LC/CAM, Yukawa-sum exchange)
      inp%tRangeSepErf = .true.
    case(xcFunctional%GGA_B97_3c)
      ! B97-3c (pure semilocal GGA)
    case(xcFunctional%HYB_B97, xcFunctional%HYB_B97_1, xcFunctional%HYB_B97_K)
      ! B97 / B97-1 / B97-K (global hybrid GGAs)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%WB97, xcFunctional%WB97X, xcFunctional%WB97X_D, xcFunctional%WB97X_D3)
      ! wB97 / wB97X / wB97X-D / wB97X-D3 (erf range-separated hybrid GGAs, Yukawa-sum exchange)
      inp%tRangeSepErf = .true.
    case(xcFunctional%HYB_O3LYP)
      ! HYB_O3LYP (global hybrid GGA)
      inp%tGlobalHybrid = .true.
    case(xcFunctional%GGA_OPBE)
      ! GGA_OPBE (pure GGA: OPTX exchange + PBE correlation)
    case default
      call error_("Unknown exchange-correlation functional!", fname, line, iline)
    end select
    inp%iXC = iXC

    if (inp%iXC == xcFunctional%HYB_B3LYP) then
      ! 20% fraction of HFX hard-coded at the moment
      inp%camAlpha = 0.2_dp
      inp%camBeta = 0.0_dp
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%iXC == xcFunctional%HYB_PBE0) then
      inp%camBeta = 0.0_dp
      call nextline_(fp, iLine, line)
      ! currently only HYB-PBE0 does support arbitrary HFX portions (HYB-B3LYP does not)
      read(line, *, iostat=iErr) inp%camAlpha
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%tGlobalHybrid) then
      ! global hybrids (B97/B97-1/B97-2/B97-3/B97-K, revPBE0, O3LYP and the hybrid meta-GGAs r2SCANh/
      ! r2SCAN0/r2SCAN50, PW6B95, TPSSh, M06/M06-2X, MN15, CF22D): fixed HFX fraction
      inp%camBeta = 0.0_dp
      if (inp%iXC == xcFunctional%HYB_B97_2) then
        inp%camAlpha = 0.21_dp
      elseif (inp%iXC == xcFunctional%HYB_B97_3) then
        inp%camAlpha = 0.269288_dp
      elseif (inp%iXC == xcFunctional%HMGGA_r2SCANh) then
        inp%camAlpha = 0.10_dp
      elseif (inp%iXC == xcFunctional%HMGGA_r2SCAN0) then
        inp%camAlpha = 0.25_dp
      elseif (inp%iXC == xcFunctional%HMGGA_PW6B95) then
        inp%camAlpha = 0.28_dp
      elseif (inp%iXC == xcFunctional%HMGGA_MN15) then
        inp%camAlpha = 0.44_dp
      elseif (inp%iXC == xcFunctional%HMGGA_M06_2X) then
        inp%camAlpha = 0.54_dp
      elseif (inp%iXC == xcFunctional%HYB_revPBE0) then
        inp%camAlpha = 0.25_dp
      elseif (inp%iXC == xcFunctional%HMGGA_TPSSh) then
        inp%camAlpha = 0.10_dp
      elseif (inp%iXC == xcFunctional%HMGGA_r2SCAN50) then
        inp%camAlpha = 0.50_dp
      elseif (inp%iXC == xcFunctional%HMGGA_M06) then
        inp%camAlpha = 0.27_dp
      elseif (inp%iXC == xcFunctional%HMGGA_CF22D) then
        inp%camAlpha = 0.462806_dp
      elseif (inp%iXC == xcFunctional%HYB_B97) then
        inp%camAlpha = 0.1943_dp
      elseif (inp%iXC == xcFunctional%HYB_B97_1) then
        inp%camAlpha = 0.21_dp
      elseif (inp%iXC == xcFunctional%HYB_B97_K) then
        inp%camAlpha = 0.42_dp
      elseif (inp%iXC == xcFunctional%HYB_O3LYP) then
        inp%camAlpha = 0.1161_dp
      end if
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%tLC) then
      inp%camAlpha = 0.0_dp
      inp%camBeta = 1.0_dp
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%omega
      if (inp%omega < 1.0e-08_dp) then
        write(stdOut,'(a)') 'Chosen omega too small!'
        stop
      end if
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%tCam) then
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%omega, inp%camAlpha, inp%camBeta
      if (inp%omega < 1.0e-08_dp) then
        write(stdOut,'(a)') 'Chosen omega too small!'
        stop
      end if
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    elseif (inp%tRangeSepErf) then
      ! erf range-separated hybrids (wB97X/wB97M): fixed camAlpha/camBeta/omega (erf exact exchange =
      ! camAlpha*K_full + camBeta*sum_i c_i K_LRYukawa(beta_i*omega)). Read M, then the Becke grid line.
      if (inp%iXC == xcFunctional%WB97X_V) then
        inp%omega = 0.30_dp; inp%camAlpha = 0.167_dp; inp%camBeta = 0.833_dp
      elseif (inp%iXC == xcFunctional%WB97M_V) then
        inp%omega = 0.30_dp; inp%camAlpha = 0.15_dp; inp%camBeta = 0.85_dp
      elseif (inp%iXC == xcFunctional%HSE06) then
        inp%omega = 0.11_dp; inp%camAlpha = 0.25_dp; inp%camBeta = -0.25_dp
      elseif (inp%iXC == xcFunctional%HSE12) then
        inp%omega = 0.0978977840165_dp; inp%camAlpha = 0.313_dp; inp%camBeta = -0.313_dp
      elseif (inp%iXC == xcFunctional%LC_WPBE) then
        inp%omega = 0.40_dp; inp%camAlpha = 0.0_dp; inp%camBeta = 1.0_dp
      elseif (inp%iXC == xcFunctional%LC_PBE) then
        inp%omega = 0.40_dp; inp%camAlpha = 0.0_dp; inp%camBeta = 1.0_dp
      elseif (inp%iXC == xcFunctional%LC_BNL) then
        inp%omega = 0.33_dp; inp%camAlpha = 0.0_dp; inp%camBeta = 1.0_dp
      elseif (inp%iXC == xcFunctional%CAM_B3LYP) then
        inp%omega = 0.33_dp; inp%camAlpha = 0.19_dp; inp%camBeta = 0.46_dp
      elseif (inp%iXC == xcFunctional%CAM_PBEH) then
        inp%omega = 0.70_dp; inp%camAlpha = 1.0_dp; inp%camBeta = -0.80_dp
      elseif (inp%iXC == xcFunctional%WHPBE0) then
        inp%omega = 0.20_dp; inp%camAlpha = 0.25_dp; inp%camBeta = 0.25_dp
      elseif (inp%iXC == xcFunctional%WB97) then
        inp%omega = 0.40_dp; inp%camAlpha = 0.0_dp; inp%camBeta = 1.0_dp
      elseif (inp%iXC == xcFunctional%WB97X) then
        inp%omega = 0.30_dp; inp%camAlpha = 0.157706_dp; inp%camBeta = 0.842294_dp
      elseif (inp%iXC == xcFunctional%WB97X_D) then
        inp%omega = 0.20_dp; inp%camAlpha = 0.222036_dp; inp%camBeta = 0.777964_dp
      elseif (inp%iXC == xcFunctional%WB97X_D3) then
        inp%omega = 0.25_dp; inp%camAlpha = 0.195728_dp; inp%camBeta = 0.804272_dp
      end if
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%mYukawa
      call checkerror_(fname, line, iLine, iErr)
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) inp%nRadial, inp%nAngular, inp%ll_max, inp%rm
      call checkerror_(fname, line, iLine, iErr)
    end if

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) inp%r0, inp%dr, inp%epsilon, inp%maxdist
    call checkerror_(fname, line, iLine, iErr)

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) inp%ninteg1, inp%ninteg2
    call checkerror_(fname, line, iLine, iErr)

    if (inp%tDensitySuperpos) then
      allocate(potcomps(2))
      potcomps = [2, 3]
    else
      allocate(potcomps(3))
      potcomps = [2, 3, 4]
    end if
    ! atom1's radial derivative is needed when the kinetic acts on atom1 (homonuclear reuses atom1
    ! for atom2) OR for the meta-GGA vtau operator 1/2 grad(phi1).grad(phi2) (needs phi1' on both
    ! centres, also for heteronuclear).
    tReadRadDerivs = (.not. inp%tHetero) .or. xcFunctional%isMGGA(inp%iXC)

    call readatom_(fname, fp, iLine, potcomps, inp%tDensitySuperpos, tReadRadDerivs,&
        & (inp%tGlobalHybrid .or. inp%tLC .or. inp%tCam .or. inp%tRangeSepErf), inp%atom1)
    if (inp%tHetero) then
      call readatom_(fname, fp, iLine, potcomps, inp%tDensitySuperpos, .true., (inp%tGlobalHybrid&
          & .or. inp%tLC .or. inp%tCam .or. inp%tRangeSepErf), inp%atom2)
    end if

    close(fp)

  end subroutine readInput


  !> Fills TAtomdata instance based on slateratom's output.
  subroutine readatom_(fname, fp, iLine, potcomps, tDensitySuperpos, tReadRadDerivs, tNonLocal,&
      & atom)

    !> filename
    character(len=*), intent(in) :: fname

    !> file identifier
    integer, intent(in) :: fp

    !> current line index
    integer, intent(inout) :: iLine

    !> potential data columns, summed up in order to receive the total atomic potential
    integer, intent(in) :: potcomps(:)

    !> true, if density superposition is requested, otherwise potential superposition is applied
    logical, intent(in) :: tDensitySuperpos

    !> true, if radial grid-orbital 1st/2nd derivative shall be read
    logical, intent(in) :: tReadRadDerivs

    !! true, there are non-local exchange contributions to calculate
    logical, intent(in) :: tNonLocal

    !> atomic properties instance
    type(TAtomdata), intent(out) :: atom

    !! character buffer
    character(maxlen) :: line, buffer

    !! temporary storage for checking radial wavefunction sign
    real(dp) :: vals(2)

    !! temporarily stores atomic wavefunction and potential
    real(dp), allocatable :: data(:,:), potval(:)

    !! error status
    integer :: iErr

    !! auxiliary variables
    integer :: ii, imax

    call nextline_(fp, iLine, line)
    if (tNonLocal) then
      read(line, *, iostat=iErr) atom%nBasis, atom%nCore
    else
      read(line, *, iostat=iErr) atom%nBasis
    end if
    call checkerror_(fname, line, iLine, iErr)

    allocate(atom%angmoms(atom%nBasis))
    allocate(atom%rad(atom%nBasis))
    if (tReadRadDerivs) then
      allocate(atom%drad(atom%nBasis))
      allocate(atom%ddrad(atom%nBasis))
    end if

    do ii = 1, atom%nBasis
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) buffer, atom%angmoms(ii)
      call checkerror_(fname, line, iLine, iErr)
      if (tReadRadDerivs) then
        call readdata_(buffer, [1, 3, 4, 5], data)
        call TGridorb2_init(atom%rad(ii), data(:, 1), data(:, 2))
        call TGridorb2_init(atom%drad(ii), data(:, 1), data(:, 3))
        call TGridorb2_init(atom%ddrad(ii), data(:, 1), data(:, 4))
      else
        call readdata_(buffer, [1, 3], data)
        call TGridorb2_init(atom%rad(ii), data(:, 1), data(:, 2))
      end if
      ! check if wave function follows the sign convention
      ! (positive where abs(r * R(r)) has its maximum)
      imax = maxloc(abs(data(:, 1) * data(:, 2)), dim=1)
      if (data(imax, 2) < 0.0_dp) then
        write(stdOut, "(A,F5.2,A)") "Wave function negative at the maximum of radial probability&
            & (r =", data(imax, 1), " Bohr)"
        write(stdOut, "(A)") "Please change the sign of the wave function (and of its derivatives)!"
        write(stdOut, "(A,A,A)") "File: '", trim(buffer), "'"
        stop
      end if
    end do
    call checkangmoms_(atom%angmoms)

    ! read core orbitals, LC functionals only
    if (tNonLocal) then
      allocate(atom%coreAngmoms(atom%nCore))
      allocate(atom%coreOcc(atom%nCore))
      allocate(atom%coreRad(atom%nCore))

      do ii = 1, atom%nCore
        call nextline_(fp, iLine, line)
        read(line, *, iostat=iErr) buffer, atom%coreAngmoms(ii), atom%coreOcc(ii)
        call checkerror_(fname, line, iLine, iErr)
        call readdata_(buffer, [1, 3], data)
        call TGridorb2_init(atom%coreRad(ii), data(:, 1), data(:, 2))
        vals = atom%coreRad(ii)%getValue([0.01_dp, 0.02_dp])
        if ((vals(2) - vals(1)) < 0.0_dp) then
          call atom%coreRad(ii)%rescale(-1.0_dp)
        end if
      end do

    end if

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) buffer
    call checkerror_(fname, line, iLine, iErr)
    call readdata_(buffer, [1, 3, 4, 5], data)
    allocate(potval(size(data, dim=1)))
    potval(:) = 0.0_dp
    do ii = 1, size(potcomps)
      potval(:) = potval + data(:, potcomps(ii))
    end do
    call TGridorb2_init(atom%pot, data(:, 1), potval)

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) buffer
    call checkerror_(fname, line, iLine, iErr)
    if (tDensitySuperpos) then
      ! density{i}.dat columns: r weight rho drho ddrho tau  (tau = col 6, meta-GGA)
      call readdata_(buffer, [1, 3, 4, 5, 6], data)
      call TGridorb2_init(atom%rho, data(:, 1), data(:, 2))
      call TGridorb2_init(atom%drho, data(:, 1), data(:, 3))
      call TGridorb2_init(atom%ddrho, data(:, 1), data(:, 4))
      call TGridorb2_init(atom%tau, data(:, 1), data(:, 5))
    else
      if (trim(line) /= "noread") then
        write(stdOut, "(A,I0,A)") "Line ", iLine, " ignored since density is not needed."
      end if
    end if

  end subroutine readatom_


  !> Reads desired colums of a data file.
  subroutine readdata_(fname, cols, data)

    !> filename
    character(len=*), intent(in) :: fname

    !> desired columns to read from file
    integer, intent(in) :: cols(:)

    !> obtained data on grid with nGrid entries
    real(dp), intent(out), allocatable :: data(:,:)

    !! temporarily stores all columns of a single line in file
    real(dp), allocatable :: tmp(:)

    !! character buffer for current line of file
    character(maxlen) :: line

    !! number of grid points stored in file
    integer :: nGrid

    !! error status
    integer :: iErr

    !! current line
    integer :: iLine

    !! file identifier
    integer :: fp

    !! auxiliary variable
    integer :: ii

    iLine = 1

    allocate(tmp(maxval(cols)))

    open(newunit=fp, file=fname, action="read", form="formatted")

    call nextline_(fp, iLine, line)
    read(line, *, iostat=iErr) nGrid
    call checkerror_(fname, line, iLine, iErr)

    allocate(data(nGrid, size(cols)))
    do ii = 1, nGrid
      call nextline_(fp, iLine, line)
      read(line, *, iostat=iErr) tmp(:)
      call checkerror_(fname, line, iLine, iErr)
      data(ii, :) = tmp(cols)
    end do

    close(fp)

  end subroutine readdata_


  !> Iterates through lines of a file, while respecting an user-def. comment string and empty lines.
  subroutine nextline_(fp, iLine, line)

    !> file identifier
    integer, intent(in) :: fp

    !> current line of the file
    integer, intent(inout) :: iLine

    !> line buffer
    character(maxlen), intent(out) :: line

    !! position of comment string in line if present, otherwise zero
    integer :: ii

    !! temporarily stores an entire line
    character(maxlen) :: buffer

    do while (.true.)
      iLine = iLine + 1
      read(fp, lineformat) buffer
      ii = index(buffer, comment)
      if (ii == 0) then
        line = adjustl(buffer)
      else
        line = adjustl(buffer(1:ii - 1))
      end if
      if (len_trim(line) > 0) exit
    end do

  end subroutine nextline_


  !> Checks range of angular momenta w.r.t. program compatibility.
  subroutine checkangmoms_(angmoms)

    !> angular momenta
    integer, intent(in) :: angmoms(:)

    if (maxval(angmoms) > 4) then
      write(stdOut,*) "Only angular momentum up to 'f' is allowed."
      stop
    end if

  end subroutine checkangmoms_


  !> Error handling.
  subroutine checkerror_(fname, line, iLine, iErr)

    !> filename
    character(len=*), intent(in) :: fname

    !> content of current line
    character(len=*), intent(in) :: line

    !> current line of parsed file
    integer, intent(in) :: iLine

    !> error status
    integer, intent(in) :: iErr

    if (iErr /= 0) then
      call error_("Bad syntax", fname, line, iLine)
    end if

  end subroutine checkerror_


  !> Throws error message.
  subroutine error_(txt, fname, line, iLine)

    !> user-specified error message
    character(len=*), intent(in) :: txt

    !> filename
    character(len=*), intent(in) :: fname

    !> content of erroneous line
    character(len=*), intent(in) :: line

    !> index of erroneous line
    integer, intent(in) :: iLine

    write(stdOut, "(A,A)") "!!! Parsing error: ", txt
    write(stdOut, "(2X,A,A)") "File: ", trim(fname)
    write(stdOut, "(2X,A,I0)") "Line number: ", iLine
    write(stdOut, "(2X,A,A,A)") "Line: '", trim(line), "'"

    flush(stdOut)
    call abortProgram()

  end subroutine error_

end module input
