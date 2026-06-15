!> Output routines for the sktwocnt code.
module output

  use common_accuracy, only : dp
  use common_globalenv, only : stdOut

  implicit none
  private

  public :: write_sktables, write_repulsive, write_onsite


contains


  !> Writes the on-site H0 correction SK block (10 SK integrals per dimer distance) to
  !! at1-at2.onsite.dat. Columns: ddsig ddpi dddel pdsig pdpi ppsig pppi sdsig spsig sssig (Hartree).
  subroutine write_onsite(eonsite, fname)

    !> on-site SK block, (10, ndist)
    real(dp), intent(in) :: eonsite(:,:)

    !> output file name (defaults to at1-at2.onsite.dat)
    character(*), intent(in), optional :: fname

    !! file identifier and distance index
    integer :: fp, ii
    character(:), allocatable :: fn

    fn = "at1-at2.onsite.dat"
    if (present(fname)) fn = fname
    if (size(eonsite, dim=2) > 0) then
      open(newunit=fp, file=fn, status="replace", action="write")
      write(fp, "(A)") "# on-site H0 correction: ddsig ddpi dddel pdsig pdpi ppsig pppi sdsig spsig sssig"
      write(fp, "(I0)") size(eonsite, dim=2)
      do ii = 1, size(eonsite, dim=2)
        write(fp, "(10ES21.12)") eonsite(:, ii)
      end do
      close(fp)
      write(stdOut, "(A)") "On-site H0 correction block written (" // fn // ")."
    else
      write(stdOut, "(A)") "Nothing to write (on-site H0 correction)."
    end if

  end subroutine write_onsite


  !> Writes the first-principles repulsive profiles (one row per dimer distance) to at1-at2.rep.dat:
  !! col 1: E_dc^xc(R)         (Hartree-free XC double counting)
  !! col 2: (E_nn - E_H)(R)    (electrostatic / charge penetration)
  !! col 3: E_rep(R) = -E_dc^xc(R) + (E_nn - E_H)(R)   (the energy added to the DFTB+ total)
  subroutine write_repulsive(edcxc, eelec)

    !> XC double-counting profile, one value per tabulated distance (Hartree)
    real(dp), intent(in) :: edcxc(:)

    !> electrostatic profile (E_nn - E_H), one value per tabulated distance (Hartree)
    real(dp), intent(in) :: eelec(:)

    !! file identifier and distance index
    integer :: fp, ii

    if (size(edcxc) > 0) then
      open(newunit=fp, file="at1-at2.rep.dat", status="replace", action="write")
      write(fp, "(A)") "# E_dc^xc(Ha)   (E_nn-E_H)(Ha)   E_rep=-E_dc^xc+(E_nn-E_H)(Ha)"
      write(fp, "(I0)") size(edcxc)
      do ii = 1, size(edcxc)
        write(fp, "(3ES21.12)") edcxc(ii), eelec(ii), -edcxc(ii) + eelec(ii)
      end do
      close(fp)
      write(stdOut, "(A)") "Repulsive profiles written (at1-at2.rep.dat)."
    else
      write(stdOut, "(A)") "Nothing to write (repulsive profiles)."
    end if

  end subroutine write_repulsive

  !> Writes tabulated Hamiltonian and overlap matrix to file.
  subroutine write_sktables(skham, skover)

    !> Hamiltonian and overlap matrix
    real(dp), intent(in) :: skham(:,:), skover(:,:)

    if (size(skham, dim=2) > 0) then
      call write_sktable_("at1-at2.ham.dat", skham)
      write(stdOut, "(A)") "SK-table (Hamiltonian) written."
    else
      write(stdOut, "(A)") "Nothing to write (Hamiltonian)."
    end if

    if (size(skover, dim=2) > 0) then
      call write_sktable_("at1-at2.over.dat", skover)
      write(stdOut, "(A)") "SK-table (overlap) written."
    else
      write(stdOut, "(A)") "Nothing to write (overlap)."
    end if

  end subroutine write_sktables


  !> Helper routine writing the SK files.
  subroutine write_sktable_(fname, sktable)

    !> file name
    character(len=*), intent(in) :: fname

    !> Slater-Koster type integrals (Hamiltonian or overlap)
    real(dp), intent(in) :: sktable(:,:)

    !! file identifier
    integer :: fp

    !! number of all nonzero two-center integrals
    integer :: ninteg

    !! number of dimer distances, i.e. lines of written file
    integer :: nline

    !! formatting string
    character(len=20) :: formstr

    ninteg = size(sktable, dim=1)
    nline = size(sktable, dim=2)
    write(formstr, "(A,I0,A)") "(", ninteg, "ES21.12)"

    open(newunit=fp, file=fname, status="replace", action="write")
    write(fp, "(I0)") nline
    write(fp, formstr) sktable
    close(fp)

  end subroutine write_sktable_

end module output
