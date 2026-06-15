#:include 'common.fypp'

!> Program to calculate two-center integrals of Slater-Koster tables.
program main

  use common_accuracy, only : dp
  use common_environment, only : TEnvironment, TEnvironment_init
  use common_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut
  use input, only : readInput
  use twocnt, only : TTwocntInp, TIntegMap, get_twocenter_integrals
  use output, only : write_sktables, write_repulsive, write_onsite
  use cmdargs, only : parse_command_arguments

  implicit none

  !> representation of parsed input for sktwocnt.
  type(TTwocntInp) :: inp

  !> specifies type for mapping integrals.
  type(TIntegMap) :: imap

  !> resulting Hamiltonian and overlap matrices
  real(dp), allocatable :: skham(:,:), skover(:,:)

  !> first-principles repulsive profiles
  real(dp), allocatable :: edcxc(:), eelec(:)

  !> on-site H0 correction SK block (10, ndist); eonsite2 = the reverse (atom 2's) block
  real(dp), allocatable :: eonsite(:,:), eonsite2(:,:)

  type(TEnvironment) :: env

  call initGlobalEnv()
  call TEnvironment_init(env)
#:if WITH_MPI
  call env%initMpi(1)
#:endif
  call parse_command_arguments()

  call readInput(inp, "sktwocnt.in")
  write(stdOut, "(A)") "Reading input done."

  call get_twocenter_integrals(env, inp, imap, skham, skover, edcxc, eelec, eonsite, eonsite2)
  write(stdOut, "(A)") "Two-center integration done."

  if (env%tGlobalLead) then
    call write_sktables(skham, skover)
    call write_repulsive(edcxc, eelec)
    call write_onsite(eonsite)
    call write_onsite(eonsite2, "at2-at1.onsite.dat")
  end if

  call env%destruct()
  call destructGlobalEnv()

end program main
