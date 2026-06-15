#:include 'common.fypp'

!> Pulay (DIIS) mixer.
module pulaymixer

  use common_accuracy, only : dp
  use lapackroutines, only : getrf, getrs
  implicit none

  private
  public :: TPulayMixer
  public :: TPulayMixer_init, TPulayMixer_reset, TPulayMixer_mix


  !> Contains data for a Pulay/DIIS mixer.
  type TPulayMixer
    private

    !> linear mixing parameter applied to the optimal residual
    real(dp) :: mixParam

    !> maximum history depth
    integer :: maxHist

    !> length of the vectors to mix
    integer :: nElem

    !> number of currently stored history vectors
    integer :: nHist

    !> residual-norm threshold below which DIIS extrapolation is engaged. 
    real(dp) :: diisStart

    !> stored input vectors (nElem, maxHist), chronological (oldest in column 1)
    real(dp), allocatable :: inpHist(:,:)

    !> stored residual vectors (nElem, maxHist)
    real(dp), allocatable :: resHist(:,:)

  end type TPulayMixer


contains

  !> Creates a Pulay mixer.
  subroutine TPulayMixer_init(this, mixParam, maxHist)

    !> Pulay mixer instance on exit
    type(TPulayMixer), intent(out) :: this

    !> linear mixing parameter
    real(dp), intent(in) :: mixParam

    !> maximum history depth (clamped to >= 2)
    integer, intent(in) :: maxHist

    this%mixParam = mixParam
    this%maxHist = max(2, maxHist)
    this%nElem = 0
    this%nHist = 0
    this%diisStart = 0.20_dp

  end subroutine TPulayMixer_init


  !> Resets the mixer for a (possibly new) vector length, clearing history.
  subroutine TPulayMixer_reset(this, nElem)

    !> Pulay mixer instance
    type(TPulayMixer), intent(inout) :: this

    !> length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    this%nElem = nElem
    this%nHist = 0
    if (allocated(this%inpHist)) deallocate(this%inpHist)
    if (allocated(this%resHist)) deallocate(this%resHist)
    allocate(this%inpHist(nElem, this%maxHist), source=0.0_dp)
    allocate(this%resHist(nElem, this%maxHist), source=0.0_dp)

  end subroutine TPulayMixer_reset


  !> Does the actual Pulay/DIIS mixing.
  subroutine TPulayMixer_mix(this, qInpResult, qDiff)

    !> Pulay mixer instance
    type(TPulayMixer), intent(inout) :: this

    !> Input vector on entry, mixed vector on exit
    real(dp), intent(inout) :: qInpResult(:)

    !> Residual (output - input); measure of lack of convergence
    real(dp), intent(in) :: qDiff(:)

    integer :: mm, ii, jj, info
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: bmat(:,:), rhs(:)
    real(dp) :: scal

    @:ASSERT(size(qInpResult) == size(qDiff))

    if (sqrt(dot_product(qDiff, qDiff)) > this%diisStart) then
      this%nHist = 0
      qInpResult(:) = qInpResult + this%mixParam * qDiff
      return
    end if

    if (this%nHist == this%maxHist) then
      this%inpHist(:, 1:this%maxHist - 1) = this%inpHist(:, 2:this%maxHist)
      this%resHist(:, 1:this%maxHist - 1) = this%resHist(:, 2:this%maxHist)
    else
      this%nHist = this%nHist + 1
    end if
    this%inpHist(:, this%nHist) = qInpResult
    this%resHist(:, this%nHist) = qDiff

    mm = this%nHist
    if (mm < 2) then
      ! too little history -> plain linear mixing
      qInpResult(:) = qInpResult + this%mixParam * qDiff
      return
    end if

    allocate(bmat(mm + 1, mm + 1), rhs(mm + 1), ipiv(mm + 1))
    do ii = 1, mm
      do jj = 1, mm
        bmat(ii, jj) = dot_product(this%resHist(:, ii), this%resHist(:, jj))
      end do
    end do
    ! normalize the residual-overlap block for conditioning
    scal = 0.0_dp
    do ii = 1, mm
      scal = max(scal, abs(bmat(ii, ii)))
    end do
    if (scal <= 0.0_dp) scal = 1.0_dp
    bmat(1:mm, 1:mm) = bmat(1:mm, 1:mm) / scal
    bmat(1:mm, mm + 1) = -1.0_dp
    bmat(mm + 1, 1:mm) = -1.0_dp
    bmat(mm + 1, mm + 1) = 0.0_dp
    rhs(:) = 0.0_dp
    rhs(mm + 1) = -1.0_dp

    call getrf(bmat, ipiv, iError=info)
    if (info /= 0) then
      qInpResult(:) = qInpResult + this%mixParam * qDiff
      return
    end if
    call getrs(bmat, ipiv, rhs, iError=info)
    if (info /= 0) then
      qInpResult(:) = qInpResult + this%mixParam * qDiff
      return
    end if

    qInpResult(:) = 0.0_dp
    do ii = 1, mm
      qInpResult(:) = qInpResult(:) + rhs(ii) * (this%inpHist(:, ii) + this%mixParam * this%resHist(:, ii))
    end do

  end subroutine TPulayMixer_mix

end module pulaymixer
