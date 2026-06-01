!> Module that contains routines for inter- and extrapolation.
module common_interpolation

  use common_accuracy, only : dp

  implicit none
  private

  public :: poly5zero, spline3_free, polyinter, get_cubic_spline, get2ndNaturalSplineDerivs
  public :: TCubicSplineEval, precompCubicSplineBrackets


  !> Precomputed cubic-spline bracket indices and evaluation weights for a fixed set of abscissas and
  !! evaluation points, shared across many splines that differ only in their ordinates. See
  !! precompCubicSplineBrackets.
  type :: TCubicSplineEval

    !> bracketing abscissa indices for each evaluation point
    integer, allocatable :: left(:), right(:)

    !> spline weights for the bracketing ordinates fct(left)/fct(right)
    real(dp), allocatable :: aa(:), bb(:)

    !> spline weights for the bracketing 2nd derivatives dds(left)/dds(right)
    real(dp), allocatable :: wl(:), wr(:)

  end type TCubicSplineEval


contains

  !> Returns the value of a polynomial of 5th degree at x.
  !! \details The polynomial is created with the following boundary conditions:
  !!   Its value, its 1st and 2nd derivatives are zero at x = 0 and agree with the provided values
  !!   at x = dx.
  pure function poly5zero(y0, y0p, y0pp, xx, dx) result(yy)

    !> value of the polynom at x = dx
    real(dp), intent(in) :: y0

    !> value of the 1st derivative at x = dx
    real(dp), intent(in) :: y0p

    !> value of the 2nd derivative at x = dx
    real(dp), intent(in) :: y0pp

    !> point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> point, where the polynomials value and first two derivatives should take the provided values
    real(dp), intent(in) :: dx

    !! value of the polynomial at xx
    real(dp) :: yy

    real(dp) :: dx1, dx2, cc, bb, aa, xr

    ! f(x) = ax^5 + bx^4 + cx^3 + dx^2 + ex + f
    ! f(0) = 0, f'(0) = 0, f''(0) = 0 --> d = e = f = 0

    dx1 = y0p * dx
    dx2 = y0pp * dx**2

    ! c * (dx)**3
    cc = 10.0_dp * y0 - 4.0_dp * dx1 + 0.5_dp * dx2

    ! b * (dx)**4
    bb = - 15.0_dp * y0 + 7.0_dp * dx1 - 1.0_dp * dx2

    ! a * (dx)**5
    aa = 6.0_dp * y0 - 3.0_dp * dx1 + 0.5_dp * dx2

    xr = xx / dx
    yy = ((aa * xr + bb) * xr + cc) * xr**3

  end function poly5zero


  !! Returns the value of a free spline at a certain point.
  !! \details The spline is created with the following boundary conditions:
  !!   Its value, 1st and 2nd derivatives agree with the provided values at
  !!   x = 0 and its value agrees with the provided value at x = dx.
  !! \note If you want the value for a derivative, you have to query them both.
  pure subroutine spline3_free(y0, y0p, y0pp, dx, ydx, xx, yy, yp, ypp)

    !> function value at x = 0
    real(dp), intent(in) :: y0

    !> first derivative at x = 0
    real(dp), intent(in) :: y0p

    !> second derivative at x = 0
    real(dp), intent(in) :: y0pp

    !> function value at dx
    real(dp), intent(in) :: ydx

    !> second fitting point
    real(dp), intent(in) :: dx

    !> point to interpolate
    real(dp), intent(in) :: xx

    !> value of the 3rd order polynomial at xx
    real(dp), intent(out), optional :: yy

    !> first derivative at xx
    real(dp), intent(out), optional :: yp

    !> second derivative at xx
    real(dp), intent(out), optional :: ypp

    !! spline coefficients
    real(dp) :: aa, bb, cc, dd

    !! reciprocal second fitting point
    real(dp) :: dx1

    ! assert(present(yp) .eqv. present(ypp))

    dx1 = 1.0_dp / dx

    aa = y0
    bb = y0p
    cc = 0.5_dp * y0pp
    dd = (((ydx - y0) * dx1 - y0p) * dx1 - 0.5_dp * y0pp) * dx1

    if (present(yy)) then
      yy = ((dd * xx + cc) * xx + bb) * xx + aa
    end if

    if (present(yp)) then
      yp = (3.0_dp * dd * xx + 2.0_dp * cc) * xx + bb
      ypp = 6.0_dp * dd * xx + 2.0_dp * cc
    end if

  end subroutine spline3_free


  !> Polynomial interpolation through given points.
  !! \note The algorithm is based on the Numerical recipes.
  pure function polyinter(xp, yp, xx) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> point, where the polynomial should be evaluated
    real(dp), intent(in) :: xx

    !! value of the polynomial
    real(dp) :: yy

    !! number of interpolation abscissas
    integer :: nn

    !! auxiliary variables
    integer :: icl, ii, mm
    !! work arrays cc/dd are kept at a fixed (compile-time) size so they live on the stack rather than
    !! being heap-allocated on every call. polyinter is invoked once per grid point via the elemental
    !! TGridorb*_getValue (i.e. millions of times), so a per-call heap allocation of these tiny arrays
    !! dominated the runtime/allocator pressure. maxInter is far above any interpolation order used in
    !! skprogs (ninter = 8); only cc(1:nn)/dd(1:nn) are ever accessed, so results are byte-identical.
    integer, parameter :: maxInter = 64
    real(dp) :: cc(maxInter), dd(maxInter)
    real(dp) :: dx, dxnew, dyy, rtmp

    nn = size(xp)

    ! assert(nn > 1)
    ! assert(size(yp) == nn)
    ! assert(nn <= maxInter)

    cc(1:nn) = yp
    dd(1:nn) = yp
    icl = 1
    dx = abs(xx - xp(icl))
    do ii = 2, nn
      dxnew = abs(xx - xp(ii))
      if (dxnew < dx) then
        icl = ii
        dx = dxnew
      end if
    end do
    yy = yp(icl)
    icl = icl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        rtmp = xp(ii) - xp(ii + mm)
        ! if (abs(rtmp) < epsilon(1.0_dp)) then
        !   write(*,*) "Polint failed"
        !   stop
        ! end if
        rtmp = (cc(ii + 1) - dd(ii)) / rtmp
        cc(ii) = (xp(ii) - xx) * rtmp
        dd(ii) = (xp(ii + mm) - xx) * rtmp
      end do
      if (2 * icl < nn - mm) then
        dyy = cc(icl + 1)
      else
        dyy = dd(icl)
        icl = icl - 1
      end if
      yy = yy + dyy
    end do

  end function polyinter


  !> Returns cubic spline value at specified point x = dx.
  !!
  !! \note: Algorithm is based on the Numerical recipes.
  subroutine get_cubic_spline(xx, fct, dds, dx, yy)

    !> abscissas
    real(dp), intent(in) :: xx(:)

    !> ordinates
    real(dp), intent(in) :: fct(:)

    !> spline's 2nd derivatives, derived from the data set via set_cubic_spline()
    real(dp), intent(in) :: dds(:)

    !> evaluation point
    real(dp), intent(in) :: dx

    !> cubic spline value at specified point dx
    real(dp), intent(out) :: yy

    !! bisection search
    integer :: left, right, middle

    !! difference of bracketing abscissas
    real(dp) :: step

    !! auxiliary variables for spline evaluation
    real(dp) :: aa, bb

    ! initialize bisection
    left = 1
    right = size(xx)

    ! bisection search
    do
      if ((right - left) <= 1) exit
      middle = (right + left) / 2
      if (dx >= xx(middle)) then
        left = middle
      else
        right = middle
      end if
    end do

    step = xx(right) - xx(left)
    aa = (xx(right) - dx) / step
    bb = (dx - xx(left)) / step

    ! calculate the spline value
    yy = aa * fct(left) + bb * fct(right) + step**2 / 6.0_dp * (aa**2 - 1.0_dp) * aa * dds(left)&
        & + step**2 / 6.0_dp * (bb**2 - 1.0_dp) * bb * dds(right)

  end subroutine get_cubic_spline


  !> Precomputes cubic-spline bracket indices and evaluation weights for a set of evaluation points
  !! against fixed abscissas. This lets a caller that evaluates many splines sharing the same
  !! abscissas (xx) and evaluation points (dxs) — but differing only in the ordinates (fct) and their
  !! 2nd derivatives (dds) — perform the bisection search and weight setup only once, then evaluate
  !!   yy = this%aa(i)*fct(this%left(i)) + this%bb(i)*fct(this%right(i))
  !!      + this%wl(i)*dds(this%left(i)) + this%wr(i)*dds(this%right(i))
  !! for each point i. The bisection and weight expressions are identical to get_cubic_spline (same
  !! operations in the same order), so the resulting spline values are byte-identical. Buffers are
  !! (re)allocated only when absent or the wrong size, so repeated calls with constant sizes are
  !! allocation-free.
  subroutine precompCubicSplineBrackets(xx, dxs, this)

    !> abscissas (shared across all splines)
    real(dp), intent(in) :: xx(:)

    !> evaluation points (shared across all splines)
    real(dp), intent(in) :: dxs(:)

    !> precomputed bracket indices and evaluation weights
    type(TCubicSplineEval), intent(inout) :: this

    !! bisection search bounds and midpoint
    integer :: lo, hi, mid

    !! evaluation-point iterator and number of evaluation points
    integer :: ii, npts

    !! difference of bracketing abscissas and local weight values
    real(dp) :: step, a, b

    npts = size(dxs)
    if (.not. allocated(this%left) .or. size(this%left) /= npts) then
      if (allocated(this%left)) then
        deallocate(this%left, this%right, this%aa, this%bb, this%wl, this%wr)
      end if
      allocate(this%left(npts), this%right(npts), this%aa(npts), this%bb(npts), this%wl(npts),&
          & this%wr(npts))
    end if

    do ii = 1, npts
      ! bisection search (identical to get_cubic_spline)
      lo = 1
      hi = size(xx)
      do
        if ((hi - lo) <= 1) exit
        mid = (hi + lo) / 2
        if (dxs(ii) >= xx(mid)) then
          lo = mid
        else
          hi = mid
        end if
      end do
      step = xx(hi) - xx(lo)
      a = (xx(hi) - dxs(ii)) / step
      b = (dxs(ii) - xx(lo)) / step
      this%left(ii) = lo
      this%right(ii) = hi
      this%aa(ii) = a
      this%bb(ii) = b
      this%wl(ii) = step**2 / 6.0_dp * (a**2 - 1.0_dp) * a
      this%wr(ii) = step**2 / 6.0_dp * (b**2 - 1.0_dp) * b
    end do

  end subroutine precompCubicSplineBrackets


  !> Finds the 2nd derivatives of natural splines.
  !!
  !! \note Works only for equidistant points!
  subroutine get2ndNaturalSplineDerivs(xx, fct, gama)

    !> abscissas
    real(dp), intent(in) :: xx(:)

    !> ordinates
    real(dp), intent(in) :: fct(:)

    !> array with 2nd derivatives
    real(dp), intent(out), allocatable :: gama(:)

    !! iterates over spline nodes/equations
    integer :: ii

    !! number of tabulated points, i.e. abscissas
    integer :: nn

    !! error status
    integer :: info

    !! tridiagonal equations
    real(dp), allocatable :: beta(:), alpha(:)

    ! recurring pre-factors
    real(dp) :: fac, step

    !! auxiliary variables for lapack solver
    integer, allocatable :: IWORK(:), IPIV(:)
    real(dp) :: RCOND, FERR(1), BERR(1)
    real(dp), allocatable :: DLF(:), DF(:), DUF(:), DU2(:), WORK(:), SOL(:)

    nn = size(xx)
    allocate(alpha(nn - 1))
    allocate(beta(nn))
    allocate(gama(nn))

    ! pre-factors
    step = 1.0_dp / real(nn + 1, dp)
    step = step**2 / 6.0_dp
    fac = step * 4.0_dp

    ! fill the matrix
    gama(1) = 0.0_dp
    gama(nn) = 0.0_dp
    beta(1) = step
    beta(nn) = step
    alpha(1) = step

    do ii = 2, nn - 1
      gama(ii) = fct(ii + 1) - 2.0_dp * fct(ii) + fct(ii - 1)
      beta(ii) = fac
      alpha(ii) = step
    end do

    ! set up lapack variables
    allocate(DLF(nn - 1))
    allocate(DUF(nn - 1))
    allocate(DF(nn))
    allocate(DU2(nn - 2))
    allocate(WORK(3 * nn))
    allocate(IWORK(nn))
    allocate(IPIV(nn))
    allocate(SOL(nn))

    ! solve the equation
    call DGTSVX('N', 'N', nn, 1, alpha, beta, alpha, DLF, DF, DUF, DU2, IPIV, gama,&
        & nn, SOL, nn, RCOND, FERR, BERR, WORK, IWORK, INFO)

    gama = SOL

  end subroutine get2ndNaturalSplineDerivs


  !> Find second derivatives of the splines at tabulated function points.
  !!
  !! \note Algorithm is based on the Numerical recipes.
  !! \details Recursion relation: u(i) = gama(i) + beta(i)u(i+1)
  subroutine set_cubic_spline(xx, fct, ds1, dsn, gama, tNatural)

    !> abscissas
    real(dp), intent(in) :: xx(:)

    !> ordinates
    real(dp), intent(in) :: fct(:)

    !> first derivative of the interpolating function at point 1
    real(dp), intent(in) :: ds1

    !> first derivative of the interpolating function at point N
    real(dp), intent(in) :: dsn

    !> 2nd derivatives of the interpolating function at the tabulated points
    real(dp), intent(out), allocatable :: gama(:)

    !> true, if natural splines are desired (2nd derivative equal to zero at boundary),
    !! otherwise lower/upper boundary condition is set to have a specified first derivative (ds1)
    logical, intent(in) :: tNatural

    !! iterates over spline nodes/equations
    integer :: ii

    !! number of tabulated points, i.e. abscissas
    integer :: nn

    !!
    real(dp) :: mu_i, xi_i

    !!
    real(dp), allocatable :: beta(:)

    allocate(beta(size(xx)))
    allocate(gama(size(xx)))

    nn = size(xx)

    ! first equation with boundary condition ds1
    if (tNatural) then
       beta(1) = 0.0_dp
       gama(1) = 0.0_dp
    else
       beta(1) = - 0.5_dp
       gama(1) = 3.0_dp / (xx(2) - xx(1)) * (fct(2) - fct(1) / (xx(2) - xx(1)) - ds1)
     end if

    ! middle equations
    do ii = 2, nn - 1
      mu_i = (xx(ii) - xx(ii - 1)) / (xx(ii + 1) - xx(ii - 1))
      xi_i = (2.0_dp + mu_i * beta(ii - 1))
      beta(ii) = (mu_i - 1.0_dp) / xi_i
      gama(ii) = ((6.0_dp / (xx(ii + 1) - xx(ii - 1))) * ((fct(ii + 1) - fct(ii))&
          & / (xx(ii + 1) - xx(ii)) - (fct(ii) - fct(ii - 1)) / (xx(ii) - xx(ii - 1)))&
          & - mu_i * gama(ii - 1) ) / xi_i
    end do

    ! last equation with boundary value dsn
    if (tNatural) then
      gama(nn) = 0.0_dp
    else
      gama(nn) = (6.0_dp / (xx(nn) - xx(nn - 1)) * (dsn - (fct(nn) - fct(nn - 1))&
          & / (xx(nn) - xx(nn - 1)) - gama(nn - 1)) / (2.0_dp + beta(nn - 1)) )
    end if

    ! backsubstitution
    do ii = nn - 1, 1, -1
      gama(ii) = gama(ii) + beta(ii) * gama(ii + 1)
    end do

  end subroutine set_cubic_spline

end module common_interpolation
