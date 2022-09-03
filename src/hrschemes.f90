module hrschemes
!>---------------------------------------------------------------------------------------------
!> This module contains a collection of high-resolution weighted essentially non-oscillatory
!> (WENO) schemes for uniform and *arbitrary* finite volume grids.
!> Source: ICASE 97-65 by Shu, 1997.
!>---------------------------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    public :: weno, calc_c, lax_friedrichs, c1, c2, c3

    integer, parameter :: rk = real64

    !> Parameter arrays for ENO and WENO methods
    real(rk), parameter :: d1(0:0) = 1._rk, &
                           d2(0:1) = [2._rk/3_rk, 1._rk/3_rk], &
                           d3(0:2) = [0.3_rk, 0.6_rk, 0.1_rk]
    real(rk), parameter :: &
        c1(0:0,-1:0) = reshape([1._rk, 1._rk], [1,2], order=[1,2]), &
        c2(0:1,-1:1) = reshape ([3._rk/2, -1._rk/2, 1._rk/2, 1._rk/2, -1._rk/2, 3._rk/2 ], &
                       [2,3], order=[1,2]), &
        c3(0:2,-1:2) = reshape ([11._rk/6, -7._rk/6, 1._rk/3, 1._rk/3, 5._rk/6, -1._rk/6, &
                       -1._rk/6, 5._rk/6, 1._rk/3, 1._rk/3, -7._rk/6, 11._rk/6], &
                       [3,4], order=[1,2])

    abstract interface
        pure function flux(u, t)
            import :: rk
            real(rk) :: flux
            real(rk), intent(in) :: u, t
        end function
    end interface

    contains

    subroutine weno(k, v, vl, vr, eps, c)
    !>-----------------------------------------------------------------------------------------
    !> This subroutine implements the (2k-1)th order WENO method for *arbitrary* finite volume
    !> grids described in ICASE 97-65 (Shu, 1997).
    !>
    !>          |---1---| ...  |---(i-1)---|-------(i)-------|---(i+1)---| .... |--nc--|
    !>                                      ^               ^
    !>                                      vl(i)           vr(i)
    !>                                    v_{i-1/2}^+     v_{i+1/2}^-
    !>
    !>
    !> ARGUMENTS:
    !> k             order of reconstruction within the cell (k = 1, 2 or 3)
    !> v             vector(1-(k-1):nc+(k-1)) with *average* cell values, including (k-1) ghost
    !>               cells on each side
    !> vl            vector(1:nc) with reconstructed value at left boundary of cell i (v_{i-1/2}^+)
    !> vr            vector(1:nc) with reconstructed value at right boundary of cell i (v_{i+1/2}^-)
    !> eps           numerical smoothing factor
    !> c(j,r,i)      optional array(0:k-1,-1:k-1,1:nc) of constants for a *non-uniform* grid
    !>               (see calc_c)
    !>
    !> INTERNAL VARIABLES:
    !> nc            number of cells
    !>-----------------------------------------------------------------------------------------
    integer, intent (in) :: k
    real(rk), intent(in) :: eps, v(2-k:)
    real(rk), intent(out) :: vl(:), vr(:)
    real(rk), intent(in), target, optional :: c(0:,-1:,:)

    real(rk), dimension(0:k-1) :: vlr, vrr, w, wtilde, alfa, alfatilde, beta
    real(rk), allocatable :: d(:), ci(:,:)
    integer :: i, r, nc
    logical :: usrgrid
    character(:), allocatable :: msg

        !> Select constant parameters according to order of the method
        select case (k)
            case(1)
                d = d1
                ci = c1
            case(2)
                d = d2
                ci = c2
            case(3)
                d = d3
                ci = c3
            case default
                msg = "Invalid input 'k' in 'weno35'. Valid range: 1 <= k <= 3."
                error stop msg
        end select

        !> Check if user supplied grid
        if (present(c)) then
            usrgrid = .true.
        else
            usrgrid = .false.
        end if

        !> Algorithm
        !> Obtain the 'k' reconstructed values vi+1/2(r) & vi-1/2(r)
        !> Todo: change to 'do concurrent'
        nc = size(vl)
        do i = 1, nc

            !> Equations 2.10, 2.51
            if (usrgrid) ci = c(:,:,i)
            do r = 0, k-1
                vrr(r) = sum(ci(:,r)*v(i-r:i-r+k-1))
                vlr(r) = sum(ci(:,r-1)*v(i-r:i-r+k-1))
            end do

            select case(k)
                
                !> Not sure about this case
                case(1)
                    beta(0) = 0._rk

                !> Equation 2.62    
                case(2)

                    beta(0) = (v(i+1) - v(i))**2
                    beta(1) = (v(i) - v(i-1))**2

                !> Equation 2.63
                case(3)

                    beta(0) = 13._rk/12*(v(i) - 2*v(i+1) + v(i+2))**2  &
                            + 1._rk/4*(3*v(i) - 4*v(i+1) + v(i+2))**2

                    beta(1) = 13._rk/12*(v(i-1) - 2*v(i) + v(i+1))**2  &
                            + 1._rk/4*(v(i-1) - v(i+1))**2

                    beta(2) = 13._rk/12*(v(i-2) - 2*v(i-1) + v(i))**2  &
                            + 1._rk/4*(v(i-2) - 4*v(i-1) + 3*v(i))**2

            end select

            !> Equations 2.58-2.59 and procedure 2.2-4
            alfa = d/(eps + beta)**2
            alfatilde = d(k-1:0:-1)/(eps + beta)**2
            w = alfa/sum(alfa)
            wtilde = alfatilde/sum(alfatilde)

            !> Procedure 2.2-5
            vr(i) = sum(w*vrr)
            vl(i) = sum(wtilde*vlr)

        end do


    end subroutine weno
    !>#########################################################################################


    pure subroutine calc_c(k, xedges, c)
    !>-----------------------------------------------------------------------------------------
    !> This subroutine computes the array of constants 'c(j,r,i)' required to use weno35 with
    !> arbitrary (i.e., non-uniform) grids.
    !> Source: ICASE 97-65 by Shu, 1997.
    !>
    !> ARGUMENTS:
    !> k             order (>=1) of reconstruction within the cell
    !> xedges(i)     vector(0:nc) of cell edges
    !                xedges(i) value of x at right boundary of cell i (x_{i+1/2})
    !                xedges(i-1) value of x at left boundary of cell i (x_{i-1/2})
    !> c(j,r,i)      array(0:k-1,-1:k-1,1:nc) of constants for a non-uniform grid
    !>-----------------------------------------------------------------------------------------
    integer, intent (in) :: k
    real(rk), intent(in) :: xedges(0:)
    real(rk), intent(out) :: c(0:,-1:,:)

    real(rk) :: prod1, prod2, sum1, sum2, dx
    real(rk), allocatable, target :: xext(:)
    real(rk), dimension(:), pointer :: xl, xr
    integer :: i, j, l, nc, ng, m, q, r
    character(:), allocatable :: msg

        !> Check input conditions
        if (k < 1) then
            msg = "Invalid input 'k' in 'calc_c'. Valid range: k >= 1."
            error stop msg
        end if

        if (k > 3) then
            msg = "Input 'k' is probably too high in 'calc_c'. Usual range: k <= 3."
            error stop msg
        end if

        if (ubound(c,1) /= k-1) then
            msg = "Invalid ubound(c,1) in 'calc_c'."
            error stop msg
        end if

        if (ubound(c,2) /= k-1) then
            msg = "Invalid ubound(cgrid,2) in 'calc_c'."
            error stop msg
        end if

        nc = ubound(xedges,1)
        if (ubound(c,3) /= nc) then
            msg = "Invalid ubound(cgrid,3) in 'calc_c'."
            error stop msg
        end if

        !> Allocate extended grid with (k+1) ghost cells on each side
        ng = k+1
        allocate(xext(0-ng:nc+ng))
        xext(0:nc) = xedges

        !> Extend to the left linearly
        dx = xext(1) - xext(0)
        do i = -1, lbound(xext,1), -1
            xext(i) = xext(i+1) - dx
        end do

        !> Extend to the right linearly
        dx = xext(nc) - xext(nc-1)
        do i = nc+1, ubound(xext,1)
            xext(i) = xext(i-1) + dx
        end do

        !> Get pointers to left and right cell boundaries
        xl(1-ng:) => xext(lbound(xext,1):ubound(xext,1)-1)
        xr(1-ng:) => xext(lbound(xext,1)+1:ubound(xext,1))

        !> Compute array of constants 'c' for each grid position
        !> Equation 2.20, page 6.
        do i = 1, nc

            do r = -1, k-1

                do j = 0, k-1

                    sum2 = 0

                    do m = j+1, k

                        prod2 = 1
                        do l = 0, k
                            if (l == m) cycle
                            prod2 = prod2*(xl(i-r+m) - xl(i-r+l))
                        end do

                        sum1 = 0
                        do l = 0, k

                            if (l == m) cycle

                                prod1 = 1
                                do q = 0, k
                                    if (q == m .or. q == l) cycle
                                    prod1 = prod1*(xr(i) - xl(i-r+q))
                                end do

                            sum1 = sum1 + prod1

                        end do

                        sum2 = sum2 + sum1/prod2

                    end do

                    c(j,r,i) = sum2*(xr(i-r+j) - xl(i-r+j))

                end do

            end do

        end do

    end subroutine calc_c
    !>#########################################################################################


    pure real(rk) function lax_friedrichs(f, vm, vp, t, alpha)
    !>-----------------------------------------------------------------------------------------
    !> Monotone Lax-Friedrichs flux.
    !> Equation 2.72, page 21.
    !>
    !> ARGUMENTS:
    !> f      flux function f(u)
    !> vm     left  side of reconstructed of v (v_{i^+1/2}^-)
    !> vp     right side of reconstructed of v (v_{i^+1/2}^+)
    !> t      time
    !> alpha  max(abs(f'(v)))
    !>-----------------------------------------------------------------------------------------
    procedure(flux) :: f
    real(rk), intent (in) :: vm, vp, t, alpha

        lax_friedrichs = (f(vm, t) + f(vp, t) - alpha*(vp - vm))/2

    end function lax_friedrichs
    !>#########################################################################################

end module hrschemes
