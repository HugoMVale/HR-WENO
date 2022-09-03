program example_hrweno
!>----------------------------------------------------------------------------------------------
!> This program illustrates the application of modules 'weno' and 'tvdode' for solving a 1D
!> hyperbolic equation (Burger's equation):
!>                    d/dt u(x,t) = - d/dx f(u(x,t))
!>                    with f(u,t) = (u^2)/2
!>
!> The spatial variable 'x' is discretized according to a finite-volume approach and the time
!> variable 't' is left continuous (method of lines), leading to:
!>                    du(i,t)/dt = -1/dx(i)*(f(i+1/2,t) - f(i-1/2,t))
!>
!>
!>                              ul=u(i-1/2)^+            ur=u(i+1/2)^-
!>               --|-----(i-1)------|---------(i)----------|------(i+1)-----|--
!>                               x(i-1/2)               x(i+1/2)
!>                                  |<---    dx(i)     --->|
!>----------------------------------------------------------------------------------------------
    use tvdode, only : rktvd123
    use weno, only : weno35
    use iso_fortran_env, only : real64, error_unit
    implicit none

    integer, parameter :: rk = real64
    integer, parameter :: nc = 20
    real(rk) :: xedges(0:nc), x(nc), dx(nc), u(nc)
    real(rk) :: dt, time, time_out, time_start, time_end, xmin, xmax
    logical :: first = .true.
    integer :: fevals, num_time_points, order, istate, itask, ii

    !> Define the spatial grid
    !> This example uses a linear grid, but you can try a geometric grid as well
    xmin = -1_rk
    xmax = 1_rk
    do ii = 0, nc
      xedges(ii) = xmin + (xmax - xmin)*ii/nc
    end do

    !> Cell width and center (general case)
    dx = xedges(1:nc) - xedges(0:nc-1)
    x = xedges(0:nc-1) + dx/2

    !> Initial condition u(x,t=0)
    u = 0_rk
    u(:nc/2) = 1_rk

    !> Open file where results will be stored
    call output(1)

    !> Call ODE time solver
    time_start = 0_rk
    time_end = 4_rk
    dt = 1e-2_rk
    time = time_start
    num_time_points = 50
    order = 2
    istate = 1
    itask = 1
    fevals = 0
    
    do ii = 0, num_time_points
        time_out = time_end*ii/num_time_points
        call rktvd123(fu, u, time, time_out, dt, order, itask, istate)
        call output(2)
    end do

    !> End of simulation
    call output(3)

    contains

    pure subroutine fu(t, v, vdot)
    !------------------------------------------------------------------------------------------
    ! This
    !------------------------------------------------------------------------------------------
    real(rk), intent(in) :: t, v(:)
    real(rk), intent(out) :: vdot(:)
    real(rk), parameter :: eps = 1e-6_rk
    integer :: i
    !//////////////////////////////////////////////////////////////////////////////////////////

    !> Allocate u arrays, including gost cells for u
    !allocate(u(1-(k-1):nc+(k-1)), ul(nc), ur(nc))
    
    !WHAT is u and upoint??

    ! if (first) then

    !     first = .FALSE.

    !     !CALL CALC_CIRJ(k,M,x(2:M+1),c) !bulshit

    !     !u(-1:0) = 0.0d0
    !     !u(nmax+1:nmax+2) = 0.0d0

    !     !H(0)  = 0.0d0

    ! end if

        do i = 1, nc
            vdot(i) = (0 - 0)/dx(i)
        end do

    end subroutine fu
    !##########################################################################################


    elemental function flux(v, t)
    !>------------------------------------------------------------------------------------------
    !> Flux function for Burger's equation.
    !>
    !> ARGUMENTS:
    !> v     function v(x,t)
    !> t     variable t
    !>------------------------------------------------------------------------------------------
    real(rk) :: flux
    real(rk), intent (in) :: v
    real(rk), intent (in), optional :: t
    
        flux = 0.5_rk*v*v
    
    end function flux
    !##########################################################################################


    subroutine output(message)
    !>------------------------------------------------------------------------------------------
    !> Auxiliary routine to save results to file.
    !>
    !> ARGUMENTS:
    !> message: parameter to select output action
    !------------------------------------------------------------------------------------------
    integer, intent (in) :: message
    integer :: i
    
        select case (message)

        !> Open files and write headers and grid
        case (1)

            print *, "Running test..."
            print *, "Start: ", fdate()

            !> Write grid
            open (unit=1, file="./output/xgrid.txt", status="replace", action="write", &
                  position="rewind")

            write (1,'(1x, a5, 2(a15))') "i", "x(i)", "dx(i)"
            do i = 1, nc 
                write (1,'(1x, i5, 2(e15.5))') i, x(i), dx(i)
            end do

            !> Write header u
            open (unit=2, file="./output/u.txt", status="replace", action="write", &
                  position="rewind")

            write (2,'(1x, a15)', advance="no") "t"
            do i = 1, nc
                write (2,'(a15)', advance="no") "u("//itoa(i)//")"
            end do
            write (2,*) ""

        !> Write values
        case (2)
            write (2,'(1x, e15.5)', advance="no") time
            do i = 1, nc
                write (2,'(e15.5)', advance="no") u(i)
            end do
            write (2,*) ""

        !> Close files
        case (3)
            close (1)
            close (2)
            print *, "End  : ", fdate()
            print '(a13, i5)', "RHS Evals.: ", fevals
            print *

        end select

    end subroutine output
    !##########################################################################################


    function itoa(i) result(res)
    !>------------------------------------------------------------------------------------------
    !> Convert integer to string.
    !>
    !> ARGUMENTS:
    !> i:    integer
    !>------------------------------------------------------------------------------------------
      character(:), allocatable :: res
      integer, intent(in) :: i
      character(range(i)+2) :: tmp
      write(tmp,'(i0)') i
      res = trim(tmp)
    end function

end program
