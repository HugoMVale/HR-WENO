program example_1d
!>----------------------------------------------------------------------------------------------
!> This program illustrates the application of modules 'hrschemes' and 'tvdode' for solving a 1D
!> hyperbolic equation (Burger's equation):
!>                    d/dt u(x,t) = - d/dx f(u(x,t))
!>                    with f(u,t) = (u^2)/2
!>
!> The spatial variable 'x' is discretized according to a finite-volume approach and the time
!> variable 't' is left continuous (method of lines), leading to:
!>                    du(i,t)/dt = -1/dx(i)*(f(i+1/2,t) - f(i-1/2,t))
!>
!>
!>                              ul=u(i-1/2)^+          ur=u(i+1/2)^-
!>               --|-----(i-1)------|---------(i)----------|------(i+1)-----|--
!>                               x(i-1/2)               x(i+1/2)
!>                                  |<---    dx(i)     --->|
!>----------------------------------------------------------------------------------------------
    use tvdode, only : rktvd
    use hrschemes, only : weno, lax_friedrichs
    use iso_fortran_env, only : real64, error_unit
    implicit none

    integer, parameter :: rk = real64
    integer, parameter :: nc = 100
    real(rk) :: xedges(0:nc), x(nc), dx(nc), u(nc)
    real(rk) :: dt, time, time_out, time_start, time_end, xmin, xmax
    integer :: num_time_points, order, istate, itask, ii

    !> Define the spatial grid
    !> In this example, we use a linear grid, but any smooth grid can be used
    xmin = -5_rk
    xmax = 5_rk
    do ii = 0, nc
      xedges(ii) = xmin + (xmax - xmin)*ii/nc
    end do

    !> Compute cell width and center (general formulas for any grid)
    dx = xedges(1:nc) - xedges(0:nc-1)
    x = xedges(0:nc-1) + dx/2

    !> Initial condition u(x,t=0)
    u = ic(x)

    !> Open file where results will be stored
    call output(1)

    !> Call ODE time solver
    time_start = 0_rk
    time_end = 4._rk
    dt = 1e-2_rk
    time = time_start
    num_time_points = 100
    order = 3
    istate = 1
    itask = 1
    
    do ii = 0, num_time_points
        time_out = time_end*ii/num_time_points
        call rktvd(fu, u, time, time_out, dt, order, itask, istate)
        call output(2)
    end do

    !> End of simulation
    call output(3)

    contains

    pure subroutine fu(t, v, vdot)
    !>------------------------------------------------------------------------------------------
    !> This subroutine computes the right hand side of:
    !>
    !>                    du(i,t)/dt = -1/dx(i)*(f(i+1/2,t) - f(i-1/2,t))
    !>
    !> INTERNAL VARIABLES:
    !> vl      vector(1:nc) with reconstructed value at left boundary of cell i (v_{i-1/2}^+)
    !> vr      vector(1:nc) with reconstructed value at right boundary of cell i (v_{i+1/2}^-)
    !>------------------------------------------------------------------------------------------
    real(rk), intent(in) :: t, v(:)
    real(rk), intent(out) :: vdot(:)
    
    integer, parameter :: k = 3
    real(rk) :: vl(nc), vr(nc), fedges(0:nc), vext(1-(k-1):nc+(k-1))
    real(rk), parameter :: eps = 1e-6_rk,  alpha = 1._rk
    integer :: i
   
        !> Populate extended vector with ghost cells
        !> There is probably a smarter way to do this
        vext(1:nc) = v
        vext(:0) = v(1)
        vext(nc+1:) = v(nc)

        !> Get reconstructed values at cell boundaries
        call weno(k, vext, vl, vr, eps)

        !> Fluxes at interior cell boundaries
        !> We use a Lax-Friedrichs method as recomended by Shu (other options possible)
        do i = 1, nc-1
            fedges(i) = lax_friedrichs(flux, vr(i), vl(i+1), t, alpha)
        end do

        !> Apply problem-specific flux constraints at domain boundaries
        fedges(0) = fedges(1)
        fedges(nc) = fedges(nc-1)

        !> Evaluate du/dt
        vdot = - (fedges(1:) - fedges(:nc-1))/dx

    end subroutine fu
    !##########################################################################################


    pure real(rk) function flux(v, t)
    !>------------------------------------------------------------------------------------------
    !> Flux function. Here we define the flux corresponding to Burger's equation.
    !>
    !> ARGUMENTS:
    !> v     function v(x,t)
    !> t     variable t
    !>------------------------------------------------------------------------------------------
    real(rk), intent (in) :: v, t
    
        flux = (v**2)/2
    
    end function flux
    !##########################################################################################


    elemental real(rk) function ic(z)
    !>------------------------------------------------------------------------------------------
    !> Initial condition. Here we used a limited linear profile.
    !>
    !> ARGUMENTS:
    !> z     spatial variable 
    !>------------------------------------------------------------------------------------------
    real(rk), intent (in) :: z
    real(rk), parameter :: z1=-4._rk, z2=2._rk, v1=1._rk, v2=-0.5_rk   
       
        ic = v1 + (v2 - v1)/(z2 - z1)*(z - z1)
        ic = max(min(ic, v1), v2)
    
    end function ic
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

            print *, "Running example1d..."
            print *, "Start: ", fdate()

            !> Write grid
            open (unit=1, file="./output/xgrid.txt", status="replace", action="write", &
                  position="rewind")

            write (1,'(a5, 2(1x, a15))') "i", "x(i)", "dx(i)"
            do i = 1, nc 
                write (1,'(i5, 2(1x, es15.5))') i, x(i), dx(i)
            end do

            !> Write header u
            open (unit=2, file="./output/u.txt", status="replace", action="write", &
                  position="rewind")

            write (2,'(a16)', advance="no") "t"
            do i = 1, nc
                write (2,'(1x, a16)', advance="no") "u("//itoa(i)//")"
            end do
            write (2,*) ""

        !> Write values
        case (2)
            write (2,'(es16.5e3)', advance="no") time
            do i = 1, nc
                write (2,'(1x, es16.5e3)', advance="no") u(i)
            end do
            write (2,*) ""

        !> Close files
        case (3)
            close (1)
            close (2)
            print *, "End  : ", fdate()
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
