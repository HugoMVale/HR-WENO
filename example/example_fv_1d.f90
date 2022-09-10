program example_fv_1d
!----------------------------------------------------------------------------------------------
!!   This program illustrates the application of modules 'weno' and 'tvdode' for solving a 1D
!! hyperbolic equation (Burger's equation):
!!```    
!!                    d/dt u(x,t) = - d/dx f(u(x,t))
!!                    with f(u,t) = (u**2)/2
!!                    and u(x,t=0) = ic(x)
!!```
!!   The spatial variable 'x' is discretized according to a finite-volume approach and the time
!! variable 't' is left continuous (method of lines), leading to:
!!```    
!!                du(i,t)/dt = -1/dx(i)*( f(u(x(i+1/2),t)) - f(u(x(i-1/2),t)) )
!!
!!                              ul=u(i-1/2)^+          ur=u(i+1/2)^-
!!               --|-----(i-1)------|---------(i)----------|------(i+1)-----|--
!!                               x(i-1/2)               x(i+1/2)
!!                                  |<---    dx(i)     --->|
!!```
!!  In this particular example, we use the 3rd order 'rktvd' ode solver (we could equally well
!! employ the 'mstvd' solver). The reconstruction is done with the 5th order WENO scheme; to 
!! try other orders, we can change the parameter 'k' in procedure 'rhs' .      
!----------------------------------------------------------------------------------------------
    use tvdode, only : rktvd
    use weno, only : wenok
    use hrutils, only : godunov, lax_friedrichs, tgrid1, grid1
    use iso_fortran_env, only : real64, error_unit
    implicit none

    integer, parameter :: rk = real64
    integer, parameter :: nc = 100
    real(rk) :: u(nc)
    real(rk) :: dt, time, time_out, time_start, time_end, xmin, xmax
    integer :: num_time_points, order, istate, itask, ii
    type(tgrid1) :: gx

    ! Define the spatial grid
    ! In this example, we use a linear grid, but any smooth grid can be used
    xmin = -5._rk
    xmax = 5._rk
    gx = grid1(xmin, xmax, nc)
    
    ! Initial condition u(x,t=0)
    u = ic(gx%c)

    ! Open file where results will be stored
    call output(1)

    ! Call ODE time solver
    time_start = 0._rk
    time_end = 12._rk
    dt = 1e-2_rk
    time = time_start
    num_time_points = 100
    order = 3
    istate = 1
    itask = 1
    
    do ii = 0, num_time_points
        time_out = time_end*ii/num_time_points
        call rktvd(rhs, u, time, time_out, dt, order, itask, istate)
        call output(2)
    end do

    ! End of simulation
    call output(3)

    contains

    pure subroutine rhs(t, v, vdot)
    !------------------------------------------------------------------------------------------
    !! This subroutine computes the *numerical approximation* to the right hand side of:
    !!```
    !!                du(i,t)/dt = -1/dx(i)*( f(u(x(i+1/2),t)) - f(u(x(i-1/2),t)) )
    !!```
    !!   There are two main steps. First, we use the WENO scheme to obtain the reconstructed 
    !! values of 'u' at the left and right cell boundaries. Note that, in general, because of
    !! discontinuities, \( u_{i+1/2}^+ \neq u_{(i+1)+1/2}^- \). Second, we use a suitable flux
    !! method (e.g., Godunov, Lax-Friedrichs) to compute the flux from the reconstructed
    !! 'u' values.
    !------------------------------------------------------------------------------------------
    real(rk), intent(in) :: t
        !! time variable
    real(rk), intent(in) :: v(:)
        !! vector(N) with v(z,t) values
    real(rk), intent(out) :: vdot(:)
        !! vector(N) with v'(z,t) values
    
    integer, parameter :: k = 3
    real(rk) :: vl(nc), vr(nc), fedges(0:nc), vext(1-(k-1):nc+(k-1))
    real(rk), parameter :: eps = 1e-6_rk,  alpha = 1._rk
    integer :: i
   
        ! Populate extended 'v' vector with ghost cells
        vext(1:nc) = v
        vext(:0) = v(1)
        vext(nc+1:) = v(nc)

        ! Get reconstructed values at cell boundaries
        call wenok(k, vext, vl, vr, eps)

        ! Fluxes at interior cell boundaries
        ! One can use the Lax-Friedrichs or the Godunov method
        do concurrent (i=1:nc-1)
            !fedges(i) = lax_friedrichs(flux, vr(i), vl(i+1), gx%r(i), t, alpha)
            fedges(i) = godunov(flux, vr(i), vl(i+1), gx%r(i), t)
        end do

        ! Apply problem-specific flux constraints at domain boundaries
        fedges(0) = fedges(1)
        fedges(nc) = fedges(nc-1)

        ! Evaluate du/dt
        vdot = - (fedges(1:) - fedges(:nc-1))/gx%d

    end subroutine rhs


    pure real(rk) function flux(v, z, t)
    !------------------------------------------------------------------------------------------
    !! Flux function. Here we define the flux corresponding to Burger's equation.
    !------------------------------------------------------------------------------------------
    real(rk), intent (in) :: v
        !! function v(z,t)
    real(rk), intent (in) :: z
        !! spatial variable
    real(rk), intent (in) :: t
        !! time variable
    
        flux = (v**2)/2
    
    end function flux


    elemental real(rk) function ic(z)
    !------------------------------------------------------------------------------------------
    !! Initial condition. Here we used a limited linear profile.
    !------------------------------------------------------------------------------------------
    real(rk), intent (in) :: z
        !! spatial variable
    real(rk), parameter :: z1=-4._rk, z2=2._rk, v1=1._rk, v2=-0.5_rk   
       
        ic = v1 + (v2 - v1)/(z2 - z1)*(z - z1)
        ic = max(min(ic, v1), v2)
    
    end function ic


    subroutine output(message)
    !------------------------------------------------------------------------------------------
    !! Auxiliary routine to save results to file.
    !------------------------------------------------------------------------------------------
    integer, intent (in) :: message
        !! parameter to select output action 
    integer :: i
    
        select case (message)

        ! Open files and write headers and grid
        case (1)

            print *, "Running example1d..."
            print *, "Start: ", fdate()

            ! Write grid
            open (unit=1, file="./output/xgrid.txt", status="replace", action="write", &
                  position="rewind")

            write (1,'(a5, 2(1x, a15))') "i", "x(i)", "dx(i)"
            do i = 1, nc 
                write (1,'(i5, 2(1x, es15.5))') i, gx%c(i), gx%d(i)
            end do

            ! Write header u
            open (unit=2, file="./output/u.txt", status="replace", action="write", &
                  position="rewind")

            write (2,'(a16)', advance="no") "t"
            do i = 1, nc
                write (2,'(1x, a16)', advance="no") "u("//itoa(i)//")"
            end do
            write (2,*) ""

        ! Write values
        case (2)
            write (2,'(es16.5e3)', advance="no") time
            do i = 1, nc
                write (2,'(1x, es16.5e3)', advance="no") u(i)
            end do
            write (2,*) ""

        ! Close files
        case (3)
            close (1)
            close (2)
            print *, "End  : ", fdate()
            print *

        end select

    end subroutine output


    function itoa(i) result(res)
    !------------------------------------------------------------------------------------------
    !! Convert integer to string.
    !------------------------------------------------------------------------------------------
    integer, intent(in) :: i
        !! Integer
    character(:), allocatable :: res
    character(range(i)+2) :: tmp
      write(tmp,'(i0)') i
      res = trim(tmp)
    end function

end program
