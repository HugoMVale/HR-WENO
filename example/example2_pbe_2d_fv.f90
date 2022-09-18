program example_pbe_2d_fv
!!   This program illustrates the application of modules 'weno' and 'mstvd' for solving a 2D
!! population balance equation (PBE):
!!```
!!                    d/dt u(x,t) + d/dx1 f1(u(x,t)) + d/dx2 f2(u(x,t)) = 0
!!                    with f1(u,t) = u
!!                    with f2(u,t) = u
!!                    and u(x,t=0) = ic(x)
!!```
!!   The internal coordinates 'x=(x1,x2)' are discretized according to a finite-volume approach
!! and the time variable 't' is left continuous (method of lines). See example1 for notation.
!!   In this particular example, we use the 3rd order 'mstvd' ode solver. The reconstruction
!! is done with the 5th order WENO scheme; to try other orders, we can change the parameter 'k'
!! in procedure 'rhs' .
   use tvdode, only: mstvd
   use weno, only: wenok
   use hrutils, only: godunov, tgrid1, grid1
   use iso_fortran_env, only: real64, stderr => error_unit, stdout => output_unit
   use stdlib_strings, only: to_string
   implicit none

   integer, parameter :: rk = real64
   integer, parameter :: nc(2) = [200, 200]
   real(rk) :: u(product(nc))
   real(rk), dimension(product(nc), 4) :: uold, udotold
   type(tgrid1) :: gx(2)
   real(rk) :: dt, time, time_out, time_start, time_end, xmin, xmax
   integer :: num_time_points, istate, ii, jj

   ! Define grids for x1 and x2
   ! In this example, we use linear grids, but any smooth grid can be used
   gx(1) = grid1(0._rk, 10._rk, nc(1))
   gx(2) = grid1(0._rk, 10._rk, nc(2))

   ! Initial condition u(x,t=0)
   do concurrent(ii=1:nc(1), jj=1:nc(2))
      u((jj - 1)*nc(1) + ii) = ic([gx(1)%center(ii), gx(2)%center(jj)])
   end do

   ! Open file where results will be stored
   call output(1)

   ! Call ODE time solver
   time_start = 0._rk
   time_end = 5._rk
   dt = 5e-3_rk

   time = time_start
   num_time_points = 100
   istate = 1
   do ii = 0, num_time_points
      time_out = time_end*ii/num_time_points
      call mstvd(rhs, u, time, time_out, dt, uold, udotold, istate)
      call output(2)
   end do

   ! End of simulation
   call output(3)

contains

   pure subroutine rhs(t, v, vdot)
   !! This subroutine computes the *numerical approximation* to the right hand side of:
   !!```
   !!     du(i,j,t)/dt = -1/dx1(i)*( f1(u(x1(i+1/2),x2(j),t)) - f1(u(x1(i-1/2),x2(j),t)) )
   !!                    -1/dx2(j)*( f2(u(x1(i),x2(j+1/2),t)) - f2(u(x1(i),x2(j-1/2),t)) )
   !!```
   !!   There are two main steps. First, we use the WENO scheme to obtain the reconstructed
   !! values of 'u' at the left and right cell boundaries along each dimension. Second, we use
   !! a suitable flux method (e.g., Godunov, Lax-Friedrichs) to compute the flux from the
   !! reconstructed 'u' values.
      real(rk), intent(in) :: t
        !! time variable
      real(rk), intent(in) :: v(:)
        !! vector(N) with v(z,t) values
      real(rk), intent(out) :: vdot(:)
        !! vector(N) with v'(z,t) values

      integer, parameter :: k = 3
      real(rk), allocatable :: vl(:), vr(:), vext(:)
      real(rk) :: varray(nc(1), nc(2)), fedges(0:nc(1), 0:nc(1), 2)
      real(rk), parameter :: eps = 1e-6_rk
      integer :: i, j

      ! Reshape v to array
      varray = reshape(v, [shape(varray)])

      ! Fluxes along x1 at interior cell boundaries
      fedges = 0
      allocate (vl(nc(1)), vr(nc(1)), vext(1 - (k - 1):nc(1) + (k - 1)))
      do concurrent(j=1:nc(2))
         ! Populate extended 'v' vector with ghost cells
         vext(1:nc(1)) = varray(:, j)
         vext(:0) = vext(1)
         vext(nc(1) + 1:) = vext(nc(1))
         !Get reconstructed values at cell boundaries
         call wenok(k, eps, vext, vl, vr)
         do concurrent(i=1:nc(1) - 1)
            fedges(i, j, 1) = godunov(flux1, vr(i), vl(i + 1), &
                                      [gx(1)%right(i), gx(2)%center(j)], t)
         end do
      end do

      ! Fluxes along x2 at interior cell boundaries
      deallocate (vl, vr, vext)
      allocate (vl(nc(2)), vr(nc(2)), vext(1 - (k - 1):nc(2) + (k - 1)))
      do concurrent(i=1:nc(1))
         ! Populate extended 'v' vector with ghost cells
         vext(1:nc(2)) = varray(i, :)
         vext(:0) = vext(1)
         vext(nc(2) + 1:) = vext(nc(2))
         !Get reconstructed values at cell boundaries
         call wenok(k, eps, vext, vl, vr)
         do concurrent(j=1:nc(2) - 1)
            fedges(i, j, 2) = godunov(flux2, vr(j), vl(j + 1), &
                                      [gx(1)%center(i), gx(2)%right(j)], t)
         end do
      end do

      ! Apply problem-specific flux constraints at domain boundaries
      fedges(0, :, 1) = 0
      fedges(nc(1), :, 1) = 0
      fedges(:, 0, 2) = 0
      fedges(:, nc(2), 2) = 0

      ! Evaluate du/dt
      do concurrent(i=1:nc(1), j=1:nc(2))
         vdot((j - 1)*nc(1) + i) = &
            -(fedges(i, j, 1) - fedges(i - 1, j, 1))/gx(1)%width(i) &
            - (fedges(i, j, 2) - fedges(i, j - 1, 2))/gx(2)%width(j)
      end do

   end subroutine rhs

   pure real(rk) function flux1(v, x, t)
   !! Flux function along x1.
      real(rk), intent(in) :: v
        !! function v(z,t)
      real(rk), intent(in) :: x(:)
        !! vector of internal coordinates
      real(rk), intent(in) :: t
        !! time

      flux1 = v !*x(1)**2

   end function flux1

   pure real(rk) function flux2(v, x, t)
   !! Flux function along x2.
      real(rk), intent(in) :: v
        !! function v(z,t)
      real(rk), intent(in) :: x(:)
        !! vector of internal coordinates
      real(rk), intent(in) :: t
        !! time

      flux2 = v !*x(1)*x(2)

   end function flux2

   pure real(rk) function ic(x)
   !! Initial condition. Here we used rectangular pulse in both coordinates.
      real(rk), intent(in) :: x(2)
         !! vector of internal coordinates

      if ((x(1) >= 1._rk .and. x(1) <= 3._rk) .and. (x(2) >= 1._rk .and. x(2) <= 3._rk)) then
         ic = 1._rk
      else
         ic = 0._rk
      end if

   end function ic

   subroutine output(message)
   !! Auxiliary routine to save results to file.
      integer, intent(in) :: message
         !! parameter to select output action

      character(*), parameter :: folder = "./output/example2/"
      real(rk) :: cpu_start = 0._rk, cpu_end = 0._rk
      integer :: i, j, funit_x(size(nc)) = 0, funit_u = 0

      select case (message)

         ! Open files and write headers and grid
      case (1)

         write (stdout, '(1x, a)') "Running example2..."
         write (stdout, '(1x, a, 1x, a)') "Start:", fdate()
         call cpu_time(cpu_start)

         ! Write grid
         do concurrent(j=1:size(nc))
            open (newunit=funit_x(j), file=folder//"x"//to_string(j)//".txt", &
                  status="replace", action="write", position="rewind")

            write (funit_x(j), '(a5, 2(1x, a15))') "i", "x(i)", "dx(i)"
            do i = 1, nc(j)
               write (funit_x(j), '(i5, 2(1x, es15.5))') i, gx(j)%center(i), gx(j)%width(i)
            end do
         end do

         ! Write header u
         open (newunit=funit_u, file=folder//"u.txt", status="replace", &
               action="write", position="rewind")

         write (funit_u, '(a16)', advance="no") "t"
         do i = 1, size(u)
            write (funit_u, '(1x, a16)', advance="no") "u("//to_string(i)//")"
         end do
         write (funit_u, *) ""

         ! Write values u(x,t)
      case (2)
         write (funit_u, '(es16.5e3)', advance="no") time
         do i = 1, size(u)
            write (funit_u, '(1x, es16.5e3)', advance="no") u(i)
         end do
         write (funit_u, *) ""

         ! Close files
      case (3)
         do concurrent(j=1:size(nc))
            close (funit_x(j))
         end do
         close (funit_u)
         write (stdout, '(1x, a, 1x, a)') "End  :", fdate()
         call cpu_time(cpu_end)
         write (stdout, '(1x, a, 1x, f7.1)') "Elaspsed time (ms) :", 1e3*(cpu_end - cpu_start)

      end select

   end subroutine output

end program