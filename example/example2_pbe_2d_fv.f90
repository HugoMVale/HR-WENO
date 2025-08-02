program example_pbe_2d_fv
!!   This program illustrates the application of modules 'hrweno' and 'mstvd' for solving a 2D
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
!! below.
   use, intrinsic :: iso_fortran_env, only: stderr => error_unit, stdout => output_unit
   use hrweno_kinds, only: rk
   use hrweno_tvdode, only: mstvd
   use hrweno_weno, only: weno
   use hrweno_fluxes, only: godunov
   use hrweno_grids, only: grid1
   use stdlib_strings, only: to_string
   ! use omp_lib
   implicit none

   integer, parameter :: nc(2) = [250, 250], k = 3
   real(rk) :: u(product(nc))
   type(grid1) :: gx(2)
   type(weno) :: myweno(2)
   type(mstvd) :: ode
   real(rk) :: dt, time, time_out, time_start, time_end
   integer :: num_time_points, ii, jj

   ! OMP settings
   ! call omp_set_num_threads(min(2, omp_get_num_procs()))

   ! Define grids for x1 and x2
   ! In this example, we use linear grids, but any smooth grid can be used
   call gx(1)%linear(xmin=0.0_rk, xmax=10.0_rk, ncells=nc(1))
   call gx(2)%linear(xmin=0.0_rk, xmax=10.0_rk, ncells=nc(2))

   ! Init weno objects
   myweno(1) = weno(ncells=nc(1), k=k, eps=1e-6_rk)
   myweno(2) = weno(ncells=nc(2), k=k, eps=1e-6_rk)

   ! Open file where results will be stored
   call output(1)

   ! Initial condition u(x,t=0)
   do concurrent(ii=1:nc(1), jj=1:nc(2))
      u((jj - 1)*nc(1) + ii) = ic([gx(1)%center(ii), gx(2)%center(jj)])
   end do

   ! Call ODE time solver
   ode = mstvd(rhs, size(u))

   time_start = 0.0_rk
   time_end = 5.0_rk
   dt = 5e-3_rk

   time = time_start
   num_time_points = 100
   do ii = 0, num_time_points
      time_out = time_end*ii/num_time_points
      call ode%integrate(u, time, time_out, dt)
      call output(2)
   end do

   ! End of simulation
   call output(3)

contains

   subroutine rhs(t, v, vdot)
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

      real(rk) :: vl1(nc(1)), vl2(nc(2)), vr1(nc(1)), vr2(nc(2))
      real(rk) :: fedges1(0:nc(1), 0:nc(2)), fedges2(0:nc(2), 0:nc(1))
      integer :: i, j

      ! Fluxes along x1 and x2 at interior cell boundaries
      !$omp parallel
      !$omp do private(vl1, vr1)
      do j = 1, nc(2)
         call myweno(1)%reconstruct(v(1+(j-1)*nc(1):j*nc(1)), vl1, vr1)
         do i = 1, nc(1) - 1
            fedges1(i, j) = godunov(flux1, vr1(i), vl1(i + 1), &
                                    [gx(1)%right(i), gx(2)%center(j)], t)
         end do
      end do
      !$omp end do nowait
      !$omp do private(vl2, vr2)
      do i = 1, nc(1)
         call myweno(2)%reconstruct(v(i:i+(nc(2)-1)*nc(1):nc(1)), vl2, vr2)
         do j = 1, nc(2) - 1
            fedges2(j, i) = godunov(flux2, vr2(j), vl2(j + 1), &
                                    [gx(1)%center(i), gx(2)%right(j)], t)
         end do
      end do
      !$omp end do
      !$omp end parallel

      ! Apply problem-specific flux constraints at domain boundaries
      fedges1(0, :) = 0
      fedges1(nc(1), :) = 0
      fedges2(0, :) = 0
      fedges2(nc(2), :) = 0

      ! Evaluate du/dt
      do concurrent(i=1:nc(1), j=1:nc(2))
         vdot((j - 1)*nc(1) + i) = &
            - (fedges1(i, j) - fedges1(i - 1, j))/gx(1)%width(i) &
            - (fedges2(j, i) - fedges2(j - 1, i))/gx(2)%width(j)
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

      if ((x(1) >= 1.0_rk .and. x(1) <= 3.0_rk) .and. (x(2) >= 1.0_rk .and. x(2) <= 3.0_rk)) then
         ic = 1.0_rk
      else
         ic = 0.0_rk
      end if

   end function ic

   subroutine output(action)
   !! Auxiliary routine to save results to file.
      integer, intent(in) :: action
         !! parameter to select output action

      character(*), parameter :: folder = ".\output\example2\"
      real(rk) :: time_elapsed
      integer, save :: funit_x(size(nc)) = 0, funit_u = 0
      integer :: i, j

      select case (action)

      ! Open files and write headers and grid
      case (1)

         write (stdout, '(1x, a)') "Running example2 ..."
         ! write (stdout, '(1x, a, i3)') "Max # threads: ", omp_get_max_threads()
         
         ! Write grid
         do concurrent(j=1:size(nc))
            open (newunit=funit_x(j), file=folder//"x"//to_string(j)//".txt", &
                  status="replace", action="write", position="rewind")

            write (funit_x(j), '(a5, 2(1x, a15))') "i", "x(i)", "dx(i)"
            do i = 1, nc(j)
               write (funit_x(j), '(i5, 2(1x, es15.5))') i, gx(j)%center(i), gx(j)%width(i)
            end do
         end do
         
         call timer()

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
         call timer(time_elapsed)
         write (stdout, '(1x, a, 1x, f5.2)') "Elapsed time (s) :", time_elapsed

      end select

   end subroutine output

   subroutine timer(res)
      !! Quick and dirty timer
      real(rk), intent(out), optional :: res
      real(rk) :: res_
      logical, save :: first_call = .true.
      integer, save :: values_previous(8)
      integer :: values_now(8), delta(8)

      call date_and_time(values=values_now)
      if (first_call) then
         values_previous = values_now
         first_call = .false.
      end if

      delta = values_now - values_previous
      res_ = 1.0_rk*delta(3)*24*60**2 + 1.0_rk*delta(5)*60**2 + 1.0_rk*delta(6)*60 &
            + 1.0_rk*delta(7) + 1.0_rk*delta(8)/1000
      
      if (present(res)) res = res_

   end subroutine

end program example_pbe_2d_fv


