module test_hrutils
!! Test for module 'hrutils' using test-drive.
   use hrutils, only: godunov, lax_friedrichs, grid1, tgrid1
   use iso_fortran_env, only: real64, stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_tests_hrutils

   integer, parameter :: rk = real64
   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_hrutils(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fluxes", test_fluxes), &
                  new_unittest("grid1", test_grid1) &
                  ]

   end subroutine

   subroutine test_fluxes(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: order, itask, istate
      real(rk) :: h, href, vm, vp, x(1), t, rtol
      integer :: i

      ! General settings
      vm = 2._rk
      x = 3._rk
      t = 5._rk
      rtol = 1e-8_rk

      ! Check h(a,a) = h(a)
      write (stderr, *) "Check h(a,a) = h(a)"
      vp = vm
      href = f(vp, x, t)
      h = godunov(f, vm, vp, x, t)
      call check(error, h, href, rel=.true., thr=rtol)
      if (allocated(error)) return
      h = lax_friedrichs(f, vm, vp, x, t, x(1)*t)
      call check(error, h, href, rel=.true., thr=rtol)
      if (allocated(error)) return

      ! Check Check ->
      write (stderr, *) "Check ->"
      vp = -2*vm
      href = f(vm, x, t)
      h = godunov(f, vm, vp, x, t)
      call check(error, h, href, rel=.true., thr=rtol)
      if (allocated(error)) return
      h = lax_friedrichs(f, vm, vp, x, t, x(1)*t)
      call check(error, h, href, rel=.true., thr=rtol)
      if (allocated(error)) return

      ! Check Check <-
      write (stderr, *) "Check <-"
      vm = -vm
      vp = -vp
      href = f(vm, x, t)
      h = godunov(f, vm, vp, x, t)
      call check(error, h, href, rel=.true., thr=rtol)
      if (allocated(error)) return
      h = lax_friedrichs(f, vm, vp, x, t, x(1)*t)
      call check(error, h, href, rel=.true., thr=rtol)
      if (allocated(error)) return

   end subroutine test_fluxes

   subroutine test_grid1(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nc
      real(rk) :: xmin, xmax
      type(tgrid1) :: gx

      ! General settings
      xmin = 0._rk
      xmax = 8._rk
      nc = 4

      ! Make grid
      gx = grid1(xmin, xmax, nc)

      ! print *, "nc", gx%nc
      ! print *, "edges", gx%edges
      ! print *, "d", gx%d
      ! print *, "c", gx%c
      ! print *, "l", gx%l
      ! print *, "r", gx%r

      ! Checks
      call check(error, gx%ncells, nc)
      if (allocated(error)) return
      call check(error, gx%edges(0), xmin)
      if (allocated(error)) return
      call check(error, gx%edges(nc), xmax)
      if (allocated(error)) return
      call check(error, gx%left, gx%edges(0:nc - 1))
      if (allocated(error)) return
      call check(error, gx%right, gx%edges(1:nc))
      if (allocated(error)) return
      call check(error, gx%width, gx%right - gx%left)
      if (allocated(error)) return
      call check(error, gx%center, (gx%left + gx%right)/2)

   end subroutine test_grid1

   !> Simple flux function to test numerical fluxes
   pure function f(u, x, t)
      real(rk) :: f
      real(rk), intent(in) :: u, x(:), t
      f = u*x(1)*t
   end function

end module test_hrutils
