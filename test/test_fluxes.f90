module test_fluxes
!! Test for module 'fluxes' using test-drive.
   use fluxes, only: godunov, lax_friedrichs
   use iso_fortran_env, only: real64, stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_tests_fluxes

   integer, parameter :: rk = real64
   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_fluxes(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fluxes", test_allfluxes) &
                  ]

   end subroutine

   subroutine test_allfluxes(error)
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

   end subroutine test_allfluxes

   !> Simple flux function to test numerical fluxes
   pure function f(u, x, t)
      real(rk) :: f
      real(rk), intent(in) :: u, x(:), t
      f = u*x(1)*t
   end function

end module test_fluxes
