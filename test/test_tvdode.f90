module test_tvdode
!! Test for module 'tvdode' using test-drive.
   use, intrinsic :: iso_fortran_env, only: stderr => error_unit
   use hrweno_kinds, only: rk
   use hrweno_tvdode, only: rktvd, mstvd
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_tests_tvdode

   logical, parameter :: verbose = .false.
   integer, parameter :: nu = 10
   integer :: ii
   real(rk), parameter :: a(nu) = [(-1._rk + real(ii - 1, rk)*4/(nu - 1), ii=1, nu)]

contains

   !> Collect all exported unit tests
   subroutine collect_tests_tvdode(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rktvd", test_rktvd), &
                  new_unittest("mstvd", test_mstvd) &
                  ]

   end subroutine

   subroutine test_rktvd(error)
      type(error_type), allocatable, intent(out) :: error
      type(rktvd) :: ode
      integer :: order
      real(rk) :: t0, t, tout, dt, u(nu), u0(nu), uref(nu), rtol(3)
      integer :: i

      ! Initial conditions and ode settings
      t0 = -1.3_rk
      tout = 1.5_rk
      u0 = 0.1_rk
      dt = ((tout - t0)/3000)

      ! Run check for each order
      rtol = [2e-2_rk, 1e-3_rk, 1e-3_rk]
      do order = 1, 3

         ! Init ode object
         ode = rktvd(fu, nu, order)

         ! Numerical solution at t=tout
         u = u0
         t = t0
         call ode%integrate(u, t, tout, dt)

         ! Check error
         uref = usol(t0, u0, t)
         call check(error, u, uref, rel=.true., thr=rtol(order))

         ! Show results if test fails (for debugging)
         if (allocated(error) .or. verbose) then
            write (stderr, '(2(a6), 2(a26))') "order", "i", "u(i)", "uref(i)"
            do i = 1, size(u)
               write (stderr, '(2(i6), 2(es26.16e3))') order, i, u(i), uref(i)
            end do
         end if

      end do

   end subroutine test_rktvd

   subroutine test_mstvd(error)
      type(error_type), allocatable, intent(out) :: error
      type(mstvd) :: ode
      real(rk) :: t0, t, tout, dt, u(nu), u0(nu), uref(nu)
      integer :: i

      ! Initial conditions and ode settings
      t0 = -1.3_rk
      tout = 1.5_rk
      u0 = 0.1_rk
      dt = 1e-3_rk

      ! Init ode object
      ode = mstvd(fu, nu)

      ! Numerical solution at t=tout
      u = u0
      t = t0
      call ode%integrate(u, t, tout, dt)

      ! Check error
      uref = usol(t0, u0, t)
      call check(error, u, uref, rel=.true., thr=1e-3_rk)

      ! Show results if test fails (for debugging)
      if (allocated(error) .or. verbose) then
         write (stderr, '(a4, 2(a26))') "i", "u(i)", "uref(i)"
         do i = 1, size(u)
            write (stderr, '(i4, 2(es26.16e3))') i, u(i), uref(i)
         end do
      end if

   end subroutine test_mstvd

   !> Simple linear u'(u) to test ode solvers
   pure subroutine fu(t, u, udot)
      real(rk), intent(in) :: t, u(:)
      real(rk), intent(out) :: udot(:)
      udot = a*u
   end subroutine

   !> Analytical solution of u'(u)=au to test ode solvers
   pure function usol(t0, u0, t)
      real(rk), intent(in) :: t0, t, u0(:)
      real(rk) :: usol(size(u0))
      usol = u0*exp(a*(t - t0))
   end function

end module test_tvdode
