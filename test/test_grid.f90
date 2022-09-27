module test_grid
    !! Test for module 'grid' using test-drive.
   use grid, only: grid1
   use iso_fortran_env, only: real64, stderr => error_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_tests_grid

   integer, parameter :: rk = real64
   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_grid(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("grid1", test_grid1) &
                  ]

   end subroutine

   subroutine test_grid1(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nc, scl
      real(rk) :: xmin, xmax, atol
      type(grid1) :: gx

      ! General settings
      xmin = 1e-1_rk
      xmax = 1e3_rk
      nc = 1000000
      atol = 1e-10_rk

      do scl = 1, 2

         ! Make grid
         call gx%new(xmin, xmax, nc, scl=scl)

         !print *, scl
         !print *, "nc", gx%ncells
         !print *, "edges", gx%edges
         ! print *, "d", gx%d
         ! print *, "c", gx%c
         ! print *, "l", gx%l
         ! print *, "r", gx%r

         ! Checks
         call check(error, gx%ncells, nc)
         if (allocated(error)) return
         call check(error, gx%edges(0), xmin, thr=atol)
         if (allocated(error)) return
         call check(error, gx%edges(nc), xmax, thr=atol)
         if (allocated(error)) return
         call check(error, gx%left, gx%edges(0:nc - 1))
         if (allocated(error)) return
         call check(error, gx%right, gx%edges(1:nc))
         if (allocated(error)) return
         call check(error, gx%width, gx%right - gx%left)
         if (allocated(error)) return
         call check(error, gx%center, (gx%left + gx%right)/2)

      end do

   end subroutine test_grid1

end module test_grid
