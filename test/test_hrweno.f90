module test_hrweno
!! Test for module 'hrweno' using test-drive.
   use, intrinsic :: iso_fortran_env, only: stderr => error_unit
   use hrweno_kinds, only: rk
   use hrweno_weno, only: weno, c1, c2, c3
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_tests_hrweno

   logical, parameter :: verbose = .false.

contains

   !> Collect all exported unit tests
   subroutine collect_tests_hrweno(testsuite)
      ! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("weno with uniform grid", test_weno_uniform), &
                  new_unittest("array 'c' for non-uniform grid", test_calc_cnu), &
                  new_unittest("weno with non-uniform grid", test_weno_nonuniform) &
                  ]

   end subroutine

   subroutine test_weno_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(weno) :: myweno
      integer, parameter :: nc = 30
      real(rk) :: atol, v(nc), vl(nc), vr(nc)
      integer :: i, k

      ! Run check for each order
      atol = 1e-6_rk
      do k = 1, 3

         ! Init weno object
         myweno = weno(nc, k, eps=1e-6_rk)

         ! Set cell average value
         ! Just a reactangular pulse __|¯¯|__
         v = 0
         v(nc/3:2*nc/3) = 1

         ! Perform reconstruction
         call myweno%reconstruct(v, vl, vr)

         ! Check error
         call check(error, v, vl, thr=atol)
         if (.not. (allocated(error))) then
            call check(error, v, vr, thr=atol)
         end if

         ! Detailed comparison for debugging
         if (allocated(error) .or. verbose) then
            write (stderr, '(2(a4),3(a26))') "k", "i", "v(i)", "vl(i)", "vr(i)"
            do i = 1, nc
               write (stderr, '(2(i4),3(es26.16e3))') k, i, v(i), vl(i), vr(i)
            end do
         end if

      end do

   end subroutine test_weno_uniform

   subroutine test_calc_cnu(error)
      type(error_type), allocatable, intent(out) :: error
      type(weno) :: myweno
      integer, parameter :: nc = 30
      real(rk) :: xedges(0:nc)
      real(rk), allocatable :: cref(:, :)
      real(rk) :: rtol, xmin, xmax
      integer :: i, k

      ! Allocate and define an abritrary uniform grid
      xmin = 0.0_rk
      xmax = 3.0_rk
      do i = 0, nc
         xedges(i) = xmin + (xmax - xmin)*i/nc
      end do

      ! Run check for each order
      rtol = 1e-5_rk
      do k = 1, 3

         ! Init weno object
         myweno = weno(nc, k, eps=1e-6_rk, xedges=xedges)

         ! Get reference solution
         if (allocated(cref)) deallocate (cref)
         select case (k)
         case (1)
            cref = c1
         case (2)
            cref = c2
         case (3)
            cref = c3
         end select

         ! Check error
         do i = 1, nc
            call check(error, reshape(myweno%cnu(:, :, i), [size(cref)]), &
                       reshape(cref, [size(cref)]), rel=.true., thr=rtol)
            if (allocated(error)) return
         end do

      end do

   end subroutine test_calc_cnu

   subroutine test_weno_nonuniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(weno) :: myweno
      integer, parameter :: nc = 30
      real(rk) :: v(nc), vl(nc), vr(nc), xedges(0:nc)
      real(rk) :: atol, xmin, xmax
      integer :: i, k

      ! Allocate and define an abritrary non-uniform grid
      xmin = 0.0_rk
      xmax = 1.0_rk
      do i = 0, nc
         xedges(i) = xmin + (xmax - xmin)*i/nc
      end do
      xedges = xedges**3

      ! Run check for each order
      atol = 1e-6_rk
      do k = 1, 3

         ! Init weno object
         myweno = weno(nc, k, xedges=xedges)

         ! Set cell average value
         ! Just a reactangular pulse __|¯¯|__
         v = 0
         v(nc/3:2*nc/3) = 1

         ! Perform reconstruction
         call myweno%reconstruct(v, vl, vr)

         ! Check error
         call check(error, vl, v, thr=atol)
         if (.not. (allocated(error))) then
            call check(error, vr, v, thr=atol)
         end if

         ! Detailed comparison for debugging
         if (allocated(error) .or. verbose) then
            write (stderr, '(2(a4),3(a26))') "k", "i", "v(i)", "vl(i)", "vr(i)"
            do i = 1, nc
               write (stderr, '(2(i4),3(es26.16e3))') k, i, v(i), vl(i), vr(i)
            end do
         end if

      end do

   end subroutine test_weno_nonuniform

end module test_hrweno
