module test_weno
!>---------------------------------------------------------------------------------------------
!> Test for module 'weno' using test-drive.
!>---------------------------------------------------------------------------------------------
    use weno, only : weno35, calc_c, c2, c3
    use iso_fortran_env, only : real64, error_unit
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_tests_weno

    integer, parameter :: rk = real64

    contains

    !> Collect all exported unit tests
    subroutine collect_tests_weno(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
        new_unittest("weno35 with uniform grid", test_weno35_uniform), &
        new_unittest("calc_c", test_calc_c), &
        new_unittest("weno35 with non-uniform grid", test_weno35_nonuniform) &
        ]

    end subroutine

    subroutine test_weno35_uniform(error)
        type(error_type), allocatable, intent(out) :: error
        real(rk), dimension(:), allocatable :: v, vl, vr
        real(rk) :: eps, atol
        integer :: i, k, nc = 30

        !> Run check for each order
        eps = 1e-6_rk
        atol = 1e-9_rk
        do k = 2, 3

          !> Allocate arrays
          if(allocated(v)) deallocate(v)
          if(allocated(vl)) deallocate(vl)
          if(allocated(vr)) deallocate(vr)
          allocate(v(1-(k-1):nc+(k-1)), vl(nc), vr(nc))

          !> Set cell average value including ghost cells
          !> Just a reactangular pulse __|¯¯|__
          v = 0_rk
          v(nc/3:2*nc/3) = 1_rk

          !> Call procedure
          call weno35(k, v, vl, vr, eps)

          !> Check error
          call check(error, v(1:nc), vl, thr=atol)
          call check(error, v(1:nc), vr, thr=atol)

          !> Detailed comparison for debugging
          if (allocated(error)) then
            write (error_unit, '(2(a4),3(a26))'), "k", "i", "v(i)", "vl(i)", "vr(i)"
            do i = 1, nc
              write (error_unit, '(2(i4),3(es26.16e3))'), k, i, v(i), vl(i), vr(i)
            end do
          end if

      end do

    end subroutine test_weno35_uniform


    subroutine test_calc_c(error)
      type(error_type), allocatable, intent(out) :: error
      real(rk), allocatable :: xedges(:), c(:,:,:), cref(:,:)
      real(rk) :: rtol, xmin, xmax
      integer :: i, k, nc

      !> Allocate and define an abritrary uniform grid
      xmin = 0_rk
      xmax = 3_rk
      nc = 30
      allocate(xedges(0:nc))
      do i = 0, nc
        xedges(i) = xmin + (xmax - xmin)*i/nc
      end do

      !> Run check for each order
      rtol = 1e-9_rk
      do k = 2, 3

        !> Allocate c array
        if(allocated(c)) deallocate(c)
        allocate(c(0:k-1,-1:k-1,1:nc))

        !> Compute c
        call calc_c(k, xedges, c)

        !> Get reference solution
        if(allocated(cref)) deallocate(cref)
        select case(k)
          case(2)
            cref = c2
          case(3)
            cref = c3
        end select

        !> Check error
        do i = 1, nc
          call check(error, reshape(c(:,:,i),[size(cref)]), &
               reshape(cref,[size(cref)]), rel=.true., thr=rtol)
          if (allocated(error)) return
        end do

      end do

    end subroutine test_calc_c


    subroutine test_weno35_nonuniform(error)
      type(error_type), allocatable, intent(out) :: error
      real(rk), allocatable :: v(:), vl(:), vr(:), xedges(:), c(:,:,:)
      real(rk) :: eps, atol, xmin, xmax
      integer :: i, k, nc

      !> Allocate and define an abritrary non-uniform grid
      xmin = 0_rk
      xmax = 1_rk
      nc = 30
      allocate(xedges(0:nc))
      do i = 0, nc
        xedges(i) = xmin + (xmax - xmin)*i/nc
      end do
      xedges = xedges**3

      !> Run check for each order
      eps = 1e-6_rk
      atol = 1e-9_rk
      do k = 2, 3

        !> Allocate 'c' array
        if(allocated(c)) deallocate(c)
        allocate(c(0:k-1,-1:k-1,1:nc))

        !> Allocate 'v' arrays
        if(allocated(v)) deallocate(v)
        if(allocated(vl)) deallocate(vl)
        if(allocated(vr)) deallocate(vr)
        allocate(v(1-(k-1):nc+(k-1)), vl(nc), vr(nc))

        !> Set cell average value including ghost cells
        !> Just a reactangular pulse __|¯¯|__
        v = 0_rk
        v(nc/3:2*nc/3) = 1_rk

        !> Call procedures
        call calc_c(k, xedges, c)
        call weno35(k, v, vl, vr, eps, c)

        !> Check error
        call check(error, v(1:nc), vl, thr=atol)
        call check(error, v(1:nc), vr, thr=atol)

        !> Detailed comparison for debugging
        if (allocated(error)) then
          write (error_unit, '(2(a4),3(a26))'), "k", "i", "v(i)", "vl(i)", "vr(i)"
          do i = 1, nc
            write (error_unit, '(2(i4),3(es26.16e3))'), k, i, v(i), vl(i), vr(i)
          end do
        end if

    end do

  end subroutine test_weno35_nonuniform

end module test_weno
