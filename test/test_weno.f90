module test_weno
!>----------------------------------------------------------------------------------------------
!> Test for module 'weno' using test-drive.
!> Hugo Vale
!>----------------------------------------------------------------------------------------------
    use weno, only : weno35, calc_cgrid
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
        new_unittest("weno35-uniform", test_weno35_uniform), &
        new_unittest("calc_cgrid", test_calc_cgrid) &
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
            write (*, '(2(a4),3(a26))'), "k", "i", "v(i)", "vl(i)", "vr(i)"
            do i = 1, nc
              write (*, '(2(i4),3(es26.16e3))'), k, i, v(i), vl(i), vr(i)
            end do
          end if

      end do

    end subroutine


    subroutine test_calc_cgrid(error)
      type(error_type), allocatable, intent(out) :: error
      real(rk), allocatable :: xedges(:), cgrid(:,:,:)
      real(rk) :: atol, xmin, xmax, xstep
      integer :: i, k, nc 
     
      !> Common grid settings
      xmin = 0_rk
      xmax = 1_rk
      nc = 4
      xstep = (xmax - xmin)/nc

      !> Run check for each order
      atol = 1e-9_rk
      do k = 2, 3
        
        !> Allocate and define an abritrary uniform grid with (k-1) ghost cells on each side
        if(allocated(xedges)) deallocate(xedges) 
        allocate(xedges(1-k:nc+(k-1)))
        xedges(0) = xmin
        do i = 1, ubound(xedges,1)
          xedges(i) = xedges(i-1) + xstep
        end do
        do i = 0, lbound(xedges,1), -1
          xedges(i) = xedges(i+1) - xstep
        end do

        print *, lbound(xedges,1)
        print *, ubound(xedges,1)
        print *, size(xedges)
        print *, xedges

        !> Allocate cgrid array
        if(allocated(cgrid)) deallocate(cgrid) 
        allocate(cgrid(0:k-1,-1:k-1,1:nc))

        !> Compute cgrid
        call calc_cgrid(k, xedges, cgrid)
               
        !> Check error
        ! do i = 1, nc
        !   call check(error, v(i), vl(i), thr=atol)
        !   if (allocated(error)) return
        !   call check(error, v(i), vr(i), thr=atol)
        !   if (allocated(error)) return
        ! end do 
      end do
    
    end subroutine

end module test_weno
