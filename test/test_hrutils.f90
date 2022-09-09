module test_hrutils
    !>---------------------------------------------------------------------------------------------
    !> Test for module 'hrutils' using test-drive.
    !>---------------------------------------------------------------------------------------------
        use hrutils, only : godunov, lax_friedrichs 
        use iso_fortran_env, only : real64, error_unit
        use testdrive, only : new_unittest, unittest_type, error_type, check
        implicit none
        private
    
        public :: collect_tests_hrutils
    
        integer, parameter :: rk = real64
        logical, parameter :: verbose = .false.
    
        contains
    
        !> Collect all exported unit tests
        subroutine collect_tests_hrutils(testsuite)
            !> Collection of tests
            type(unittest_type), allocatable, intent(out) :: testsuite(:)
    
            testsuite = [ &
            new_unittest("fluxes", test_fluxes) &
            !> new_unittest("mstvd", test_mstvd) &
            ]
    
        end subroutine
    
        subroutine test_fluxes(error)
            type(error_type), allocatable, intent(out) :: error
            integer :: order, itask, istate
            real(rk) :: h, href, vm, vp, x, t, rtol
            integer :: i
    
            !> General settings
            vm = 2._rk
            x = 3._rk
            t = 5._rk
            rtol = 1e-8_rk

            !> Check h(a,a) = h(a)
            write(error_unit, *) "Check h(a,a) = h(a)"  
            vp = vm
            href = vm*x*t 
            h = godunov(f, vm, vp, x, t)
            call check(error, h, href, rel=.true., thr=rtol)
            if (allocated(error)) return
            h = lax_friedrichs(f, vm, vp, x, t, x*t)
            call check(error, h, href, rel=.true., thr=rtol)
            if (allocated(error)) return

            !> Check wind blows from left to right, h(a,b) =~ h(a)
            write(error_unit, *) "h(a,b) =~ h(a)"
            vp = -2*vm
            href = vm*x*t 
            h = godunov(f, vm, vp, x, t)
            call check(error, h, href, rel=.true., thr=rtol)
            if (allocated(error)) return
            h = lax_friedrichs(f, vm, vp, x, t, x*t)
            write(error_unit, *) href, h
            call check(error, h, href, rel=.true., thr=rtol)
            if (allocated(error)) return
            stop

            !> Check wind blows from right to left, h(a,b) =~ h(b)
            vm = - vm
            vp = - vp
            write(error_unit, *) "h(a,b) =~ h(b)"  
            href = vp*x*t 
            h = godunov(f, vm, vp, x, t)
            write(error_unit, *) href, h
            call check(error, h, href, rel=.true., thr=rtol)
            if (allocated(error)) return
            !h = lax_friedrichs(f, vm, vp, x, t, x*t)
            write(error_unit, *) href, h
            call check(error, h, href, rel=.true., thr=rtol)
            if (allocated(error)) return
    
        end subroutine test_fluxes
       
        pure function f(u, x, t)
            !> Simple flux function to test numerical fluxes
            real(rk) :: f
            real(rk), intent(in) :: u, x, t
            f = u*x*t
        end function
        
    end module test_hrutils
    