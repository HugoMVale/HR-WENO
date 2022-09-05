module hrutils
    !>---------------------------------------------------------------------------------------------
    !>   This module contains a collection of high-resolution weighted essentially non-oscillatory
    !> (WENO) schemes for *arbitrary* (uniform or non-uniform) finite volume/difference methods.
    !>   Source: ICASE 97-65 by Shu, 1997.
    !>---------------------------------------------------------------------------------------------
        use, intrinsic :: iso_fortran_env, only : real64
        implicit none
        private
    
        public :: lax_friedrichs, godunov
    
        integer, parameter :: rk = real64
     
        abstract interface
            pure function flux(u, t)
                import :: rk
                real(rk) :: flux
                real(rk), intent(in) :: u, t
            end function
        end interface
    
        contains
      
        pure real(rk) function lax_friedrichs(f, vm, vp, t, alpha)
        !>-----------------------------------------------------------------------------------------
        !>   Monotone Lax-Friedrichs flux. It is more dissipative than the Godunov method, but 
        !> computationally less demanding.
        !>   Source: Equation 2.72, page 21.
        !>
        !> ARGUMENTS:
        !> f      flux function f(v, t)
        !> vm     left (minus) reconstruction v_{i^+1/2}^-
        !> vp     right (plus) reconstruction v_{i^+1/2}^+ = v_{(i+1)^+1/2}^-
        !> t      time
        !> alpha  max(abs(f'(v))) in the domain on the problem
        !>
        !> NOTES:
        !> - Although potentially useful, this procedure cannot be defined as *elemental*,
        !>   because it has a dummy procedure as input argument.
        !>-----------------------------------------------------------------------------------------
        procedure(flux) :: f
        real(rk), intent (in) :: vm, vp, t, alpha
    
            lax_friedrichs = (f(vm, t) + f(vp, t) - alpha*(vp - vm))/2
    
        end function lax_friedrichs
        
    
        pure real(rk) function godunov(f, vm, vp, t)
        !>-----------------------------------------------------------------------------------------
        !>   Monotone Godunov flux. It is less dissipative than the Lax-Friedrichs method, but 
        !> computationally more demanding because of the if constructs.
        !>   Source: Equation 2.70, page 21.
        !>
        !> ARGUMENTS:
        !> f      flux function f(v, t)
        !> vm     left (minus) reconstruction v_{i^+1/2}^-
        !> vp     right (plus) reconstruction v_{i^+1/2}^+ = v_{(i+1)^+1/2}^-
        !> t      time
        !>
        !> NOTES:
        !> - See note about *elemental* in 'lax_friedrichs'.
        !>-----------------------------------------------------------------------------------------
        procedure(flux) :: f
        real(rk), intent (in) :: vm, vp, t
        real(rk) :: fm, fp
    
            fm = f(vm, t)
            fp = f(vp, t)
            
            if (vm <= vp) then
                godunov = min(fm, fp)
            else
                godunov = max(fm, fp)
            end if
    
        end function godunov
    
    end module hrutils
    