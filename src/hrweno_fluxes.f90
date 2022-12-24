module hrweno_fluxes
!!   This module contains basic flux schemes for *scalar* problems. They are mostly
!! intented to help test the other modules. The WENO schemes themselves are applicable to
!! scalar and multiple component problems.
!!   Source: ICASE 97-65 by Shu, 1997.
   use hrweno_kinds, only: rk
   implicit none
   private

   public :: lax_friedrichs, godunov

   abstract interface
      pure function flux(u, x, t)
         import :: rk
         real(rk) :: flux
         real(rk), intent(in) :: u, x(:), t
      end function
   end interface

contains

   pure real(rk) function lax_friedrichs(f, vm, vp, x, t, alpha)
    !!   Monotone Lax-Friedrichs flux. It is more dissipative than the Godunov method, but
    !! computationally less demanding.
    !!   Source: Equation 2.72, page 21.
    !!
    !! @note
    !!   Although it might be useful, this procedure cannot be defined as *elemental*,
    !! because it has a dummy procedure as argument.
      procedure(flux) :: f
        !! flux function, \( f(v, x, t) \)
      real(rk), intent(in) :: vm
        !! left (minus) reconstruction, \( v_{i+1/2}^- \)
      real(rk), intent(in) :: vp
        !! right (plus) reconstruction, \( v_{i+1/2}^+ = v_{(i+1)-1/2}^- \)
      real(rk), intent(in) :: x(:)
        !! \(x\) at flux interface, \( x_{i+1/2} \)
      real(rk), intent(in) :: t
        !! time, \( t \)
      real(rk), intent(in) :: alpha
        !! \( \max(|f'(v)|) \) in the domain on the problem

      lax_friedrichs = (f(vm, x, t) + f(vp, x, t) - alpha*(vp - vm))/2

   end function lax_friedrichs

   pure real(rk) function godunov(f, vm, vp, x, t)
    !!   Monotone Godunov flux. It is less dissipative than the Lax-Friedrichs method, but
    !! computationally more demanding because of the if constructs.
    !!   Source: Equation 2.70, page 21.
    !!
    !! @note
    !!   See note about *elemental* in 'lax_friedrichs'.
      procedure(flux) :: f
        !! flux function, \( f(v, x, t) \)
      real(rk), intent(in) :: vm
        !! left (minus) reconstruction, \( v_{i+1/2}^- \)
      real(rk), intent(in) :: vp
        !! right (plus) reconstruction, \( v_{i+1/2}^+ = v_{(i+1)-1/2}^- \)
      real(rk), intent(in) :: x(:)
        !! \(x\) at flux interface, \( x_{i+1/2} \)
      real(rk), intent(in) :: t
        !! time, \( t \)

      real(rk) :: fm, fp

      fm = f(vm, x, t)
      fp = f(vp, x, t)

      if (vm <= vp) then
         godunov = min(fm, fp)
      else
         godunov = max(fm, fp)
      end if

   end function godunov

end module hrweno_fluxes
