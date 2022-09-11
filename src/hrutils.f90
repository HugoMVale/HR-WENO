module hrutils
!!   This module contains basic flux schemes for *scalar* problems. They are mostly
!! intented to help test the other modules. The WENO schemes themselves are applicable to
!! scalar and multiple component problems.
!!   Source: ICASE 97-65 by Shu, 1997.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   public :: lax_friedrichs, godunov, grid1, tgrid1

   integer, parameter :: rk = real64

   abstract interface
      pure function flux(u, x, t)
         import :: rk
         real(rk) :: flux
         real(rk), intent(in) :: u, x, t
      end function
   end interface

   type :: tgrid1
    !! 1D grid
      real(rk), allocatable :: edges(:)
        !! vector(0:nc) of cell edges
      real(rk), allocatable :: c(:)
        !! vector(nc) of cell centers, \( x_i \)
      real(rk), allocatable :: d(:)
        !! vector(nc) of cell widths,  \( x_{i+1/2} - x_{i-1/2} \)
      real(rk), allocatable :: l(:)
        !! vector(nc) of left cell boundaries, \( x_{i-1/2} \)
      real(rk), allocatable :: r(:)
        !! vector(nc) of right cell boundaries, , \( x_{i+1/2} \)
      integer :: nc
        !! number of cells
   end type

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
        !! flux function f(v, x, t)
      real(rk), intent(in) :: vm
        !! left (minus) reconstruction \( v_{i+1/2}^- \)
      real(rk), intent(in) :: vp
        !! right (plus) reconstruction \( v_{i+1/2}^+ = v_{(i+1)+1/2}^- \)
      real(rk), intent(in) :: x
        !! x at flux interface, \( x_{i+1/2} \)
      real(rk), intent(in) :: t
        !! time
      real(rk), intent(in) :: alpha
        !! max(abs(f'(v))) in the domain on the problem

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
        !! flux function f(v, x, t)
      real(rk), intent(in) :: vm
        !! left (minus) reconstruction \( v_{i+1/2}^- \)
      real(rk), intent(in) :: vp
        !! right (plus) reconstruction \( v_{i+1/2}^+ = v_{(i+1)+1/2}^- \)
      real(rk), intent(in) :: x
        !! x at flux interface, \( x_{i+1/2} \)
      real(rk), intent(in) :: t
        !! time

      real(rk) :: fm, fp

      fm = f(vm, x, t)
      fp = f(vp, x, t)

      if (vm <= vp) then
         godunov = min(fm, fp)
      else
         godunov = max(fm, fp)
      end if

   end function godunov

   pure type(tgrid1) function grid1(xmin, xmax, nc)
    !!   Function to generate a 1D linear grid.
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: nc
        !! number of grid cells

      real(rk) :: xedges(0:nc), rx
      integer :: i

      ! Compute linear mesh
      rx = (xmax - xmin)/nc
      do concurrent(i=0:nc)
         xedges(i) = xmin + rx*i
      end do

      ! Map values to grid object
      grid1%nc = nc
      grid1%edges = xedges
      grid1%l = xedges(0:nc - 1)
      grid1%r = xedges(1:nc)
      grid1%c = (grid1%l + grid1%r)/2
      grid1%d = grid1%r - grid1%l

   end function grid1

end module hrutils
