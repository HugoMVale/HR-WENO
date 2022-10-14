module hrweno
!!   This module contains a collection of high-resolution weighted essentially non-oscillatory
!! (WENO) schemes for *arbitrary* (uniform or non-uniform) finite volume/difference methods.
!!   Source: ICASE 97-65 by Shu, 1997.
   use hrweno_kinds, only: rk
   implicit none
   private

   public :: weno, c1, c2, c3

   !> Parameter arrays for WENO methods
   real(rk), parameter :: d1(0:0) = 1._rk, &
                          d2(0:1) = [2._rk/3, 1._rk/3], &
                          d3(0:2) = [0.3_rk, 0.6_rk, 0.1_rk]
   real(rk), parameter :: &
      c1(0:0, -1:0) = reshape([1._rk, 1._rk], [1, 2], order=[1, 2]), &
      c2(0:1, -1:1) = reshape([3._rk/2, -1._rk/2, 1._rk/2, 1._rk/2, -1._rk/2, 3._rk/2], &
                              [2, 3], order=[1, 2]), &
      c3(0:2, -1:2) = reshape([11._rk/6, -7._rk/6, 1._rk/3, 1._rk/3, 5._rk/6, -1._rk/6, &
                               -1._rk/6, 5._rk/6, 1._rk/3, 1._rk/3, -7._rk/6, 11._rk/6], &
                              [3, 4], order=[1, 2])

   type :: weno
   !! WENO class.
      character(:), allocatable :: msg
         !! error message
      integer :: ierr = 0
         !! error status
      integer :: k
         !! order of reconstruction within the cell (k = 1, 2 or 3)
      integer :: ncells
         !! number of cells
      real(rk) :: eps = 1e-6_rk
         !! numerical smoothing factor
      real(rk), allocatable, private :: d(:)
         !! array of constants
      real(rk), allocatable, private :: c(:, :)
         !! array of constants for uniform grid
      logical, private :: uniform_grid
         !! flag for uniform grid
      real(rk), allocatable :: cnu(:, :, :)
         !! array of constants for a non-uniform grid
   contains
      procedure, pass(self) :: reconstruct => weno_reconstruct
      procedure, pass(self), private :: calc_cnu => weno_calc_cnu
   end type weno

   interface weno
      module procedure :: weno_init
   end interface weno

contains

   pure type(weno) function weno_init(ncells, k, eps, xedges) result(self)
   !! Initialize 'weno' object.
   !!
   !! @note
   !! If the grid is not uniform, the user must supply the edges of the grid.
      integer, intent(in) :: ncells
         !! number of cells
      integer, intent(in) :: k
         !! order of reconstruction within the cell (1 <= k <= 3)
      real(rk), intent(in), optional :: eps
         !! numerical smoothing factor
      real(rk), intent(in), optional :: xedges(0:)
         !! vector(0:ncells) of cell edges;
         !! xedges(i) is the value of x at right boundary of cell i, \( x_{i+1/2} \);
         !! xedges(i-1) is the value of x at left boundary of cell i, \( x_{i-1/2} \).

      ! Check inputs
      if (ncells > 0) then
         self%ncells = ncells
      else
         self%msg = "Invalid input 'ncells'. Valid range: ncells > 0."
         self%ierr = 1
         error stop self%msg
      end if

      if ((k >= 1) .and. (k <= 3)) then
         self%k = k
      else
         self%msg = "Invalid input 'k'. Valid range: 1 <= k <= 3."
         self%ierr = 1
         error stop self%msg
      end if

      if (present(eps)) then
         if (eps > 0._rk) then
            self%eps = eps
         else
            self%msg = "Invalid input 'eps'. Valid range: eps > 0."
            self%ierr = 1
            error stop self%msg
         end if
      end if

      if (present(xedges)) then
         if (size(xedges) == ncells + 1) then
            self%uniform_grid = .false.
            allocate (self%cnu(0:k - 1, -1:k - 1, 1:ncells))
            call self%calc_cnu(xedges)
         else
            self%msg = "Invalid input 'xedges': size(xedges) /= ncells + 1."
            self%ierr = 1
            error stop self%msg
         end if
      else
         self%uniform_grid = .true.
      end if

      ! Select constant parameters according to order of the method
      select case (k)
      case (1)
         self%d = d1
         self%c = c1
      case (2)
         self%d = d2
         self%c = c2
      case (3)
         self%d = d3
         self%c = c3
      end select

   end function weno_init

   pure subroutine weno_reconstruct(self, v, vl, vr)
   !!   This subroutine implements the (2k-1)th order WENO method for *arbitrary* (uniform or
   !! non-uniform) finite volume/difference schemes described in ICASE 97-65 (Shu, 1997).
   !!   The method is applicable to scalar as well as multicomponent problems. In the later
   !! case, the reconstruction is applied in a component by component fashion.
   !!   The scheme below depics a generic finite *volume* discretization and the notation
   !! used (see Arguments).
   !!```
   !! --0--|--1--| ...  |--(i-1)--|-------(i)-------|--(i+1)--| .... |--nc--|--(ncells+1)---
   !!                              ^               ^
   !!                              vl(i)           vr(i)
   !!                              v_{i-1/2}^+     v_{i+1/2}^-
   !!```
   !!   The procedure can equally be used for finite difference methods. In that case, 'v'
   !! is not the average cell value, but rather the flux! See section 2.3.2, page 22.
   !!
   !! @note
   !!   For a scalar 1D problem, this procedure is called once per time step. In contrast,
   !! for a scalar 2D problem, it is called (ncells1+ncells2) times per step. So, efficiency
   !! is very important. The current implementation is rather general in terms of order and
   !! grid type, but at the cost of a number of 'select case' constructs. I wonder if it
   !! would be wise to make a specific version for k=2 (3rd order) and uniform grids to get
   !! maximum performance for multi-dimensional problems.
      class(weno), intent(in) :: self
         !! object
      real(rk), intent(in) :: v(:)
         !! vector(ncells) of *average* cell values (if finite volume)
      real(rk), intent(out) :: vl(:)
         !! vector(ncells) of reconstructed v at left boundary of cell i, \(v_{i-1/2}^+\)
      real(rk), intent(out) :: vr(:)
         !! vector(ncells) of reconstructed v at right boundary of cell i, \(v_{i+1/2}^-\)

      real(rk), dimension(0:self%k - 1) :: vlr, vrr, w, wtilde, alfa, alfatilde, beta
      real(rk), dimension(0:self%k - 1, -1:self%k - 1) :: ci
      real(rk), dimension(1 - (self%k - 1):self%ncells + (self%k - 1)) :: vext
      integer :: i, r

      ! Algorithm
      ! Obtain the 'k' reconstructed values vi+1/2(r) & vi-1/2(r)
      associate (k => self%k, nc => self%ncells, eps => self%eps, d => self%d, c => self%c, &
                 cnu => self%cnu)
         ci = c
         vext(1:nc) = v
         vext(:0) = v(1)
         vext(nc + 1:) = v(nc)
         do concurrent(i=1:nc)

            ! Equations 2.10, 2.51
            if (.not. (self%uniform_grid)) ci = cnu(:, :, i)
            do concurrent(r=0:k - 1)
               vrr(r) = sum(ci(:, r)*vext(i - r:i - r + k - 1))
               vlr(r) = sum(ci(:, r - 1)*vext(i - r:i - r + k - 1))
            end do

            select case (k)

            case (1)
               beta(0) = 0._rk

               ! Equation 2.62
            case (2)
               beta(0) = (vext(i + 1) - vext(i))**2
               beta(1) = (vext(i) - vext(i - 1))**2

               ! Equation 2.63
            case (3)
               beta(0) = 13._rk/12*(vext(i) - 2*vext(i + 1) + vext(i + 2))**2 &
                         + 1._rk/4*(3*vext(i) - 4*vext(i + 1) + vext(i + 2))**2

               beta(1) = 13._rk/12*(vext(i - 1) - 2*vext(i) + vext(i + 1))**2 &
                         + 1._rk/4*(vext(i - 1) - vext(i + 1))**2

               beta(2) = 13._rk/12*(vext(i - 2) - 2*vext(i - 1) + vext(i))**2 &
                         + 1._rk/4*(vext(i - 2) - 4*vext(i - 1) + 3*vext(i))**2

            end select

            ! Equations 2.58-2.59 and procedure 2.2-4
            alfa = d/(eps + beta)**2
            alfatilde = d(k - 1:0:-1)/(eps + beta)**2
            w = alfa/sum(alfa)
            wtilde = alfatilde/sum(alfatilde)

            ! Procedure 2.2-5
            vr(i) = sum(w*vrr)
            vl(i) = sum(wtilde*vlr)

         end do
      end associate

   end subroutine weno_reconstruct

   pure subroutine weno_calc_cnu(self, xedges)
   !!   This subroutine computes the array of constants 'c(j,r,i)' required to use 'weno'
   !! with non-uniform grids.
   !!
   !! @note
   !!   This procedure is only called a very small of times (as many as the number of spatial
   !! dimensions) at the *start* of the simulation. So, there is no point in doing complex
   !! optimizations.
      class(weno), intent(inout) :: self
         !! object
      real(rk), intent(in) :: xedges(0:)
         !! vector(0:nc) of cell edges

      real(rk) :: prod1, prod2, sum1, sum2, dx
      real(rk), target :: xext(0 - (self%k + 1):self%ncells + (self%k + 1))
      real(rk), dimension(:), pointer :: xl, xr
      integer :: i, j, l, ng, m, q, r

      associate (k => self%k, nc => self%ncells, cnu => self%cnu)

         ! Allocate extended grid with (k+1) ghost cells on each side
         ng = k + 1
         xext(0:nc) = xedges

         ! Extend to the left linearly
         dx = xext(1) - xext(0)
         do i = -1, lbound(xext, 1), -1
            xext(i) = xext(i + 1) - dx
         end do

         ! Extend to the right linearly
         dx = xext(nc) - xext(nc - 1)
         do i = nc + 1, ubound(xext, 1)
            xext(i) = xext(i - 1) + dx
         end do

         ! Get pointers to left and right cell boundaries
         xl(1 - ng:) => xext(lbound(xext, 1):ubound(xext, 1) - 1)
         xr(1 - ng:) => xext(lbound(xext, 1) + 1:ubound(xext, 1))

         ! Compute array of constants 'c' for each grid position
         ! Equation 2.20, page 6.
         do concurrent(i=1:nc, r=-1:k - 1, j=0:k - 1)

            sum2 = 0._rk
            do m = j + 1, k

               prod2 = 1._rk
               do l = 0, k
                  if (l == m) cycle
                  prod2 = prod2*(xl(i - r + m) - xl(i - r + l))
               end do

               sum1 = 0._rk
               do l = 0, k

                  if (l == m) cycle

                  prod1 = 1._rk
                  do q = 0, k
                     if (q == m .or. q == l) cycle
                     prod1 = prod1*(xr(i) - xl(i - r + q))
                  end do

                  sum1 = sum1 + prod1

               end do

               sum2 = sum2 + sum1/prod2

            end do

            cnu(j, r, i) = sum2*(xr(i - r + j) - xl(i - r + j))

         end do
      end associate
   end subroutine weno_calc_cnu

end module hrweno
