module grid
!!   This module implements a 1D grid class.
   use, intrinsic :: iso_fortran_env, only: real64
   use stdlib_optval, only: optval
   implicit none
   private

   public :: grid1

   integer, parameter :: rk = real64

   type :: grid1
    !! 1D grid
      character(20) :: name = "x [-]"
        !! variable name
      real(rk), allocatable :: edges(:)
        !! vector(0:nc) of cell edges
      real(rk), allocatable :: center(:)
        !! vector(nc) of cell centers, \( x_i \)
      real(rk), allocatable :: width(:)
        !! vector(nc) of cell widths,  \( x_{i+1/2} - x_{i-1/2} \)
      real(rk), dimension(:), pointer :: left
        !! vector(nc) of left cell boundaries, \( x_{i-1/2} \)
      real(rk), dimension(:), pointer :: right
        !! vector(nc) of right cell boundaries, , \( x_{i+1/2} \)
      integer :: ncells
        !! number of cells
   contains
      procedure, pass(self) :: new
   end type grid1

contains

   pure subroutine new(self, xmin, xmax, nc, scl)
      !! Constructor grid1
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: nc
        !! number of grid cells
      integer, intent(in), optional :: scl
        !! grid scale: linear (1), logarithmic (2)

      real(rk) :: xedges(0:nc), rx
      integer :: i

      ! Compute mesh
      select case (optval(scl, 1))
      case (1)
         rx = (xmax - xmin)/nc
         do concurrent(i=0:nc)
            xedges(i) = xmin + rx*i
         end do
      case (2)
         rx = (log(xmax/xmin))/nc
         do concurrent(i=0:nc)
            xedges(i) = log(xmin) + rx*i
         end do
         xedges = exp(xedges)
      end select

      ! Map values to grid object
      self%ncells = nc
      self%edges = xedges
      self%left => self%edges(0:nc - 1)
      self%right => self%edges(1:nc)
      self%center = (self%left + self%right)/2
      self%width = self%right - self%left

   end subroutine new

end module grid
