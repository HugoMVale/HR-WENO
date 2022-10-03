module grid
!!   This module implements a 1D grid class.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   public :: grid1

   integer, parameter :: rk = real64

   type :: grid1
    !! 1D grid
      character(:), allocatable :: name
        !! variable name
      real(rk), allocatable :: edges(:)
        !! vector(0:nc) of cell edges
      real(rk), allocatable :: center(:)
        !! vector(nc) of cell centers, \( x_i \)
      real(rk), allocatable :: width(:)
        !! vector(nc) of cell widths,  \( x_{i+1/2} - x_{i-1/2} \)
      real(rk), dimension(:), pointer :: left => null()
        !! vector(nc) of left cell boundaries, \( x_{i-1/2} \)
      real(rk), dimension(:), pointer :: right => null()
        !! vector(nc) of right cell boundaries, , \( x_{i+1/2} \)
      integer :: ncells = 0
        !! number of cells
      integer :: scl = 0
        !! scale (1: linear, 2: bilinear, 3: log, 4: geometric)
   contains
      procedure, pass(self) :: linear => grid1_linear
      procedure, pass(self) :: bilinear => grid1_bilinear
      procedure, pass(self) :: log => grid1_log
      procedure, pass(self) :: geometric => grid1_geometric
      procedure, pass(self), private :: clear => grid1_clear
      procedure, pass(self), private :: compute => grid1_compute
   end type grid1

contains

   pure subroutine grid1_linear(self, xmin, xmax, nc)
    !! Constructor linear grid
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: nc
        !! number of grid cells

      real(rk) :: xedges(0:nc), rx
      integer :: i

      ! Check input
      if (.not. (xmax > xmin)) then
         error stop "Invalid input 'xmin', 'xmax'. Valid range: xmax > xmin."
      end if
      if (nc < 1) then
         error stop "Invalid input 'nc'. Valid range: nc > 1."
      end if

      ! Clear grid, if we are trying to reallocate
      call self%clear

      ! Compute mesh
      rx = (xmax - xmin)/nc
      do concurrent(i=0:nc)
         xedges(i) = xmin + rx*i
      end do

      self%scl = 1
      call self%compute(xedges)

   end subroutine grid1_linear

   pure subroutine grid1_bilinear(self, xmin, xcross, xmax, nc1, nc2)
   !! Constructor bilinear grid
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xcross
        !! cross-over boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: nc1
        !! number of grid cells in range [xmin, xcross]
      integer, intent(in) :: nc2
        !! number of grid cells in range [xcross, xmax]

      real(rk) :: xedges(0:nc1 + nc2), rx
      integer :: i

      if (.not. (xcross > xmin)) then
         error stop "Invalid input 'xmin', 'xcross'. Valid range: xcross > xmin."
      end if
      if (.not. (xmax > xcross)) then
         error stop "Invalid input 'xcross', 'xmax'. Valid range: xmax > xcross."
      end if
      if (nc1 < 1) then
         error stop "Invalid input 'nc1'. Valid range: nc1 > 1."
      end if
      if (nc2 < 1) then
         error stop "Invalid input 'nc2'. Valid range: nc2 > 1."
      end if

      ! Compute mesh [xmin, xcross]
      rx = (xcross - xmin)/nc1
      do concurrent(i=0:nc1)
         xedges(i) = xmin + rx*i
      end do

      ! Compute mesh [xcross, xmax]
      rx = (xmax - xcross)/nc2
      do concurrent(i=1:nc2)
         xedges(nc1 + i) = xcross + rx*i
      end do

      self%scl = 2
      call self%compute(xedges)

   end subroutine grid1_bilinear

   pure subroutine grid1_log(self, xmin, xmax, nc)
   !! Constructor log grid
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: nc
        !! number of grid cells

      real(rk) :: xedges(0:nc), rx
      integer :: i

      ! Check input
      if (xmin <= 0._rk) then
         error stop "Invalid input 'xmin'. Valid range: xmin > 0."
      end if
      if (.not. (xmax > xmin)) then
         error stop "Invalid input 'xmin', 'xmax'. Valid range: xmax > xmin."
      end if
      if (nc < 1) then
         error stop "Invalid input 'nc'. Valid range: nc > 1."
      end if

      ! Clear grid, if we are trying to reallocate
      call self%clear

      ! Compute mesh
      rx = (log(xmax/xmin))/nc
      do concurrent(i=0:nc)
         xedges(i) = log(xmin) + rx*i
      end do
      xedges = exp(xedges)

      self%scl = 3
      call self%compute(xedges)

   end subroutine grid1_log

   pure subroutine grid1_geometric(self, xmin, xmax, ratio, nc)
   !! Constructor geometric grid
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      real(rk), intent(in) :: ratio
        !! constant ratio of geometric grid
      integer, intent(in) :: nc
        !! number of grid cells

      real(rk) :: xedges(0:nc), sum_ratio, a1
      integer :: i

      ! Check input
      if (.not. (xmax > xmin)) then
         error stop "Invalid input 'xmin', 'xmax'. Valid range: xmax > xmin"
      end if
      if (ratio <= 0._rk) then
         error stop "Invalid input 'ratio'. Valid range: ratio > 0"
      end if
      if (nc < 1) then
         error stop "Invalid input 'nc'. Valid range: nc > 1"
      end if

      ! Clear grid, if we are trying to reallocate
      call self%clear

      ! Compute mesh
      sum_ratio = 0._rk
      do i = 1, nc
         sum_ratio = sum_ratio + ratio**(i - 1)
      end do
      a1 = (xmax - xmin)/sum_ratio
      xedges(0) = xmin
      do i = 1, nc
         xedges(i) = xedges(i - 1) + a1*ratio**(i - 1)
      end do

      self%scl = 4
      call self%compute(xedges)

   end subroutine grid1_geometric

   pure subroutine grid1_compute(self, xedges)
   !! Aux procedure to compute grid features
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xedges(0:)
        !! grid edges

      ! Map values to grid object
      self%ncells = ubound(xedges, 1)
      self%edges = xedges
      self%left => self%edges(0:self%ncells - 1)
      self%right => self%edges(1:self%ncells)
      self%center = (self%left + self%right)/2
      self%width = self%right - self%left

   end subroutine grid1_compute

   pure subroutine grid1_clear(self)
   !! Clear grid
      class(grid1), intent(inout), target :: self
        !! object

      ! Map values to grid object
      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%edges)) deallocate (self%edges)
      if (allocated(self%center)) deallocate (self%center)
      if (allocated(self%width)) deallocate (self%width)
      if (associated(self%left)) nullify (self%left)
      if (associated(self%right)) nullify (self%right)
      self%ncells = 0
      self%scl = 0

   end subroutine grid1_clear

end module grid
