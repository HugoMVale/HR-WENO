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

   pure subroutine grid1_linear(self, xmin, xmax, ncells)
    !! Constructor linear grid
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: ncells
        !! number of grid cells

      real(rk) :: xedges(0:ncells), rx
      integer :: i

      ! Check input
      if (.not. (xmax > xmin)) then
         error stop "Invalid input 'xmin', 'xmax'. Valid range: xmax > xmin."
      end if
      if (ncells < 1) then
         error stop "Invalid input 'ncells'. Valid range: ncells > 1."
      end if

      ! Clear grid, if we are trying to reallocate
      call self%clear

      ! Compute mesh
      rx = (xmax - xmin)/ncells
      do concurrent(i=0:ncells)
         xedges(i) = xmin + rx*i
      end do

      self%scl = 1
      call self%compute(xedges)

   end subroutine grid1_linear

   pure subroutine grid1_bilinear(self, xmin, xcross, xmax, ncells)
   !! Constructor bilinear grid
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xcross
        !! cross-over boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: ncells(2)
        !! number of grid cells in range [xmin, xcross] and [xcross, xmax]

      real(rk) :: xedges(0:sum(ncells)), rx
      integer :: i

      if (.not. (xcross > xmin)) then
         error stop "Invalid input 'xmin', 'xcross'. Valid range: xcross > xmin."
      end if
      if (.not. (xmax > xcross)) then
         error stop "Invalid input 'xcross', 'xmax'. Valid range: xmax > xcross."
      end if
      if (any(ncells < 1)) then
         error stop "Invalid input 'ncells'. Valid range: ncells(i) >= 1."
      end if

      ! Compute mesh [xmin, xcross]
      rx = (xcross - xmin)/ncells(1)
      do concurrent(i=0:ncells(1))
         xedges(i) = xmin + rx*i
      end do

      ! Compute mesh [xcross, xmax]
      rx = (xmax - xcross)/ncells(2)
      do concurrent(i=1:ncells(2))
         xedges(ncells(1) + i) = xcross + rx*i
      end do

      self%scl = 2
      call self%compute(xedges)

   end subroutine grid1_bilinear

   pure subroutine grid1_log(self, xmin, xmax, ncells)
   !! Constructor log grid
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: ncells
        !! number of grid cells

      real(rk) :: xedges(0:ncells), rx
      integer :: i

      ! Check input
      if (xmin <= 0._rk) then
         error stop "Invalid input 'xmin'. Valid range: xmin > 0."
      end if
      if (.not. (xmax > xmin)) then
         error stop "Invalid input 'xmin', 'xmax'. Valid range: xmax > xmin."
      end if
      if (ncells < 1) then
         error stop "Invalid input 'ncells'. Valid range: ncells > 1."
      end if

      ! Clear grid, if we are trying to reallocate
      call self%clear

      ! Compute mesh
      rx = (log(xmax/xmin))/ncells
      do concurrent(i=0:ncells)
         xedges(i) = log(xmin) + rx*i
      end do
      xedges = exp(xedges)

      self%scl = 3
      call self%compute(xedges)

   end subroutine grid1_log

   pure subroutine grid1_geometric(self, xmin, xmax, ratio, ncells)
   !! Constructor geometric grid
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      real(rk), intent(in) :: ratio
        !! constant ratio of geometric grid
      integer, intent(in) :: ncells
        !! number of grid cells

      real(rk) :: xedges(0:ncells), sum_ratio, a1
      integer :: i

      ! Check input
      if (.not. (xmax > xmin)) then
         error stop "Invalid input 'xmin', 'xmax'. Valid range: xmax > xmin"
      end if
      if (ratio <= 0._rk) then
         error stop "Invalid input 'ratio'. Valid range: ratio > 0"
      end if
      if (ncells < 1) then
         error stop "Invalid input 'ncells'. Valid range: ncells > 1"
      end if

      ! Clear grid, if we are trying to reallocate
      call self%clear

      ! Compute mesh
      sum_ratio = 0._rk
      do i = 1, ncells
         sum_ratio = sum_ratio + ratio**(i - 1)
      end do
      a1 = (xmax - xmin)/sum_ratio
      xedges(0) = xmin
      do i = 1, ncells
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
