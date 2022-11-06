module grids
!!   This module implements a convenience 1D grid class. The module is mostly indented to help
!! write the examples. The WENO schemes themselves are completely independent from this module,
!! so that they can be used anywhere.
   use hrweno_kinds, only: rk
   use stdlib_optval, only: optval
   implicit none
   private

   public :: grid1

   type :: grid1
    !! 1D grid class.
      character(:), allocatable :: name
        !! variable name
      character(:), allocatable :: scl
        !! scale type
      integer :: ncells
        !! number of cells
      real(rk), allocatable :: edges(:)
        !! vector(0:ncells) of cell edges
      real(rk), allocatable :: center(:)
        !! vector(ncells) of cell centers, \( x_i \)
      real(rk), allocatable :: width(:)
        !! vector(ncells) of cell widths,  \( x_{i+1/2} - x_{i-1/2} \)
      real(rk), dimension(:), pointer :: left => null()
        !! vector(ncells) of left cell boundaries, \( x_{i-1/2} \)
      real(rk), dimension(:), pointer :: right => null()
        !! vector(ncells) of right cell boundaries, \( x_{i+1/2} \)
   contains
      procedure, pass(self) :: linear => grid1_linear
      procedure, pass(self) :: bilinear => grid1_bilinear
      procedure, pass(self) :: log => grid1_log
      procedure, pass(self) :: geometric => grid1_geometric
      procedure, pass(self), private :: clear => grid1_clear
      procedure, pass(self), private :: compute => grid1_compute
   end type grid1

contains

   pure subroutine grid1_linear(self, xmin, xmax, ncells, name)
    !! Initialize *linear* grid. <br>
    !! Constant width: \( width(i) = width(i+1) \)
    !!
    !!```
    !!     |  ...  |------(i)------|------(i+1)------|  ...  |
    !!    xmin      <-  width(i) -> <- width(i+1)  ->       xmax
    !!```
    !!
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: ncells
        !! number of grid cells
      character(*), intent(in), optional :: name
        !! grid name

      real(rk) :: xedges(0:ncells), rx
      integer :: i

      ! Check input
      if (xmax <= xmin) then
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

      self%scl = "linear"
      call self%compute(xedges, optval(name, ""))

   end subroutine grid1_linear

   pure subroutine grid1_bilinear(self, xmin, xcross, xmax, ncells, name)
    !! Initialize *bilinear* grid.
    !! Equivalent to 2 linear grids in series.
    !!
    !!```
    !!     |  ...  |--|--|--| ... | ... |---|---|---|  ...  |
    !!    xmin                  xcross                     xmax
    !!```
    !!
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xcross
        !! cross-over boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: ncells(2)
        !! number of grid cells in range [xmin, xcross] and [xcross, xmax]
      character(*), intent(in), optional :: name
        !! grid name

      real(rk) :: xedges(0:sum(ncells)), rx
      integer :: i

      if (xcross <= xmin) then
         error stop "Invalid input 'xmin', 'xcross'. Valid range: xcross > xmin."
      end if
      if (xmax <= xcross) then
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

      self%scl = "bilinear"
      call self%compute(xedges, optval(name, ""))

   end subroutine grid1_bilinear

   pure subroutine grid1_log(self, xmin, xmax, ncells, name)
   !! Initialize *logarithmic* grid. <br>
   !! Equivalent to a linear grid in terms of \( y=\log(x) \).
   !!
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain (xmin>0)
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      integer, intent(in) :: ncells
        !! number of grid cells
      character(*), intent(in), optional :: name
        !! grid name

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

      self%scl = "log"
      call self%compute(xedges, optval(name, ""))

   end subroutine grid1_log

   subroutine grid1_geometric(self, xmin, xmax, ratio, ncells, name)
    !! Initialize *geometric* grid. <br>
    !! Geometrically increasing/decreasing width: \( width(i+1) = R \; width(i) \)
    !!
    !!```
    !!     |  ...  |------(i)------|------(i+1)------|  ...  |
    !!    xmin      <-  width(i) -> <- width(i+1)  ->       xmax
    !!```
    !!
      class(grid1), intent(inout) :: self
        !! object
      real(rk), intent(in) :: xmin
        !! lower boundary of grid domain
      real(rk), intent(in) :: xmax
        !! upper boundary of grid domain
      real(rk), intent(in) :: ratio
        !! constant ratio \(R\) of geometric grid (R>0)
      integer, intent(in) :: ncells
        !! number of grid cells
      character(*), intent(in), optional :: name
        !! grid name

      real(rk) :: xedges(0:ncells), a
      integer :: i

      ! Check input
      if (xmax <= xmin) then
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
      a = (xmax - xmin)/(ratio**ncells - 1._rk)
      do concurrent(i=0:ncells)
         xedges(i) = xmin + a*(ratio**i - 1)
      end do

      self%scl = "geometric"
      call self%compute(xedges, optval(name, ""))

   end subroutine grid1_geometric

   pure subroutine grid1_compute(self, xedges, name)
   !! Auxiliar procedure to compute/assign grid components.
      class(grid1), intent(inout), target :: self
        !! object
      real(rk), intent(in) :: xedges(0:)
        !! grid edges
      character(*), intent(in), optional :: name
        !! grid name

      ! Map values to grid object
      self%ncells = ubound(xedges, 1)
      self%edges = xedges
      self%left => self%edges(0:self%ncells - 1)
      self%right => self%edges(1:self%ncells)
      self%center = (self%left + self%right)/2
      self%width = self%right - self%left
      self%name = name

   end subroutine grid1_compute

   pure subroutine grid1_clear(self)
   !! Clear grid object.
      class(grid1), intent(inout), target :: self
        !! object

      ! Reset all variables
      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%edges)) deallocate (self%edges)
      if (allocated(self%center)) deallocate (self%center)
      if (allocated(self%width)) deallocate (self%width)
      if (associated(self%left)) nullify (self%left)
      if (associated(self%right)) nullify (self%right)
      self%ncells = 0
      self%scl = ""

   end subroutine grid1_clear

end module grids
