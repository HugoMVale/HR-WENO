module tvdode
!!   This module contains two total variation diminishing (TVD) high-order schemes for solving
!! initial value problems. It is very important to use TVD schemes for time integration. Even
!! with a TVD spacial discretization, if the time discretization is done by a non-TVD method,
!! the result may be oscillatory.
!!   Source: ICASE 97-65 by Shu, 1997.
   use, intrinsic :: iso_fortran_env, only: real64
   use stdlib_optval, only: optval
   implicit none
   private

   public :: tvdode_class, rktvd, mstvd

   integer, parameter :: rk = real64

   type :: tvdode_class
   !! Abstract class for tvode-like types
      procedure(integrand), pointer, private :: fu => null()
         !! subroutine with the derivative u'(t)
      integer, private :: neq = 0
         !! number of equations
      integer, private :: order = 0
         !! order of the method
      integer :: istate = 0
         !! flag indicating the state of the integration:
         !! 1 first call for a problem,
         !! 2 subsequent call for a problem.
      character(:), allocatable:: msg
      !! error message
      real(rk), allocatable, private :: ui(:)
      real(rk), allocatable, private :: udot(:)
   end type

   type, extends(tvdode_class) :: rktvd
   !!  'rktvd' class
   contains
      procedure, pass(self) :: init => rktvd_init
      procedure, pass(self) :: integrate => rktvd_integrate
   end type

   type, extends(tvdode_class) :: mstvd
   !!  'mstvd' class
      real(rk), allocatable, private :: uold(:, :)
      real(rk), allocatable, private :: udotold(:, :)
   contains
      procedure, pass(self) :: init => mstvd_init
      procedure, pass(self) :: integrate => mstvd_integrate
   end type

   abstract interface
      pure subroutine integrand(self, t, u, udot)
      !!  Integrand for 'tvdode_class'
         import :: rk, tvdode_class
         class(tvdode_class), intent(inout) :: self
         real(rk), intent(in) :: t, u(:)
         real(rk), intent(out) :: udot(:)
      end subroutine
   end interface

contains

   pure subroutine rktvd_init(self, fu, neq, order)
   !! Initialize 'rktvd' object.
      class(rktvd), intent(inout) :: self
         !! object
      procedure(integrand) :: fu
         !! subroutine with the derivative u'(t)
      integer, intent(in) :: neq
         !! number of equations
      integer, intent(in) :: order
         !! order of the method (1, 2 or 3)

      ! Clear object if required
      if (allocated(self%msg)) deallocate (self%msg)
      if (allocated(self%ui)) deallocate (self%ui)
      if (allocated(self%udot)) deallocate (self%udot)

      self%fu => fu

      if (neq > 0) then
         self%neq = neq
      else
         self%msg = "Invalid input 'neq'. Valid range: neq >= 1."
         self%istate = -1
         error stop self%msg
      end if

      if ((order >= 1) .and. (order <= 3)) then
         self%order = order
      else
         self%msg = "Invalid input 'order' in 'rktvd'. Valid range: 1 <= k <= 3."
         self%istate = -1
         error stop self%msg
      end if

      allocate (self%ui(self%neq), self%udot(self%neq))
      self%istate = 1

   end subroutine rktvd_init

   pure subroutine rktvd_integrate(self, u, t, tout, dt, itask)
   !!   This subroutine implements the optimal 1st, 2nd and 3rd order TVD RK methods described
   !! in ICASE 97-65 (Shu, 1997).
   !!   The routine was built to work similarly to LSODE.
   !!
   !! @note
   !!   There are also 4th and 5th order methods, but they have lower CFL coeffiecients and
   !! are more difficult to implement. See Equation 4.15, page 44.
   !!
   !! @todo
   !! - Adjust dt in final step to avoid overshoting tout by some fraction of dt.
   !! - Maybe include an optional work array that could be transfered to fu.
      class(rktvd), intent(inout) :: self
         !! object
      real(rk), intent(inout) :: u(:)
         !! vector(neq) with the variables to integrate u(t)
      real(rk), intent(inout) :: t
         !! time; on return it will be the current value of t (close to tout)
      real(rk), intent(in) :: tout
         !! time where next output is desired
      real(rk), intent(in) :: dt
         !! time step
      integer, intent(in), optional :: itask
         !! flag indicating the task to be performed:
         !! 1 normal integration until tout;
         !! 2 single dt step.

      integer :: itask_

      ! Check input conditions
      if (self%istate < 1) return
      if (is_done(t, tout, dt)) return
      itask_ = optval(itask, 1)

      ! Algorthm selection
      associate (ui => self%ui, udot => self%udot)
         select case (self%order)

            ! ------------------------------- 1st-order RK (Euler) -------------------------------
            ! Equation (4.10), page 43.
         case (1)
            do
               call self%fu(t, u, udot)
               u = u + dt*udot
               t = t + dt
               if (is_done(t, tout, dt) .or. itask_ == 2) exit
            end do

            ! --------------------------------- 2nd-order RK -------------------------------------
            ! Equation (4.10), page 43.
         case (2)
            do
               call self%fu(t, u, udot)
               ui = u + dt*udot
               call self%fu(t, ui, udot)
               u = (u + ui + dt*udot)/2
               t = t + dt
               if (is_done(t, tout, dt) .or. itask_ == 2) exit
            end do

            ! --------------------------------- 3rd-order RK -------------------------------------
            ! Equation (4.11), page 43.
         case (3)
            do
               call self%fu(t, u, udot)
               ui = u + dt*udot
               call self%fu(t, ui, udot)
               ui = (3*u + ui + dt*udot)/4
               call self%fu(t, ui, udot)
               u = (u + 2*ui + 2*dt*udot)/3
               t = t + dt
               if (is_done(t, tout, dt) .or. itask_ == 2) exit
            end do

         end select
      end associate

      if (self%istate == 1) self%istate = 2

   end subroutine rktvd_integrate

   pure subroutine mstvd_init(self, fu, neq)
   !! Initialize 'mstvd' object.
      class(mstvd), intent(inout) :: self
         !! object
      procedure(integrand) :: fu
         !! subroutine with the derivative u'(t)
      integer, intent(in) :: neq
         !! number of equations

      ! Clear object if required
      if (allocated(self%msg)) deallocate (self%msg)
      if (allocated(self%ui)) deallocate (self%ui)
      if (allocated(self%udot)) deallocate (self%udot)
      if (allocated(self%uold)) deallocate (self%uold)
      if (allocated(self%udotold)) deallocate (self%udotold)

      self%fu => fu

      if (neq > 0) then
         self%neq = neq
      else
         self%msg = "Invalid input 'neq'. Valid range: neq >= 1."
         self%istate = -1
         error stop self%msg
      end if

      self%order = 3

      allocate (self%ui(self%neq), self%udot(self%neq), &
                self%uold(self%neq, self%order + 1), self%udotold(self%neq, self%order + 1))
      self%istate = 1

   end subroutine mstvd_init

   pure subroutine mstvd_integrate(self, u, t, tout, dt)
   !!   This subroutine implements a 5-step, 3rd order TVD multi-step method described
   !! in ICASE 97-65 (Shu, 1997). In theory, this method should have an efficiency 1.5 times
   !! higher than the RK method of the same order. However, in practice they appear to be
   !! almost identical.
   !!   The routine was built to work similarly to LSODE.
   !!
   !! @note
   !!   There is a 2nd order multi-step method, but the corresponding CFL value is half that
   !! of the 2nd order RK method. Thus, thre is no reason to implement it.
   !!
   !! @todo
   !! - Maybe include an optional work array that could be transfered to 'fu'.
      class(mstvd), intent(inout) :: self
         !! object
      real(rk), intent(inout) :: u(:)
         !! vector(neq) with the variables to integrate u(t)
      real(rk), intent(inout) :: t
         !! time; on return it will be the current value of t (close to tout)
      real(rk), intent(in) :: tout
         !! time where next output is desired
      real(rk), intent(in) :: dt
         !! time step

      type(rktvd) :: ode_start
      integer :: i

      ! Check input conditions
      if (self%istate < 1) return
      if (is_done(t, tout, dt)) return

      ! The first starting values must be computed with a single-step method: we chose
      ! the RK method of the same order.
      ! The factor 2 in 't+2*dt' is not important, it just needs to be larger than 1.0
      ! so that one full 'dt' step can be computed.
      associate (ui => self%ui, udot => self%udot, uold => self%uold, udotold => self%udotold)
      if (self%istate == 1) then

         call ode_start%init(self%fu, self%neq, self%order)

         do i = (self%order + 1), 1, -1
            uold(:, i) = u
            call self%fu(t, u, udotold(:, i))
            call ode_start%integrate(u, t, t + 2*dt, dt, itask=2)
         end do

         self%istate = 2

      end if

      ! Equation (4.26), page 48.
      do

         if (is_done(t, tout, dt)) exit

         call self%fu(t, u, udot)
         ui = (25*u + 50*dt*udot + 7*uold(:, 4) + 10*dt*udotold(:, 4))/32
         t = t + dt

         ! Shift u and udot values one step into the past
         udotold = eoshift(udotold, shift=-1, dim=2)
         uold = eoshift(uold, shift=-1, dim=2)
         udotold(:, 1) = udot
         uold(:, 1) = u
         u = ui

      end do
      end associate

   end subroutine mstvd_integrate

   pure logical function is_done(t, tout, dt)
    !! Aux function to check if the integration is finished.
      real(rk), intent(in) :: t
        !! current time
      real(rk), intent(in) :: tout
        !! time where next output is desired
      real(rk), intent(in) :: dt
        !! time step

      is_done = (t - tout)*sign(1._rk, dt) > 0._rk

   end function is_done

end module tvdode
