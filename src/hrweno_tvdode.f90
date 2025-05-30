module hrweno_tvdode
!!   This module contains two total variation diminishing (TVD) high-order schemes for solving
!! initial value problems. It is very important to use TVD schemes for time integration. Even
!! with a TVD spacial discretization, if the time discretization is done by a non-TVD method,
!! the result may be oscillatory.
!!   Source: ICASE 97-65 by Shu, 1997.
   use hrweno_kinds, only: rk
   use stdlib_optval, only: optval
   implicit none
   private

   public :: rktvd, mstvd

   type, abstract :: tvdode
   !! Abstract class for TVD ODE solvers.
      procedure(integrand), pointer, nopass, private :: fu => null()
         !! subroutine with the derivative \( u'(t,u) \)
      integer :: neq
         !! number of equations
      integer :: order
         !! order of the method
      integer :: fevals = 0
         !! number of function evaluations
      integer :: istate = 0
         !! flag indicating the state of the integration:
         !! 1 first call for a problem,
         !! 2 subsequent call for a problem.
      character(:), allocatable:: msg
         !! error message
      real(rk), allocatable, private :: ui(:)
      real(rk), allocatable, private :: udot(:)
   contains
      procedure, pass(self) :: error_msg
   end type tvdode

   type, extends(tvdode) :: rktvd
   !! Runge-Kutta TVD ODE solver class.
   contains
      procedure, pass(self) :: integrate => rktvd_integrate
   end type rktvd

   type, extends(tvdode) :: mstvd
   !! Multi-step TVD ODE solver class.
      real(rk), allocatable, private :: uold(:, :)
      real(rk), allocatable, private :: udotold(:, :)
   contains
      procedure, pass(self) :: integrate => mstvd_integrate
   end type mstvd

   abstract interface
      subroutine integrand(t, u, udot)
      !!  Integrand for `tvdode` class
         import :: rk
         real(rk), intent(in) :: t, u(:)
         real(rk), intent(out) :: udot(:)
      end subroutine
   end interface

   interface rktvd
      module procedure :: rktvd_init
   end interface rktvd

   interface mstvd
      module procedure :: mstvd_init
   end interface mstvd

contains

   type(rktvd) function rktvd_init(fu, neq, order) result(self)
   !! Initialize `rktvd` object.
      procedure(integrand) :: fu
         !! subroutine with the derivative \( u'(t,u) \)
      integer, intent(in) :: neq
         !! number of equations
      integer, intent(in) :: order
         !! order of the method (1, 2 or 3)

      self%fu => fu

      if (neq > 0) then
         self%neq = neq
      else
         call self%error_msg("Invalid input 'neq'. Valid range: neq >= 1.")
      end if

      if ((order >= 1) .and. (order <= 3)) then
         self%order = order
      else
         call self%error_msg("Invalid input 'order' in 'rktvd'. Valid range: 1 <= k <= 3.")
      end if

      allocate (self%ui(self%neq), self%udot(self%neq))
      self%istate = 1

   end function rktvd_init

   subroutine rktvd_integrate(self, u, t, tout, dt, itask)
   !!   This subroutine implements the optimal 1st, 2nd and 3rd order TVD RK methods described
   !! in ICASE 97-65 (Shu, 1997).
   !!   The routine was built to work similarly to LSODE.
   !!
   !! @note
   !!   There are also 4th and 5th order methods, but they have lower CFL coefficients and
   !! are more difficult to implement. See Equation 4.15, page 44.
   !!
   !! @todo
   !! - Adjust dt in final step to avoid overshoting tout by some fraction of dt.
      class(rktvd), intent(inout) :: self
         !! object
      real(rk), intent(inout) :: u(:)
         !! vector(neq) with the variables to integrate \( u(t) \)
      real(rk), intent(inout) :: t
         !! time; on return it will be the current value of \(t\) (close to tout)
      real(rk), intent(in) :: tout
         !! time where next output is desired
      real(rk), intent(in) :: dt
         !! time step
      integer, intent(in), optional :: itask
         !! flag indicating the task to be performed:
         !! 1 normal integration until tout;
         !! 2 single `dt` step.

      integer :: itask_

      ! Check input conditions
      if (self%istate < 1) return
      if (is_done(t, tout, dt)) return
      itask_ = optval(itask, 1)

      ! Algorthm selection
      associate (ui => self%ui, udot => self%udot)
         select case (self%order)

            ! ------------------------------- 1st-order RK (Euler) ----------------------------
            ! Equation (4.10), page 43.
         case (1)
            do
               call self%fu(t, u, udot)
               u = u + dt*udot
               t = t + dt
               self%fevals = self%fevals + 1
               if (is_done(t, tout, dt) .or. itask_ == 2) exit
            end do

            ! --------------------------------- 2nd-order RK ----------------------------------
            ! Equation (4.10), page 43.
         case (2)
            do
               call self%fu(t, u, udot)
               ui = u + dt*udot
               call self%fu(t + dt, ui, udot)
               u = (u + ui + dt*udot)/2
               t = t + dt
               self%fevals = self%fevals + 2
               if (is_done(t, tout, dt) .or. itask_ == 2) exit
            end do

            ! --------------------------------- 3rd-order RK -------------------------------------
            ! Equation (4.11), page 43.
         case (3)
            do
               call self%fu(t, u, udot)
               ui = u + dt*udot
               call self%fu(t + dt, ui, udot)
               ui = (3*u + ui + dt*udot)/4
               call self%fu(t + dt/2, ui, udot)
               u = (u + 2*ui + 2*dt*udot)/3
               t = t + dt
               self%fevals = self%fevals + 3
               if (is_done(t, tout, dt) .or. itask_ == 2) exit
            end do

         end select
      end associate

      if (self%istate == 1) self%istate = 2

   end subroutine rktvd_integrate

   type(mstvd) function mstvd_init(fu, neq) result(self)
   !! Initialize `mstvd` object.
      procedure(integrand) :: fu
         !! subroutine with the derivative \( u'(t,u) \)
      integer, intent(in) :: neq
         !! number of equations

      self%fu => fu

      if (neq > 0) then
         self%neq = neq
      else
         call self%error_msg("Invalid input 'neq'. Valid range: neq >= 1.")
      end if

      self%order = 3

      allocate (self%ui(self%neq), self%udot(self%neq), &
                self%uold(self%neq, self%order + 1), self%udotold(self%neq, self%order + 1))
      self%istate = 1

   end function mstvd_init

   subroutine mstvd_integrate(self, u, t, tout, dt)
   !!   This subroutine implements a 5-step, 3rd order TVD multi-step method described
   !! in ICASE 97-65 (Shu, 1997). In theory, this method should have an efficiency 1.5 times
   !! higher than the RK method of the same order. However, in practice they appear to be
   !! almost identical.
   !!   The routine was built to work similarly to LSODE.
   !!
   !! @note
   !!   There is a 2nd order multi-step method, but the corresponding CFL value is half that
   !! of the 2nd order RK method. Thus, there is no reason to implement it.
      class(mstvd), intent(inout) :: self
         !! object
      real(rk), intent(inout) :: u(:)
         !! vector(neq) with the variables to integrate \( u(t) \)
      real(rk), intent(inout) :: t
         !! time; on return it will be the current value of \(t\) (close to tout)
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

         ode_start = rktvd(self%fu, self%neq, self%order)

         do i = (self%order + 1), 1, -1
            uold(:, i) = u
            call self%fu(t, u, udotold(:, i))
            call ode_start%integrate(u, t, t + 2*dt, dt, itask=2)
         end do

         self%fevals = ode_start%fevals
         self%istate = 2

      end if

      ! Equation (4.26), page 48.
      do

         if (is_done(t, tout, dt)) exit

         call self%fu(t, u, udot)
         ui = (25*u + 50*dt*udot + 7*uold(:, 4) + 10*dt*udotold(:, 4))/32
         t = t + dt
         self%fevals = self%fevals + 1

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

      is_done = (t - tout)*sign(1.0_rk, dt) > 0.0_rk

   end function is_done

   pure subroutine error_msg(self, msg)
   !! Error method.
      class(tvdode), intent(inout) :: self
         !! object
      character(*), intent(in) :: msg
         !! message

      self%msg = msg
      self%istate = -1
      error stop self%msg

   end subroutine

end module hrweno_tvdode
