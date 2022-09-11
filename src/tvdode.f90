module tvdode
!!   This module contains two total variation diminishing (TVD) high-order schemes for solving
!! initial value problems. It is very important to use TVD schemes for time integration. Even
!! with a TVD spacial discretization, if the time discretization is done by a non-TVD method,
!! the result may be oscillatory.
!!   Source: ICASE 97-65 by Shu, 1997.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   public :: rktvd, mstvd

   integer, parameter :: rk = real64

   abstract interface
      pure subroutine integrand(t, u, udot)
         import :: rk
         real(rk), intent(in) :: t, u(:)
         real(rk), intent(out) :: udot(:)
      end subroutine
   end interface

contains

   pure subroutine rktvd(fu, u, t, tout, dt, order, itask, istate)
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
      procedure(integrand) :: fu
        !! subroutine with the derivative u'(t)
      real(rk), intent(inout) :: u(:)
        !! vector(N) with the variables to integrate u(t)
      real(rk), intent(inout) :: t
        !! time; on return it will be the current value of t (close to tout)
      real(rk), intent(in) :: tout
        !! time where next output is desired
      real(rk), intent(in) :: dt
        !! time step
      integer, intent(in) :: order
        !! order of the method (1, 2 or 3)
      integer, intent(in) :: itask
        !! flag indicating the task to be performed:
        !! 1 normal integration until tout;
        !! 2 single dt step.
      integer, intent(inout) :: istate
        !! flag indicating the state of the integration:
        !! 1 first call for a problem,
        !! 2 subsequent call for a problem.
      real(rk), dimension(size(u)) :: ui, udot
      character(:), allocatable :: errmsg

      ! Check input conditions
      if (is_done(t, tout, dt)) return

      if (istate == 1) then
         if (order < 1 .or. order > 3) then
            errmsg = "Invalid input 'order' in 'rktvd'. Valid range: 1 <= k <= 3."
            error stop errmsg
         end if
         if (itask < 1 .or. itask > 2) then
            errmsg = "Invalid input 'itask' in 'rktvd'. Valid set: {1, 2}."
            error stop errmsg
         end if
         istate = 2
      else if (istate < 1 .or. istate > 2) then
         errmsg = "Invalid value 'istate' in 'rktvd'. Valid set: {1, 2}."
         error stop errmsg
      end if

      ! Algorthm selection
      select case (order)
         ! ------------------------------- 1st-order RK (Euler) -------------------------------
         ! Equation (4.10), page 43.
      case (1)
         do
            call fu(t, u, udot)
            u = u + dt*udot
            t = t + dt
            if (is_done(t, tout, dt) .or. itask == 2) exit
         end do

         ! --------------------------------- 2nd-order RK -------------------------------------
         ! Equation (4.10), page 43.
      case (2)
         do
            call fu(t, u, udot)
            ui = u + dt*udot
            call fu(t, ui, udot)
            u = (u + ui + dt*udot)/2
            t = t + dt
            if (is_done(t, tout, dt) .or. itask == 2) exit
         end do

         ! --------------------------------- 3rd-order RK -------------------------------------
         ! Equation (4.11), page 43.
      case (3)
         do
            call fu(t, u, udot)
            ui = u + dt*udot
            call fu(t, ui, udot)
            ui = (3*u + ui + dt*udot)/4
            call fu(t, ui, udot)
            u = (u + 2*ui + 2*dt*udot)/3
            t = t + dt
            if (is_done(t, tout, dt) .or. itask == 2) exit
         end do

      end select

   end subroutine rktvd

   pure subroutine mstvd(fu, u, t, tout, dt, uold, udotold, istate)
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
      procedure(integrand) :: fu
        !! subroutine with the derivative u'(t)
      real(rk), intent(inout) :: u(:)
        !! vector(N) with the variables to integrate u(t)
      real(rk), intent(inout) :: t
        !! time; on return it will be the current value of t (close to tout)
      real(rk), intent(in) :: tout
        !! time where next output is desired
      real(rk), intent(in) :: dt
        !! time step
      real(rk), intent(inout) :: uold(:, :)
        !! array(N,4) with the 4 previous values of u(t)
      real(rk), intent(inout) :: udotold(:, :)
        !! array(N,4) with the 4 previous values of utot(t) = u'(t)
      integer, intent(inout) :: istate
        !! flag indicating the state of the integration:
        !! 1 first call for a problem,
        !! 2 subsequent call for a problem.

      integer, parameter :: order = 3
      real(rk), dimension(size(u)) :: ui, udot
      character(:), allocatable :: errmsg
      integer :: itask_rktvd, istate_rktvd, i

      ! Check input conditions
      if (is_done(t, tout, dt)) return

      if (istate == 1) then
         if (size(uold, 2) /= 4 .or. size(udotold, 2) /= 4) then
            errmsg = "Invalid dimensions of arrays 'uold' or 'udotold' in 'mstvd'."
            error stop errmsg
         end if
      else if (istate < 1 .or. istate > 2) then
         errmsg = "Invalid input 'istate' in 'mstvd'. Valid set: {1, 2}."
         error stop errmsg
      end if

      ! The first starting values must be computed with a single-step method: we chose
      ! the RK method of the same order.
      ! The factor 2 in 't+2*dt' is not important, it just needs to be larger than 1.0
      ! so that one full 'dt' step can be computed.
      if (istate == 1) then

         itask_rktvd = 2
         istate_rktvd = 1

         do i = (order + 1), 1, -1
            uold(:, i) = u
            call fu(t, u, udotold(:, i))
            call rktvd(fu, u, t, t + 2*dt, dt, order, itask_rktvd, istate_rktvd)
         end do

         istate = 2

      end if

      ! Equation (4.26), page 48.
      do

         if (is_done(t, tout, dt)) exit

         call fu(t, u, udot)
         ui = (25*u + 50*dt*udot + 7*uold(:, 4) + 10*dt*udotold(:, 4))/32
         t = t + dt

         ! Shift u and udot values one step into the past
         do i = order, 1, -1
            udotold(:, i + 1) = udotold(:, i)
            uold(:, i + 1) = uold(:, i)
         end do
         udotold(:, 1) = udot
         uold(:, 1) = u
         u = ui

      end do

   end subroutine mstvd

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
