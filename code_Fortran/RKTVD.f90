SUBROUTINE RKTVD(DU,neq,dt,t,tout,U,order,iopt)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine implements the optimal 1st, 2nd and 3rd order TVD RK methods (ICASE 97-65 by
!Shu, 1997).
!It is very important to use TVD schemes for time discretization. Even with a TVD spacial
!discretization, if the time discretization is by a non-TVD method, the result may be
!oscillatory.
!
!The routine was built to work like LSODE.
!
!VARIABLES:
!
!DU		function with the derivatives
!neq	number of equations
!dt		time step
!t		time; on return it will be the current value of t (close to tout)
!tout	time where next output is desired
!U		dependent variables
!Unew	auxiliary vector
!L		derivative of U
!order  order of the method (1,2 or 3)
!iopt	parameter for a single step
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: neq, order, iopt
DOUBLE PRECISION, INTENT(IN) :: dt, tout
DOUBLE PRECISION, INTENT(INOUT) :: t, U(neq)
EXTERNAL DU

DOUBLE PRECISION :: Unew(neq), L(neq)
INTEGER :: i

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!	*** CHECK INPUT CONDITIONS ***

IF (tout-t<dt) RETURN

IF (order<1 .OR. order>3) STOP



!	*** CALCULATIONS ***

SELECT CASE (order)

! ------------------------------ 1st order RK (Euler) -----------------------------------------
! Equation (4.10), page 43.

    CASE (1)

    DO

        CALL DU(neq,t,U,L)

        DO i=1,neq

            U(i) = U(i) + dt*L(i)

        END DO


        t=t+dt

        IF (t>=tout .OR. iopt==1) EXIT

        END DO



! -------------------------------- 2nd order RK -----------------------------------------------
! Equation (4.10), page 43.

    CASE (2)

    DO

        CALL DU(neq,t,U,L)

        DO i=1,neq

            Unew(i) = U(i) + dt*L(i)

        END DO



        CALL DU(neq,t,Unew,L)

        DO i=1,neq

            U(i) = 0.5d0*U(i) + 0.5d0*Unew(i) + 0.5d0*dt*L(i)

        END DO


        t=t+dt

        IF (t>=tout .OR. iopt==1) EXIT

    END DO


! -------------------------------- 3rd order RK ------------------------------------------------
! Equation (4.11), page 43.

    CASE (3)

    DO

        CALL DU(neq,t,U,L)

        DO i=1,neq

            Unew(i) = U(i) + dt*L(i)

        END DO



        CALL DU(neq,t,Unew,L)

        DO i=1,neq

            Unew(i) = (3.0d0/4.0d0)*U(i) + 0.25d0*Unew(i) + 0.25d0*dt*L(i)

        END DO



        CALL DU(neq,t,Unew,L)

        DO i=1,neq

            U(i) = (1.0d0/3.0d0)*U(i) + (2.0d0/3.0d0)*Unew(i) + (2.0d0/3.0d0)*dt*L(i)

        END DO


        t=t+dt

        IF (t>=tout .OR. iopt==1) EXIT

    END DO


END SELECT


END SUBROUTINE RKTVD
!##############################################################################################
