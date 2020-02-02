SUBROUTINE MSTVD3(DU,neq,dt,t,tout,U,istate,Uold,Lold)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine implements a 5-step, 3rd order TVD multi-step method (ICASE 97-65 by Shu, 1997).
!It is very important to use TVD schemes for time discretization. Even with a TVD spacial
!discretization, if the time discretization is by a non-TVD method, the result may be
!oscillatory.
!In theory, this method should have an efficiency 1.5 times higher than the RK method of the
!same order. However, in practice they appear to be almost identical.
!
!The routine was built to work like LSODE.
!
!VARIABLES:
!
!DU         function with the derivatives
!neq        number of equations
!dt         time step
!t			time; on return it will be the current value of t (close to tout)
!tout       time where next output is desired
!U          dependent variables at t(n)
!Unew       value of U(n+1)
!L          derivative of U at t(n)
!istate     index to specify the state of the calculation (1 for the first call)
!Uold       vector with the 4 previous values of U
!Lold       vector with the 4 previous values of L
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: neq
INTEGER, INTENT (INOUT) :: istate
DOUBLE PRECISION, INTENT(IN) :: dt, tout
DOUBLE PRECISION, INTENT(INOUT) :: t, U(neq), Uold(4,neq), Lold(4,neq)
EXTERNAL DU

DOUBLE PRECISION :: Unew(neq), L(neq)
INTEGER :: i

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

IF (tout-t<dt) RETURN


IF (istate==1) THEN

    istate=2

    Uold(4,:) = U
    CALL DU(neq,t,U,Lold(4,:))

    CALL RKTVD3(DU,neq,dt,t,t+2.0d0*dt,U,1)
    Uold(3,:) = U
    CALL DU(neq,t,U,Lold(3,:))

    CALL RKTVD3(DU,neq,dt,t,t+2.0d0*dt,U,1)
    Uold(2,:) = U
    CALL DU(neq,t,U,Lold(2,:))

    CALL RKTVD3(DU,neq,dt,t,t+2.0d0*dt,U,1)
    Uold(1,:) = U
    CALL DU(neq,t,U,Lold(1,:))

    CALL RKTVD3(DU,neq,dt,t,t+2.0d0*dt,U,1)

END IF


! Equation (4.26), page 48.
DO

    IF (t>=tout) EXIT

    CALL DU(neq,t,U,L)

    DO i=1,neq

        Unew(i) = (25.0d0/32.0d0)*U(i) + (25.0d0/16.0d0)*dt*L(i) + (7.0d0/32.0d0)*Uold(4,i) &

                + (5.0d0/16.0d0)*dt*Lold(4,i)

    END DO

    Uold(4,:) = Uold(3,:)
    Uold(3,:) = Uold(2,:)
    Uold(2,:) = Uold(1,:)
    Uold(1,:) = U
    U  = Unew

    Lold(4,:) = Lold(3,:)
    Lold(3,:) = Lold(2,:)
    Lold(2,:) = Lold(1,:)
    Lold(1,:) = L

    t=t+dt

END DO


END SUBROUTINE MSTVD3
!##############################################################################################
