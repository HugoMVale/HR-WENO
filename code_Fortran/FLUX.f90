FUNCTION GODUNOV(uL,uR,G)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!"Classical" Godunov upwind flux, for a growth term.
!
!VARIABLES:
!
!G      Growth rate (Flux = G.u)
!uL     left  side of reconstructed value of u (u_{i^+1/2}^-)
!uR     right side of reconstructed value of u (u_{i^+1/2}^+)
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

DOUBLE PRECISION, INTENT (IN) :: uL, uR, G
DOUBLE PRECISION :: GODUNOV

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

IF(G>=0.0d0) THEN

    GODUNOV = G*uL

ELSE

    GODUNOV = G*uR

END IF


END FUNCTION GODUNOV
!##############################################################################################









FUNCTION LaxF(uL,uR,a,G)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Lax-Friedrichs flux, for a growth term.
!
!VARIABLES:
!
!a      constant equal to max(G) for all u.
!G      Growth rate (Flux = G.u)
!uL     left  side of reconstructed of u (u_{i^+1/2}^-)
!uR     right side of reconstructed of u (u_{i^+1/2}^+)
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

DOUBLE PRECISION, INTENT (IN) :: uL, uR, a, G
DOUBLE PRECISION :: LaxF

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


LaxF = 0.5d0*(G*uL + G*uR - a*(uR - uL))


END FUNCTION LaxF
!##############################################################################################
