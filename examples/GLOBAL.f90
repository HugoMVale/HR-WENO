MODULE GLOBAL
!----------------------------------------------------------------------------------------------
!ADD COMMENTS
!
!----------------------------------------------------------------------------------------------
IMPLICIT NONE

    INTEGER, PARAMETER :: M=4
    !DOUBLE PRECISION,  ::
    DOUBLE PRECISION :: xright(0:M), x(M), dx(M), u(M), H(0:M)
    DOUBLE PRECISION :: time
    INTEGER :: fevals=0
    logical :: flag_output, first=.TRUE.

END MODULE GLOBAL
