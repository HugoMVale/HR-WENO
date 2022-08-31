PROGRAM Test_WENO
!----------------------------------------------------------------------------------------------
! This program illustrates the application of modules WENO and TDVODE for solving a 1D
! hyperbolic equation:
!                    du(x,t)/dt = - d(f(u(x,t)))/dx
!                    with f(u,t)) = (u^2)/2
!
! The spatial variable 'x' is discretized according to a finite-volume approach and the time
! variable 't' is left continuous (method of lines), leading to:
!                    du(i,t)/dt = -1/dx(i)*(f(i+1/2,t) - f(i-1/2,t))
!
!
!                               u(i-1/2)^+            u(i+1/2)^-
!               --|-----(i-1)------|---------(i)----------|------(i+1)-----|--
!                               x(i-1/2)               x(i+1/2)
!                                  |<---    dx(i)     --->|
!----------------------------------------------------------------------------------------------
    USE TVDODE
    USE GLOBAL
    IMPLICIT NONE
    DOUBLE PRECISION :: dt, time_out, time_start, time_end, Rx, xmin, xmax
    INTEGER :: num_time_points, order, i
    !//////////////////////////////////////////////////////////////////////////////////////////

    ! Define grid limits
    xmin = -1.0d0
    xmax = 1.0d0

    ! Build a mesh (you can try geometric as well)
    Rx = (xmax - xmin)/DBLE(M)
    DO i=0,M
        xright(i) = xmin + Rx*DBLE(i)
    END DO

    ! Initial condition u(x,t=0), e.g. a simple rectangular pulse
    u = 0.0d0
    u(1:M/2) = 1.0d0

    ! Time settings
    time_start = 0.0d0
    time_end = 4.0d0
    dt = 1.0d-2
    order = 2

    ! Cell size and center (general case)
    dx = xright(1:M) - xright(0:M-1)
    x = xright(0:M-1) + 0.5d0*dx

    ! Open file where results will be stored
    CALL OUTPUT(1)

    ! Call ODE time solver
    num_time_points = 50
    time = time_start
    DO i=0,num_time_points

        time_out = DBLE(i)/DBLE(num_time_points)*time_end

        CALL rktvd(RHS,M,dt,time,time_out,u,order,0)

        flag_output = .TRUE.

        CALL RHS(M,time,u)

        flag_output = .FALSE.

    END DO

    ! End of simulation
    CALL OUTPUT(3)

    CONTAINS

    FUNCTION FLUX(u,t)
    !------------------------------------------------------------------------------------------
    ! Flux function.
    ! ARGUMENTS:
    ! u     function u(x,t)
    ! t     variable t
    !------------------------------------------------------------------------------------------
        IMPLICIT NONE
        DOUBLE PRECISION, dimension(:), INTENT (IN) :: u
        DOUBLE PRECISION, dimension(size(u)) :: FLUX
        DOUBLE PRECISION, INTENT (IN), OPTIONAL :: t
    !//////////////////////////////////////////////////////////////////////////////////////////

        FLUX = 0.5d0*u*u

    END FUNCTION FLUX
    !##########################################################################################


    SUBROUTINE RHS(M,time,u,du)
    !------------------------------------------------------------------------------------------
    ! This
    !------------------------------------------------------------------------------------------
        USE WENO
        USE GLOBAL, ONLY : flag_output, first, fevals, dx, x
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: M
        DOUBLE PRECISION, INTENT (IN) :: time, u(M)
        DOUBLE PRECISION, INTENT (OUT), OPTIONAL :: du(M)
        INTEGER, PARAMETER :: k=3
        DOUBLE PRECISION, PARAMETER :: EPSILON=1.0d-30
        !DOUBLE PRECISION :: u(-1:nmax+2), upoint(nmax,2)
        !DOUBLE PRECISION, SAVE :: c(1:M,-1:k-1,0:k-1)=0

        INTEGER :: i
    !//////////////////////////////////////////////////////////////////////////////////////////

        !WHAT is u and upoint??

        IF (first) THEN

            first = .FALSE.

            !CALL CALC_CIRJ(k,M,x(2:M+1),c) !bulshit

            !u(-1:0) = 0.0d0
            !u(nmax+1:nmax+2) = 0.0d0

            !H(0)  = 0.0d0

        END IF

        IF (.NOT.(flag_output)) THEN

            fevals = fevals + 1

            DO i=1,M

                du(i) = (0 - 0)/dx(i)

            END DO

        END IF

        !   OUTPUT
        IF (flag_output) THEN

            CALL OUTPUT(2)

        END IF


    END SUBROUTINE RHS
    !##########################################################################################


    SUBROUTINE OUTPUT(message)
    !------------------------------------------------------------------------------------------
    ! Auxiliary routine to save results to file.
    ! ARGUMENTS:
    ! message: parameter to select output action
    !------------------------------------------------------------------------------------------
        USE GLOBAL, ONLY : M, u, x, time, fevals
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: message
        INTEGER :: i
    !//////////////////////////////////////////////////////////////////////////////////////////

        SELECT CASE (message)

        ! Open files and write headers and grid
        CASE (1)

            PRINT *, "Running test..."
            PRINT *, "Start: ", FDATE()

            ! Write grid
            OPEN (UNIT=1, FILE="./output/xgrid.txt", STATUS="REPLACE", ACTION="WRITE", &
                  POSITION="REWIND")

            WRITE (1,'(1X, A5, 2(A15))') "i", "x(i)", "dx(i)"
            DO i=1,M
                WRITE (1,'(1X, I5, 2(E15.5))') i, x(i), dx(i)
            END DO

            !Write header u
            OPEN (UNIT=2, FILE="./output/u.txt", STATUS="REPLACE", ACTION="WRITE", &
                  POSITION="REWIND")

            WRITE (2,'(1X, A15)', ADVANCE="NO") "t"
            DO i=1,M
                WRITE (2,'(A15)', ADVANCE="NO") "u("//itoa(i)//")"
            END DO
            WRITE (2,*) ""

        ! Write values
        CASE (2)
            WRITE (2,'(1X, E15.5)', ADVANCE="NO") time
            DO i=1,M
                WRITE (2,'(E15.5)', ADVANCE="NO") u(i)
            END DO
            WRITE (2,*) ""

        ! Close files
        CASE (3)
            CLOSE (1)
            CLOSE (2)
            PRINT *, "End  : ", FDATE()
            PRINT '(A13,I5)', "RHS Evals.: ", fevals
            PRINT *

        END SELECT

    END SUBROUTINE OUTPUT
    !##########################################################################################


    FUNCTION itoa(i) RESULT(res)
    !------------------------------------------------------------------------------------------
    ! Convert integer to string.
    ! ARGUMENTS:
    ! 1:    integer
    !------------------------------------------------------------------------------------------
      CHARACTER(:), ALLOCATABLE :: res
      INTEGER, INTENT(IN) :: i
      CHARACTER(RANGE(i)+2) :: tmp
      WRITE(tmp,'(i0)') i
      res = TRIM(tmp)
    END FUNCTION

END PROGRAM
