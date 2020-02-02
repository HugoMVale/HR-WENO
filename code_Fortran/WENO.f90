SUBROUTINE WENO(k,M,eps,v,c,vl,vr)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine implements the (2k-1)th order WENO method for ARBITRARY finite volume grids.
!k = 2 or 3.
!Source: ICASE 97-65 by Shu, 1997.
!
!VARIABLES:
!
!eps        numerical smoothing factor
!M          number of cells
!k          order of reconstruction within the cell
!v(i)       average value of cell i
!c(i,r,j)   array of constants for the grid
!vl(i)      reconstructed value at left  boundary of cell i (v_{i-1/2}^+)
!vr(i)      reconstructed value at right boundary of cell i (v_{i+1/2}^-)
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: k, M
DOUBLE PRECISION, INTENT(IN) :: eps, v(1-(k-1):M+(k-1)), c(1:M,-1:k-1,0:k-1)
DOUBLE PRECISION, INTENT(OUT) :: vl(1:M), vr(1:M)

DOUBLE PRECISION :: vlr(0:k-1), vrr(0:k-1), w(0:k-1), wtil(0:k-1), alfa(0:k-1), &
                    alfatil(0:k-1), beta(0:k-1), d(0:k-1)
INTEGER :: i, j, r

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


! SELECT PARAMETERS ACCORDING TO ORDER OF THE METHOD

IF (k==2) THEN

    d = (/2.0d0/3.0d0, 1.0d0/3.0d0/)

ELSE IF (k==3) THEN

    d = (/0.3d0, 0.6d0, 0.1d0/)

END IF



! ALGORITHM

vl = 0.0d0
vr = 0.0d0


DO i=1,M


	!	OBTAIN THE k RECONSTRUCTED VALUES vi+1/2(r) & vi-1/2(r)

    vrr = 0.0d0
    vlr = 0.0d0

    DO r=0,k-1

        DO j=0,k-1

            vrr(r) = vrr(r) + c(i,r,j)*v(i-r+j)

            vlr(r) = vlr(r) + c(i,r-1,j)*v(i-r+j)

        END DO

    END DO


    !	FORM THE WEIGHTS W AND W~

    IF (k==2) THEN

        beta(0) = (v(i+1) - v(i))**2

        beta(1) = (v(i) - v(i-1))**2


    ELSE IF (k==3) THEN

        beta(0) = 13.0d0/12.0d0*(v(i) - 2.0d0*v(i+1) + v(i+2))**2 &

                  + 1.0d0/4.0d0*(3.0d0*v(i) - 4.0d0*v(i+1) + v(i+2))**2

        beta(1) = 13.0d0/12.0d0*(v(i-1) - 2.0d0*v(i) + v(i+1))**2 &

                  + 1.0d0/4.0d0*(v(i-1) - v(i+1))**2

        beta(2) = 13.0d0/12.0d0*(v(i-2) - 2.0d0*v(i-1) + v(i))**2 &

                  + 1.0d0/4.0d0*(v(i-2) - 4.0d0*v(i-1) + 3.0d0*v(i))**2

    END IF


    DO r=0,k-1

        alfa(r) = d(r)/(eps + beta(r))**2

        alfatil(r) = d(k-1-r)/(eps + beta(r))**2

    END DO


    DO r=0,k-1

        w(r) = alfa(r)/SUM(alfa)

        wtil(r) = alfatil(r)/SUM(alfatil)

    END DO




    !	FIND THE (2k-1)-th ORDER RECONSTRUCTION

    DO r=0,k-1

        vr(i) = vr(i) + w(r)*vrr(r)

        vl(i) = vl(i) + wtil(r)*vlr(r)

    END DO


END DO


END SUBROUTINE WENO
!##############################################################################################

















SUBROUTINE CALC_CIRJ(k,M,xright,c)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine computes the array for constants to required as input to WENO.
!Source: ICASE 97-65 by Shu, 1997.
!Equation (2.20), page 6.
!
!VARIABLES:
!
!M          number of cells
!k          order of reconstruction within the cell
!xright(i)  value of x at right boundary of cell i (x_{i+1/2})
!c(i,r,j)   array of constants for the grid
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: k, M
DOUBLE PRECISION, INTENT(IN) :: xright(0:M)
DOUBLE PRECISION, INTENT(OUT) :: c(1:M,-1:k-1,0:k-1)

DOUBLE PRECISION :: test, prod1, prod2, sum1, sum2, xl(-2:M+3), xr(-2:M+3)
INTEGER :: i, j, r, l, n, q

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


xr(0:M) = xright(0:M)
xl(1:M) = xright(0:M-1)


!	cirj VALUES FOR UNIFORM GRID

test = DABS((xr(1) - xl(1))/(xr(M) - xl(M)) - 1.0d0)


IF (test>1.0d-9) THEN   !cuidado <

    IF (k==2) THEN

        DO i=1,M

            c(i,-1,:) = (/3.0d0/2.0d0,-1.0d0/2.0d0/)
            c(i,0,:)  = (/1.0d0/2.0d0,1.0d0/2.0d0/)
            c(i,1,:)  = (/-1.0d0/2.0d0,3.0d0/2.0d0/)

        END DO

    ELSE IF (k==3) THEN

        DO i=1,M

            c(i,-1,:) = (/11.0d0/6.0d0,-7.0d0/6.0d0,1.0d0/3.0d0/)
            c(i,0,:)  = (/1.0d0/3.0d0,5.0d0/6.0d0,-1.0d0/6.0d0/)
            c(i,1,:)  = (/-1.0d0/6.0d0,5.0d0/6.0d0,1.0d0/3.0d0/)
            c(i,2,:)  = (/1.0d0/3.0d0,-7.0d0/6.0d0,11.0d0/6.0d0/)

        END DO

    END IF

    RETURN

END IF



!	ALGORITHM FOR NON-UNIFORM GRID


DO i=1,M

    DO r=-1,2

        DO j=0,2

            sum2 = 0.0d0

            DO n=j+1,k

                prod2 = 1.0d0

                DO l=0,k

                    IF (l==n) CYCLE

                    prod2 = prod2*(xl(i-r+n) - xl(i-r+l))

                END DO

                sum1=0.0d0

                DO l=0,k

                    IF (l==n) CYCLE

                        prod1 = 1.0d0

                        DO q=0,k

                            IF (q==n .OR. q==l) CYCLE

                            prod1 = prod1*(xr(i) - xl(i-r+q))

                        END DO

                    sum1 = sum1 + prod1

                END DO

                sum2 = sum2 + sum1/prod2

            END DO

            c(i,r,j) = sum2*(xr(i-r+j) - xl(i-r+j))

        END DO

    END DO

END DO


END SUBROUTINE CALC_CIRJ
!##############################################################################################
