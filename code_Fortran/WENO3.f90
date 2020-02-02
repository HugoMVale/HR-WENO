SUBROUTINE WENO3(M,eps,v,vl,vr)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine implements th 3rd order WENO method for UNIFORM finite volume grids.
!k = 2.
!Source: ICASE 97-65 by Shu, 1997.
!
!VARIABLES:
!
!eps    numerical smoothing factor
!M      number of cells
!v(i)   average value of cell i
!vl(i)  reconstructed value at left  boundary of cell i (v_{i-1/2}^+)
!vr(i)  reconstructed value at right boundary of cell i (v_{i+1/2}^-)
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: M
DOUBLE PRECISION, INTENT(IN) :: eps, v(-1:M+2)
DOUBLE PRECISION, INTENT(OUT) :: vl(1:M), vr(1:M)

INTEGER, PARAMETER :: k=2
DOUBLE PRECISION, DIMENSION(0:k-1), PARAMETER :: d=(/2.0d0/3.0d0, 1.0d0/3.0d0/)
DOUBLE PRECISION, DIMENSION(-1:k-1,0:k-1), PARAMETER :: &
c = RESHAPE ((/3.0d0/2.0d0,-1.0d0/2.0d0,1.0d0/2.0d0,1.0d0/2.0d0,-1.0d0/2.0d0,3.0d0/2.0d0/), &
(/k+1,k/), ORDER=(/2,1/))

DOUBLE PRECISION :: vlr(0:k-1), vrr(0:k-1), w(0:k-1), wtil(0:k-1), alfa(0:k-1), &
                    alfatil(0:k-1), beta(0:k-1)
INTEGER :: i, j, r

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


vl = 0.0d0
vr = 0.0d0


DO i=1,M


    !	OBTAIN THE k RECONSTRUCTED VALUES vi+1/2(r) & vi-1/2(r)

    vrr = 0.0d0
    vlr = 0.0d0

    DO r=0,k-1

        DO j=0,k-1

            vrr(r) = vrr(r) + c(r,j)*v(i-r+j)

            vlr(r) = vlr(r) + c(r-1,j)*v(i-r+j)

        END DO

    END DO


    !	FORM THE WEIGHTS W AND W~

    beta(0) = (v(i+1) - v(i))**2

    beta(1) = (v(i) - v(i-1))**2



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


END SUBROUTINE WENO3
!##############################################################################################
