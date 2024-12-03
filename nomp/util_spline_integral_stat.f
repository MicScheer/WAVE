*CMZ :          16/08/2024  14.55.31  by  Michael Scheer
*CMZ :  4.01/03 16/05/2023  19.38.31  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_INTEGRAL_STAT(X,Y,N,RESULT
     &                                 ,COEF,WORK1,WORK2,WORK3,WORK4,ISTAT)

C---  CALCULATES INTERGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N,ISTAT
      REAL*8 X(N),Y(N),RESULT
      REAL*8 COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF_STATus(X,Y,N,0.0d0,0.0d0,COEF,WORK1,WORK2,WORK3,WORK4,ISTAT)

C--- INTEGRATION

      RESULT=0.0D0
      DO I=1,N-1

      RESULT=RESULT
     &          +(X(I+1)-X(I))*0.5D0
     &          *(Y(I)+Y(I+1))
     &          -(X(I+1)-X(I))**3/24.D0
     &          *(COEF(I)+COEF(I+1))

      ENDDO

      RETURN
      END
