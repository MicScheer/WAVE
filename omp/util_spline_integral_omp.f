*CMZ :  4.01/03 09/06/2023  12.33.19  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.03.55  by  Michael Scheer
*CMZ :  2.63/03 21/05/2008  13.37.36  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_INTEGRAL_omp(N,result)
*KEEP,gplhint.
*KEND.

C---  CALCULATES INTEGRAL OF Y(X) VIA SPLINES

      use wobsvmod

      IMPLICIT NONE

      INTEGER I,N
      !REAL*8 X(N),Y(N),RESULT
      real*8 RESULT
      !REAL*8 COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF_omp(N,-9999.0d0,-9999.0d0)

C--- INTEGRATION

      RESULT=0.0D0
      DO I=1,N-1

      RESULT=RESULT
     &          +(x_th(I+1)-x_th(I))*0.5D0
     &          *(wobsv1_th(I)+wobsv1_th(I+1))
     &          -(x_th(I+1)-x_th(I))**3/24.D0
     &          *(wobsv3_th(I)+wobsv3_th(I+1))

      ENDDO

      RETURN
      END
