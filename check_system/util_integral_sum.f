*CMZ : 00.00/16 18/03/2014  17.02.27  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.21.42  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.26.06  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_INTEGRAL_SUM(X,Y,N,resultat)

C---  CALCULATES INTERGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N
      REAL*8 X(N),Y(N),resultat


C--- INTEGRATION

      resultat=0.0D0
      DO I=1,N-1
         resultat=resultat+(Y(I)+Y(I+1))/2.*(X(I+1)-X(I))
      ENDDO

      RETURN
      END
