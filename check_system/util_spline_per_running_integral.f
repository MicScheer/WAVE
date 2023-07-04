*CMZ : 00.00/16 18/03/2014  17.02.27  by  Michael Scheer
*CMZ : 00.00/11 26/05/2011  12.31.14  by  Michael Scheer
*-- Author :    Michael Scheer   11/02/2011
      SUBROUTINE util_spline_periodic_running_integral(X,Y,N,resultat
     &  ,COEF,WORK1,WORK2,WORK3,WORK4,WORK5,IFAIL)

C---  CALCULATES RUNNING INTERGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N,IFAIL
      REAL*8 X(N),Y(N),resultat(n)
      REAL*8 COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N),WORK5(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF_PERIODIC(X,Y,N,COEF,
     &  WORK1,WORK2,WORK3,WORK4,WORK5,IFAIL)

C--- INTEGRATION

      resultat(1)=0.0D0
      DO I=1,N-1

        resultat(i+1)=resultat(i)
     &          +(X(I+1)-X(I))*0.5D0
     &          *(Y(I)+Y(I+1))
     &          -(X(I+1)-X(I))**3/24.D0
     &          *(COEF(I)+COEF(I+1))

      ENDDO

      RETURN
      END
