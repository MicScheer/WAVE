*CMZ : 00.00/02 19/03/99  14.44.17  by  Michael Scheer
*-- Author :    Michael Scheer   15/03/99

      SUBROUTINE UTIL_POLINT(XA,YA,N,X,Y,DY)

C POLYNOMIAL INTERPOLATION OF YA(XA) AT X. DY IS ESTIMATED ERROR OF Y(X)
C LITERATUR: NUMERICAL RECIPIES

      IMPLICIT NONE

      INTEGER N,NMAX,NS,I,M
      PARAMETER (NMAX=10)

      REAL*8 XA(N),YA(N),X,Y,DY,HP,HO,DEN
      REAL*8 C(NMAX),D(NMAX),DIF,W,DIFT

      NS=1
      DIF=ABS(X-XA(1))

      DO I=1,N

          DIFT=ABS(X-XA(I))
          IF (DIFT.LT.DIF) THEN
         NS=I
         DIF=DIFT
          ENDIF

          C(I)=YA(I)
          D(I)=YA(I)

      ENDDO !N

      Y=YA(NS)
      NS=NS-1

      DO M=1,N-1

          DO I=1,N-M
         HO=XA(I)-X
         HP=XA(I+M)-X
         W=C(I+1)-D(I)
         DEN=HO-HP
         IF (DEN.EQ.0.)
     &          STOP '*** ERROR UTIL_POLINT: BAD INTERVALS IN X  ***'
         DEN=W/DEN
         D(I)=HP*DEN
         C(I)=HO*DEN
          ENDDO   !I

          IF (2*NS.LT.N-M) THEN
         DY=C(NS+1)
          ELSE
         DY=D(NS)
         NS=NS-1
          ENDIF !(2*NS.LT.N-M)

          Y=Y+DY

      ENDDO !M

      RETURN
      END
