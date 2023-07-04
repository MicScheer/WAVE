*CMZ : 00.00/07 28/01/2008  16.49.49  by  Michael Scheer
*CMZ : 00.00/02 27/06/97  10.48.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_PARABEL_DIFF(X,Y,A,YP,XOPT,AOPT,IFAIL
     &                            ,PLEVEL,XLEVEL)

C--- CALCULATES A(1),A(2),A(3), THE DERIVATIVES YP(X(1)),YP(X(2)),YP(X(3)),
C    AND THE EXTREMUM (XOPT,A(XOPT)) OF PARABOLA Y=A1+A2*X+A3*X**2
C    FROM COORDINATES OF THE THREE POINTS (X(1),Y(1)),(X(2),Y(2)),(X(3),Y(3))
C
C     IN ADDITION X1,2=XLEVEL(1:2) IS RETURNED WITH
C         Y(XLEVEL)=PLEVEL*Y(2) OR
C         Y(XLEVEL)=Y(2)/PLEVEL
C     (DEPENDING ON SECOND DERIVATIVE)
C
C
C     IFAIL=-1    FAILURE IN UTIL_PARABEL
C     IFAIL=1      NO SOLUTION FOUND
C     IFAIL=2      XLEVEL(1)=XLEVEL(2)
C     IFAIL=3      A(3).EQ.0
C
C
C
      IMPLICIT NONE

      INTEGER IFAIL

      REAL*8 A(3),X(3),Y(3)
      REAL*8 YP(3),XOPT,AOPT
      REAL*8 XLEVEL(2),PLEVEL,PPLEVEL

      CALL UTIL_PARABEL(X,Y,A,YP,XOPT,AOPT,IFAIL)

      IF (IFAIL.NE.0) THEN
        IFAIL=-1
        RETURN
      ENDIF

      IF (A(3).NE.0.D0) THEN

        XOPT=-A(2)/(2.D0*A(3))
        AOPT=A(1)+A(2)*XOPT+A(3)*XOPT**2

        IF (A(3).LT.0.D0) THEN
          PPLEVEL=1.D0/PLEVEL
        ELSE
          PPLEVEL=PLEVEL
        ENDIF

        XLEVEL(2)=
     &      (-4.D0*A(1)*A(3)+A(2)*A(2)+4.D0*A(3)*PPLEVEL*Y(2))
        IF(XLEVEL(2).LT.0.D0) THEN
          XLEVEL(1)=0.D0
          XLEVEL(2)=0.D0
          IFAIL=2
          RETURN
        ELSE
          XLEVEL(2)=DSQRT(XLEVEL(2))
          XLEVEL(1)=( XLEVEL(2)-A(2))/(2.D0*A(3))
          XLEVEL(2)=(-XLEVEL(2)-A(2))/(2.D0*A(3))
        ENDIF

      ELSE

        IFAIL=3
        XLEVEL(1)=0.D0
        XLEVEL(2)=0.D0
        AOPT=-1.D30

        IF (Y(1).GT.AOPT) THEN
          AOPT=Y(1)
          XOPT=X(1)
        ENDIF
        IF (Y(2).GT.AOPT) THEN
          AOPT=Y(2)
          XOPT=X(2)
        ENDIF
        IF (Y(3).GT.AOPT) THEN
          AOPT=Y(3)
          XOPT=X(3)
        ENDIF

      ENDIF


      RETURN
      END
