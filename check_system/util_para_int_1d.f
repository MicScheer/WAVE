*CMZ : 00.00/14 23/09/2011  10.51.30  by  Michael Scheer
*CMZ : 00.00/02 07/08/97  12.21.55  by  Michael Scheer
*-- Author :    Michael Scheer   07/08/97

      SUBROUTINE UTIL_PARA_INT_1D(X3,Y3,X,Y,IFAIL)

      IMPLICIT NONE

      INTEGER IFAIL
      REAL*8 X3(3),Y3(3),X,Y,XOPT,AOPT,YP(3),A(3)

      IFAIL=0

      CALL UTIL_PARABEL(X3,Y3,A,YP,XOPT,AOPT,IFAIL)

      Y=A(1)+A(2)*X+A(3)*X*X

      RETURN
      END
