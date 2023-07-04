*CMZ : 00.00/14 23/09/2011  10.51.55  by  Michael Scheer
*CMZ : 00.00/02 07/08/97  12.22.07  by  Michael Scheer
*-- Author :    Michael Scheer   07/08/97

      SUBROUTINE UTIL_PARA_INT_2D(X3,Y3,F33,X,Y,F,JFAIL)

      IMPLICIT NONE

      INTEGER JFAIL,IFAIL,IX,IY
      REAL*8 X,Y,F,X3(3),Y3(3)
      REAL*8 F3(3),FF3(3),F33(3,3)

      JFAIL=0

      DO IY=1,3

         DO IX=1,3
             F3(IX)=F33(IX,IY)
         ENDDO !IX

         CALL UTIL_PARABEL_INTER_1D(X3,F3,X,FF3(IY),IFAIL)
         IF (IFAIL.NE.0) JFAIL=IFAIL

      ENDDO !IY

      CALL UTIL_PARABEL_INTER_1D(Y3,FF3,Y,F,IFAIL)
      IF (IFAIL.NE.0) JFAIL=IFAIL

      RETURN
      END
