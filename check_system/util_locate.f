*CMZ : 00.00/02 22/03/99  11.54.48  by  Michael Scheer
*-- Author :    Michael Scheer   19/03/99

      SUBROUTINE UTIL_LOCATE(NDIM,XA,X,IOLD,I,IFAIL)

      IMPLICIT NONE

      INTEGER IOLD,NDIM,IFAIL,I1,I2,IM,ID,I

      REAL*8 XA(NDIM),X

      IFAIL=0

      IF (X.LT.XA(1).OR.X.GT.XA(NDIM)) THEN
          IFAIL=1
          RETURN
      ENDIF

      IF (XA(1).GE.XA(NDIM)) THEN
          IFAIL=2
          RETURN
      ENDIF

      IF (IOLD.GT.NDIM.OR.IOLD.LT.1) IOLD=1

C HUNTING

      I1=IOLD
      IF (X.GE.XA(I1)) THEN
          ID=1
      ELSE
          ID=-1
      ENDIF

10        I2=I1+ID
          IF (I2.GE.NDIM) THEN
         I2=NDIM
         GOTO 100
          ELSEIF (I2.LE.1) THEN
         I2=I1
         I1=1
         GOTO 100
          ELSEIF (ID.GT.0.AND.X.GE.XA(I2).OR.ID.LT.0.AND.X.LT.XA(I2)) THEN
         I1=I2
         ID=2*ID
         GOTO 10
          ENDIF

100   CONTINUE

      IF (I2.LT.I1) THEN
          IM=I1
          I1=I2
          I2=IM
      ENDIF

      IF (X.LT.XA(I1).OR.X.GT.XA(I2)) THEN
          IFAIL=3
          RETURN
      ENDIF

C BISECTIONING

200   IF (I2-I1.LE.1) THEN
          IF (X.EQ.XA(I2)) THEN
         I=I2
          ELSE
         I=I1
          ENDIF
          IOLD=I
          IF (X.LT.XA(I1).OR.X.GT.XA(I2)) THEN
         IFAIL=4
         RETURN
          ENDIF
          RETURN
      ENDIF

      IM=(I2+I1)/2
      IF (X.GT.XA(IM)) THEN
          I1=IM
      ELSE
          I2=IM
      ENDIF

      GOTO 200

      END
