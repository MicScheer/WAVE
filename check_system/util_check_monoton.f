*CMZ : 00.00/02 26/10/2004  11.55.46  by  Michael Scheer
*-- Author :    Michael Scheer   26/10/2004
      SUBROUTINE UTIL_CHECK_MONOTON(N,X,IMONO)

C     Returns IMONO: 2 for strictly raising,
C                    1 for raising,
C                    0 not monoton or const.
C                   -1 for falling,
C                   -2 for strictly falling,


      IMPLICIT NONE

      INTEGER I,N,IMONO,IMONO0,IUP
      DOUBLE PRECISION X(N)

      IUP=0
      DO I=1,N
        IF (X(1).LT.X(I)) THEN
          IUP=1
          GOTO 1
        ELSE IF (X(1).GT.X(I)) THEN
          IUP=-1
          GOTO 1
        ENDIF
      ENDDO

1     IMONO0=0
      IMONO=0

      IF (IUP.EQ.1) THEN

        DO I=1,N-1
          IF (X(I).EQ.X(I+1)) IMONO0=1
          IF (X(I).GT.X(I+1)) RETURN
        ENDDO

        IF (IMONO0.EQ.0) THEN
          IMONO=2
        ELSE
          IMONO=1
        ENDIF

      ELSE IF (IUP.EQ.-1) THEN

        DO I=1,N-1
          IF (X(I).EQ.X(I+1)) IMONO0=1
          IF (X(I).LT.X(I+1)) RETURN
        ENDDO

        IF (IMONO0.EQ.0) THEN
          IMONO=-2
        ELSE
          IMONO=-1
        ENDIF

      ELSE

        IMONO=0

      ENDIF

      RETURN
      END
