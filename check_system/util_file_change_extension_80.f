*CMZ : 00.00/02 26/01/2004  09.37.25  by  Michael Scheer
*CMZ : 00.00/01 29/01/97  12.00.35  by  Michael Scheer
*-- Author :    Michael Scheer   29/01/97

      SUBROUTINE UTIL_FILE_CHANGE_EXTENSION_80(FILE,EXT)

      IMPLICIT NONE

      CHARACTER(80) FILE,EXT,WS
      INTEGER I,ILAST,IDOT,ILENEXT

      DO I=1,80
          WS(I:I)=' '
      ENDDO

      IDOT=-9
      DO I=1,80
          IF(FILE(I:I).NE.'.'.AND.FILE(I:I).NE.' ') THEN
         WS(I:I)=FILE(I:I)
         ILAST=I
          ELSE
         IDOT=I
         WS(IDOT:IDOT)='.'
         GOTO 99
          ENDIF
      ENDDO

99    CONTINUE

      IF (IDOT.EQ.-9.AND.ILAST.LT.80) THEN
          IDOT=ILAST+1
          WS(IDOT:IDOT)='.'
      ENDIF

      ILENEXT=0
      DO I=1,80
          IF(EXT(I:I).NE.'.') GOTO 999
          IF(EXT(I:I).NE.' ') ILENEXT=ILENEXT+1
      ENDDO

999   IF (IDOT+1+ILENEXT.GT.80) THEN
          WRITE(6,*)'*** ERROR IN UTIL_FILE_CHANGE_EXTENSION_80:'
          WRITE(6,*)
     &'FILENAME AND EXTENSION REQUIRES MORE THAN 80 CHARACTERS'
          STOP '*** PROGRAM ABORTED DUE TO ERROR ***'
      ENDIF

      DO I=1,ILENEXT
          IF (IDOT+I.LE.80) THEN
            IF (EXT(I:I).NE.' '.AND.EXT(I:I).NE.'.') THEN
         WS(IDOT+I:IDOT+I)=EXT(I:I)
            ELSE
         WS(IDOT+I:IDOT+I)=' '
            ENDIF
          ENDIF
      ENDDO

      FILE=WS

      RETURN
      END
