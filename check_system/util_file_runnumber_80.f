*CMZ : 00.00/02 26/01/2004  09.37.25  by  Michael Scheer
*CMZ : 00.00/01 29/01/97  12.14.02  by  Michael Scheer
*-- Author :    Michael Scheer   29/01/97

      SUBROUTINE UTIL_FILE_RUNNUMBER_80(FILE,IRUN)

      IMPLICIT NONE

      CHARACTER(80) FILE,WS,CRUN
      INTEGER I,ILAST,ILENCRUN,IRUN,ILEN

      CALL IZITOC(IRUN,CRUN)

      ILEN=0
      DO I=1,80
          WS(I:I)=' '
          IF (FILE(ILEN:ILEN).NE.' ') ILEN=ILEN+1
      ENDDO

      ILENCRUN=0
      DO I=1,80
          IF(CRUN(I:I).NE.' ') ILENCRUN=ILENCRUN+1
      ENDDO

      IF (ILENCRUN+ILEN.GT.80) THEN
          WRITE(6,*)'*** ERROR IN UTIL_FILE_RUNNUMBER_80:'
          WRITE(6,*)
     &'FILENAME AND RUNNUMBER REQUIRES MORE THAN 80 CHARACTERS'
          STOP '*** PROGRAM ABORTED DUE TO ERROR ***'
      ENDIF

      DO I=1,80
          IF(FILE(I:I).NE.'.') THEN
         WS(I:I)=FILE(I:I)
         ILAST=I
          ELSE
         GOTO 99
          ENDIF
      ENDDO
99    CONTINUE

      DO I=1,ILENCRUN
         WS(ILAST+I:ILAST+I)=CRUN(I:I)
      ENDDO
      ILAST=ILAST+ILENCRUN

      DO I=ILAST+1,80
          WS(I:I)=FILE(I:I)
      ENDDO

      FILE=WS

      RETURN
      END
