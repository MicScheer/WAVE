*CMZ : 00.00/15 24/09/2013  09.51.24  by  Michael Scheer
*CMZ : 00.00/07 21/07/2009  14.01.57  by  Michael Scheer
*CMZ : 00.00/02 26/01/2004  09.37.25  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.28.07  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_WARTE(IBATCH)

      INTEGER IBATCH
      CHARACTER(5) ANS

      WRITE(6,*)
      WRITE(6,*)'STOP CAUSED BY ROUTINE UTIL_WARTE'
      WRITE(6,*)
      READ(5,'(1A5)',ERR=99) ANS

      IF (IBATCH.NE.0) RETURN

      IF (ANS.EQ.'BLITZ'.OR.ANS.EQ.'blitz') THEN
        CALL sleep( 1)
      ELSE IF (ANS.EQ.'SHORT'.OR.ANS.EQ.'short') THEN
        CALL sleep( 15)
      ELSE IF (ANS.EQ.'LONG'.OR.ANS.EQ.'long') THEN
        CALL sleep( 300)
      ELSE IF (ANS.EQ.'PAUSE'.OR.ANS.EQ.'pause') THEN
        CALL sleep( 3600)

      ENDIF
99    RETURN
      END
