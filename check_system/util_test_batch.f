*CMZ : 00.00/02 26/01/2004  09.37.25  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.51  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_TEST_BATCH(IBATCH)

      CHARACTER(11) TESTBATCH
      INTEGER IBATCH

      CALL UTIL_GET_MODE(TESTBATCH)

C----------- OLD -------------------------------------------
C     CHARACTER(11) TESTBATCH
C     CALL LIB$SPAWN('@UTIL:UTIL_TEST_BATCH.COM')
C     OPEN(UNIT=10,FILE='UTIL:TEST_BATCH.DAT',STATUS='OLD')
C     READ(10,'(1A11)')TESTBATCH
C----------- OLD -------------------------------------------

      IF (TESTBATCH.EQ.'BATCH') THEN
        IBATCH=1
      ELSE
        IBATCH=0
      ENDIF

      RETURN
      END
