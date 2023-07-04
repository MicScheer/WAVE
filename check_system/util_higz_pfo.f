*CMZ : 00.00/07 21/07/2009  14.25.11  by  Michael Scheer
*CMZ : 00.00/02 26/01/2004  09.37.26  by  Michael Scheer
*-- Author :    Michael Scheer   07/01/98

      SUBROUTINE UTIL_HIGZ_PFO(LUN,FILE)

      IMPLICIT NONE

      INTEGER LUN
      CHARACTER(80) FILE

      OPEN(UNIT=LUN,FILE=FILE,STATUS='NEW')
      CALL IGMETA(LUN,-111)

      RETURN
      END
