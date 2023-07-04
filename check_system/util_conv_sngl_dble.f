*CMZ :          22/06/2017  12.43.46  by  Michael Scheer
*CMZ : 00.00/02 26/01/2004  09.37.26  by  Michael Scheer
*-- Author :    Michael Scheer   21/01/98

      SUBROUTINE UTIL_CONV_SNGL_DBLE(X4,X8)

      REAL*8 X8
      REAL*4 X4
      CHARACTER(17) C17

        WRITE(C17,'(E15.7e3)')X4
      READ(C17,*)X8

      RETURN
      END
