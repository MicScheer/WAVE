*CMZ : 00.00/01 03/01/97  11.52.32  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.17  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_CHECK_GFLOAT_1(G,IGFLOAT)

      IMPLICIT NONE

      REAL*4 G,GG
      INTEGER J,JCOMP,IGFLOAT

      EQUIVALENCE(GG,J)

      DATA JCOMP/16512/

      GG=G

      IF (J.EQ.JCOMP) THEN
         IGFLOAT=0
      ELSE
         IGFLOAT=1
      ENDIF

      RETURN
      END
