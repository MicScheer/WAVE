*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ : 00.00/05 22/12/2006  19.08.21  by  Michael Scheer
*CMZ : 00.00/02 13/09/2002  13.37.13  by  Michael Scheer
*CMZ :  1.19/07 22/08/2002  15.44.21  by  Michael Scheer
*CMZ :  1.00/00 22/11/2001  16.14.10  by  Michael Scheer
*CMZ :  0.01/02 20/11/2001  15.36.03  by  Michael Scheer
*CMZ :  0.01/01 19/11/2001  15.17.00  by  Michael Scheer
*CMZ :  0.01/00 16/11/2001  17.38.53  by  Michael Scheer
*-- Author :    Michael Scheer   09/11/2001

      SUBROUTINE UTIL_GET_NUMBER_OF_WORDS(LEN,CLINE,NWORDS)

      IMPLICIT NONE

      INTEGER LEN,I,NWORDS,IN,INO,IC
      CHARACTER(*) CLINE
      CHARACTER C1

      EQUIVALENCE (IC,C1)
        ic=0

        NWORDS=0
        IN=0
        INO=0
        DO I=1,LEN
          C1=CLINE(I:I)
            IF (C1.NE.' '.AND.IC.NE.9.AND.C1.NE.',') THEN
              IN=1
            ELSE
              IN=0
            ENDIF
            IF (NWORDS.EQ.0.AND.IN.EQ.1) THEN
              NWORDS=1
              INO=1
            ENDIF
            IF (IN.EQ.1.AND.INO.EQ.0) THEN
              NWORDS=NWORDS+1
            ENDIF
            INO=IN
      ENDDO

      RETURN
      END
