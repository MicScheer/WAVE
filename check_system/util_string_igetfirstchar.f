*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ : 00.00/06 07/01/2008  14.32.30  by  Michael Scheer
*CMZ :  1.19/07 22/08/2002  15.44.21  by  Michael Scheer
*-- Author :    Michael Scheer   09/11/2001
      SUBROUTINE util_string_igetfirstchar(N1,CLINE,C)

      IMPLICIT NONE

      INTEGER N1,I,IC
      CHARACTER(*) CLINE
      CHARACTER C,C1

      EQUIVALENCE (IC,C1)
      ic=0
      N1=-1

      DO I=1,LEN(CLINE)
        C1=CLINE(I:I)
        IF (C1.NE.' '.AND.IC.NE.9) THEN !NO BLANKS, NO TABS
          n1=I
          C=CLINE(I:I)
          RETURN
        ENDIF
      ENDDO

      RETURN
      END
