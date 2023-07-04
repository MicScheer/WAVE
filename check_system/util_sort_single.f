*CMZ : 00.00/15 20/04/2012  09.12.02  by  Michael Scheer
*CMZ : 00.00/07 12/10/2009  13.32.06  by  Michael Scheer
*CMZ : 00.00/00 11/01/95  11.42.12  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SORT_SINGLE(N,RA)

C--- HEAPSORT ROUTINE; SEE NUMERICAL RECIPES 8.2 S 231

      IMPLICIT NONE

      INTEGER N,L,IR,I,J
      REAL*4 RA(N),RRA

      if (n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
