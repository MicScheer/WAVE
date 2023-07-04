*CMZ : 00.00/16 03/06/2014  08.39.05  by  Michael Scheer
*CMZ : 00.00/15 05/01/2012  13.52.39  by  Michael Scheer
*CMZ : 00.00/00 11/01/95  11.41.04  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SORT_FUNC_INTEGER(N,IA,IYA)

C--- HEAPSORT ROUTINE; SEE NUMERICAL RECIPES 8.2 S 231
C--- ARRAY YA IS FUNCTION OF IA AND SORTED ACCORDINGLY

      IMPLICIT NONE

      INTEGER N,L,IR,I,J

      integer IA(N),IRA
      integer IYA(N),IYYA

      if (n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          IRA=IA(L)
          IYYA=IYA(L)
        ELSE
          IRA=IA(IR)
          IYYA=IYA(IR)
          IA(IR)=IA(1)
          IYA(IR)=IYA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            IA(1)=IRA
            IYA(1)=IYYA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(IA(J).LT.IA(J+1))J=J+1
          ENDIF
          IF(IRA.LT.IA(J))THEN
            IA(I)=IA(J)
            IYA(I)=IYA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        IA(I)=IRA
        IYA(I)=IYYA
      GO TO 10
      END
