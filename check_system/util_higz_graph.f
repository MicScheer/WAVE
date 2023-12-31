*CMZ : 00.00/02 26/01/2004  09.37.26  by  Michael Scheer
*-- Author :    Michael Scheer   07/01/98

      SUBROUTINE UTIL_HIGZ_GRAPH(NPLOT,XPL,YPL,CHOPT,MTYP)

      IMPLICIT NONE

      CHARACTER(4) CHOPT
        INTEGER NPLOT,NPLOTP,MTYP,I
      PARAMETER (NPLOTP=1000)

        REAL XPL(NPLOTP),YPL(NPLOTP),XPL1,XPL2,YPL1,YPL2
     &      ,XPLMAX,XPLMIN,YPLMAX,YPLMIN,DXPL,DYPL,RMTYP

      IF (NPLOT.GT.NPLOTP) THEN
          NPLOT=NPLOTP
      ENDIF

      RMTYP=MTYP
      CALL IGSET('MTYP',RMTYP)

      XPLMAX=-1.E30
      XPLMIN=+1.E30
      YPLMAX=-1.E30
      YPLMIN=+1.E30

      DO I=1,NPLOT
             IF (XPL(I).LT.XPLMIN) XPLMIN=XPL(I)
             IF (XPL(I).GT.XPLMAX) XPLMAX=XPL(I)
             IF (YPL(I).LT.YPLMIN) YPLMIN=YPL(I)
             IF (YPL(I).GT.YPLMAX) YPLMAX=YPL(I)
      ENDDO

             DXPL=XPLMAX-XPLMIN
             DYPL=YPLMAX-YPLMIN
           IF (DXPL.EQ.0) DXPL=1.
           IF (DYPL.EQ.0) DYPL=1.
             XPL2=XPLMAX+0.1*DXPL
             XPL1=XPLMIN-0.1*DXPL
             YPL2=YPLMAX+0.1*DYPL
             YPL1=YPLMIN-0.1*DYPL
              CALL HPLFRA(XPL1,XPL2,YPL1,YPL2,' ')
              CALL IGRAPH(NPLOT,XPL,YPL,CHOPT)
              CALL IUWK(0,1)

      RETURN
      END
