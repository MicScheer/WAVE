*CMZ :  4.01/03 28/06/2023  13.05.56  by  Michael Scheer
*CMZ :  4.00/11 11/06/2021  12.58.32  by  Michael Scheer
*CMZ :  4.00/07 29/04/2020  12.10.47  by  Michael Scheer
*CMZ :  4.00/04 05/08/2019  11.36.57  by  Michael Scheer
*CMZ :  3.02/03 04/11/2014  13.04.08  by  Michael Scheer
*CMZ :  3.02/00 10/09/2014  11.51.49  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.66/09 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/07 25/02/2010  12.59.17  by  Michael Scheer
*CMZ :  2.66/06 24/11/2009  14.22.00  by  Michael Scheer
*CMZ :  2.65/03 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.64/01 19/08/2009  15.11.55  by  Michael Scheer
*CMZ :  2.64/00 14/08/2009  14.31.12  by  Michael Scheer
*CMZ :  2.63/05 14/08/2009  13.06.13  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.37/04 03/12/2001  20.47.32  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.10  by  Michael Scheer
*CMZ :  2.37/01 14/11/2001  11.04.48  by  Michael Scheer
*CMZ :  2.33/00 27/04/2001  15.15.42  by  Michael Scheer
*CMZ :  2.31/01 25/04/2001  15.49.16  by  Michael Scheer
*CMZ :  2.30/02 12/04/2001  19.03.14  by  Michael Scheer
*CMZ :  2.20/03 23/02/2001  11.33.57  by  Michael Scheer
*CMZ :  2.20/01 24/11/2000  15.17.32  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.49.06  by  Michael Scheer
*CMZ :  2.16/06 27/08/2000  21.54.15  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/00 30/11/99  14.54.00  by  Michael Scheer
*CMZ :  2.10/01 17/03/99  17.04.42  by  Michael Scheer
*CMZ :  2.02/00 12/02/99  18.05.25  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  18.25.54  by  Michael Scheer
*CMZ :  1.03/00 16/01/98  14.52.08  by  Michael Scheer
*CMZ :  1.00/04 21/10/97  13.41.08  by  Michael Scheer
*CMZ :  1.00/01 20/10/97  14.55.47  by  Michael Scheer
*CMZ : 00.01/05 01/02/95  10.49.58  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  17.35.35  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.19.37  by  Michael Scheer
*CMZ : 00.01/01 23/06/94  10.04.05  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.31  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.36  by  Michael Scheer
*-- Author :  Michael Scheer
      SUBROUTINE WFILL0
*KEEP,gplhint.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.

      use clustermod
      use bunchmod

C--- REFERENCE ORBIT IS SCANNED AND SOURCE POINT ARE DETECTED
C    SOURCE POINTS ARE POINTS OF WHICH LIGHT PASSES COLLIMATOR

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,specdip.
      include 'specdip.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      INTEGER INSIDE,IPOI,IOVER,IBFCNT,IL,JLIGHT,I,J,ISA,ISE,ICEN,IBSIGN
      INTEGER IBUFFP,lun,ins
      PARAMETER (IBUFFP=1000)
      INTEGER IBUFFT(IBUFFP,2)

      DOUBLE PRECISION FLEFT,FRIGHT,FUP,FDOWN,BLEFT,BRIGHT,BUP,BDOWN
      DOUBLE PRECISION BMAG2,B0CUT2,BCUT2,BHYS2,TOPANG,OPANG
      DOUBLE PRECISION DIRVER0,DIRVERU,DIRVERD
      DOUBLE PRECISION DIRHOR0,DIRHORL,DIRHORR
      DOUBLE PRECISION DVERUF,DVERDF,DHORLF,DHORRF,DVERUB,DVERDB,DHORLB,DHORRB
      DOUBLE PRECISION DHOLD,DVOLD,DXF,dtf,DXB,X0,Y0,Z0
      DOUBLE PRECISION ANGMXHF,ANGMXVF,ANGMXHB,ANGMXVB
      DOUBLE PRECISION BY,BYOLD,DYU,DYD,DZR,DZL
      DOUBLE PRECISION X1,X2,Y1,Y2,Z1,Z2
      DOUBLE PRECISION VX1,VY1,VZ1
      double precision dxanf,dxend

      integer IBUFF,ixanf,ixend,jxanf,jxend

      integer ISTAT
      double precision, dimension (:), allocatable :: xbuff,ybuff

      IF (IUNDULATOR.EQ.0
     &    .AND.
     &    IAMPLI.EQ.0
     &    .AND.
     &    icluster.EQ.0
     &    .AND
     &    .IBL0CUT.GT.-2) THEN

        INSIDE=0  !IF CURRENT POINT OF CO INSIDE SOURCE, INSIDE=1
        IBFCNT=0  !BUFFER COUNTER

C---  LIMITS OF THE TWO PINHOLES

        IF(WID1*WID2*HIG1*HIG2.LE.0. .OR. CX1.GT.CX2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WFILL0 ***'
          WRITE(LUNGFO,*)'COLLIMATOR HAS MEANINGLESS DIMENSIONS'
          WRITE(LUNGFO,*)'CHECK VARIABLES CX1,CX2,WID1,WID2,HIG1,HIG2'
          WRITE(LUNGFO,*)'IN NAMELIST PINHOLE'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WFILL0 ***'
          WRITE(6,*)'COLLIMATOR HAS MEANINGLESS DIMENSIONS'
          WRITE(6,*)'CHECK VARIABLES CX1,CX2,WID1,WID2,HIG1,HIG2'
          WRITE(6,*)'IN NAMELIST PINHOLE'
          WRITE(6,*)
          STOP
        ENDIF

        FLEFT=CZ1-WID1/2.
        FRIGHT=CZ1+WID1/2.
        FUP=CY1+HIG1/2.
        FDOWN=CY1-HIG1/2.

        BLEFT=CZ2-WID2/2.
        BRIGHT=CZ2+WID2/2.
        BUP=CY2+HIG2/2.
        BDOWN=CY2-HIG2/2.

C--- LOOP OVER ALL POINTS OF THE REFERENCE ORBIT, THE DIRECTION OF THE
C     ELECTRON IS CALCULATED, THE OPENING ANGLE OF THE LIGHT IS ADDED,
C     THIS NEW DIRECTION IS EXTRAPOLATED TO EACH PINHOLE, THE TRAJECTORY
C     POINT IS CHECKED, IF IT BELONGS TO THE SOURCE OR NOT

        BHYS2=WBL0HYS**2

        IF (WBL0HYS.LE.0.D0) BHYS2=1.D0

        B0CUT2=WBL0CUT**2
        BCUT2=B0CUT2*BHYS2 !CUT ON MAGNETIC FIELD
        IF (IBL0CUT.LT.0) BCUT2=BCUT2-1.0D-30
        OPANG=DABS(WGWINFC/DMYGAMMA) !OPENING ANGLE

        IF (OPANG.GT.PI1/4.D0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** WARNING IN WFILL0: WGWINFC CORRESPONDS TO ANGLE GREATER THAN 90 DEGREE'
C         WRITE(LUNGFO,*)'SET TO 45 DEGREE'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)
     &      '*** WARNING IN WFILL0: WGWINFC CORRESPONDS TO ANGLE GREATER THAN 90 DEGREE'
C         WRITE(6,*)'SET TO 45 DEGREE'
          WRITE(6,*)
C         WGWINFC=PI1/4.D0*DMYGAMMA
C         OPANG=PI1/4.D0
        ENDIF !(OPANG.GT.PI1/4.D0)
        TOPANG=DTAN(OPANG) !TANGENS OF OPENING ANGLE

        DHOLD=WTRA(3,2,1)/WTRA(1,2,1)
        DVOLD=WTRA(2,2,1)/WTRA(1,2,1)

        BYOLD=WTRA(2,3,1)

        DO IPOI=1,NCO

C- DIRECTIONS

C           DIRHOR0=VZ/VX, HOR. DIRECTION

          DIRHOR0=WTRA(3,2,IPOI)/WTRA(1,2,IPOI)
          DIRHORR=(DIRHOR0+TOPANG)/(1.D0-DIRHOR0*TOPANG) !TAN(X+Y)=...
          DIRHORL=(DIRHOR0-TOPANG)/(1.D0+DIRHOR0*TOPANG) !TAN(X-Y)=...

C        DIRVER0=VY/VX, VER. DIRECTION

          DIRVER0=WTRA(2,2,IPOI)/WTRA(1,2,IPOI)
          DIRVERU=(DIRVER0+TOPANG)/(1.D0-DIRVER0*TOPANG) !TAN(X+Y)=...
          DIRVERD=(DIRVER0-TOPANG)/(1.D0+DIRVER0*TOPANG) !TAN(X-Y)=...

C- POSITION AND DISTANCES

          X0=WTRA(1,1,IPOI)
          Y0=WTRA(2,1,IPOI)
          Z0=WTRA(3,1,IPOI)

          BY=WTRA(2,3,IPOI)
          IF (IBL0CUT.GT.0.AND.(BY.GE.0.D0 .AND. BYOLD.LT.0.D0
     &        .OR.BY.LE.0.D0 .AND. BYOLD.GT.0.D0)) THEN
            IBSIGN=1
          ELSE
            IBSIGN=0
          ENDIF

          DXF=CX1-X0 !POSITIVE, IF ELECTRON IN FRONT OF THE PINHOLE
          DXB=CX2-X0

          DYD=Y0+DIRVERD*DXF
          DYU=Y0+DIRVERU*DXF
          DZR=Z0+DIRHORR*DXF
          DZL=Z0+DIRHORL*DXF

          DVERUF=FUP-DYD
          DVERDF=FDOWN-DYU
          DHORLF=FLEFT-DZR
          DHORRF=FRIGHT-DZL

C--- IF LIGHT CONE HITS EDGE OF PINHOLE CUT CONE

          IF  (DXF.GT.0.0 .AND.
     &        FUP.LT.DYU
     &        .AND.
     &        FUP.GT.DYD
     &        )
     &        THEN
            DIRVERU=(FUP-Y0)/DXF
          ENDIF

          IF  (DXF.GT.0.0 .AND.
     &        FDOWN.LT.DYU
     &        .AND.
     &        FDOWN.GT.DYD
     &        )
     &        THEN
            DIRVERD=(FDOWN-Y0)/DXF
          ENDIF

          IF  (DXF.GT.0.0 .AND.
     &        FRIGHT.LT.DZR
     &        .AND.
     &        FRIGHT.GT.DZL
     &        )
     &        THEN
            DIRHORR=(FRIGHT-Z0)/DXF
          ENDIF

          IF  (DXF.GT.0.0 .AND.
     &        FLEFT.LT.DZR
     &        .AND.
     &        FLEFT.GT.DZL
     &        )
     &        THEN
            DIRHORL=(FLEFT-Z0)/DXF
          ENDIF

          DVERUB=BUP-(Y0+DIRVERD*DXB)
          DVERDB=BDOWN-(Y0+DIRVERU*DXB)
          DHORLB=BLEFT-(Z0+DIRHORR*DXB)
          DHORRB=BRIGHT-(Z0+DIRHORL*DXB)

          IF  (
     &      DVERUF.GE.0. .AND. DVERDF.LE.0.    .AND.
     &      DHORLF.LE.0. .AND. DHORRF.GE.0.
     &      .AND.
     &      DXF.GE.0
     &      .AND.
     &      DVERUB.GE.0. .AND. DVERDB.LE.0.    .AND.
     &      DHORLB.LE.0. .AND. DHORRB.GE.0.
     &      .AND.
     &      DXB.GE.0
     &        ) THEN

            IOVER=1

          ELSE

            IOVER=0

          ENDIF

C- CHECK SPACING OF POINTS
C  DELTA OF ANGLE MUST BE LOWER THAN THE MAXIMUM OF THE LIGHT CONE
C  AND OPENING ANGLE OF THE PINHOLE

          IF(DXF.NE.0.0) THEN
            ANGMXHF=DMAX1(DABS(WID1/DXF),DABS(TOPANG))
            ANGMXVF=DMAX1(DABS(HIG1/DXF),DABS(TOPANG))
          ENDIF
          IF(DXB.NE.0.0) THEN
            ANGMXHB=DMAX1(DABS(WID2/DXB),DABS(TOPANG))
            ANGMXVB=DMAX1(DABS(HIG2/DXB),DABS(TOPANG))
          ENDIF

          IF(
     &        DABS((DIRHOR0-DHOLD)/(1.D0+DIRHOR0*DHOLD)).GE.ANGMXHF
     &        .OR.
     &        DABS((DIRVER0-DVOLD)/(1.D0+DIRVER0*DVOLD)).GE.ANGMXVF
     &        .OR.
     &        DABS((DIRHOR0-DHOLD)/(1.D0+DIRHOR0*DHOLD)).GE.ANGMXHB
     &        .OR.
     &        DABS((DIRVER0-DVOLD)/(1.D0+DIRVER0*DVOLD)).GE.ANGMXVB
     &        )
     &        THEN

            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WFILL0 ***'
            WRITE(LUNGFO,*)'SPACING OF POINTS IN REFERENCE ORBIT TOO LARGE'
            WRITE(LUNGFO,*)'INCREASE NUMBER OF POINTS OR OPENING ANGLE'
            WRITE(LUNGFO,*)'(WGWINFC IN NAMELIST COLLIN)'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WFILL0 ***'
            WRITE(6,*)'SPACING OF POINTS IN REFERENCE ORBIT TOO LARGE'
            WRITE(6,*)'INCREASE NUMBER OF POINTS OR OPENING ANGLE'
            WRITE(6,*)'(WGWINFC IN NAMELIST COLLIN)'
            WRITE(6,*)
            STOP

          ENDIF

          DHOLD=DIRHOR0
          DVOLD=DIRVER0

          BMAG2=WTRA(1,3,IPOI)*WTRA(1,3,IPOI)+    !SQUARE OF MAG.-FIELD
     &      WTRA(2,3,IPOI)*WTRA(2,3,IPOI)+
     &      WTRA(3,3,IPOI)*WTRA(3,3,IPOI)

          IF (
     &        BMAG2.GT.BCUT2
     &        .AND.
     &        IOVER.EQ.1 .AND. INSIDE.EQ.0
     &        ) THEN

            IBFCNT=IBFCNT+1

            IF (IBFCNT.GT.IBUFFP) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*) '*** ERROR IN WFILL0 ***'
              WRITE(LUNGFO,*) '    TOO MANY SOURCES'
              WRITE(LUNGFO,*) '    INCREASE IBUFFP IN SUBROUTINE WFILL0'
              WRITE(LUNGFO,*) '    OR REDUCED NUMBER OF SOURCES TO BE TREATED'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'     Number of sources:',IBFCNT
              WRITE(LUNGFO,*)'     Begin and end of sources (x,y,z):'
              DO JLIGHT=1,IBFCNT-1
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'     '
     &            ,JLIGHT,(SNGL(WTRA(I,1,IBUFFT(JLIGHT,1))),I=1,3)
                WRITE(LUNGFO,*)'     ',
     &            JLIGHT,(SNGL(WTRA(I,1,IBUFFT(JLIGHT,2))),I=1,3)
              ENDDO
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*) '*** ERROR IN WFILL0 ***'
              WRITE(6,*) '    TOO MANY SOURCES'
              WRITE(6,*) '    INCREASE IBUFFP IN SUBROUTINE WFILL0'
              WRITE(6,*) '    OR REDUCED NUMBER OF SOURCES TO BE TREATED'
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*)'     Number of sources:',IBFCNT
              WRITE(6,*)'     Begin and end of sources (x,y,z):'
              DO JLIGHT=1,IBFCNT-1
                WRITE(6,*)
                WRITE(6,*)'     '
     &            ,JLIGHT,(SNGL(WTRA(I,1,IBUFFT(JLIGHT,1))),I=1,3)
                WRITE(6,*)'     ',
     &            JLIGHT,(SNGL(WTRA(I,1,IBUFFT(JLIGHT,2))),I=1,3)
              ENDDO
              WRITE(6,*)
              STOP
            ENDIF

            IBUFFT(IBFCNT,1)=IPOI
            INSIDE=1
            BCUT2=B0CUT2/BHYS2
            IF (IBL0CUT.LT.0) BCUT2=BCUT2-1.0D-30

          END IF

          IF (
     &        INSIDE.EQ.1 .AND.
     &        (BMAG2.LT.BCUT2.OR.IBSIGN.EQ.1.OR.
     &        IOVER.NE.1 .OR. IPOI.EQ.NCO)
     &        ) THEN

            IBUFFT(IBFCNT,2)=IPOI
            INSIDE=0
            BCUT2=B0CUT2*BHYS2
            IF (IBL0CUT.LT.0) BCUT2=BCUT2-1.0D-30

          END IF

          BYOLD=BY

        END DO !LOOP OVER ALL CO-POINTS

C--- EXTEND SOURCE A LITTLE BIT

        DO IL=1,IBFCNT
          IF (IBUFFT(IL,1).EQ.IBUFFT(IL,2)) THEN
            IBUFFT(IL,1)=MAX(1,IBUFFT(IL,1)-1)
            IBUFFT(IL,2)=MIN(NCO,IBUFFT(IL,2)+1)
          ENDIF
        ENDDO  !IBFCNT

        IF (ISOUREXT.GT.0) THEN
          DO IL=1,IBFCNT
            IBUFFT(IL,1)=MAX(1,IBUFFT(IL,1)-ISOUREXT)
            IBUFFT(IL,2)=MIN(NCO,IBUFFT(IL,2)+ISOUREXT)
          END DO
        ENDIF

C--- CHECK ON DOUBLE COUNTING

        DO IL=1,IBFCNT-1
          IF (IBUFFT(IL,2).GE.IBUFFT(IL+1,1)) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WFILL0 ***'
            WRITE(LUNGFO,*)'SOURCES HAVE COMMON POINTS'
            WRITE(LUNGFO,*)
     &        'PROBABLY THRESHOLD WBL0CUT OR WGWINFC (NAMELIST COLLIN)'
            WRITE(LUNGFO,*)'MUST BE INCREASED'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WFILL0 ***'
            WRITE(6,*)'SOURCES HAVE COMMON POINTS'
            WRITE(6,*)
     &        'PROBABLY THRESHOLD WBL0CUT OR WGWINFC (NAMELIST COLLIN)'
            WRITE(6,*)'MUST BE INCREASED'
            WRITE(6,*)
            IF (ISOUREXT.LT.0) THEN
              WRITE(LUNGFO,*)'*** ERROR IGNORED DUE TO FLAG ISOUREXT ***'
              WRITE(6,*)'*** ERROR IGNORED DUE TO FLAG ISOUREXT ***'
            ELSE
              STOP
            ENDIF
          END IF
        END DO

      ELSE  !IUNDULATOR

        IBFCNT=1
        IBUFFT(1,1)=1
        IBUFFT(1,2)=NCO

      ENDIF !IUNDULATOR

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR WFILL0:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Parameters:'
      WRITE(LUNGFO,*)'     WGWINFC:',WGWINFC
      WRITE(LUNGFO,*)'     WBL0CUT:',WBL0CUT
      WRITE(LUNGFO,*)'     WBL0HYS:',WBL0HYS
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Pinhole of collimators:'
      WRITE(LUNGFO,*)'     (X,Y,Z), width, height:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,'(6H       ,5(1PE12.4))')
     &  SNGL(CX1),SNGL(CY1),SNGL(CZ1),SNGL(WID1),SNGL(HIG1)
      WRITE(LUNGFO,'(6H       ,5(1PE12.4))')
     &  SNGL(CX2),SNGL(CY2),SNGL(CZ2),SNGL(WID2),SNGL(HIG2)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Number of sources:',IBFCNT

      IF (IBFCNT.EQ.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN WFILL0: NO SOURCE FOUND ***'
        WRITE(LUNGFO,*)'CHECK COLLIMATOR, MAG. FIELD...'
        WRITE(LUNGFO,*)'*** PROGRAMM WAVE TERMINATED ***'

        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN WFILL0: NO SOURCE FOUND ***'
        WRITE(6,*)'CHECK COLLIMATOR, MAG. FIELD...'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAMM WAVE TERMINATED ***'

        STOP

      ENDIF !(IBFCNT.GT.0)

      WRITE(LUNGFO,*)'     Begin and end of sources (x,y,z):'
      DO JLIGHT=1,IBFCNT
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     '
     &    ,JLIGHT,(SNGL(WTRA(I,1,IBUFFT(JLIGHT,1))),I=1,3)
        WRITE(LUNGFO,*)'     ',
     &    JLIGHT,(SNGL(WTRA(I,1,IBUFFT(JLIGHT,2))),I=1,3)
      ENDDO
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     Horiz., vert. slopes and mag. field'
      WRITE(LUNGFO,*)
     &  '     at begin and end of sources:'
      DO JLIGHT=1,IBFCNT
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     ',JLIGHT,
     &    SNGL(WTRA(3,2,IBUFFT(JLIGHT,1))/WTRA(1,2,IBUFFT(JLIGHT,1)))
     &    ,SNGL(WTRA(2,2,IBUFFT(JLIGHT,1))/WTRA(1,2,IBUFFT(JLIGHT,1)))
     &    ,SNGL(DSQRT(
     &    WTRA(1,3,IBUFFT(JLIGHT,1))**2+
     &    WTRA(2,3,IBUFFT(JLIGHT,1))**2+
     &    WTRA(3,3,IBUFFT(JLIGHT,1))**2))
        WRITE(LUNGFO,*)'     ',JLIGHT,
     &    SNGL(WTRA(3,2,IBUFFT(JLIGHT,2))/WTRA(1,2,IBUFFT(JLIGHT,2)))
     &    ,SNGL(WTRA(2,2,IBUFFT(JLIGHT,2))/WTRA(1,2,IBUFFT(JLIGHT,2)))
     &    ,SNGL(DSQRT(
     &    WTRA(1,3,IBUFFT(JLIGHT,2))**2+
     &    WTRA(2,3,IBUFFT(JLIGHT,2))**2+
     &    WTRA(3,3,IBUFFT(JLIGHT,2))**2))
      ENDDO
      WRITE(LUNGFO,*)

C--- WRITE FILE FILEL0

      IF(IWFILL0.NE.0) THEN
        IF (LUNL0.NE.6) THEN
          OPEN (UNIT=LUNL0,FILE=FILEL0,FORM ='FORMATTED'
     &      ,STATUS = 'UNKNOWN')
        ENDIF

        WRITE(LUNL0,*)ICODE
        WRITE(LUNL0,*)DMYGAMMA
        WRITE(LUNL0,*)IBFCNT
        WRITE(LUNL0,*)WGWINFC,WBL0CUT,WBL0HYS
        WRITE(LUNL0,*)CX1,CY1,CZ1,WID1,HIG1
        WRITE(LUNL0,*)CX2,CY2,CZ2,WID2,HIG2

        DO JLIGHT=1,IBFCNT
          WRITE(LUNL0,*)IBUFFT(JLIGHT,1),IBUFFT(JLIGHT,2)
          WRITE(LUNL0,*)((WTRA(I,J,IBUFFT(JLIGHT,1)),I=1,3),J=1,4)
          WRITE(LUNL0,*)((WTRA(I,J,IBUFFT(JLIGHT,2)),I=1,3),J=1,4)
        ENDDO

        IF (LUNL0.NE.6) then
          flush(lunl0)
          CLOSE(LUNL0)
        endif

C     WRITE(6,*)'*** SR WFILL0: SOURCE WRITTEN TO FILE ***'
C     WRITE(6,*)'FILE:',FILEL0
C     WRITE(LUNGFO,*)
C     WRITE(LUNGFO,*)
C     WRITE(LUNGFO,*)'*** SR WFILL0: SOURCE WRITTEN TO FILE ***'
C     WRITE(LUNGFO,*)'FILE:',FILEL0
C     WRITE(LUNGFO,*)
      ENDIF

C- COPY SOURCES TO COMMON SOURCE

      NSOURCE=IBFCNT

      IF (ISPECDIP.GT.0.AND.NDIP.GT.NSOURCE) THEN
        ALLOCATE(SOURCEG(3,2,NDIP))
        ALLOCATE(SOURCEA(3,4,NDIP))
        ALLOCATE(SOURCEE(3,4,NDIP))
        ALLOCATE(SOURCEAO(3,4,NDIP))
        ALLOCATE(SOURCEEO(3,4,NDIP))
        ALLOCATE(SOURCEN(3,4,NDIP))
        ALLOCATE(SOURCET(3,NDIP))
        ALLOCATE(ISOURAE(2,NDIP))
        ALLOCATE(ISOURCEN(NDIP))
        ALLOCATE(WTRA2IS(NDIP))
      ELSE  !(ISPECDIP.NE.0.AND.NDIP.GT.NSOURCE)
        ALLOCATE(SOURCEG(3,2,NSOURCE))
        ALLOCATE(SOURCEA(3,4,NSOURCE))
        ALLOCATE(SOURCEE(3,4,NSOURCE))
        ALLOCATE(SOURCEAO(3,4,NSOURCE))
        ALLOCATE(SOURCEEO(3,4,NSOURCE))
        ALLOCATE(SOURCEN(3,4,NSOURCE))
        ALLOCATE(SOURCET(3,NSOURCE))
        ALLOCATE(ISOURAE(2,NSOURCE))
        ALLOCATE(ISOURCEN(NSOURCE))
        ALLOCATE(WTRA2IS(NSOURCE))
      ENDIF !(ISPECDIP.NE.0.AND.NDIP.GT.NSOURCE)

      sourceg=0.0d0
      sourcea=0.0d0
      sourcee=0.0d0
      sourceao=0.0d0
      sourceeo=0.0d0
      sourcet=0.0d0
      wtra2is=0.0d0
      isourae=0
      isourcen=0

      DO JLIGHT=1,NSOURCE

      ISA=IBUFFT(JLIGHT,1)
      ISE=IBUFFT(JLIGHT,2)
      ICEN=ISA+(ISE-ISA)/2

        SOURCEG(1,1,JLIGHT)=WTRA(1,5,ISA)
        SOURCEG(2,1,JLIGHT)=WTRA(2,5,ISA)
        SOURCEG(3,1,JLIGHT)=WTRA(3,5,ISA)

        SOURCEG(1,2,JLIGHT)=WTRA(1,5,ISE)
        SOURCEG(2,2,JLIGHT)=WTRA(2,5,ISE)
        SOURCEG(3,2,JLIGHT)=WTRA(3,5,ISE)

        DO I=1,3

          DO J=1,2
            SOURCEA(I,J,JLIGHT)=WTRA(I,J,ISA)
            SOURCEE(I,J,JLIGHT)=WTRA(I,J,ISE)
            ISOURAE(1,JLIGHT)=ISA
            ISOURAE(2,JLIGHT)=ISE
            ISOURCEN(JLIGHT)=ICEN
            SOURCEN(I,J,JLIGHT)=WTRA(I,J,ICEN)
          ENDDO   !J=1,2

          SOURCEA(I,4,JLIGHT)=WTRA(I,3,ISA)
          SOURCEE(I,4,JLIGHT)=WTRA(I,3,ISE)
          ISOURCEN(JLIGHT)=ICEN
          SOURCEN(I,4,JLIGHT)=WTRA(I,3,ICEN)

          SOURCET(1,JLIGHT)=WTIM0(ISA)
          SOURCET(2,JLIGHT)=WTIM0(ISE)
          SOURCET(3,JLIGHT)=WTIM0(ICEN)

      ENDDO !I=1,3

c Wegen IAMPLI, damit MYINUM keinen Einfluss hat. 12.8.2009,
c Nachtrag: Hat es aber trotzdem, da wtra2is von MYINUM abhaengt.

        DXF=XSTOP-SOURCEE(1,1,JLIGHT)
        DTF=DXF/WTRA(1,2,ISE)

        IF (DXF.LT.1.0D-6) THEN
          SOURCEE(1,1,JLIGHT)=XSTOP
          SOURCEE(2,1,JLIGHT)=SOURCEE(2,1,JLIGHT)
     &      +SOURCEE(2,2,JLIGHT)*DTF
          SOURCEE(3,1,JLIGHT)=SOURCEE(3,1,JLIGHT)
     &      +SOURCEE(3,2,JLIGHT)*DTF
          SOURCET(2,JLIGHT)=SOURCET(2,JLIGHT)+DTF
        ENDIF

        WTRA2IS(JLIGHT)=0.0D0

        IBUFF=0
        DXANF=1.0D30
        DXEND=1.0D30
        ixanf=-1
        ixend=-1

        allocate(xbuff(ise-isa+1))
        allocate(ybuff(ise-isa+1))

        DO IPOI=ISA,ISE

          X1=WTRA(1,1,IPOI)
          Y1=WTRA(2,1,IPOI)
          Z1=WTRA(3,1,IPOI)

          if (ipoi.lt.nco) then
            VX1=(WTRA(1,2,IPOI)+wtra(1,2,ipoi+1))/2.0d0
            Vy1=(WTRA(2,2,IPOI)+wtra(2,2,ipoi+1))/2.0d0
            Vz1=(WTRA(3,2,IPOI)+wtra(3,2,ipoi+1))/2.0d0
          else
            VX1=(WTRA(1,2,IPOI)+wtra(1,2,ipoi-1))/2.0d0
            Vy1=(WTRA(2,2,IPOI)+wtra(2,2,ipoi-1))/2.0d0
            Vz1=(WTRA(3,2,IPOI)+wtra(3,2,ipoi-1))/2.0d0
          endif

          dxf=abs(x1-xianf)
          IF (dxf.LT.DXANF) then
            DXANF=dxf
            jxanf=ipoi
          endif
          if (dxf.lt.1.0d0-9) ixanf=ipoi

          dxf=abs(x1-xiend)
          IF (dxf.LT.DXend) then
            DXend=dxf
            jxend=ipoi
          endif
          if (dxf.lt.1.0d0-9) ixend=ipoi

          IF (X1.GE.XIANF-1.0d-9.AND.X1.LE.XIEND+1.0d-9) THEN
            IBUFF=IBUFF+1
            xbuff(ibuff)=wtim0(ipoi)
            ybuff(ibuff)=0.5d0*vx1*((vy1/vx1)**2+(vz1/vx1)**2)
          endif

        enddo

        if (xianf.ge.xstart.and.ixanf.lt.0) then
          write(lungfo,*)
          write(lungfo,*)
     &      '*** Warning in WFILL0: XIANF NOT POINT OF THE TRAJECTORY ARRAY'
          write(lungfo,*)'closest x:',wtra(1,1,jxanf)
          write(lungfo,*)
          write(lungfo,*)'ROIS may be used to force XIANF onto trajectory'
          write(6,*)
          write(6,*)
     &      '*** Warning in WFILL0: XIANF NOT POINT OF THE TRAJECTORY ARRAY'
          write(6,*)'closest x:',wtra(1,1,jxanf)
          write(6,*)
          write(6,*)'ROIS may be used to force XIANF onto trajectory'
          write(6,*)
        endif

        if (xiend.le.xstop.and.ixend.lt.0) then
          write(lungfo,*)
          write(lungfo,*)
     &      '*** Warning in WFILL0: XIEND NOT POINT OF THE TRAJECTORY ARRAY'
          write(lungfo,*)'closest x:',wtra(1,1,jxend)
          write(lungfo,*)
          write(lungfo,*)'ROIS may be used to force XIEND onto trajectory'
          write(lungfo,*)
          write(6,*)
          write(6,*)
     &      '*** Warning in WFILL0: XIANF NOT POINT OF THE TRAJECTORY ARRAY'
          write(6,*)'closest x:',wtra(1,1,jxend)
          write(6,*)
          write(6,*)
          write(6,*)'ROIS may be used to force XEND onto trajectory'
        endif

        call util_higher_simpson_equidist_integral(ibuff,xbuff,ybuff
     &    ,dxf,istat)

        WTRA2IS(JLIGHT)=0.0D0
      DO IPOI=ISA,ISE-1
          X1=WTRA(1,1,IPOI)
          Y1=WTRA(2,1,IPOI)
          Z1=WTRA(3,1,IPOI)
          X2=WTRA(1,1,IPOI+1)
          Y2=WTRA(2,1,IPOI+1)
          Z2=WTRA(3,1,IPOI+1)
          IF (X2.NE.X1
     &     .AND.X1.GE.XIANF.AND.X1.LE.XIEND
     &     .AND.X2.GE.XIANF.AND.X2.LE.XIEND
     &      ) THEN
         WTRA2IS(JLIGHT)=WTRA2IS(JLIGHT)
     &           +0.5D0*(X2-X1)*(((Y2-Y1)/(X2-X1))**2+((Z2-Z1)/(X2-X1))**2)
          ENDIF    !(X2.NE.X1)
      ENDDO !ISA,ISE

        wtra2is(jlight)=(wtra2is(jlight)+dxf)/2.0d0 ! Mittelung scheint besser zu sein.

        deallocate(xbuff)
        deallocate(ybuff)

      ENDDO   !NSOURCE

      IF (ISPECDIP.GT.0.AND.NDIP.GT.NSOURCE) NSOURCE=NDIP

      if (icluster.lt.0.and.ibunch.ne.0.and.iubunch.eq.3) then
        open(newunit=lun,file='wave.ins')
        read(lun,*) ins
        close(lun)
        if (ins.eq.0) then
          open(newunit=lun,file='wave_source.clu')
          write(lun,*) sourcea,sourcee
          flush(lun)
          close(lun)
          open(newunit=lun,file='wave.status')
          write(lun,*) '0'
          flush(lun)
          close(lun)
          stop '--- WAVE terminated in WFILL0 due to ICLUSTER ---'
        else
          open(newunit=lun,file='wave_source.clu')
          read(lun,*) sourceaclu,sourceeclu
          close(lun)
        endif
      endif

      RETURN
      END
