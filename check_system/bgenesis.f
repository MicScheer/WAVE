*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  12.21.08  by  Michael Scheer
*CMZ :  2.67/00 13/02/2012  10.58.17  by  Michael Scheer
*CMZ :  2.47/12 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.47/10 30/05/2003  08.16.14  by  Michael Scheer
*CMZ :  2.47/09 26/05/2003  16.00.38  by  Michael Scheer
*CMZ :  2.47/08 16/05/2003  15.52.55  by  Michael Scheer
*CMZ :  2.36/00 05/11/2001  11.14.30  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/10 21/08/96  12.30.33  by  Michael Scheer
*CMZ : 00.01/04 30/11/94  14.09.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.18  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.57  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BGENESIS(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT
     &                               ,AXOUT,AYOUT,AZOUT)

*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,genesis.
      include 'genesis.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION
     &  XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  XKX,YKY,ZKZ,DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ,
     &  BXH,BYH,BZH,AXH,AYH,AZH,XKHALBA,YKHALBA,ZKHALBA,B0HALBA,
     &  AW,QF,CONSTAW,FIELD,GSTYPE,GSTYPEO,CURLEN,BX,BY,BZ,
     &  QYCEN,QZCEN,SL,CX,CY,XLHALBA,YLHALBA,ROOT2

      INTEGER LC,NTUPP
      PARAMETER (LC=80,NTUPP=13)

      double precision FILL(NTUPP)

      INTEGER ICAL,I,IH,IL,NREP,NDIST,IGST,ICHA

      CHARACTER(LC) CHEADER,CVERSION,CUNITLENGTH,CLINE
      CHARACTER(2) CHTAGS(NTUPP)
      CHARACTER(80) CHANGE,COMGEN

      DATA ICAL/0/,IL/1/,IH/2/

      data chtags/'x','y','z','bx','by','bz','aw','dp',
     &'qf','ad','sl','cx','cy'/

      IF (ICAL.EQ.0) THEN

        ROOT2=SQRT(2.D0)

        CALL GENESISIN

        OPEN(UNIT=99,FILE='wave_genesis.bnd',STATUS='NEW',FORM='FORMATTED')
        OPEN(UNIT=LUNGENL,FILE=FILEGENL,STATUS='OLD',FORM='FORMATTED',ERR=9)

        READ(LUNGENL,'(A)',ERR=99) CHEADER

        IF (CHEADER.NE.'# Header is included') THEN
          WRITE(6,*)'*** ERROR IN BGENESIS: BAD HEADERLINE ON LATTICE FILE'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: BAD HEADERLINE ON LATTICE FILE'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF

        READ(LUNGENL,'(A)',ERR=99) CVERSION

        READ(LUNGENL,'(A)',ERR=99) CUNITLENGTH

        IF (CUNITLENGTH(1:12).NE.'? UNITLENGTH') THEN
          WRITE(6,*)'*** ERROR IN BGENESIS: LINE FOR UNITLENGTH FILE'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: LINE FOR UNITLENGTH FILE'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF

        READ(CUNITLENGTH(14:21),*) GSUNILEN

        GSLATT(1,1)=GSUNILEN
        DO IH=2,10
          GSLATT(IH,1)=-9999.
        ENDDO

        DO I=2,NGSDIMP
          GSLATT(1,I)=GSLATT(1,I-1)+GSUNILEN
          DO IH=2,10
            GSLATT(IH,I)=-9999.
          ENDDO
        ENDDO

        NGSELEM=0
        GSTOTLEN=0.D0
        GSTYPEO=-1.

11      READ(LUNGENL,'(A)',END=88,ERR=99) CLINE

        IF (CLINE(1:2).EQ.'AW') THEN
          GSTYPE=2.
        ELSE IF (CLINE(1:2).EQ.'DP') THEN
          GSTYPE=3.
        ELSE IF (CLINE(1:2).EQ.'QF') THEN
          GSTYPE=4. ! quadrupol, gradient [T/m],
          !QF > 0 is horizontally focusing
          !QF < 0 is vertically focusing
        ELSE IF (CLINE(1:2).EQ.'QX') THEN
          GSTYPE=5. ! quadrupol offset in x
        ELSE IF (CLINE(1:2).EQ.'QY') THEN
          GSTYPE=6. ! quadrupol offset in y
        ELSE IF (CLINE(1:2).EQ.'AD') THEN
          GSTYPE=7.  ! for the time being (26.5.2003) just a drift
        ELSE IF (CLINE(1:2).EQ.'SL') THEN
          GSTYPE=8. ! solenoid, field strength [T]?
        ELSE IF (CLINE(1:2).EQ.'CX') THEN
          GSTYPE=9. ! pur vertical field, angle [mrad]
        ELSE IF (CLINE(1:2).EQ.'CY') THEN
          GSTYPE=10.  ! pur horizontal field, angle [mrad]
        ELSE
          WRITE(6,*)' '
          WRITE(6,*)
     &      '*** ERROR IN BGENESIS: BAD MAGNET TYPE ON LATTICE FILE'
          WRITE(6,*)CLINE
          WRITE(6,*)' '
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BGENESIS: BAD MAGNET TYPE ON LATTICE FILE'
          WRITE(LUNGFO,*)CLINE
          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
        ENDIF !(CLINE(1:2).EQ.'AW')

        IGST=GSTYPE

        READ(CLINE(3:LC),*) FIELD,NREP,NDIST

        IF (GSTYPE.NE.GSTYPEO) THEN
          IL=0
          CURLEN=0.D0
        ENDIF  ! GSTYPE

        DO I=1,NDIST
          IL=IL+1
          CURLEN=CURLEN+GSUNILEN
          IF (CURLEN.GT.GSTOTLEN) THEN
            GSTOTLEN=CURLEN
            IF (IL.GT.NGSDIMP-1) GOTO 999
            NGSELEM=IL
          ENDIF
        ENDDO !NDIST

        DO I=1,NREP
          IL=IL+1
          CURLEN=CURLEN+GSUNILEN
          IF (CURLEN.GT.GSTOTLEN) THEN
            GSTOTLEN=CURLEN
            IF (IL.GT.NGSDIMP-1) GOTO 999
            NGSELEM=IL
          ENDIF
          IF (GSLATT(IGST,IL).NE.-9999.) THEN
            WRITE(6,*)
     &        '*** ERROR IN BGENESIS: COLLIDING ELEMENTS ON LATTICE FILE'
            WRITE(6,*)IL,IGST,GSLATT(IGST,IL),GSLATT(1,IL)
            WRITE(6,*)CLINE
            WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
            WRITE(LUNGFO,*)
     &        '*** ERROR IN BGENESIS: COLLIDING ELEMENTS ON LATTICE FILE'
            WRITE(LUNGFO,*)IL,IGST,GSLATT(IGST,IL),GSLATT(1,IL)
            WRITE(LUNGFO,*)CLINE
            WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
            STOP
          ENDIF

          GSLATT(IGST,IL)=FIELD

        ENDDO ! NREP

        GSTYPEO=GSTYPE

        GOTO 11

88      CLOSE(LUNGENL)

        IF (ROIX(1).GT.0..OR.ROIX(2).LT.GSTOTLEN) THEN
            WRITE(6,*)' '
            WRITE(6,*)'*** WARNING IN BGENESIS: INCOMPATIBLE ROIS'
            WRITE(6,*)'ROIS ARE SET OR OVERWRITTEN BY BGENESIS'
            WRITE(6,*)'CHECK NAMELIST $ROIN IN INPUT FILE WAVE.IN'
            WRITE(6,*)' '
            WRITE(LUNGFO,*)' '
            WRITE(LUNGFO,*)'*** WARNING IN BGENESIS: INCOMPATIBLE ROIS'
            WRITE(LUNGFO,*)'ROIS ARE SET OR OVERWRITTEN BY BGENESIS'
            WRITE(LUNGFO,*)'CHECK NAMELIST $ROI IN INPUT FILE WAVE.IN'
            WRITE(LUNGFO,*)' '
        ENDIF

      IF (NGSELEM.LE.2) THEN
        WRITE(6,*)'*** ERROR IN BGENESIS:'
        WRITE(6,*)'TOO FEW ELEMENTS ON LATTICE FILE.'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
        WRITE(LUNGFO,*)'*** ERROR IN BGENESIS:'
        WRITE(LUNGFO,*)'TOO FEW ELEMENTS ON LATTICE FILE.'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
        STOP
      ENDIF

      NROI=2
      ROIX(1)=GSLATT(1,1)-GSUNILEN
      ROIX(2)=GSLATT(1,NGSELEM)

      DO I=1,NGSELEM-1

        CHANGE=''
        ICHA=1

        IF(
     &      GSLATT(2,I).NE.GSLATT(2,I+1) .OR. !AW
     &      GSLATT(3,I).NE.GSLATT(3,I+1) .OR. !DP
     &      GSLATT(4,I).NE.GSLATT(4,I+1) .OR. !QF
     &      GSLATT(5,I).NE.GSLATT(5,I+1) .OR. !QX
     &      GSLATT(6,I).NE.GSLATT(6,I+1) .OR. !QY
     &      GSLATT(7,I).NE.GSLATT(7,I+1) .OR. !AD
     &      GSLATT(8,I).NE.GSLATT(8,I+1) .OR. !SL
     &      GSLATT(9,I).NE.GSLATT(9,I+1) .OR. !CX
     &      GSLATT(10,I).NE.GSLATT(10,I+1)
     &      ) THEN

          ROIX(NROI)=GSLATT(1,I)

          IF(GSLATT(2,I).NE.GSLATT(2,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' AW'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(3,I).NE.GSLATT(3,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' DP'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(4,I).NE.GSLATT(4,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' QF'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(5,I).NE.GSLATT(5,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' QX'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(6,I).NE.GSLATT(6,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' QY'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(7,I).NE.GSLATT(7,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' AD'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(8,I).NE.GSLATT(8,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' SL'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(9,I).NE.GSLATT(9,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' CX'
            ICHA=ICHA+3
          ENDIF

          IF(GSLATT(10,I).NE.GSLATT(10,I+1)) THEN
            CHANGE(ICHA:ICHA+2)=' CY'
            ICHA=ICHA+3
          ENDIF

          WRITE(99,*)'N, ROIX(N):',NROI,ROIX(NROI),CHANGE(1:ICHA)

          NROI=NROI+1

          IF (NROI.GT.NROIP) THEN
            WRITE(6,*)' '
            WRITE(6,*)
     &        '*** ERROR IN BGENESIS: DIMENSION NROIP EXCEEDED'
            WRITE(6,*)' '
            WRITE(LUNGFO,*)' '
            WRITE(LUNGFO,*)
     &        '*** ERROR IN BGENESIS: DIMENSION NROIP EXCEEDED'
            WRITE(LUNGFO,*)' '
            STOP
          ENDIF

          ROIX(NROI)=GSLATT(1,NGSELEM)

        ENDIF

      ENDDO !NGSELEM

      CLOSE(99)

      GSZLHAL=GSXLAMD

      IF (GSZLHAL.NE.0.D0) THEN
        ZKHALBA=2.D0*PI1/GSZLHAL
      ELSE
        ZKHALBA=0.D0
      ENDIF

      IF (IGSWITYP.EQ.0) THEN ! planar undulator

        IF (GSXKX.LE.0.0) THEN
          XKHALBA=SQRT(ABS(GSXKX))*ZKHALBA
        ELSE !(GSXKX.LT.0)
          WRITE(6,*)' '
          WRITE(6,*)
     &      '*** WARNING IN BGENESIS: horizontal focusing of undulator > 0'
          WRITE(6,*)' '
          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)
     &      '*** WARNING IN BGENESIS: horizontal focusing of undulator > 0'
          WRITE(LUNGFO,*)' '
        ENDIF !(GSXKX.LT.0)

        IF (GSXKY.GE.0.0) THEN
          YKHALBA=SQRT(GSXKY)*ZKHALBA
        ELSE !(GSXKY.LT.0)
          WRITE(6,*)' '
          WRITE(6,*)
     &      '*** ERROR IN BGENESIS: vertical focusing of undulator < 0'
          WRITE(6,*)' '
          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BGENESIS: vertical focusing of undulator < 0'
          WRITE(LUNGFO,*)' '
          STOP
        ENDIF !(GSXKY.LT.0)

          ELSE !IGSWITYP

            WRITE(6,*)' '
            WRITE(6,*)'*** ERROR IN BGENESIS: ONLY PLANAR UNDULATOR AVAILABLE'
            WRITE(6,*)' '
            WRITE(LUNGFO,*)' '
            WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: ONLY PLANAR UNDULATOR AVAILABLE'
            WRITE(LUNGFO,*)' '
            STOP

          ENDIF !IGSWITYP

          IF (XKHALBA.NE.0.D0) THEN
            XLHALBA=2.D0*PI1/XKHALBA
          ELSE
            XLHALBA=1.D30
          ENDIF

          IF (YKHALBA.NE.0.D0) THEN
            YLHALBA=2.D0*PI1/YKHALBA
          ELSE
            YLHALBA=1.D30
          ENDIF


          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)'First call to BGENESIS:'
          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)'      name of GENESIS input file: ',FILEGENI
          OPEN(UNIT=99,FILE='genesis.out',STATUS='OLD')
            READ(99,'(A)')COMGEN
            READ(99,'(A)')COMGEN
            READ(99,'(A)')COMGEN
            WRITE(LUNGFO,*)COMGEN
            READ(99,'(A)')COMGEN
            READ(99,'(A)')COMGEN
            READ(99,'(A)')COMGEN
            WRITE(LUNGFO,*)COMGEN
          CLOSE(99)
          WRITE(LUNGFO,*)'      AW0, XLAMD:',
     &      SNGL(GSAW0),SNGL(GSXLAMD)
          WRITE(LUNGFO,*)'      XKX, XKY (normalized to ZKZ):',
     &      SNGL(GSXKX),SNGL(GSXKY)
          WRITE(LUNGFO,*)'      corresponding XLHALBA,YLHALBA [m]:',
     &      SNGL(XLHALBA),SNGL(YLHALBA)
          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)'      name of lattice file: ',FILEGENL
          WRITE(LUNGFO,*)'      unitlength, periodlength:',
     &      SNGL(GSUNILEN),SNGL(GSZLHAL)
          WRITE(LUNGFO,*)'      number of items on lattice file:',NGSELEM
          WRITE(LUNGFO,*)'      total length in unitlengths and meter:',
     &      SNGL(GSTOTLEN/GSUNILEN),SNGL(GSTOTLEN)

          IF (ABS(GSZLHAL/GSUNILEN-1.D0) .GT. 1.D-6 .AND.
     &        ABS(GSZLHAL/GSUNILEN-2.D0) .GT. 1.D-6) THEN
            WRITE(6,*)'*** WARNING IN BGENESIS: STRANGE PERIODLENGTH (GSZLHAL)'
            WRITE(6,*)'GENESIS NAMELIST- AND LATTICE-FILES'
            WRITE(LUNGFO,*)'*** WARNING IN BGENESIS: STRANGE PERIODLENGTH (GSZLHAL)'
            WRITE(LUNGFO,*)'GENESIS NAMELIST- AND LATTICE-FILES'
          ENDIF

          CONSTAW=EMASSKG1*CLIGHT1/ECHARGE1

          WRITE(6,*)
          WRITE(LUNGFO,*)

          ICAL=1

        ENDIF !ICAL

        IF (ICAL.EQ.3.AND.IHISINI_C.EQ.0) THEN
          CALL HISINI
          IF (KBGENESIS.GT.0) THEN
            NIDGENESIS=13
            CALL hbookm(NIDGENESIS,'GENESIS$',NTUPP,'//WAVE',1024,CHTAGS)
          ENDIF !KBGENESIS
        ELSE IF (IHISINI_C.EQ.0) THEN
          ICAL=ICAL+1
        ENDIF !ICAL


        IF (XIN.LT.0.D0.OR.XIN.GE.GSLATT(1,NGSELEM)) THEN
          BXOUT=0.D0
          BYOUT=0.D0
          BZOUT=0.D0
          AXOUT=0.D0
          AYOUT=0.D0
          AZOUT=0.D0
          GOTO 8888 ! RETURN
        ENDIF

C get position and index

          IL=XIN/GSUNILEN+1

          BXOUT=0.D0
          BYOUT=0.D0
          BZOUT=0.D0

          AXOUT=0.D0
          AYOUT=0.D0
          AZOUT=0.D0

        IF (GSLATT(2,IL).NE.-9999.) THEN

          IF(IGSWITYP.EQ.0) THEN ! planar device

            AW=GSLATT(2,IL)
            B0HALBA=AW*ROOT2*ZKHALBA*CONSTAW

            XKX=XKHALBA*(-ZIN)
            YKY=YKHALBA*YIN
            ZKZ=ZKHALBA*XIN

            IF (XKX.LT.0.D0) THEN
              DSNXKX=DSIN(XKX)
              DCSXKX=DCOS(XKX)
            ELSE
              DSNXKX=DSINH(XKX)
              DCSXKX=DCOSH(XKX)
            ENDIF

            DSHYKY=DSINH(YKY)
            DCHYKY=DSQRT(1.D0+DSHYKY*DSHYKY)
            DSNZKZ=DSIN(ZKZ)
            DCSZKZ=DCOS(ZKZ)

            BXH=-XKHALBA/YKHALBA*B0HALBA*DSNXKX*DSHYKY*DCSZKZ
            BYH=                 B0HALBA*DCSXKX*DCHYKY*DCSZKZ
            BZH=-ZKHALBA/YKHALBA*B0HALBA*DCSXKX*DSHYKY*DSNZKZ

            AXH=B0HALBA/ZKHALBA*                DCSXKX*DCHYKY*DSNZKZ
            AYH=B0HALBA/ZKHALBA*XKHALBA/YKHALBA*DSNXKX*DSHYKY*DSNZKZ
            AZH=0.

            BZOUT=BZOUT-BXH
            BYOUT=BYOUT+BYH
            BXOUT=BXOUT+BZH

            AZOUT=AZOUT-AXH
            AYOUT=AYOUT+AYH
            AXOUT=AXOUT+AZH

          ELSE !IGSWITYP

            WRITE(6,*)' '
            WRITE(6,*)'*** ERROR IN BGENESIS: ONLY PLANAR UNDULATOR AVAILABLE'
            WRITE(6,*)' '
            WRITE(LUNGFO,*)' '
            WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: ONLY PLANAR UNDULATOR AVAILABLE'
            WRITE(LUNGFO,*)' '
            STOP

          END IF !IGSWITYP

        ENDIF !AW

C        IF (GSLATT(3,IL).NE.-9999.) THEN
C
C          BX=0.D0
C          BY=0.D0
C          BZ=0.D0
C
C          BXOUT=BXOUT+BX
C          BYOUT=BYOUT+BY
C          BZOUT=BZOUT+BZ
C
C        ENDIF !DP

        IF (GSLATT(4,IL).NE.-9999.) THEN

          QF=GSLATT(4,IL) ! QUADRUPOL GRADIENT [T]

          IF (GSLATT(5,IL).NE.-9999.) THEN
            QZCEN=GSLATT(5,IL)
          ELSE
            QZCEN=GSLATT(5,IL)
          ENDIF !QX

          IF (GSLATT(6,IL).NE.-9999.) THEN
            QYCEN=GSLATT(6,IL)
          ELSE
            QYCEN=GSLATT(6,IL)
          ENDIF !QY

          BX=0.D0
          BY=QF*(ZIN-QZCEN)
          BZ=QF*(YIN-QYCEN)

          BXOUT=BXOUT+BX
          BYOUT=BYOUT+BY
          BZOUT=BZOUT+BZ

        ENDIF !QF


C        IF (GSLATT(7,IL).NE.-9999.) THEN
C
C          BX=0.D0
C          BY=0.D0
C          BZ=0.D0
C
C          BXOUT=BXOUT+BX
C          BYOUT=BYOUT+BY
C          BZOUT=BZOUT+BZ
C
C        ENDIF !AD

        IF (GSLATT(8,IL).NE.-9999.) THEN

          SL=GSLATT(8,IL)

          BX=SL
          BY=0.D0
          BZ=0.D0

          BXOUT=BXOUT+BX
          BYOUT=BYOUT+BY
          BZOUT=BZOUT+BZ

        ENDIF !SL

        IF (GSLATT(9,IL).NE.-9999.) THEN

          CX=GSLATT(9,IL)/1000.D0 ! mrad to rad

          BX=0.D0
          BY=CX*DBRHO/GSUNILEN
          BZ=0.D0

          BXOUT=BXOUT+BX
          BYOUT=BYOUT+BY
          BZOUT=BZOUT+BZ

        ENDIF !CX

        IF (GSLATT(10,IL).NE.-9999.) THEN

          CY=GSLATT(10,IL)/1000.D0

          BX=0.D0
          BY=0.D0
          BZ=CY*DBRHO/GSUNILEN

          BXOUT=BXOUT+BX
          BYOUT=BYOUT+BY
          BZOUT=BZOUT+BZ

        ENDIF !CY

8888  CONTINUE

      IF (KBGENESIS.GT.0.AND.IHISINI_C.NE.0) THEN
        FILL(1)=XIN
        FILL(2)=YIN
        FILL(3)=ZIN
        FILL(4)=BXOUT
        FILL(5)=BYOUT
        FILL(6)=BZOUT
        FILL(7)=GSLATT(2,IL) !AW
        FILL(8)=GSLATT(3,IL) !DP
        FILL(9)=GSLATT(4,IL) !QF
        FILL(10)=GSLATT(7,IL) !AD
        FILL(11)=GSLATT(8,IL) !SL
        FILL(12)=GSLATT(9,IL) !CX
        FILL(13)=GSLATT(10,IL) !CY
        CALL hfm(NIDGENESIS,FILL)
      ENDIF !KBGENESIS

      RETURN

9       WRITE(6,*)'*** ERROR IN BGENESIS: LATTICE FILE NOT FOUND'
        WRITE(6,*)'FILE: ',FILEGENL
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'

        WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: LATTICE FILE NOT FOUND'
        WRITE(LUNGFO,*)'FILE: ',FILEGENL
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'

        STOP

99      WRITE(6,*)'*** ERROR IN BGENESIS: BAD LINE IN LATTICE FILE'
        WRITE(6,*)CLINE
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'

        WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: BAD LINE IN LATTICE FILE'
        WRITE(LUNGFO,*)CLINE
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'

        STOP

999     WRITE(6,*)'*** ERROR IN BGENESIS: DIMENSION EXCEEDED'
        WRITE(6,*)'TOO MANY ELEMENTS ON LATTICE FILE.'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'

        WRITE(LUNGFO,*)'*** ERROR IN BGENESIS: DIMENSION EXCEEDED'
        WRITE(LUNGFO,*)'TOO MANY ELEMENTS ON LATTICE FILE.'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'

        STOP

      END
