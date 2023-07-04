*CMZ :  4.00/14 30/12/2021  11.43.05  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  18.47.30  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.02/05 22/03/2015  19.39.25  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  15.19.42  by  Michael Scheer
*CMZ :  3.01/06 18/06/2014  09.18.20  by  Michael Scheer
*CMZ :  3.01/05 12/06/2014  08.52.10  by  Michael Scheer
*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  15.31.59  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  11.51.13  by  Michael Scheer
*CMZ :  2.68/02 15/06/2012  15.18.04  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.47/23 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.47/22 03/12/2003  13.17.28  by  Michael Scheer
*CMZ :  2.44/01 04/12/2002  17.24.27  by  Michael Scheer
*CMZ :  2.44/00 08/11/2002  12.15.59  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.34/08 17/09/2001  19.11.02  by  Michael Scheer
*CMZ :  2.34/05 23/08/2001  17.35.07  by  Michael Scheer
*CMZ :  2.34/01 02/07/2001  17.22.17  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.25.47  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  18.22.03  by  Michael Scheer
*CMZ :  1.00/00 29/07/97  16.33.37  by  Michael Scheer
*CMZ : 00.02/05 17/04/97  17.59.48  by  Michael Scheer
*CMZ : 00.02/04 12/02/97  14.52.46  by  Michael Scheer
*CMZ : 00.02/00 28/11/96  17.04.35  by  Michael Scheer
*-- Author :    Michael Scheer   11/11/96
      SUBROUTINE BMESSINI
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

*KEEP,bmessf90u.
      include 'bmessf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,bmessf90.
      include 'bmessf90.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      CHARACTER(65) COMMENT

      INTEGER IREAD,IX,IY,IZ,ICAL,I,ISTAT,IERR
      INTEGER IBMESSX,IBMESSY,IBMESSZ
      INTEGER NOENT,NTUP,ICODEBM

      REAL*4 X,Y,Z,BX,BY,BZ,XO,YO,ZO,TUP(6)
      DOUBLE PRECISION DUM,DUMO

      INTEGER NTUP_P,ICYCLE
      PARAMETER (NTUP_P=6)
      REAL*8 TUP_D(NTUP_P)
      CHARACTER(3) CHTAGS_D(NTUP_P)
      data chtags_d/'x','y','z','bx','by','bz'/

      DATA ICAL/0/
      DATA XO,YO,ZO/3*0.D0/
      DATA NTUP/599/

      IF (ICAL.NE.0) RETURN

      NBMESSX=9999
      NBMESSY=9999
      NBMESSZ=9999

      OPEN(UNIT=99,STATUS='SCRATCH')

      ICAL=1

      IF (IHLIMIT_C.EQ.0) THEN
        CALL hlimitm(mhbookp)
        IHLIMIT_C=1
      ENDIF

      BMESSEPS=0.00001

      BMXMIN=1.D30
      BMXMAX=-1.D30
      BMYMIN=1.D30
      BMYMAX=-1.D30
      BMZMIN=1.D30
      BMZMAX=-1.D30

      BMBXMIN=1.D30
      BMBXMAX=-1.D30
      BMBYMIN=1.D30
      BMBYMAX=-1.D30
      BMBZMIN=1.D30
      BMBZMAX=-1.D30

      BMESSDX=9999.
      BMESSDY=9999.
      BMESSDZ=9999.

      IF (NTUPGRID.LE.0) THEN

        IF (NTUPGRID.LT.0) THEN
          CALL HISINI
          CALL hbookm(NTUP,'WBMAP',NTUP_P,'//WAVE',
     &      nbmessx*nbmessy*nbmessz,CHTAGS_D)
        ENDIF

        OPEN(UNIT=LUNB0,FILE=FILEB0,STATUS='OLD',FORM='UNFORMATTED',ERR=9999)
        REWIND(LUNB0)

        READ(LUNB0,ERR=901)ICODEBM,COMMENT
        GOTO 902
901     WRITE(6,*)ICODEBM,COMMENT
902     CONTINUE
        READ(LUNB0,ERR=903)NBMESSX,BMXMIN,BMXMAX
        GOTO 904
903     WRITE(6,*)NBMESSX,BMXMIN,BMXMAX
904     CONTINUE
        READ(LUNB0)NBMESSY,BMYMIN,BMYMAX
        READ(LUNB0)NBMESSZ,BMZMIN,BMZMAX

        ALLOCATE(BDATA(3,NBMESSZ,NBMESSY,NBMESSX))
        DO IX=1,NBMESSX
          DO IY=1,NBMESSY
            DO IZ=1,NBMESSZ
              DO I=1,3
                BDATA(I,IZ,IY,IX)=-9999.
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        BMESSDX=(BMXMAX-BMXMIN)/(NBMESSX-1)
        BMESSDY=(BMYMAX-BMYMIN)/(NBMESSY-1)
        BMESSDZ=(BMZMAX-BMZMIN)/(NBMESSZ-1)

        DO IBMESSX=1,NBMESSX
          DO IBMESSY=1,NBMESSY
            DO IBMESSZ=1,NBMESSZ

              READ(LUNB0)BX,BY,BZ

C*** ATTENTION LAST INDEX IS IN LONGITUDINAL DIRECTION (I.E. X)

              BDATA(1,IBMESSZ,IBMESSY,IBMESSX)=BX
              BDATA(2,IBMESSZ,IBMESSY,IBMESSX)=BY
              BDATA(3,IBMESSZ,IBMESSY,IBMESSX)=BZ

              IF (BX.LT.BMBXMIN) BMBXMIN=BX
              IF (BX.GT.BMBXMAX) BMBXMAX=BX
              IF (BY.LT.BMBYMIN) BMBYMIN=BY
              IF (BY.GT.BMBYMAX) BMBYMAX=BY
              IF (BZ.LT.BMBZMIN) BMBZMIN=BZ
              IF (BZ.GT.BMBZMAX) BMBZMAX=BZ

            ENDDO
          ENDDO
        ENDDO

        IF (NTUPGRID.LT.0) THEN

          DO IBMESSX=1,NBMESSX
            DO IBMESSY=1,NBMESSY
              DO IBMESSZ=1,NBMESSZ

                TUP_D(1)=BMXMIN+(IBMESSX-1)*BMESSDX
                TUP_D(2)=BMYMIN+(IBMESSY-1)*BMESSDY
                TUP_D(3)=BMZMIN+(IBMESSZ-1)*BMESSDZ
                TUP_D(4)=BDATA(1,IBMESSZ,IBMESSY,IBMESSX)
                TUP_D(5)=BDATA(2,IBMESSZ,IBMESSY,IBMESSX)
                TUP_D(6)=BDATA(3,IBMESSZ,IBMESSY,IBMESSX)
                CALL hfm(NTUP,TUP_D)

              ENDDO
            ENDDO
          ENDDO

          CALL MHROUT(NTUP,ICYCLE,' ')
          CALL hdeletm(NTUP)

        ENDIF   !NTUPGRID

        CLOSE(LUNB0)

        IREAD=NBMESSX*NBMESSY*NBMESSZ

      ELSE  !NTUPGRID

        CALL hropenm(LUNB0,'BMAP',FILEB0,' ',1024,ISTAT)

        IF (ISTAT.NE.0) THEN
          WRITE(6,*)
     &      '*** SR BMESSINI: ERROR OPENING NTUPLE FILE FILEB0'
          STOP
        ENDIF

        CALL hgnparm(NTUP,'BMESSINI')
        CALL hnoentm(NTUP,NOENT)

        DO IREAD=1,NOENT

          CALL hgnfm(NTUP,IREAD,TUP,IERR)
          IF (IERR.NE.0) THEN
            WRITE(6,*)'*** SR BMESSINI: ERROR READING NTUPLE'
            STOP
          ENDIF

          X=TUP(1)
          Y=TUP(2)
          Z=TUP(3)
          BX=TUP(4)
          BY=TUP(5)
          BZ=TUP(6)

          IF (IREAD.EQ.1) THEN
            XO=X
            YO=Y
            ZO=Z
          ENDIF

          IF (X.NE.XO) THEN
            IF (BMESSDX.EQ.9999.) THEN
              WRITE(99,*)X,XO
              BACKSPACE(99)
              READ(99,*)DUM,DUMO
              BMESSDX=DABS(DUM-DUMO)
            ENDIF
          ENDIF

          IF (Y.NE.YO) THEN
            IF (BMESSDY.EQ.9999.) THEN
              WRITE(99,*)Y,YO
              BACKSPACE(99)
              READ(99,*)DUM,DUMO
              BMESSDY=DABS(DUM-DUMO)
            ENDIF
          ENDIF

          IF (Z.NE.ZO) THEN
            IF (BMESSDZ.EQ.9999.) THEN
              WRITE(99,*)Z,ZO
              BACKSPACE(99)
              READ(99,*)DUM,DUMO
              BMESSDZ=DABS(DUM-DUMO)
            ENDIF
          ENDIF

          IF (X.LT.BMXMIN) BMXMIN=X
          IF (X.GT.BMXMAX) BMXMAX=X
          IF (Y.LT.BMYMIN) BMYMIN=Y
          IF (Y.GT.BMYMAX) BMYMAX=Y
          IF (Z.LT.BMZMIN) BMZMIN=Z
          IF (Z.GT.BMZMAX) BMZMAX=Z

          IF (BX.LT.BMBXMIN) BMBXMIN=BX
          IF (BX.GT.BMBXMAX) BMBXMAX=BX
          IF (BY.LT.BMBYMIN) BMBYMIN=BY
          IF (BY.GT.BMBYMAX) BMBYMAX=BY
          IF (BZ.LT.BMBZMIN) BMBZMIN=BZ
          IF (BZ.GT.BMBZMAX) BMBZMAX=BZ

          IBMESSX=NINT((X-BMXMIN)/BMESSDX)
          IF (ABS(X-(BMXMIN+IBMESSX*BMESSDX)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI:'
            WRITE(6,*)'BAD X-VALUE:',X
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
            WRITE(LUNGFO,*)'BAD X-VALUE:',X
            STOP
          ENDIF
          IBMESSX=IBMESSX+1

          IBMESSY=NINT((Y-BMYMIN)/BMESSDY)
          IF (ABS(Y-(BMYMIN+IBMESSY*BMESSDY)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI:'
            WRITE(6,*)'BAD Y-VALUE:',Y
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
            WRITE(LUNGFO,*)'BAD Y-VALUE:',Y
            STOP
          ENDIF
          IBMESSY=IBMESSY+1

          IBMESSZ=NINT((Z-BMZMIN)/BMESSDZ)
          IF (ABS(Z-(BMZMIN+IBMESSZ*BMESSDZ)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI:'
            WRITE(6,*)'BAD Z-VALUE:',Z
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
            WRITE(LUNGFO,*)'BAD Z-VALUE:',Z
            STOP
          ENDIF
          IBMESSZ=IBMESSZ+1

        ENDDO  !IREAD

        IREAD=NOENT
        NBMESSX=NINT((BMXMAX-BMXMIN)/BMESSDX)+1
        NBMESSY=NINT((BMYMAX-BMYMIN)/BMESSDY)+1
        NBMESSZ=NINT((BMZMAX-BMZMIN)/BMESSDZ)+1

        if (iroottrees.ge.0) then
          CALL hrendm('BMAP')
        endif
        CLOSE(LUNB0)

        ALLOCATE(BDATA(3,NBMESSZ,NBMESSY,NBMESSX))
        DO IX=1,NBMESSX
          DO IY=1,NBMESSY
            DO IZ=1,NBMESSZ
              DO I=1,3
                BDATA(I,IZ,IY,IX)=-9999.
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        CALL hropenm(LUNB0,'BMAP',FILEB0,' ',1024,ISTAT)

        IF (ISTAT.NE.0) THEN
          WRITE(6,*)
     &      '*** SR BMESSINI: ERROR OPENING NTUPLE FILE FILEB0'
          STOP
        ENDIF

        CALL hgnparm(NTUP,'BMESSINI')
        CALL hnoentm(NTUP,NOENT)

        DO IREAD=1,NOENT

          CALL hgnfm(NTUP,IREAD,TUP,IERR)
          IF (IERR.NE.0) THEN
            WRITE(6,*)'*** SR BMESSINI: ERROR READING NTUPLE'
            STOP
          ENDIF

          X=TUP(1)
          Y=TUP(2)
          Z=TUP(3)
          BX=TUP(4)
          BY=TUP(5)
          BZ=TUP(6)

          IF (IREAD.EQ.1) THEN
            XO=X
            YO=Y
            ZO=Z
          ENDIF

          IF (X.NE.XO) THEN
            IF (BMESSDX.EQ.9999.) THEN
              WRITE(99,*)X,XO
              BACKSPACE(99)
              READ(99,*)DUM,DUMO
              BMESSDX=DABS(DUM-DUMO)
            ENDIF
          ENDIF

          IF (Y.NE.YO) THEN
            IF (BMESSDY.EQ.9999.) THEN
              WRITE(99,*)Y,YO
              BACKSPACE(99)
              READ(99,*)DUM,DUMO
              BMESSDY=DABS(DUM-DUMO)
            ENDIF
          ENDIF

          IF (Z.NE.ZO) THEN
            IF (BMESSDZ.EQ.9999.) THEN
              WRITE(99,*)Z,ZO
              BACKSPACE(99)
              READ(99,*)DUM,DUMO
              BMESSDZ=DABS(DUM-DUMO)
            ENDIF
          ENDIF

          IF (X.LT.BMXMIN) BMXMIN=X
          IF (X.GT.BMXMAX) BMXMAX=X
          IF (Y.LT.BMYMIN) BMYMIN=Y
          IF (Y.GT.BMYMAX) BMYMAX=Y
          IF (Z.LT.BMZMIN) BMZMIN=Z
          IF (Z.GT.BMZMAX) BMZMAX=Z

          IF (BX.LT.BMBXMIN) BMBXMIN=BX
          IF (BX.GT.BMBXMAX) BMBXMAX=BX
          IF (BY.LT.BMBYMIN) BMBYMIN=BY
          IF (BY.GT.BMBYMAX) BMBYMAX=BY
          IF (BZ.LT.BMBZMIN) BMBZMIN=BZ
          IF (BZ.GT.BMBZMAX) BMBZMAX=BZ

          IBMESSX=NINT((X-BMXMIN)/BMESSDX)
          IF (ABS(X-(BMXMIN+IBMESSX*BMESSDX)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI:'
            WRITE(6,*)'BAD X-VALUE:',X
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
            WRITE(LUNGFO,*)'BAD X-VALUE:',X
            STOP
          ENDIF
          IBMESSX=IBMESSX+1

          IBMESSY=NINT((Y-BMYMIN)/BMESSDY)
          IF (ABS(Y-(BMYMIN+IBMESSY*BMESSDY)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI:'
            WRITE(6,*)'BAD Y-VALUE:',Y
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
            WRITE(LUNGFO,*)'BAD Y-VALUE:',Y
            STOP
          ENDIF
          IBMESSY=IBMESSY+1

          IBMESSZ=NINT((Z-BMZMIN)/BMESSDZ)
          IF (ABS(Z-(BMZMIN+IBMESSZ*BMESSDZ)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI:'
            WRITE(6,*)'BAD Z-VALUE:',Z
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
            WRITE(LUNGFO,*)'BAD Z-VALUE:',Z
            STOP
          ENDIF
          IBMESSZ=IBMESSZ+1

C*** ATTENTION LAST INDEX IS IN LONGITUDINAL DIRECTION (I.E. X)

          BDATA(1,IBMESSZ,IBMESSY,IBMESSX)=BX
          BDATA(2,IBMESSZ,IBMESSY,IBMESSX)=BY
          BDATA(3,IBMESSZ,IBMESSY,IBMESSX)=BZ

        ENDDO  !IREAD

        IREAD=NOENT

        if (iroottrees.ge.0) then
          CALL hrendm('BMAP')
        endif
        CLOSE(LUNB0)

      ENDIF !NTUPGRID

      X=BMXMAX
      IBMESSX=NINT((X-BMXMIN)/BMESSDX)

      IF (ABS(X-(BMXMIN+IBMESSX*BMESSDX)).GT.BMESSEPS) THEN
        WRITE(6,*)'*** ERROR IN BMESSINI:'
        WRITE(6,*)'BAD X-VALUE:',X
        WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
        WRITE(LUNGFO,*)'BAD X-VALUE:',X
        STOP
      ENDIF

      IF (NBMESSX.EQ.9999) THEN
        NBMESSX=IBMESSX+1
      ENDIF

      Y=BMYMAX
      IBMESSY=NINT((Y-BMYMIN)/BMESSDY)

      IF (ABS(Y-(BMYMIN+IBMESSY*BMESSDY)).GT.BMESSEPS) THEN
        WRITE(6,*)'*** ERROR IN BMESSINI:'
        WRITE(6,*)'BAD Y-VALUE:',Y
        WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
        WRITE(LUNGFO,*)'BAD Y-VALUE:',Y
        STOP
      ENDIF

      IF (NBMESSY.EQ.9999) THEN
        NBMESSY=IBMESSY+1
      ENDIF

      Z=BMZMAX
      IBMESSZ=NINT((Z-BMZMIN)/BMESSDZ)

      IF (ABS(Z-(BMZMIN+IBMESSZ*BMESSDZ)).GT.BMESSEPS) THEN
        WRITE(6,*)'*** ERROR IN BMESSINI:'
        WRITE(6,*)'BAD Z-VALUE:',Z
        WRITE(LUNGFO,*)'*** ERROR IN BMESSINI:'
        WRITE(LUNGFO,*)'BAD Z-VALUE:',Z
        STOP
      ENDIF

      IF (NBMESSZ.EQ.9999) THEN
        NBMESSZ=IBMESSZ+1
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE BMESSINI:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     DATA FILE:'
      WRITE(LUNGFO,*)'     ',FILEB0

      IF (NTUPGRID.EQ.0) THEN
        WRITE(LUNGFO,*)'     ',ICODEBM,' ',COMMENT
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     NUMBER OF DATA POINTS READ: ',IREAD
      WRITE(LUNGFO,*)'     NBMESSX, NBMESSY, NBMESSZ:  '
     &  ,NBMESSX,NBMESSY,NBMESSZ
      WRITE(LUNGFO,*)'     XMIN, XMAX:   ',BMXMIN,BMXMAX
      WRITE(LUNGFO,*)'     YMIN, YMAX:   ',BMYMIN,BMYMAX
      WRITE(LUNGFO,*)'     ZMIN, ZMAX:   ',BMZMIN,BMZMAX
      WRITE(LUNGFO,*)'     BXMIN, BXMAX: ',BMBXMIN,BMBXMAX
      WRITE(LUNGFO,*)'     BYMIN, BYMAX: ',BMBYMIN,BMBYMAX
      WRITE(LUNGFO,*)'     BZMIN, BZMAX: ',BMBZMIN,BMBZMAX
      WRITE(LUNGFO,*)'     BMESSDX, BMESSDY, BMESSDZ:  '
      WRITE(LUNGFO,*)'     ',BMESSDX,BMESSDY,BMESSDZ
      WRITE(LUNGFO,*)'     EPSILON:                    ',BMESSEPS

      IF (IRFILB0.GT.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     MORD3DG:     ',MORD3DG
        WRITE(LUNGFO,*)'     NBMDATX,NBMDATY,NBMDATZ:',NBMDATX,NBMDATY,NBMDATZ

        IF (NBMESSX.LT.NBMDATX) THEN
          WRITE(6,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE NBMDATX DIFFERENT X-VALUES'
          WRITE(6,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE NBMDATX DIFFERENT X-VALUES'
          WRITE(LUNGFO,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
          STOP
        ENDIF

        IF (NBMESSY.LT.NBMDATY) THEN
          WRITE(6,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE NBMDATY DIFFERENT Y-VALUES'
          WRITE(6,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE NBMDATY DIFFERENT Y-VALUES'
          WRITE(LUNGFO,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          STOP
        ENDIF

        IF (NBMESSZ.LT.NBMDATZ) THEN
          WRITE(6,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE NBMDATZ DIFFERENT Z-VALUES'
          WRITE(6,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE NBMDATZ DIFFERENT Z-VALUES'
          WRITE(LUNGFO,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          STOP
        ENDIF

      ELSE !IRFILB0.GT.0

        IF (NBMESSX.LT.2) THEN
          WRITE(6,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE 2 DIFFERENT X-VALUES'
          WRITE(6,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE 2 DIFFERENT X-VALUES'
          WRITE(LUNGFO,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
          STOP
        ENDIF

        IF (NBMESSY.LT.2) THEN
          WRITE(6,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE 2 DIFFERENT Y-VALUES'
          WRITE(6,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE 2 DIFFERENT Y-VALUES'
          WRITE(LUNGFO,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          STOP
        ENDIF

        IF (NBMESSZ.LT.2) THEN
          WRITE(6,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE 2 DIFFERENT Z-VALUES'
          WRITE(6,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          WRITE(LUNGFO,*)
     &      '*** ERROR IN BMESSINI: MUST HAVE 2 DIFFERENT Z-VALUES'
          WRITE(LUNGFO,*)'CHECK DATA FILE FILEB0 AND NAMELIST BGRIDN'
          STOP
        ENDIF

      ENDIF !IRFILB0.GT.0

      DO IX=1,NBMESSX
        DO IY=1,NBMESSY
          DO IZ=1,NBMESSZ
            IF (
     &          BDATA(1,IZ,IY,IX).EQ.-9999. .OR.
     &          BDATA(2,IZ,IY,IX).EQ.-9999. .OR.
     &          BDATA(3,IZ,IY,IX).EQ.-9999.) THEN
              WRITE(6,*)'*** ERROR IN BMESSINI: STRANGE INPUT FILE'
              WRITE(LUNGFO,*)'*** ERROR IN BMESSINI: STRANGE INPUT FILE'
              STOP
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (XSTART.EQ.9999.) THEN
        IF (IPERIODG.NE.0.AND.SIGNG.LT.0.D0) THEN
          XSTART=2.D0*BMXMIN
        ELSE
          XSTART=BMXMIN
        ENDIF
      ENDIF

      IF (XSTOP.EQ.9999.) THEN
        IF (IPERIODG.NE.0.AND.SIGNG.LT.0.D0) THEN
          XSTOP=2.D0*BMXMAX
        ELSE
          XSTOP=BMXMAX
        ENDIF
      ENDIF

      CLOSE(99)

      RETURN

9999  WRITE(6,*) '*** Error in BMESSINI: Could not open file'
      WRITE(6,*) FILEB0

      STOP
      END
