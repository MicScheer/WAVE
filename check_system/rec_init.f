*CMZ :  4.00/15 27/04/2022  08.35.57  by  Michael Scheer
*CMZ :  4.00/13 10/11/2021  16.56.28  by  Michael Scheer
*CMZ :  4.00/07 04/06/2020  16.23.19  by  Michael Scheer
*CMZ :  4.00/04 17/05/2019  14.22.20  by  Michael Scheer
*CMZ :  3.06/00 18/02/2019  17.12.04  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.03/04 29/11/2017  11.56.17  by  Michael Scheer
*CMZ :  3.03/02 04/03/2016  18.03.53  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.08.57  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  07.54.21  by  Michael Scheer
*CMZ :  2.63/04 22/07/2009  07.39.04  by  Michael Scheer
*CMZ :  2.63/03 02/06/2009  16.18.47  by  Michael Scheer
*CMZ :  2.62/02 06/06/2007  11.45.09  by  Michael Scheer
*CMZ :  2.61/05 12/04/2007  09.30.49  by  Michael Scheer
*CMZ :  2.61/04 29/03/2007  16.05.58  by  Michael Scheer
*CMZ :  2.58/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.57/05 14/12/2006  10.21.21  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.42/00 03/09/2002  15.19.21  by  Michael Scheer
*CMZ :  2.41/13 03/09/2002  14.27.45  by  Michael Scheer
*CMZ :  2.41/12 22/08/2002  10.56.13  by  Michael Scheer
*CMZ :  2.41/11 21/08/2002  10.55.08  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.41/08 14/08/2002  17.00.53  by  Michael Scheer
*CMZ :  2.39/00 14/01/2002  14.35.25  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.20/10 04/04/2001  18.24.16  by  Michael Scheer
*CMZ :  2.20/09 03/04/2001  14.30.31  by  Michael Scheer
*CMZ :  2.15/01 28/03/2001  13.23.50  by  Michael Scheer
*CMZ :  2.20/05 13/03/2001  13.41.14  by  Michael Scheer
*CMZ :  2.16/07 22/09/2000  10.44.05  by  Michael Scheer
*CMZ :  2.16/04 20/06/2000  17.04.46  by  Michael Scheer
*CMZ :  2.15/00 19/05/2000  11.01.39  by  Michael Scheer
*CMZ :  2.14/02 25/04/2000  17.04.16  by  Michael Scheer
*CMZ :  2.14/01 19/04/2000  13.47.47  by  Michael Scheer
*CMZ :  2.14/00 18/04/2000  18.12.12  by  Michael Scheer
*CMZ :  2.13/11 19/03/2000  12.48.34  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.26.25  by  Michael Scheer
*CMZ :  1.03/06 05/08/98  17.54.46  by  Michael Scheer
*CMZ :  1.02/03 14/01/98  10.03.25  by  Michael Scheer
*CMZ :  1.02/00 18/12/97  11.55.11  by  Michael Scheer
*CMZ :  1.01/00 28/10/97  18.46.15  by  Michael Scheer
*CMZ :  1.00/00 10/07/97  17.39.21  by  Michael Scheer
*CMZ : 00.01/10 02/06/96  12.33.40  by  Michael Scheer
*CMZ : 00.01/08 29/06/95  10.57.36  by  Michael Scheer
*CMZ : 00.01/07 22/03/95  15.54.56  by  Michael Scheer
*CMZ : 00.00/01 03/03/95  15.56.37  by  Michael Scheer
*CMZ : 00.00/00 03/03/95  15.50.11  by  Johannes Bahrdt
*-- Author : Michael Scheer
      SUBROUTINE REC_INIT
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

C*****************************************************************
c       Berechnung des Magnetfeldes eines Arrays von
c       mehreren Permanentmagnetkloetzen
c*****************************************************************

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,modulator.
      include 'modulator.cmn'
*KEND.

      INTEGER IANZX,IANZZ,I,IMAG,IIMOD,II,ICAL,IUREC

C        DOUBLE PRECISION bxfield(501,41),byfield(501,41),bzfield(501,41)
      DOUBLE PRECISION XMUE0,PI,GRARAD,XSHIFT,YSHIFT,ZSHIFT
      DOUBLE PRECISION XMIN,XMAX,ZMIN,ZMAX,DELTAX,DELTAZ
      DOUBLE PRECISION OFFX,OFFX_I,OFFX_E

      INTEGER IMAGUP_I,IMAGLO_I,IMAGNET_I,IMOD_I
      DOUBLE PRECISION XLEN_I,YLEN_I,ZLEN_I,SHIFT_IU,SHIFT_ID
      DOUBLE PRECISION DUP_I,DLO_I
      DOUBLE PRECISION DX0_I(NKLOTZ),DY0_I(NKLOTZ),DZ0_I(NKLOTZ)
      DOUBLE PRECISION THETA0_I(NKLOTZ),PHI0_I(NKLOTZ),BC0_I(NKLOTZ)
      DOUBLE PRECISION DDX_I(NKLOTZ),DDY_I(NKLOTZ),DDZ_I(NKLOTZ)

      INTEGER IMAGUP_E,IMAGLO_E,IMAGNET_E,IMOD_E
      DOUBLE PRECISION XLEN_E,YLEN_E,ZLEN_E,SHIFT_EU,SHIFT_ED
      DOUBLE PRECISION DUP_E,DLO_E
      DOUBLE PRECISION DX0_E(NKLOTZ),DY0_E(NKLOTZ),DZ0_E(NKLOTZ)
      DOUBLE PRECISION THETA0_E(NKLOTZ),PHI0_E(NKLOTZ),BC0_E(NKLOTZ)
      DOUBLE PRECISION DDX_E(NKLOTZ),DDY_E(NKLOTZ),DDZ_E(NKLOTZ)

      INTEGER IMAGTOTS
      DOUBLE PRECISION XLENS(NKLOTZ),YLENS(NKLOTZ),ZLENS(NKLOTZ)
      DOUBLE PRECISION DXS(NKLOTZ),DYS(NKLOTZ),DZS(NKLOTZ)
      DOUBLE PRECISION THETAS(NKLOTZ),PHIS(NKLOTZ),BCS(NKLOTZ)
      DOUBLE PRECISION FTX(NKLOTZ),FTY(NKLOTZ),FT2(NKLOTZ),XSCL,TSCL,DYFT
      DOUBLE PRECISION F1(NKLOTZ),F2(NKLOTZ),F3(NKLOTZ),F4(NKLOTZ)

      DOUBLE PRECISION XLEN_A,YLEN_A,ZLEN_A,COSPHI,SINPHI,SINTHE,COSTHE

      REAL xran(1),xranO,VBX,VBY,VBZ,B0,rr(2)
      INTEGER NTOTIN,NTOT2IN,NFUNTAP

      CHARACTER(2) CSTAR
      CHARACTER(80) COMTAP

      DATA ICAL/0/
      DATA PI/3.141592653589793D0/

      IF (ICAL.NE.0) RETURN

      irecsolve=0

      ICAL=1

      IF (BCSTART.EQ.9999.) BCSTART=-1.D10
      IF (BCEND.EQ.9999.) BCEND=1.D10
      IF (BCRANSIG.LE.0.0) BCRANSIG=1.E10

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      SR REC_INIT: REC-structure read from file'
      WRITE(LUNGFO,*)'      ',FILEREC
      WRITE(LUNGFO,*)'      WINREC,RANGREC:'
      WRITE(LUNGFO,*)'      ',WINREC,RANGREC
      WRITE(LUNGFO,*)'      BCRAN, BCRANSIG:'
      WRITE(LUNGFO,*)'      ',BCRAN,BCRANSIG
      WRITE(LUNGFO,*)'      K90270:'
      WRITE(LUNGFO,*)'      ',K90270
      WRITE(LUNGFO,*)'      IRECSEED, NURANMOD:'
      WRITE(LUNGFO,*)'      ',IRECSEED,NURANMOD
      WRITE(LUNGFO,*)'      BCSTART, BCEND:'
      WRITE(LUNGFO,*)'      ',BCSTART, BCEND
      WRITE(LUNGFO,*)'      SCALKL,SCALADD:'
      WRITE(LUNGFO,*)'      ',SCALKL,SCALADD
      WRITE(LUNGFO,*)'      RECGAP:',RECGAP
      WRITE(LUNGFO,*)'      USHIFT,DSHIFT:'
      WRITE(LUNGFO,*)'      ',USHIFT,DSHIFT
      WRITE(LUNGFO,*)

      IF (IRECSEED.NE.0) CALL RMARIN(IRECSEED,NTOTIN,NTOT2IN) !CERN V113

      IF (IRECU.NE.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      Parameter for additional undulators:'
        WRITE(LUNGFO,*)

        DO IUREC=1,IRECU

          IF (ubansig(iurec).LE.0.0) ubansig(iurec)=1.0E10

          IF (USIGOFFY(IRECU).EQ.0.D0) THEN
            USIGOFFY(IRECU)=1.D30
          ENDIF

          if (urecmupar(irecu).eq.0.0d0) urecmupar(irecu)=1.0d0
          if (urecmuper(irecu).eq.0.0d0) urecmuper(irecu)=1.0d0

c          if (urecmupar(irecu).ne.1.0d0) irecsolve=1
c          if (urecmuper(irecu).ne.1.0d0) irecsolve=1

          WRITE(LUNGFO,*)'      KRECPER:',KRECPER(IRECU)
          WRITE(LUNGFO,*)'      URECLX:',URECLX(IRECU)
          WRITE(LUNGFO,*)'      URECLY:',URECLY(IRECU)
          WRITE(LUNGFO,*)'      URECLZ:',URECLZ(IRECU)
          WRITE(LUNGFO,*)'      URECGAP:',URECGAP(IRECU)
          WRITE(LUNGFO,*)'      UTAPER:',UTAPER(IRECU)
          WRITE(LUNGFO,*)'      URECCX:',URECCX(IRECU)
          WRITE(LUNGFO,*)'      URECCZ:',URECCZ(IRECU)
          WRITE(LUNGFO,*)'      URECBC:',URECBC(IRECU)
c         WRITE(LUNGFO,*)'      URECMUPAR:',URECmupar(IRECU)
c          WRITE(LUNGFO,*)'      URECMUPER:',URECmuper(IRECU)
          WRITE(LUNGFO,*)'      UBANGERR:',UBANGERR(IRECU)
          WRITE(LUNGFO,*)'      UBANSIG:',UBANSIG(IRECU)
          WRITE(LUNGFO,*)'      USIGOFFY:',USIGOFFY(IRECU)
          WRITE(LUNGFO,*)'      IUHELI:',IUHELI(IRECU)
          WRITE(LUNGFO,*)'      URSPLIT:',URSPLIT(IRECU)
          WRITE(LUNGFO,*)'      URSHIFT:',URSHIFT(IRECU)
          WRITE(LUNGFO,*)'      URSHADD:',URSHADD(IRECU)
        ENDDO

      ENDIF

      xmue0=1.2566d-6
      grarad=pi/180.d0

      IF (IKRESTOR.LE.0) THEN

        IF (SCALADD.NE.0.D0.OR.SCALKL.NE.0.D0) THEN

          open(unit=LUNREC,file=FILEREC,status='old')

          call util_skip_comment(lunrec)
          read(LUNREC,*)xmin,xmax,ianzx
          call util_skip_comment(lunrec)
          read(LUNREC,*)zmin,zmax,ianzz
          if(ianzx.eq.1)then
            deltax=0.
          else
            deltax=(xmax-xmin)/FLOAT(ianzx-1)
          endif
          if(ianzz.eq.1)then
            deltaz=0.
          else
            deltaz=(zmax-zmin)/FLOAT(ianzz-1)
          endif
          call util_skip_comment(lunrec)
          read(LUNREC,*)xlen(1),ylen(1),zlen(1)
          call util_skip_comment(lunrec)
          read(LUNREC,*)shift_U
          IF (SHIFT_U.EQ.9999.) THEN
            SHIFT_U=USHIFT
          ELSE IF (SHIFT_U.EQ.-9999.) THEN
            SHIFT_U=-USHIFT
          ENDIF
          call util_skip_comment(lunrec)
          read(LUNREC,*)imagup
          IF(IMAGUP.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(LUNREC,*)dup
          IF (DUP.EQ.9999.) DUP=RECGAP
          IF (DUP.EQ.-9999.) DUP=-RECGAP
          do i=1,imagup
            call util_skip_comment(lunrec)
            read(LUNREC,*)theta0(i),phi0(i),bc0(i)
            BC0(I)=BC0(I)*SCALKL
            theta0(i)=theta0(i)*grarad
            phi0(i)=phi0(i)*grarad
            call util_skip_comment(lunrec)
            read(LUNREC,*)dx0(i),dy0(i),dz0(i)
            dy0(i)=dy0(i)+(dup+0.5*ylen(1))/ylen(1)
          enddo

          call util_skip_comment(lunrec)
          read(LUNREC,*)shift_D
          IF (SHIFT_D.EQ.9999.) THEN
            SHIFT_D=DSHIFT
          ELSE IF (SHIFT_D.EQ.-9999.) THEN
            SHIFT_D=-DSHIFT
          ELSE IF (SHIFT_D.EQ.8888.) THEN
            SHIFT_D=-2.D0*SHIFT_U
          ENDIF
          call util_skip_comment(lunrec)
          read(LUNREC,*)imaglo
          IF(IMAGUP+IMAGLO.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(LUNREC,*)dlo
          IF (DLO.EQ.9999.) DLO=RECGAP
          IF (DLO.EQ.-9999.) DLO=-RECGAP
          do i=1,imaglo
            call util_skip_comment(lunrec)
            read(LUNREC,*)theta0(i+imagup),phi0(i+imagup),bc0(i+imagup)
            BC0(I+IMAGUP)=BC0(I+IMAGUP)*SCALKL
            theta0(i+imagup)=theta0(i+imagup)*grarad
            phi0(i+imagup)=phi0(i+imagup)*grarad
            call util_skip_comment(lunrec)
            read(LUNREC,*)dx0(i+imagup),dy0(i+imagup),dz0(i+imagup)
            dy0(i+imagup)=dy0(i+imagup)+(dlo-0.5*ylen(1))/ylen(1)
          enddo

          do i=1,imagup/2
            dx0(i)=dx0(i)+shift_U/XLEN(1)
          enddo

          do i=imagup+imaglo/2+1,imagup+imaglo
            dx0(i)=dx0(i)+shift_U/XLEN(1)+SHIFT_D/XLEN(1)
          enddo

          imagnet=imagup+imaglo

          call util_skip_comment(lunrec)
          read(LUNREC,*)imod,OFFX
          IF(IMAGNET*IMOD.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          do i=1,imod
            call util_skip_comment(lunrec)
            read(LUNREC,*)ddx(i),ddy(i),ddz(i)
            DDX(I)=DDX(I)+OFFX/XLEN(1)
          enddo

c************************************************************

c End poles

          IMAGTOT=IMOD*IMAGNET

          call util_skip_comment(lunrec)
          READ(LUNREC,*)XLEN_I,YLEN_I,ZLEN_I
          call util_skip_comment(lunrec)
          READ(LUNREC,*)SHIFT_IU
          IF (SHIFT_IU.EQ.9999.) THEN
            SHIFT_IU=USHIFT
          ELSE IF (SHIFT_IU.EQ.-9999.) THEN
            SHIFT_IU=-USHIFT
          ENDIF

          call util_skip_comment(lunrec)
          READ(LUNREC,*)IMAGUP_I
          IF(IMAGTOT+IMAGUP_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          READ(LUNREC,*)DUP_I
          IF (DUP_I.EQ.9999.) DUP_I=RECGAP
          IF (DUP_I.EQ.-9999.) DUP_I=-RECGAP
          DO I=1,IMAGUP_I
            call util_skip_comment(lunrec)
            READ(LUNREC,*)THETA0_I(I),PHI0_I(I),BC0_I(I)
            BC0_I(I)=BC0_I(I)*SCALKL
            THETA0_I(I)=THETA0_I(I)*GRARAD
            PHI0_I(I)=PHI0_I(I)*GRARAD
            call util_skip_comment(lunrec)
            READ(LUNREC,*)DX0_I(I),DY0_I(I),DZ0_I(I)
            DY0_I(I)=DY0_I(I)+(DUP_I+0.5D0*YLEN_I)/YLEN_I
          ENDDO   !IMAGUP_I

          call util_skip_comment(lunrec)
          READ(LUNREC,*)SHIFT_ID
          IF (SHIFT_ID.EQ.9999.) THEN
            SHIFT_ID=DSHIFT
          ELSE IF (SHIFT_ID.EQ.-9999.) THEN
            SHIFT_ID=-DSHIFT
          ELSE IF (SHIFT_ID.EQ.8888.) THEN
            SHIFT_ID=-2.D0*SHIFT_U
          ENDIF
          call util_skip_comment(lunrec)
          READ(LUNREC,*)IMAGLO_I
          IF(IMAGTOT+IMAGUP_I+IMAGLO_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          READ(LUNREC,*)DLO_I
          IF (DLO_I.EQ.9999.) DLO_I=RECGAP
          IF (DLO_I.EQ.-9999.) DLO_I=-RECGAP
          DO I=1,IMAGLO_I
            call util_skip_comment(lunrec)
            READ(LUNREC,*)THETA0_I(I+IMAGUP_I),PHI0_I(I+IMAGUP_I),BC0_I(I+IMAGUP_I)
            BC0_I(I+IMAGUP_I)=BC0_I(I+IMAGUP_I)*SCALKL
            THETA0_I(I+IMAGUP_I)=THETA0_I(I+IMAGUP_I)*GRARAD
            PHI0_I(I+IMAGUP_I)=PHI0_I(I+IMAGUP_I)*GRARAD
            call util_skip_comment(lunrec)
            READ(LUNREC,*)DX0_I(I+IMAGUP_I),DY0_I(I+IMAGUP_I),DZ0_I(I+IMAGUP_I)
            DY0_I(I+IMAGUP_I)=DY0_I(I+IMAGUP_I)+(DLO_I-0.5D0*YLEN_I)/YLEN_I
          ENDDO   !IMAGLO_I

          DO I=1,IMAGUP_I/2
            DX0_I(I)=DX0_I(I)+SHIFT_IU/XLEN_I
          ENDDO

          DO I=IMAGUP_I+IMAGLO_I/2+1,IMAGUP_I+IMAGLO_I
            DX0_I(I)=DX0_I(I)+SHIFT_IU/XLEN_I+SHIFT_ID/XLEN_I
          ENDDO

          IMAGNET_I=IMAGUP_I+IMAGLO_I

          call util_skip_comment(lunrec)
          READ(LUNREC,*)IMOD_I,OFFX_I
          IF(IMAGTOT+IMAGNET_I*IMOD_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          DO I=1,IMOD_I
            call util_skip_comment(lunrec)
            READ(LUNREC,*)DDX_I(I),DDY_I(I),DDZ_I(I)
            DDX_I(I)=DDX_I(I)+OFFX_I/XLEN_I
          ENDDO

          call util_skip_comment(lunrec)
          READ(LUNREC,*)XLEN_E,YLEN_E,ZLEN_E
          call util_skip_comment(lunrec)
          READ(LUNREC,*)SHIFT_EU
          IF (SHIFT_EU.EQ.9999.) THEN
            SHIFT_EU=USHIFT
          ELSE IF (SHIFT_EU.EQ.-9999.) THEN
            SHIFT_EU=-USHIFT
          ENDIF

          call util_skip_comment(lunrec)
          READ(LUNREC,*)IMAGUP_E
          IF(IMAGTOT+IMAGNET_I*IMOD_I+IMAGUP_E.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          READ(LUNREC,*)DUP_E
          IF (DUP_E.EQ.9999.) DUP_E=RECGAP
          IF (DUP_E.EQ.-9999.) DUP_E=-RECGAP
          DO I=1,IMAGUP_E
            call util_skip_comment(lunrec)
            READ(LUNREC,*)THETA0_E(I),PHI0_E(I),BC0_E(I)
            BC0_E(I)=BC0_E(I)*SCALKL
            THETA0_E(I)=THETA0_E(I)*GRARAD
            PHI0_E(I)=PHI0_E(I)*GRARAD
            call util_skip_comment(lunrec)
            READ(LUNREC,*)DX0_E(I),DY0_E(I),DZ0_E(I)
            DY0_E(I)=DY0_E(I)+(DUP_E+0.5D0*YLEN_E)/YLEN_E
          ENDDO   !IMAGUP_E

          call util_skip_comment(lunrec)
          READ(LUNREC,*)SHIFT_ED
          IF (SHIFT_ED.EQ.9999.) THEN
            SHIFT_ED=DSHIFT
          ELSE IF (SHIFT_ED.EQ.-9999.) THEN
            SHIFT_ED=-DSHIFT
          ELSE IF (SHIFT_ED.EQ.8888.) THEN
            SHIFT_ED=-2.D0*SHIFT_U
          ENDIF
          call util_skip_comment(lunrec)
          READ(LUNREC,*)IMAGLO_E
          IF(IMAGTOT+IMAGNET_I*IMOD_I+IMAGUP_E+IMAGLO_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          READ(LUNREC,*)DLO_E
          IF (DLO_E.EQ.9999.) DLO_E=RECGAP
          IF (DLO_E.EQ.-9999.) DLO_E=-RECGAP
          DO I=1,IMAGLO_E
            call util_skip_comment(lunrec)
            READ(LUNREC,*)THETA0_E(I+IMAGUP_E),PHI0_E(I+IMAGUP_E),BC0_E(I+IMAGUP_E)
            BC0_E(I+IMAGUP_E)=BC0_E(I+IMAGUP_E)*SCALKL
            THETA0_E(I+IMAGUP_E)=THETA0_E(I+IMAGUP_E)*GRARAD
            PHI0_E(I+IMAGUP_E)=PHI0_E(I+IMAGUP_E)*GRARAD
            call util_skip_comment(lunrec)
            READ(LUNREC,*)DX0_E(I+IMAGUP_E),DY0_E(I+IMAGUP_E),DZ0_E(I+IMAGUP_E)
            DY0_E(I+IMAGUP_E)=DY0_E(I+IMAGUP_E)+(DLO_E-0.5D0*YLEN_E)/YLEN_E
          ENDDO   !IMAGLO_E

          DO I=1,IMAGUP_E/2
            DX0_E(I)=DX0_E(I)+SHIFT_EU/XLEN_E
          ENDDO

          DO I=IMAGUP_E+IMAGLO_E/2+1,IMAGUP_E+IMAGLO_E
            DX0_E(I)=DX0_E(I)+SHIFT_EU/XLEN_E+SHIFT_ED/XLEN_E
          ENDDO

          IMAGNET_E=IMAGUP_E+IMAGLO_E

          call util_skip_comment(lunrec)
          READ(LUNREC,*)IMOD_E,OFFX_E
          IF(IMAGTOT+IMAGNET_I*IMOD_I+IMAGNET_E*IMOD_E.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          DO I=1,IMOD_E
            call util_skip_comment(lunrec)
            READ(LUNREC,*)DDX_E(I),DDY_E(I),DDZ_E(I)
            DDX_E(I)=DDX_E(I)+OFFX_E/XLEN_E
          ENDDO

c************ Loop over modules ******************************

          do iimod=1,imod

            xshift=ddx(iimod)
            yshift=ddy(iimod)
            zshift=ddz(iimod)

c************ Loop over magnets ******************************

            do imag=1,imagnet

              XLEN((IIMOD-1)*IMAGNET+IMAG)=XLEN(1)
              YLEN((IIMOD-1)*IMAGNET+IMAG)=YLEN(1)
              ZLEN((IIMOD-1)*IMAGNET+IMAG)=ZLEN(1)

              bc((IIMOD-1)*IMAGNET+IMAG)=bc0(imag)
              dx((IIMOD-1)*IMAGNET+IMAG)=dx0(imag)*xlen(1)+xshift*xlen(1)
              dy((IIMOD-1)*IMAGNET+IMAG)=dy0(imag)*ylen(1)+yshift*ylen(1)
              dz((IIMOD-1)*IMAGNET+IMAG)=dz0(imag)*zlen(1)+zshift*zlen(1)
              theta((IIMOD-1)*IMAGNET+IMAG)=theta0(imag)
              phi((IIMOD-1)*IMAGNET+IMAG)=phi0(imag)


            ENDDO   !MAGNETS

          ENDDO   !MODULES


c************************************************************

c End poles

          do iimod=1,imod_I

            xshift=ddx_I(iimod)
            yshift=ddy_I(iimod)
            zshift=ddz_I(iimod)

            do imag=1,imagnet_I

              II=IMAGTOT+(IIMOD-1)*IMAGNET_I+IMAG

              XLEN(II)=XLEN_I
              YLEN(II)=YLEN_I
              ZLEN(II)=ZLEN_I

              bc(II)=bc0_I(imag)
              dx(II)=dx0_I(imag)*xlen_I+xshift*xlen_I
              dy(II)=dy0_I(imag)*ylen_I+yshift*ylen_I
              dz(II)=dz0_I(imag)*zlen_I+zshift*zlen_I
              theta(II)=theta0_I(imag)
              phi(II)=phi0_I(imag)

            ENDDO   !MAGNETS

          ENDDO   !MODULES

          IMAGTOT=IMAGTOT+IMAGNET_I*IMOD_I

          do iimod=1,imod_E

            xshift=ddx_E(iimod)
            yshift=ddy_E(iimod)
            zshift=ddz_E(iimod)

            do imag=1,imagnet_E

              II=IMAGTOT+(IIMOD-1)*IMAGNET_E+IMAG

              XLEN(II)=XLEN_E
              YLEN(II)=YLEN_E
              ZLEN(II)=ZLEN_E

              bc(II)=bc0_E(imag)
              dx(II)=dx0_E(imag)*xlen_E+xshift*xlen_E
              dy(II)=dy0_E(imag)*ylen_E+yshift*ylen_E
              dz(II)=dz0_E(imag)*zlen_E+zshift*zlen_E
              theta(II)=theta0_E(imag)
              phi(II)=phi0_E(imag)

            ENDDO   !MAGNETS

          ENDDO   !MODULES

          IMAGTOT=IMAGTOT+IMAGNET_E*IMOD_E

C COPY ARRAY IN ORDER TO SAVE THEM{

          DO I=1,IMAGTOT
            BCS(I)=BC(I)
            XLENS(I)=XLEN(I)
            YLENS(I)=YLEN(I)
            ZLENS(I)=ZLEN(I)
            THETAS(I)=THETA(I)
            PHIS(I)=PHI(I)
            DXS(I)=DX(I)
            DYS(I)=DY(I)
            DZS(I)=DZ(I)
          ENDDO

          IMAGTOTS=IMAGTOT
          IMAGTOT=0

C COPY ARRAY IN ORDER TO SAVE THEM}

C****************************************************
C { WHOLE STORY AGAIN FOR SECOND DEVICE

          call util_skip_comment(lunrec)
          read(lunrec,*)xlen(1),ylen(1),zlen(1)
          call util_skip_comment(lunrec)
          read(lunrec,*)shift_U
          IF (SHIFT_U.EQ.9999.) THEN
            SHIFT_U=USHIFT
          ELSE IF (SHIFT_U.EQ.-9999.) THEN
            SHIFT_U=-USHIFT
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)imagup
          IF(IMAGUP.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)dup
          IF (DUP.EQ.9999.) DUP=RECGAP
          IF (DUP.EQ.-9999.) DUP=-RECGAP
          do i=1,imagup
            call util_skip_comment(lunrec)
            read(lunrec,*)theta0(i),phi0(i),bc0(i)
            BC0(I)=BC0(I)*SCALKL
            theta0(i)=theta0(i)*grarad
            phi0(i)=phi0(i)*grarad
            call util_skip_comment(lunrec)
            read(lunrec,*)dx0(i),dy0(i),dz0(i)
C03JUN97    DUMZ=DABS(DZ0(I))-1.D0
C03JUN97    DZ0(I)=DSIGN (ZLENS(1)/ZLEN(1)+DUMZ,DZ0(I))
            dy0(i)=dy0(i)+(dup+0.5*ylen(1))/ylen(1)
          enddo

          call util_skip_comment(lunrec)
          read(lunrec,*)shift_D
          IF (SHIFT_D.EQ.9999.) THEN
            SHIFT_D=DSHIFT
          ELSE IF (SHIFT_D.EQ.-9999.) THEN
            SHIFT_D=-DSHIFT
          ELSE IF (SHIFT_D.EQ.8888.) THEN
            SHIFT_D=-2.D0*SHIFT_U
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)imaglo
          IF(IMAGUP+IMAGLO.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)dlo
          IF (DLO.EQ.9999.) DLO=RECGAP
          IF (DLO.EQ.-9999.) DLO=-RECGAP
          do i=1,imaglo
            call util_skip_comment(lunrec)
            read(lunrec,*)theta0(i+imagup),phi0(i+imagup),bc0(i+imagup)
            BC0(I+imagup)=BC0(I+imagup)*SCALKL
            theta0(i+imagup)=theta0(i+imagup)*grarad
            phi0(i+imagup)=phi0(i+imagup)*grarad
            call util_skip_comment(lunrec)
            read(lunrec,*)dx0(i+imagup),dy0(i+imagup),dz0(i+imagup)
C03JUN97    DUMZ=DABS(DZ0(I+imagup))-1.D0
C03JUN97    DZ0(I+imagup)=DSIGN (ZLENS(1)/ZLEN(1)+DUMZ,DZ0(I+imagup))
            dy0(i+imagup)=dy0(i+imagup)+(dlo-0.5*ylen(1))/ylen(1)
          enddo

          do i=1,imagup/2
            dx0(i)=dx0(i)+shift_U/XLEN(1)
          enddo

          do i=imagup+imaglo/2+1,imagup+imaglo
            dx0(i)=dx0(i)+shift_U/XLEN(1)+SHIFT_D/XLEN(1)
          enddo

          imagnet=imagup+imaglo

          call util_skip_comment(lunrec)
          read(lunrec,*)imod,OFFX
          IF(IMAGNET*IMOD.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          do i=1,imod
            call util_skip_comment(lunrec)
            read(lunrec,*)ddx(i),ddy(i),ddz(i)
            DDX(I)=DDX(I)+OFFX/XLEN(1)
          enddo

c************************************************************

c End poles

          IMAGTOT=IMOD*IMAGNET

          call util_skip_comment(lunrec)
          read(lunrec,*)XLEN_I,YLEN_I,ZLEN_I
          call util_skip_comment(lunrec)
          read(lunrec,*)SHIFT_IU
          IF (SHIFT_IU.EQ.9999.) THEN
            SHIFT_IU=USHIFT
          ELSE IF (SHIFT_IU.EQ.-9999.) THEN
            SHIFT_IU=-USHIFT
          ENDIF

          call util_skip_comment(lunrec)
          read(lunrec,*)IMAGUP_I
          IF(IMAGTOT+IMAGUP_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)DUP_I
          IF (DUP_I.EQ.9999.) DUP_I=RECGAP
          IF (DUP_I.EQ.-9999.) DUP_I=-RECGAP
          DO I=1,IMAGUP_I
            call util_skip_comment(lunrec)
            read(lunrec,*)THETA0_I(I),PHI0_I(I),BC0_I(I)
            BC0_I(I)=BC0_I(I)*SCALKL
            THETA0_I(I)=THETA0_I(I)*GRARAD
            PHI0_I(I)=PHI0_I(I)*GRARAD
            call util_skip_comment(lunrec)
            read(lunrec,*)DX0_I(I),DY0_I(I),DZ0_I(I)
C03JUN97    DUMZ=DABS(DZ0_I(I))-1.D0
C03JUN97    DZ0_I(I)=DSIGN (ZLENS(1)/ZLEN_I+DUMZ,DZ0_I(I))
            DY0_I(I)=DY0_I(I)+(DUP_I+0.5D0*YLEN_I)/YLEN_I
          ENDDO   !IMAGUP_I

          call util_skip_comment(lunrec)
          read(lunrec,*)SHIFT_ID
          IF (SHIFT_ID.EQ.9999.) THEN
            SHIFT_ID=DSHIFT
          ELSE IF (SHIFT_ID.EQ.-9999.) THEN
            SHIFT_ID=-DSHIFT
          ELSE IF (SHIFT_ID.EQ.8888.) THEN
            SHIFT_ID=-2.D0*SHIFT_U
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)IMAGLO_I
          IF(IMAGTOT+IMAGUP_I+IMAGLO_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)DLO_I
          IF (DLO_I.EQ.9999.) DLO_I=RECGAP
          IF (DLO_I.EQ.-9999.) DLO_I=-RECGAP
          DO I=1,IMAGLO_I
            call util_skip_comment(lunrec)
            read(lunrec,*)THETA0_I(I+IMAGUP_I),PHI0_I(I+IMAGUP_I),BC0_I(I+IMAGUP_I)
            BC0_I(I+IMAGUP_I)=BC0_I(I+IMAGUP_I)*SCALKL
            THETA0_I(I+IMAGUP_I)=THETA0_I(I+IMAGUP_I)*GRARAD
            PHI0_I(I+IMAGUP_I)=PHI0_I(I+IMAGUP_I)*GRARAD
            call util_skip_comment(lunrec)
            read(lunrec,*)DX0_I(I+IMAGUP_I),DY0_I(I+IMAGUP_I),DZ0_I(I+IMAGUP_I)
C03JUN97    DUMZ=DABS(DZ0_I(I+IMAGUP_I))-1.D0
C03JUN97    DZ0_I(I+IMAGUP_I)=DSIGN (ZLENS(1)/ZLEN_I+DUMZ,DZ0_I(I+IMAGUP_I))
            DY0_I(I+IMAGUP_I)=DY0_I(I+IMAGUP_I)+(DLO_I-0.5D0*YLEN_I)/YLEN_I
          ENDDO   !IMAGLO_I

          DO I=1,IMAGUP_I/2
            DX0_I(I)=DX0_I(I)+SHIFT_IU/XLEN_I
          ENDDO

          DO I=IMAGUP_I+IMAGLO_I/2+1,IMAGUP_I+IMAGLO_I
            DX0_I(I)=DX0_I(I)+SHIFT_IU/XLEN_I+SHIFT_ID/XLEN_I
          ENDDO

          IMAGNET_I=IMAGUP_I+IMAGLO_I

          call util_skip_comment(lunrec)
          read(lunrec,*)IMOD_I,OFFX_I
          IF(IMAGTOT+IMAGNET_I*IMOD_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          DO I=1,IMOD_I
            call util_skip_comment(lunrec)
            read(lunrec,*)DDX_I(I),DDY_I(I),DDZ_I(I)
            DDX_I(I)=DDX_I(I)+OFFX_I/XLEN_I
          ENDDO

          call util_skip_comment(lunrec)
          read(lunrec,*)XLEN_E,YLEN_E,ZLEN_E
          call util_skip_comment(lunrec)
          read(lunrec,*)SHIFT_EU
          IF (SHIFT_EU.EQ.9999.) THEN
            SHIFT_EU=USHIFT
          ELSE IF (SHIFT_EU.EQ.-9999.) THEN
            SHIFT_EU=-USHIFT
          ENDIF

          call util_skip_comment(lunrec)
          read(lunrec,*)IMAGUP_E
          IF(IMAGTOT+IMAGNET_I*IMOD_I+IMAGUP_E.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)DUP_E
          IF (DUP_E.EQ.9999.) DUP_E=RECGAP
          IF (DUP_E.EQ.-9999.) DUP_E=-RECGAP
          DO I=1,IMAGUP_E
            call util_skip_comment(lunrec)
            read(lunrec,*)THETA0_E(I),PHI0_E(I),BC0_E(I)
            BC0_E(I)=BC0_E(I)*SCALKL
            THETA0_E(I)=THETA0_E(I)*GRARAD
            PHI0_E(I)=PHI0_E(I)*GRARAD
            call util_skip_comment(lunrec)
            read(lunrec,*)DX0_E(I),DY0_E(I),DZ0_E(I)
C03JUN97    DUMZ=DABS(DZ0_E(I))-1.D0
C03JUN97    DZ0_E(I)=DSIGN (ZLENS(1)/ZLEN_E+DUMZ,DZ0_E(I))
            DY0_E(I)=DY0_E(I)+(DUP_E+0.5D0*YLEN_E)/YLEN_E
          ENDDO   !IMAGUP_E

          call util_skip_comment(lunrec)
          read(lunrec,*)SHIFT_ED
          IF (SHIFT_ED.EQ.9999.) THEN
            SHIFT_ED=DSHIFT
          ELSE IF (SHIFT_ED.EQ.-9999.) THEN
            SHIFT_ED=-DSHIFT
          ELSE IF (SHIFT_ED.EQ.8888.) THEN
            SHIFT_ED=-2.D0*SHIFT_U
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)IMAGLO_E
          IF(IMAGTOT+IMAGNET_I*IMOD_I+IMAGUP_E+IMAGLO_I.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          call util_skip_comment(lunrec)
          read(lunrec,*)DLO_E
          IF (DLO_E.EQ.9999.) DLO_E=RECGAP
          IF (DLO_E.EQ.-9999.) DLO_E=-RECGAP
          DO I=1,IMAGLO_E
            call util_skip_comment(lunrec)
            read(lunrec,*)THETA0_E(I+IMAGUP_E),PHI0_E(I+IMAGUP_E),BC0_E(I+IMAGUP_E)
            BC0_E(I+IMAGUP_E)=BC0_E(I+IMAGUP_E)*SCALKL
            THETA0_E(I+IMAGUP_E)=THETA0_E(I+IMAGUP_E)*GRARAD
            PHI0_E(I+IMAGUP_E)=PHI0_E(I+IMAGUP_E)*GRARAD
            call util_skip_comment(lunrec)
            read(lunrec,*)DX0_E(I+IMAGUP_E),DY0_E(I+IMAGUP_E),DZ0_E(I+IMAGUP_E)
C03JUN97    DUMZ=DABS(DZ0_E(I+IMAGUP_E))-1.D0
C03JUN97    DZ0_E(I+IMAGUP_E)=DSIGN (ZLENS(1)/ZLEN_E+DUMZ,DZ0_E(I+IMAGUP_E))
            DY0_E(I+IMAGUP_E)=DY0_E(I+IMAGUP_E)+(DLO_E-0.5D0*YLEN_E)/YLEN_E
          ENDDO   !IMAGLO_E

          DO I=1,IMAGUP_E/2
            DX0_E(I)=DX0_E(I)+SHIFT_EU/XLEN_E
          ENDDO

          DO I=IMAGUP_E+IMAGLO_E/2+1,IMAGUP_E+IMAGLO_E
            DX0_E(I)=DX0_E(I)+SHIFT_EU/XLEN_E+SHIFT_ED/XLEN_E
          ENDDO

          IMAGNET_E=IMAGUP_E+IMAGLO_E

          call util_skip_comment(lunrec)
          read(lunrec,*)IMOD_E,OFFX_E
          IF(IMAGTOT+IMAGNET_I*IMOD_I+IMAGNET_E*IMOD_E.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF
          DO I=1,IMOD_E
            call util_skip_comment(lunrec)
            read(lunrec,*)DDX_E(I),DDY_E(I),DDZ_E(I)
            DDX_E(I)=DDX_E(I)+OFFX_E/XLEN_E
          ENDDO

c************ Loop over modules ******************************

          do iimod=1,imod

            xshift=ddx(iimod)
            yshift=ddy(iimod)
            zshift=ddz(iimod)

c************ Loop over magnets ******************************

            do imag=1,imagnet

              XLEN((IIMOD-1)*IMAGNET+IMAG)=XLEN(1)
              YLEN((IIMOD-1)*IMAGNET+IMAG)=YLEN(1)
              ZLEN((IIMOD-1)*IMAGNET+IMAG)=ZLEN(1)

              bc((IIMOD-1)*IMAGNET+IMAG)=bc0(imag)
              dx((IIMOD-1)*IMAGNET+IMAG)=dx0(imag)*xlen(1)+xshift*xlen(1)
              dy((IIMOD-1)*IMAGNET+IMAG)=dy0(imag)*ylen(1)+yshift*ylen(1)
              dz((IIMOD-1)*IMAGNET+IMAG)=dz0(imag)*zlen(1)+zshift*zlen(1)
              theta((IIMOD-1)*IMAGNET+IMAG)=theta0(imag)
              phi((IIMOD-1)*IMAGNET+IMAG)=phi0(imag)

            ENDDO   !MAGNETS

          ENDDO   !MODULES

c************************************************************

c End poles

          do iimod=1,imod_I

            xshift=ddx_I(iimod)
            yshift=ddy_I(iimod)
            zshift=ddz_I(iimod)

            do imag=1,imagnet_I

              II=IMAGTOT+(IIMOD-1)*IMAGNET_I+IMAG

              XLEN(II)=XLEN_I
              YLEN(II)=YLEN_I
              ZLEN(II)=ZLEN_I

              bc(II)=bc0_I(imag)
              dx(II)=dx0_I(imag)*xlen_I+xshift*xlen_I
              dy(II)=dy0_I(imag)*ylen_I+yshift*ylen_I
              dz(II)=dz0_I(imag)*zlen_I+zshift*zlen_I
              theta(II)=theta0_I(imag)
              phi(II)=phi0_I(imag)

            ENDDO   !MAGNETS

          ENDDO   !MODULES

          IMAGTOT=IMAGTOT+IMAGNET_I*IMOD_I

          do iimod=1,imod_E

            xshift=ddx_E(iimod)
            yshift=ddy_E(iimod)
            zshift=ddz_E(iimod)

            do imag=1,imagnet_E

              II=IMAGTOT+(IIMOD-1)*IMAGNET_E+IMAG

              XLEN(II)=XLEN_E
              YLEN(II)=YLEN_E
              ZLEN(II)=ZLEN_E

              bc(II)=bc0_E(imag)
              dx(II)=dx0_E(imag)*xlen_E+xshift*xlen_E
              dy(II)=dy0_E(imag)*ylen_E+yshift*ylen_E
              dz(II)=dz0_E(imag)*zlen_E+zshift*zlen_E
              theta(II)=theta0_E(imag)
              phi(II)=phi0_E(imag)

            ENDDO   !MAGNETS

          ENDDO   !MODULES

          IMAGTOT=IMAGTOT+IMAGNET_E*IMOD_E

C WHOLE STORY AGAIN }
C****************************************************

C MERGE FIRST AND SECOND DEVICE

          IF(IMAGTOTS+IMAGTOT.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF

          DO IMAG=1,IMAGTOTS
            BC(IMAGTOT+IMAG)=BCS(IMAG)
            XLEN(IMAGTOT+IMAG)=XLENS(IMAG)
            YLEN(IMAGTOT+IMAG)=YLENS(IMAG)
            ZLEN(IMAGTOT+IMAG)=ZLENS(IMAG)
            THETA(IMAGTOT+IMAG)=THETAS(IMAG)
            PHI(IMAGTOT+IMAG)=PHIS(IMAG)
            DX(IMAGTOT+IMAG)=DXS(IMAG)
            DY(IMAGTOT+IMAG)=DYS(IMAG)
            DZ(IMAGTOT+IMAG)=DZS(IMAG)
          ENDDO   !IMAG

          IMAGTOT=IMAGTOTS+IMAGTOT

          IF (SCALKL.EQ.0.D0) IMAGTOT=0

C****************************************************
C ADDITIONAL MAGNETS

          call util_skip_comment(lunrec)
          read(lunrec,*)IMAGADD

          IF(IMAGTOT+IMAGADD.GT.NKLOTZ) THEN
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
            WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(6,*)
            STOP
          ENDIF

          DO IMAG=1,IMAGADD
            call util_skip_comment(lunrec)
            read(lunrec,*)XLEN_A,YLEN_A,ZLEN_A
            XLEN(IMAGTOT+IMAG)=XLEN_A
            YLEN(IMAGTOT+IMAG)=YLEN_A
            ZLEN(IMAGTOT+IMAG)=ZLEN_A
            call util_skip_comment(lunrec)
            read(lunrec,*)THETA(IMAGTOT+IMAG),PHI(IMAGTOT+IMAG),BC(IMAGTOT+IMAG)
            call util_skip_comment(lunrec)
            read(lunrec,*)DX(IMAGTOT+IMAG),DY(IMAGTOT+IMAG),DZ(IMAGTOT+IMAG)
            BC(IMAGTOT+IMAG)=BC(IMAGTOT+IMAG)*SCALADD
          ENDDO

          DO IMAG=IMAGTOT+1,IMAGTOT+IMAGADD
C03JUN97  DX(IMAG)=XLEN_A*DX(IMAG)
C03JUN97  DY(IMAG)=YLEN_A*DY(IMAG)
C03JUN97  DZ(IMAG)=ZLEN_A*DZ(IMAG)
            THETA(IMAG)=THETA(IMAG)*GRARAD
            PHI(IMAG)=PHI(IMAG)*GRARAD
          ENDDO

          IF (SCALADD.EQ.0.D0) IMAGADD=0
          IMAGTOT=IMAGTOT+IMAGADD

C****************************************************

          close(LUNREC)
        ENDIF !(SCALADD.NE.0.D0.OR.SCALKL.NE.0.D0)

C--- GENERATE FIELD ERRORS{

        IF (BCRAN.NE.0.D0) THEN

          DO I=1,IMAGTOT
            IF (DX(I).GE.BCSTART.AND.DX(I).LE.BCEND) THEN
              IF (K90270.NE.0) THEN
                IF (ABS(ABS(THETA(I)/GRARAD)-90.D0).LT.1.
     &              .OR.ABS(ABS(THETA(I)/GRARAD)-270.D0).LT.1.)  THEN
                  IF(BCRAN.GT.0.D0) THEN
                    CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
1                   IF (ABS(xran(1)).GT.BCRANSIG) THEN
                      CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                      GOTO 1
                    ENDIF
                    BC(I)=BC(I)*(1.+BCRAN*xran(1))
                  ELSE IF(BCRAN.LT.0.D0) THEN
                    BC(I)=BC(I)*(1.-BCRAN)
                  ENDIF
                ENDIF
              ELSE
                IF(BCRAN.GT.0.D0) THEN
                  CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
11                IF (ABS(xran(1)).GT.BCRANSIG) THEN
                    CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                    GOTO 11
                  ENDIF
                  BC(I)=BC(I)*(1.+BCRAN*xran(1))
                ELSE IF(BCRAN.LT.0.D0) THEN
                  BC(I)=BC(I)*(1.-BCRAN)
                ENDIF
              ENDIF
            ENDIF
          ENDDO

        ENDIF  !BCRAN.NE.0.D0

C--- GENERATE FIELD ERRORS}

C--- SIMPLE PLANAR UNDULATOR{

        IF (IRECU.GT.IURECP) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
          WRITE(6,*)'INCREASE PARAMETER IURECP IN FILE KLOTZ.CMN'
          WRITE(6,*)
          STOP
        ENDIF

        IF (IRECU.NE.0) THEN

          DO IUREC=1,IRECU

            IF (IUHELI(IUREC).EQ.0) THEN

              CALL REC_PLAN(IUREC,GRARAD)

            ELSE    !UHELI=0

              IF(IMAGTOT.GT.NKLOTZ) THEN
                WRITE(6,*)
                WRITE(6,*)'*** ERROR IN REC_INIT: DIMENSION EXCEEDED ***'
                WRITE(6,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
                WRITE(6,*)
                STOP
              ENDIF

              IF (IUHELI(IUREC).GT.0) THEN

                IMAGTOTS=IMAGTOT
                CALL REC_HELI(IUREC,GRARAD,1,0,0)

                IF (BCRAN.NE.0.D0) THEN
                  IF (NURANMOD.EQ.0) THEN
                    xranO=xran(1)
                    DO I=IMAGTOTS+1,IMAGTOT
                      IF (DX(I).GE.BCSTART.AND.DX(I).LE.BCEND) THEN
                        IF (K90270.NE.0) THEN
                          IF (ABS(ABS(THETA(I)/GRARAD)-90.D0).LT.1.
     &                        .OR.ABS(ABS(THETA(I)/GRARAD)-270.D0).LT.1.)  THEN
                            IF(BCRAN.GT.0.D0) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
111                           IF (ABS(xran(1)).GT.BCRANSIG) THEN
                                CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                                GOTO 111
                              ENDIF
                              BC(I)=BC(I)*(1.+BCRAN*xran(1))
                            ELSE IF(BCRAN.LT.0.D0) THEN
                              BC(I)=BC(I)*(1.-BCRAN)
                            ENDIF
                          ENDIF
                        ELSE
                          IF(BCRAN.GT.0.D0) THEN
                            CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
1111                        IF (ABS(xran(1)).GT.BCRANSIG) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                              GOTO 1111
                            ENDIF
                            BC(I)=BC(I)*(1.+BCRAN*xran(1))
                          ELSE IF(BCRAN.LT.0.D0) THEN
                            BC(I)=BC(I)*(1.-BCRAN)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF   !NURANMOD
                ENDIF   !BCRAN

                IMAGTOTS=IMAGTOT
                CALL REC_HELI(IUREC,GRARAD,-1,-1,-1)

                IF (BCRAN.NE.0.D0) THEN
                  IF (NURANMOD.EQ.0) THEN
                    xran(1)=xranO
                    DO I=IMAGTOTS+1,IMAGTOT
                      IF (DX(I).GE.BCSTART.AND.DX(I).LE.BCEND) THEN
                        IF (K90270.NE.0) THEN
                          IF (ABS(ABS(THETA(I)/GRARAD)-90.D0).LT.1.
     &                        .OR.ABS(ABS(THETA(I)/GRARAD)-270.D0).LT.1.)  THEN
                            IF(BCRAN.GT.0.D0) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
2                             IF (ABS(xran(1)).GT.BCRANSIG) THEN
                                CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                                GOTO 2
                              ENDIF
                              BC(I)=BC(I)*(1.+BCRAN*xran(1))
                            ELSE IF(BCRAN.LT.0.D0) THEN
                              BC(I)=BC(I)*(1.-BCRAN)
                            ENDIF
                          ENDIF
                        ELSE
                          IF(BCRAN.GT.0.D0) THEN
                            CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
22                          IF (ABS(xran(1)).GT.BCRANSIG) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                              GOTO 22
                            ENDIF
                            BC(I)=BC(I)*(1.+BCRAN*xran(1))
                          ELSE IF(BCRAN.LT.0.D0) THEN
                            BC(I)=BC(I)*(1.-BCRAN)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF   !NURANMOD
                ENDIF   !BCRAN

              ELSE    !UHELI.GT.0

                IMAGTOTS=IMAGTOT

                CALL REC_HELI(IUREC,GRARAD,-1,0,0)

                IF (BCRAN.NE.0.D0) THEN
                  IF (NURANMOD.EQ.0) THEN
                    xranO=xran(1)
                    DO I=IMAGTOTS+1,IMAGTOT
                      IF (DX(I).GE.BCSTART.AND.DX(I).LE.BCEND) THEN
                        IF (K90270.NE.0) THEN
                          IF (ABS(ABS(THETA(I)/GRARAD)-90.D0).LT.1.
     &                        .OR.ABS(ABS(THETA(I)/GRARAD)-270.D0).LT.1.)  THEN
                            IF(BCRAN.GT.0.D0) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
222                           IF (ABS(xran(1)).GT.BCRANSIG) THEN
                                CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                                GOTO 222
                              ENDIF
                              BC(I)=BC(I)*(1.+BCRAN*xran(1))
                            ELSE IF(BCRAN.LT.0.D0) THEN
                              BC(I)=BC(I)*(1.-BCRAN)
                            ENDIF
                          ENDIF
                        ELSE
                          IF(BCRAN.GT.0.D0) THEN
                            CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
21                          IF (ABS(xran(1)).GT.BCRANSIG) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                              GOTO 21
                            ENDIF
                            BC(I)=BC(I)*(1.+BCRAN*xran(1))
                          ELSE IF(BCRAN.LT.0.D0) THEN
                            BC(I)=BC(I)*(1.-BCRAN)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF   !NURANMOD
                ENDIF   !BCRAN

                CALL REC_HELI(IUREC,GRARAD,1,-1,1)

                IF (BCRAN.NE.0.D0) THEN
                  IF (NURANMOD.EQ.0) THEN
                    xran(1)=xranO
                    DO I=IMAGTOTS+1,IMAGTOT
                      IF (DX(I).GE.BCSTART.AND.DX(I).LE.BCEND) THEN
                        IF (K90270.NE.0) THEN
                          IF (ABS(ABS(THETA(I)/GRARAD)-90.D0).LT.1.
     &                        .OR.ABS(ABS(THETA(I)/GRARAD)-270.D0).LT.1.)  THEN
                            IF(BCRAN.GT.0.D0) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
211                           IF (ABS(xran(1)).GT.BCRANSIG) THEN
                                CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                                GOTO 211
                              ENDIF
                              BC(I)=BC(I)*(1.+BCRAN*xran(1))
                            ELSE IF(BCRAN.LT.0.D0) THEN
                              BC(I)=BC(I)*(1.-BCRAN)
                            ENDIF
                          ENDIF
                        ELSE
                          IF(BCRAN.GT.0.D0) THEN
                            CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
221                         IF (ABS(xran(1)).GT.BCRANSIG) THEN
                              CALL RNORML(xran,1,rr) !RANDOM NOISE OFF FIELD
                              GOTO 221
                            ENDIF
                            BC(I)=BC(I)*(1.+BCRAN*xran(1))
                          ELSE IF(BCRAN.LT.0.D0) THEN
                            BC(I)=BC(I)*(1.-BCRAN)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF   !NURANMOD
                ENDIF   !(BCRAN.NE.0.D0)

              ENDIF

            ENDIF   !UHELI

          ENDDO   !IRECU

        ENDIF  !(IRECU.NE.0)

C--- SIMPLE PLANAR UNDULATOR}

        IF (IRECMODU.NE.0) THEN


cerror? 2jun09          IF (ITHEMSYM.EQ.0) ITHEMSYM=1
cerror? 2jun09          IF (ITHESYML.EQ.0) ITHESYML=1
cerror? 2jun09          IF (ITHESYMD.EQ.0) ITHESYMD=1

          CALL MODULATOR(GRARAD)

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     SR MODULATOR:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     NMAGMOD, NSLICE: ',NMAGMOD,NSLICE
          WRITE(LUNGFO,*)'     SCALMOD, SCALRAD, SCALTHE:'
     &      ,SNGL(SCALMOD),SNGL(SCALRAD),SNGL(SCALTHE)
          WRITE(LUNGFO,*)'     ITHEMSYM,ITHESYML,ITHESYML:',
     &      ITHEMSYM,ITHESYML,ITHESYML
          WRITE(LUNGFO,*)'     THEGROTU,THEGROTL:'
     &      ,SNGL(SCALMOD),SNGL(THEGROTU),SNGL(THEGROTL)
          WRITE(LUNGFO,*)
     &      '     RADIMOD, ZLENMOD, CENMODX, CENMODY, CENMODZ, THEROT, BCMOD:'
          DO I=1,NMAGMOD
            WRITE(LUNGFO,*)
     &        SNGL(RADIMOD(I)),SNGL(ZLENMOD(I)),SNGL(CENMODX(I))
     &        ,SNGL(CENMODY(I))
     &        ,SNGL(CENMODZ(I))
     &        ,SNGL(THEROT(I)),SNGL(BCMOD(I))
          ENDDO
        ENDIF

C APPLY TAPERFUNCTION {

        IF (IHTAPER.NE.0) THEN

          OPEN(unit=LUNFT,file=FILEFTH,status='old')
          READ(LUNFT,'(A80)')COMTAP
          READ(LUNFT,*)XSCL,TSCL
          I=1
100       IF (I.GT.NKLOTZ) THEN
            WRITE(6,*)
     &        '*** ERROR IN REC_INIT: DIMENSION NKLOTZ EXCEEDED  ***'
            WRITE(6,*)
     &        '*** (OPTION IHTAPER)'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(LUNGFO,*)
     &        '*** ERROR IN REC_INIT: DIMENSION NKLOTZ EXCEEDED  ***'
            WRITE(LUNGFO,*)
     &        '*** (OPTION IHTAPER)'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            STOP
          ENDIF
          READ(LUNFT,*,END=900)FTX(I),FTY(I)
          FTX(I)=FTX(I)*XSCL
          FTY(I)=FTY(I)*TSCL
          I=I+1
          GOTO 100
900       NFUNTAP=I-1
          CLOSE(LUNFT)

          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)'      taper function applied form file:'
          WRITE(LUNGFO,*)'      ',FILEFTH
          WRITE(LUNGFO,*)'      user comment on file:'
          WRITE(LUNGFO,*)'      ',COMTAP
          WRITE(LUNGFO,*)' '

          CALL UTIL_SPLINE_COEF(FTX,FTY,NFUNTAP,-9999.0d0,-9999.0d0,FT2,F1,F2,F3,F4)

          DO I=1,IMAGTOT

            CALL UTIL_SPLINE_INTER(FTX,FTY,FT2,NFUNTAP,DX(I),DYFT,0)
c REALLY ASKING FOR DY!!
            IF (DY(I).GE.0.D0) DZ(I)=DZ(I)+DYFT
            IF (DY(I).LT.0.D0) DZ(I)=DZ(I)-DYFT

          ENDDO

        ENDIF  !IHTAPER

C APPLY TAPERFUNCTION }

C APPLY TAPERFUNCTION {

        IF (IVTAPER.NE.0) THEN

          OPEN(unit=LUNFT,file=FILEFTV,status='old')
          READ(LUNFT,'(A80)')COMTAP
          READ(LUNFT,*)XSCL,TSCL
          I=1
101       IF (I.GT.NKLOTZ) THEN
            WRITE(6,*)
     &        '*** ERROR IN REC_INIT: DIMENSION NKLOTZ EXCEEDED  ***'
            WRITE(6,*)
     &        '*** (OPTION IVTAPER)'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            WRITE(LUNGFO,*)
     &        '*** ERROR IN REC_INIT: DIMENSION NKLOTZ EXCEEDED  ***'
            WRITE(LUNGFO,*)
     &        '*** (OPTION IVTAPER)'
            WRITE(LUNGFO,*)'INCREASE PARAMETER NKLOTZ IN FILE KLOTZ.CMN'
            STOP
          ENDIF
          READ(LUNFT,*,END=901)FTX(I),FTY(I)
          FTX(I)=FTX(I)*XSCL
          FTY(I)=FTY(I)*TSCL
          I=I+1
          GOTO 101
901       NFUNTAP=I-1
          CLOSE(LUNFT)

          WRITE(LUNGFO,*)' '
          WRITE(LUNGFO,*)'      taper function applied form file:'
          WRITE(LUNGFO,*)'      ',FILEFTV
          WRITE(LUNGFO,*)'      user comment on file:'
          WRITE(LUNGFO,*)'      ',COMTAP
          WRITE(LUNGFO,*)' '

          CALL UTIL_SPLINE_COEF(FTX,FTY,NFUNTAP,-9999.0d0,-9999.0d0,FT2,F1,F2,F3,F4)

          DO I=1,IMAGTOT

            CALL UTIL_SPLINE_INTER(FTX,FTY,FT2,NFUNTAP,DX(I),DYFT,0)
            IF (DY(I).GE.0.D0) DY(I)=DY(I)+DYFT
            IF (DY(I).LT.0.D0) DY(I)=DY(I)-DYFT

          ENDDO

        ENDIF  !IVTAPER

C APPLY TAPERFUNCTION }


C SORT MAGNETS {

        DO I=1,IMAGTOT
          DXS(I)=DX(I)
        ENDDO

        CALL UTIL_SORT_FUNC(IMAGTOT,DX,DY)

        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO
        CALL UTIL_SORT_FUNC(IMAGTOT,DX,DZ)

        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO
        CALL UTIL_SORT_FUNC(IMAGTOT,DX,THETA)

        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO
        CALL UTIL_SORT_FUNC(IMAGTOT,DX,PHI)

        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO

        CALL UTIL_SORT_FUNC(IMAGTOT,DX,BC)
        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO

        CALL UTIL_SORT_FUNC(IMAGTOT,DX,XLEN)
        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO

        CALL UTIL_SORT_FUNC(IMAGTOT,DX,YLEN)
        DO I=1,IMAGTOT
          DX(I)=DXS(I)
        ENDDO

        CALL UTIL_SORT_FUNC(IMAGTOT,DX,ZLEN)

      ENDIF   !IKRESTOR

      IF (IKRESTOR.LT.0) THEN

        open(unit=99,FILE='rec.store',STATUS='NEW')

        WRITE(99,'(I10,A80)')ICODE,CODE
        WRITE(99,*)IMAGTOT
        DO IMAG=1,IMAGTOT
          WRITE(99,*)XLEN(IMAG),YLEN(IMAG),ZLEN(IMAG)
          WRITE(99,*)THETA(IMAG),PHI(IMAG),BC(IMAG)
          WRITE(99,*)DX(IMAG),DY(IMAG),DZ(IMAG)
        ENDDO

        CLOSE(99)

      ELSE IF (IKRESTOR.GT.0) THEN

        open(unit=99,FILE='rec.restore',STATUS='OLD')

        READ(99,*)IMAGTOT
        READ(99,*)IMAGTOT
        DO IMAG=1,IMAGTOT
          READ(99,*)XLEN(IMAG),YLEN(IMAG),ZLEN(IMAG)
          READ(99,*)THETA(IMAG),PHI(IMAG),BC(IMAG)
          READ(99,*)DX(IMAG),DY(IMAG),DZ(IMAG)
        ENDDO

        CLOSE(99)

      ENDIF !IKRESTOR

      IF (KBPOLYMAG.LT.0) THEN

        CSTAR='* '

        open(unit=99,FILE='wave_to_polymag.dat',STATUS='NEW',RECL=128)

        WRITE(99,*)'*** ------------------------------------------------'
        WRITE(99,*)'*** Begin of wave_to_polymag.dat'

        WRITE(99,*)CSTAR,ICODE,CODE(1:64)
        WRITE(99,*)IMAGTOT,' ! number of magnets'

        DO IMAG=1,IMAGTOT

          COSPHI=DCOS(PHI(IMAG))
          SINPHI=DSIN(PHI(IMAG))
          COSTHE=DCOS(THETA(IMAG))
          SINTHE=DSIN(THETA(IMAG))

          IF (ABS(COSTHE).LT.1.0E-12) COSTHE=0.0D0
          IF (ABS(COSPHI).LT.1.0E-12) COSPHI=0.0D0
          IF (ABS(SINTHE).LT.1.0E-12) SINTHE=0.0D0
          IF (ABS(SINPHI).LT.1.0E-12) SINPHI=0.0D0

          IF (ABS(COSTHE-1.0D0).LT.1.0E-12) COSTHE=1.0D0
          IF (ABS(SINTHE-1.0D0).LT.1.0E-12) SINTHE=1.0D0
          IF (ABS(COSPHI-1.0D0).LT.1.0E-12) COSPHI=1.0D0
          IF (ABS(SINPHI-1.0D0).LT.1.0E-12) SINPHI=1.0D0

          IF (ABS(COSTHE+1.0D0).LT.1.0E-12) COSTHE=-1.0D0
          IF (ABS(SINTHE+1.0D0).LT.1.0E-12) SINTHE=-1.0D0
          IF (ABS(COSPHI+1.0D0).LT.1.0E-12) COSPHI=-1.0D0
          IF (ABS(SINPHI+1.0D0).LT.1.0E-12) SINPHI=-1.0D0

          B0=BC(IMAG)

          VBX=COSPHI*SINTHE
          VBY=       COSTHE
          VBZ=SINPHI*SINTHE

          IF (B0.LT.0) THEN
            VBX=-VBX
            VBY=-VBY
            VBZ=-VBZ
            B0=-B0
          ENDIF

          write(99,*)'* magnet ',imag
          WRITE(99,*)sngl(DX(IMAG)),sngl(DY(IMAG)),sngl(DZ(IMAG)),
     &      ' !position of magnet'
          WRITE(99,*)B0,VBX,VBY,VBZ,
     &      ' !length bc and components of mag. vector (after rotation)'
          WRITE(99,*)'-6 1'
          WRITE(99,*)SNGL(ABS(XLEN(IMAG))),SNGL(ABS(YLEN(IMAG))),
     &      SNGL(ABS(ZLEN(IMAG))),
     &      ' !x,y,z dimensions of rectangular magnet'
          write(99,*)'*'
        ENDDO

        WRITE(99,*)'*** End of wave_to_polymag.dat'
        WRITE(99,*)'*** ------------------------------------------------'

        CLOSE(99)

      ENDIF !IKBPOLYMAG
C DEFAULT TRACKING RANGE

      DO I=1,IMAGTOT

        IF (XSTART.EQ.9999..AND.BC(I).NE.0.0)
     &    XSTART=DX(I)/1000.D0-RANGREC
        IF (XSTOP.EQ.9999..AND.BC(IMAGTOT+1-I).NE.0.0)
     &    XSTOP= DX(IMAGTOT+1-I)/1000.D0+RANGREC

      ENDDO

C SORT MAGNETS }

      IF (IPLREC.NE.0) CALL REC_PLOTM
      IF (IPLREC.EQ.9999.OR.IPLREC.EQ.-9999)
     &  STOP '*** WAVE TERMINATED DUE TO IPLREC***'

      RETURN
      end
