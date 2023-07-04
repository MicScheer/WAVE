*CMZ :  3.06/00 15/02/2019  14.50.44  by  Michael Scheer
*CMZ :  3.03/02 08/12/2015  13.47.21  by  Michael Scheer
*CMZ :  3.02/00 10/10/2014  12.26.46  by  Michael Scheer
*CMZ :  3.01/09 12/08/2014  15.15.26  by  Michael Scheer
*CMZ :  3.01/08 11/08/2014  21.00.00  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.57/05 05/12/2006  10.21.08  by  Michael Scheer
*CMZ :  2.52/09 29/10/2004  11.41.44  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.14/00 17/04/2000  15.08.07  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ :  1.00/00 02/06/97  10.57.45  by  Michael Scheer
*CMZ : 00.01/08 21/06/95  10.07.32  by  Michael Scheer
*CMZ : 00.01/07 09/03/95  14.40.13  by  Michael Scheer
*CMZ : 00.00/01 03/03/95  15.54.17  by  Michael Scheer
*-- Author : Michael Scheer
C*****************************************************************
      SUBROUTINE REC_PLOTK(IMAG,XLEFT,XRIGHT,YDOWN,YUP,
     &  THETA,PHI,MODE,ISIGN,PSCAL)
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

      IMPLICIT NONE

      INTEGER MODE,ISIGN,IMAG

      real xmin_ps,xmax_ps,ymin_ps,ymax_ps,rmsiz

      REAL*4 XLEFT,XRIGHT,YDOWN,YUP,THETA,PHI
      REAL*4 DX,DY,XP(2),YP(2),dxp
      REAL*4 XLEF,XRIG,YDOW,YU
      REAL*4 DOT,CIRC
      REAL*4 DOT0,CIRC0,rmtyp20,rmtyp24,rmtyp31
      REAL*4 PSCAL

      DATA CIRC0/2./
      DATA DOT0/1./

      data dot0/25./
      data circ0/5./
      data rmtyp20/20./
      data rmtyp24/-9999./
      data rmtyp31/31./

      call mshplt_get_frame(xmin_ps,xmax_ps,ymin_ps,ymax_ps)

      IF (THETA.NE.0.
     &    .AND.THETA.NE.90.
     &    .AND.THETA.NE.180.
     &    .AND.THETA.NE.270.) THEN
        CALL REC_PLOTTILT(IMAG,MODE)
C01JUN97 RETURN
C22MAY97         STOP '*** ERROR IN REC_PLOTK: WRONG ANGLE THETA (MUST BE 0 OR 90) ***'
        WRITE(6,*)
     &    '*** ERROR IN REC_PLOTK: WRONG ANGLE THETA (MUST BE 0, 90, 180 or 270) ***'
        STOP
      ENDIF

      IF (PHI.NE.0.AND.PHI.NE.90) THEN
        WRITE(6,*) '*** WARNING SR REC_PLOTK: WRONG ANGLE PHI (MUST BE 0 OR 90) ***'
C01JUN97   RETURN
        STOP
      ENDIF

      DX=0.0
      DY=0.0

      DOT=DOT0*PSCAL
      CIRC=CIRC0*PSCAL

      XLEF=XLEFT+DX
      YDOW=YDOWN+DY
      XRIG=XRIGHT+DX
      YU=YUP+DY

      XP(1)=XLEF
      XP(2)=XRIG
      YP(1)=YDOW
      YP(2)=YDOW

      dxp=(xmax_ps-xmin_ps)*0.02 !clip

      xp(1)=min(xp(1),xmax_ps)
      xp(2)=min(xp(2),xmax_ps)
      xp(1)=max(xp(1),xmin_ps)
      xp(2)=max(xp(2),xmin_ps)
      yp(1)=min(yp(1),ymax_ps)
      yp(2)=min(yp(2),ymax_ps)
      yp(1)=max(yp(1),ymin_ps)
      yp(2)=max(yp(2),ymin_ps)
      call mshplt_pline(2,xp,yp)

      XP(1)=XRIG
      XP(2)=XRIG
      YP(1)=YDOW
      YP(2)=YU

      if (xp(1).ge.xmin_ps+dxp.and.xp(1).le.xmax_ps-dxp) then
        yp(1)=min(yp(1),ymax_ps)
        yp(2)=min(yp(2),ymax_ps)
        yp(1)=max(yp(1),ymin_ps)
        yp(2)=max(yp(2),ymin_ps)
        call mshplt_pline(2,xp,yp)
      endif

      XP(1)=XRIG
      XP(2)=XLEF
      YP(1)=YU
      YP(2)=YU

      xp(1)=min(xp(1),xmax_ps)
      xp(2)=min(xp(2),xmax_ps)
      xp(1)=max(xp(1),xmin_ps)
      xp(2)=max(xp(2),xmin_ps)
      yp(1)=min(yp(1),ymax_ps)
      yp(2)=min(yp(2),ymax_ps)
      yp(1)=max(yp(1),ymin_ps)
      yp(2)=max(yp(2),ymin_ps)
      call mshplt_pline(2,xp,yp)

      XP(1)=XLEF
      XP(2)=XLEF
      YP(1)=YU
      YP(2)=YDOW

      if (xp(1).ge.xmin_ps+dxp.and.xp(1).le.xmax_ps-dxp) then
        yp(1)=min(yp(1),ymax_ps)
        yp(2)=min(yp(2),ymax_ps)
        yp(1)=max(yp(1),ymin_ps)
        yp(2)=max(yp(2),ymin_ps)
        call mshplt_pline(2,xp,yp)
      endif

      XP(1)=XLEF+(XRIG-XLEF)/2.
      YP(1)=YDOW+(YU-YDOW)/2.

      call mshplt_get_marker_size(rmsiz)

      IF(MODE.EQ.0) THEN

c        print*,imag,phi,theta,isign

        IF (PHI.EQ.0.0) THEN

          IF(THETA.EQ.0.0.AND.ISIGN.EQ.1.OR.THETA.EQ.180.AND.ISIGN.EQ.-1) THEN

            if (xp(1).ge.xmin_ps+dxp.and.xp(1).le.xmax_ps-dxp.and.
     &        (yp(1).gt.ymin_ps+rmsiz.and.yp(2).lt.ymax_ps+rmsiz.or.
     &        yp(2).gt.ymin_ps+rmsiz.and.yp(1).lt.ymax_ps+rmsiz)
     &          ) then
c              call mshplt_set_marker_type(1)
c              call mshplt_marker(1,xp,yp)
c              call mshplt_set_marker_type(9)
              call mgset('MTYP',rmtyp20)
              call mgset('MSCF',pscal)
              call mpm(1,xp,yp)
            endif

          ELSE IF(THETA.EQ.180.AND.ISIGN.EQ.1.OR.THETA.EQ.0.AND.ISIGN.EQ.-1) THEN

            if (xp(1).ge.xmin_ps+dxp.and.xp(1).le.xmax_ps-dxp.and.
     &        (yp(1).gt.ymin_ps+rmsiz.and.yp(2).lt.ymax_ps+rmsiz.or.
     &        yp(2).gt.ymin_ps+rmsiz.and.yp(1).lt.ymax_ps+rmsiz)
     &          ) then
              yp(1)=min(yp(1),ymax_ps)
              yp(2)=min(yp(2),ymax_ps)
              yp(1)=max(yp(1),ymin_ps)
              yp(2)=max(yp(2),ymin_ps)
c              call mshplt_scale_marker_size(1.5)
c              call mshplt_set_marker_type(6)
c              call mshplt_marker(1,xp,yp)
c              call mshplt_scale_marker_size(1./1.5)
c              call mshplt_set_marker_type(9)
              call mgset('MTYP',rmtyp31)
              call mgset('MSCF',pscal)
              call mpm(1,xp,yp)
            endif
          ELSE !THETA
            CALL ARROW(XLEF,XRIG,YU,YDOW,THETA-ISIGN*90.)
          ENDIF     !THETA
        ELSE IF(PHI.EQ.90.) THEN
          IF (THETA.EQ.90. .AND.ISIGN.GT.0 .OR.
     &        THETA.EQ.270..AND.ISIGN.LT.0) THEN
            CALL ARROW(XLEF,XRIG,YU,YDOW,270.)
          ELSE IF(THETA.EQ.270..AND.ISIGN.GT.0 .OR.
     &        THETA.EQ.90. .AND.ISIGN.LT.0) THEN
            CALL ARROW(XLEF,XRIG,YU,YDOW,90.)
          ELSE
            WRITE(6,*) '*** ERROR IN REC_PLOTK: WRONG ANGLE THETA FOR PHI=90. ***'
            RETURN
C             STOP
          ENDIF !THETA
        ENDIF    !PHI

      ENDIF   !MODE

      RETURN
      END
