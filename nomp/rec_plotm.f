*CMZ :  3.06/00 15/02/2019  13.33.23  by  Michael Scheer
*CMZ :  3.03/02 08/12/2015  15.33.15  by  Michael Scheer
*CMZ :  3.02/00 10/10/2014  13.45.06  by  Michael Scheer
*CMZ :  3.01/10 19/08/2014  11.16.18  by  Michael Scheer
*CMZ :  3.01/09 12/08/2014  17.18.45  by  Michael Scheer
*CMZ :  3.01/08 11/08/2014  19.17.31  by  Michael Scheer
*CMZ :  3.01/05 12/06/2014  09.31.15  by  Michael Scheer
*CMZ :  3.01/02 28/01/2014  17.00.03  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/11 20/02/2013  15.08.59  by  Michael Scheer
*CMZ :  2.57/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/05 11/02/2005  15.24.43  by  Michael Scheer
*CMZ :  2.52/09 29/10/2004  11.38.57  by  Michael Scheer
*CMZ :  2.52/05 13/08/2004  08.50.14  by  Michael Scheer
*CMZ :  2.48/04 16/03/2004  10.48.47  by  Michael Scheer
*CMZ :  2.48/03 03/03/2004  12.49.39  by  Michael Scheer
*CMZ :  2.48/00 26/02/2004  14.13.15  by  Michael Scheer
*CMZ :  2.47/19 01/12/2003  08.23.40  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.40/00 11/03/2002  17.29.02  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.20/05 13/03/2001  13.41.14  by  Michael Scheer
*CMZ :  2.20/01 10/11/2000  11.27.21  by  Michael Scheer
*CMZ :  2.16/07 22/09/2000  10.44.05  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.46  by  Michael Scheer
*CMZ :  2.14/00 16/04/2000  14.58.50  by  Michael Scheer
*CMZ :  2.13/11 19/03/2000  11.30.34  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.17.28  by  Michael Scheer
*CMZ :  1.03/06 30/06/98  12.42.25  by  Michael Scheer
*CMZ :  1.02/03 14/01/98  10.02.58  by  Michael Scheer
*CMZ :  1.01/00 27/11/97  15.36.37  by  Michael Scheer
*CMZ :  1.00/00 23/07/97  14.06.18  by  Michael Scheer
*CMZ : 00.02/00 25/11/96  14.10.09  by  Michael Scheer
*CMZ : 00.01/08 23/06/95  20.54.17  by  Michael Scheer
*CMZ : 00.01/07 22/03/95  12.25.18  by  Michael Scheer
*CMZ : 00.00/01 03/03/95  15.54.17  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE REC_PLOTM
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

C--- TO PLOT MAGNET CONFIGURATION

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,bforce.
      include 'bforce.cmn'
*KEEP,mplot.
      include 'mplot.cmn'
*KEND.

      INTEGER NXZON,NYZON,IFIRST,I,ISIGN,lundum

      CHARACTER(60) CHOPT
      CHARACTER(1) GANS
      CHARACTER(60) XTIT,YTIT
      CHARACTER(80) TITLE
      character(16) c16
      character c1
      byte ic1

      REAL*4 XP(5),YP(5),X1,X2,Y1,Y2,Z1,Z2,XMIN,XMAX,ZMIN,ZMAX
      REAL*4 XLEFT,XRIGHT,YUP,YDOWN,ZUP,ZDOWN,THE,SPHI,PI

      REAL*4 PSCAL,hpic,hgap,hbase

      real xleft_ps,ybottom_ps,xright_ps,ytop_ps,xsiz_ps,ysiz_ps
      real xlefto_ps,ybottomo_ps,xrighto_ps,ytopo_ps,xsizo_ps,ysizo_ps

      equivalence(ic1,c1)

      DATA GANS/'"'/

      DATA PI/3.141592653589793D0/

      if (iplrec.eq.0) return

      write(c16,*)icode
      do i=1,16
        c1=c16(i:i)
        if (ic1.ne.32) goto 1
      enddo

1     continue

      title='Run '//c16(i:len_trim(c16))//' '//CODE(1:min(70,len_trim(code)))

      lundum=223344
      call util_get_free_lun(lundum)

      fileeps_mshplt='rec_plotm.eps'
      viewer_mshplt=''

      call mshplt_set_box(0)
      call mshplt_set_scale(20.)
      call mshplt_init(lundum,-20.,-20.,0,0,575,575,'rec_plotm.eps','','',0.0)
      call mshplt_set_title_offset(-3.,1.5)
      call mshplt_hplset('GSIZ',0.4)
      call mshplt_title(title(1:len_trim(title)))

      hpic=1.9
      hbase=1.
      hgap=1.

      call mshplt_zone(1,3,1,'')

      X1=1.E30
      X2=-1.E30
      Y1=1.E30
      Y2=-1.E30
      Z1=1.E30
      Z2=-1.E30

      DO I=1,IMAGTOT
        IF (BC(I).NE.0.0) THEN
          XLEFT=DX(I)-XLEN(I)/2.D0
          XRIGHT=DX(I)+XLEN(I)/2.D0
          YDOWN=DY(I)-YLEN(I)/2.D0
          YUP=DY(I)+YLEN(I)/2.D0
          ZDOWN=DZ(I)-ZLEN(I)/2.D0
          ZUP=DZ(I)+ZLEN(I)/2.D0
          IF (XLEFT.LT.X1) X1=XLEFT
          IF (XRIGHT.GT.X2) X2=XRIGHT
          IF (YDOWN.LT.Y1) Y1=YDOWN
          IF (YUP.GT.Y2) Y2=YUP
          IF (ZDOWN.LT.Z1) Z1=ZDOWN
          IF (ZUP.GT.Z2) Z2=ZUP
        ENDIF
      ENDDO   !IMAGTOT

      XMIN=X1-0.1*(X2-X1)
      XMAX=X2+0.1*(X2-X1)
      ZMIN=Z1-0.15*(Z2-Z1)
      ZMAX=Z2+0.15*(Z2-Z1)

      IF (RPLXMN.EQ.-9999.) THEN
        XMIN=XSTART*1000.
      ELSE IF (RPLXMN.NE.9999.) THEN
        XMIN=RPLXMN
      ENDIF
      IF (RPLXMX.EQ.-9999.) THEN
        XMAX=XSTOP*1000.
      ELSE IF (RPLXMX.NE.9999.) THEN
        XMAX=RPLXMX
      ENDIF

      IF (RPLZMN.NE.9999.) ZMIN=RPLZMN
      IF (RPLZMX.NE.9999.) ZMAX=RPLZMX

      PSCAL=500./(XMAX-XMIN)

      IF (PSCAL.GT.20.) PSCAL=20.

c      call mshplt_set_character_height(0.23)
c      call mshplt_set_tic_size(0.065)

      call mshplt_frame(xmin,xmax,zmin,zmax,'','','LBrtc')

      XP(1)=XMIN
      XP(2)=XMAX
      YP(1)=0.
      YP(2)=0.

      call mshplt_pline(2,xp,yp)
      call mshplt_set_line_color(2,0,0,0)

      DO I=1,IMAGTOT
        XLEFT=DX(I)-XLEN(I)/2.D0
        XRIGHT=DX(I)+XLEN(I)/2.D0
        YDOWN=DY(I)-YLEN(I)/2.D0
        YUP=DY(I)+YLEN(I)/2.D0
        ZDOWN=DZ(I)-ZLEN(I)/2.D0
        ZUP=DZ(I)+ZLEN(I)/2.D0
        THE=180.D0/PI*THETA(I)
        SPHI=180.D0/PI*PHI(I)
        IF (BC(I).GT.0) THEN
          ISIGN=1
        ELSE
          ISIGN=-1
        ENDIF
        IF(YDOWN.GE.0.D0) THEN
          IF (BC(I).NE.0.0)
     &      CALL REC_PLOTK(I,XLEFT,XRIGHT,ZDOWN,ZUP,THE,SPHI,0,ISIGN,PSCAL)
        ELSE IF (YUP.GT.0.0.AND.YDOWN.LT.0.0) THEN
          WRITE(16,*)
          WRITE(16,*)'*** WARNING SR REC_PLOTM: STRANGE INPUT ON FILE REC.PAR***'
          WRITE(16,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR REC_PLOTM: STRANGE INPUT ON FILE REC.PAR***'
          WRITE(6,*)
        ENDIF
      ENDDO !IMAGTOT

      call mshplt_set_line_color(1,0,0,0)

      IF (IBFORCE.NE.0) THEN
        XP(1)=(BFCENX-BFLENX/2.)*1000.
        YP(1)=(BFCENZ-BFLENZ/2.)*1000.
        XP(3)=(BFCENX+BFLENX/2.)*1000.
        YP(3)=(BFCENZ+BFLENZ/2.)*1000.
        XP(2)=XP(3)
        YP(2)=YP(1)
        XP(4)=XP(1)
        YP(4)=YP(3)
        XP(5)=XP(1)
        YP(5)=YP(1)
        call mshplt_set_line_color(3,0,0,0)
        call mshplt_pline(5,xp,yp)
        call mshplt_set_line_color(1,0,0,0)
      ENDIF !IBFORCE

      call mshplt_reset_clipping
      call mshplt_set_character_height(0.35)
      call mshplt_text_NDC(0.8,-0.3,'long. coord. [mm]')
      call mshplt_set_text_angle(90.)
      call mshplt_text_NDC(-0.08,0.1,'trans. coord. [mm]')
      call mshplt_set_text_angle(0.)
      call mshplt_set_character_height(0.4)
      call mshplt_text_NDC(0.4,-0.3,'upper magnets')

      call mshplt_frame(xmin,xmax,zmin,zmax,'','','LBrtc')

      XP(1)=XMIN
      XP(2)=XMAX
      YP(1)=0.
      YP(2)=0.

      call mshplt_pline(2,xp,yp)
      call mshplt_set_line_color(4,0,0,0)

      DO I=1,IMAGTOT
        XLEFT=DX(I)-XLEN(I)/2.D0
        XRIGHT=DX(I)+XLEN(I)/2.D0
        YDOWN=DY(I)-YLEN(I)/2.D0
        YUP=DY(I)+YLEN(I)/2.D0
        ZDOWN=DZ(I)-ZLEN(I)/2.D0
        ZUP=DZ(I)+ZLEN(I)/2.D0
        THE=180.D0/PI*THETA(I)
        SPHI=180.D0/PI*PHI(I)
        IF (BC(I).GT.0) THEN
          ISIGN=1
        ELSE
          ISIGN=-1
        ENDIF
        IF(YUP.LE.0.D0) THEN
          IF (BC(I).NE.0.0)
     &      CALL REC_PLOTK(I,XLEFT,XRIGHT,ZDOWN,ZUP,THE,SPHI,0,ISIGN,PSCAL)
        ENDIF
      ENDDO !IMAGTOT
      call mshplt_set_line_color(1,0,0,0)

      IF (IBFORCE.NE.0) THEN
        XP(1)=(BFCENX-BFLENX/2.)*1000.
        YP(1)=(BFCENZ-BFLENZ/2.)*1000.
        XP(3)=(BFCENX+BFLENX/2.)*1000.
        YP(3)=(BFCENZ+BFLENZ/2.)*1000.
        XP(2)=XP(3)
        YP(2)=YP(1)
        XP(4)=XP(1)
        YP(4)=YP(3)
        XP(5)=XP(1)
        YP(5)=YP(1)
        call mshplt_set_line_color(3,0,0,0)
        call mshplt_pline(5,xp,yp)
        call mshplt_set_line_color(1,0,0,0)
      ENDIF !IBFORCE

      YP(2)=AMAX1(ABS(Y1),ABS(Y2))*1.3
      YP(1)=-YP(2)

      IF (RPLYMN.NE.9999.) YP(1)=RPLYMN
      IF (RPLYMX.NE.9999.) YP(2)=RPLYMX

      call mshplt_reset_clipping
      call mshplt_set_character_height(0.35)
      call mshplt_text_NDC(0.8,-0.3,'long. coord. [mm]')
      call mshplt_set_text_angle(90.)
      call mshplt_text_NDC(-0.08,0.1,'trans. coord. [mm]')
      call mshplt_set_text_angle(0.)
      call mshplt_set_character_height(0.4)
      call mshplt_text_NDC(0.4,-0.3,'lower magnets')

      call mshplt_frame(xmin,xmax,yp(1),yp(2),'','','LBrtc')

      XP(1)=XMIN
      XP(2)=XMAX
      YP(1)=0.
      YP(2)=0.
      call mshplt_pline(2,xp,yp)

      DO I=1,IMAGTOT
        XLEFT=DX(I)-XLEN(I)/2.D0
        XRIGHT=DX(I)+XLEN(I)/2.D0
        YDOWN=DY(I)-YLEN(I)/2.D0
        YUP=DY(I)+YLEN(I)/2.D0
        ZDOWN=DZ(I)-ZLEN(I)/2.D0
        ZUP=DZ(I)+ZLEN(I)/2.D0
        THE=180.D0/PI*THETA(I)
        SPHI=180.D0/PI*PHI(I)
        IF(DZ(I).LT.0.0D0) THEN
          call mshplt_set_line_style(1)
        ELSE IF(DZ(I).GT.0.D0) THEN
          call mshplt_set_line_style(3)
        ENDIF
        IF (BC(I).NE.0.0) then
          if (ydown.gt.0) then
            call mshplt_set_line_color(2,0,0,0)
          else
            call mshplt_set_line_color(4,0,0,0)
          endif
          CALL REC_PLOTK(I,XLEFT,XRIGHT,YDOWN,YUP,THE,SPHI,1,0,PSCAL)
        endif
      ENDDO !IMAGTOT
      call mshplt_set_line_color(1,0,0,0)

      IF (IBFORCE.NE.0) THEN
        XP(1)=(BFCENX-BFLENX/2.)*1000.
        YP(1)=(BFCENY-BFLENY/2.)*1000.
        XP(3)=(BFCENX+BFLENX/2.)*1000.
        YP(3)=(BFCENY+BFLENY/2.)*1000.
        XP(2)=XP(3)
        YP(2)=YP(1)
        XP(4)=XP(1)
        YP(4)=YP(3)
        XP(5)=XP(1)
        YP(5)=YP(1)
        call mshplt_set_line_color(3,0,0,0)
        call mshplt_pline(5,xp,yp)
        call mshplt_set_line_color(1,0,0,0)
      ENDIF !IBFORCE

      call mshplt_reset_clipping
      call mshplt_set_character_height(0.35)
      call mshplt_text_NDC(0.8,-0.3,'long. coord. [mm]')
      call mshplt_set_text_angle(90.)
      call mshplt_text_NDC(-0.08,0.1,'vert. coord. [mm]')
      call mshplt_set_text_angle(0.)
      call mshplt_set_character_height(0.4)
      call mshplt_text_NDC(0.45,-0.3,'side view')

C- TERMINATE PLOTTING

9999  call mshplt_end

      RETURN
      END
