*CMZ :  3.01/09 12/08/2014  15.00.42  by  Michael Scheer
*CMZ :  3.01/08 12/08/2014  11.36.16  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.13/09 08/03/2000  17.22.30  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  14.25.26  by  Michael Scheer
*-- Author :    Michael Scheer   21/06/95
      SUBROUTINE ARROW(XL,XR,YU,YD,ANG)

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

      REAL*4 X(2),Y(2),XC,YC,XL,XR,YU,YD,Y1,Y2,Y3,Y4,X1,X2,X3,X4,ANG,DX,DY
      real xmin_ps,xmax_ps,ymin_ps,ymax_ps

      call mshplt_get_frame(xmin_ps,xmax_ps,ymin_ps,ymax_ps)

      XC=(XR+XL)/2.
      YC=(YU+YD)/2.
      DX=XR-XL
      DY=YU-YD

      if (yc.lt.ymin_ps) return

      IF (ANG.EQ.270..OR.ANG.EQ.-90.) THEN
        X1=XC+DX/10.
        X2=XC
        X3=XC-DX/10.
        X4=X2
        Y1=YC+DY/10.
        Y2=YU-DY/10.
        Y3=Y1
        Y4=YD+DY/10.
        X(1)=X2
        Y(1)=Y2
        X(2)=X1
        Y(2)=Y1

        call mshplt_pline(2,x,y)
        X(2)=X3
        Y(2)=Y3
        call mshplt_pline(2,x,y)
        X(2)=X4
        Y(2)=Y4
        call mshplt_pline(2,x,y)
      ELSE IF (ANG.EQ.90..OR.ANG.EQ.-270.) THEN
        X1=XC+DX/10.
        X2=XC
        X3=XC-DX/10.
        X4=X2
        Y1=YC-DY/10.
        Y2=YU-DY/10.
        Y3=Y1
        Y4=YD+DY/10.
        X(1)=X4
        Y(1)=Y4
        X(2)=X1
        Y(2)=Y1
        call mshplt_pline(2,x,y)
        X(2)=X3
        Y(2)=Y3
        call mshplt_pline(2,x,y)
        X(2)=X2
        Y(2)=Y2
        call mshplt_pline(2,x,y)
      ELSE IF (ANG.EQ.0.0 .OR. ANG.EQ.360.) THEN
        X1=XL+DX/5.
        X2=XR-DX*2./5.
        X3=XR-DX/5.
        X4=X2
        Y1=YC
        Y2=YC-DY/5.
        Y3=Y1
        Y4=YC+DY/5.
        X(1)=X3
        Y(1)=Y3
        X(2)=X1
        Y(2)=Y1
        call mshplt_pline(2,x,y)
        X(2)=X2
        Y(2)=Y2
        call mshplt_pline(2,x,y)
        X(2)=X4
        Y(2)=Y4
        call mshplt_pline(2,x,y)
      ELSE IF (ANG.EQ.180..OR.ANG.EQ.-180.) THEN
        X1=XL+DX/5.
        X2=XL+DX*2./5.
        X3=XR-DX/5.
        X4=X2
        Y1=YC
        Y2=YC-DY/5.
        Y3=Y1
        Y4=YC+DY/5.
        X(1)=X1
        Y(1)=Y1
        X(2)=X3
        Y(2)=Y3
        call mshplt_pline(2,x,y)
        X(2)=X2
        Y(2)=Y2
        call mshplt_pline(2,x,y)
        X(2)=X4
        Y(2)=Y4
        call mshplt_pline(2,x,y)
      ELSE
        WRITE(6,*)'*** WARNING SR ARROW: STRANGE ANGLE',ANG,'***'
        RETURN
      ENDIF

      RETURN
      END
