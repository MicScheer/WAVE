*CMZ :  2.70/05 02/01/2013  12.50.20  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.03.55  by  Michael Scheer
*CMZ :  2.63/03 15/05/2009  12.38.54  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_spline_integral_window(x,y,n,xmin,xmax,result
     &  ,coef,work1,work2,work3,work4,mode,istat)
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

c---  calculates integral of y(x) in interval [xmin,xmax] via splines

      implicit none

      real*8 x(*),y(*),result,xmin,xmax,ymin,ymax,cmin,cmax
      real*8 coef(*),work1(*),work2(*),work3(*),work4(*)

      integer i,n,mode,imn,imx,i1,i2,istat

      result=0.0d0
      istat=-1

      if (x(1).gt.x(n)) return

c---  spline-coefficients

      if (mode.lt.0) then
        call util_spline_coef(x,y,n,-9999.0d0,-9999.0d0,coef,work1,work2,work3,work4)
      endif

      if (xmin.lt.x(1)) then
        istat=-2
        return
      endif
      if (xmax.gt.x(n)) then
        istat=-3
        return
      endif

c intervals

      call util_bisec(n,x,xmin,imn)
      call util_bisec(n,x,xmax,imx)

      if (imn.eq.imx) then

c from xmin to xmax

        call util_spline_inter(x,y,coef,n,xmin,ymin,0)
        call util_spline_inter(x,y,coef,n,xmax,ymax,0)

        cmin=coef(imn)+(coef(imn+1)-coef(imn))/(x(imn+1)-x(imn))*(xmin-x(imn))
        cmax=coef(imn)+(coef(imn+1)-coef(imn))/(x(imn+1)-x(imn))*(xmax-x(imn))

        result=result+
     &    (xmax-xmin)*0.5d0
     &    *(ymin+ymax)-
     &    (xmax-xmin)**3/24.d0
     &    *(cmin+cmax)

        return

      endif

      i1=imn+1
      i2=imx-1

c--- integration

c from xmin to x(imn+1)
      call util_spline_inter(x,y,coef,n,xmin,ymin,0)
      cmin=coef(imn)+(coef(imn+1)-coef(imn))/(x(imn+1)-x(imn))*(xmin-x(imn))
      result=result+
     &  (x(imn+1)-xmin)*0.5d0
     &  *(ymin+y(imn+1))-
     &  (x(imn+1)-xmin)**3/24.d0
     &  *(cmin+coef(imn+1))

c from x(imx) to xmax
      call util_spline_inter(x,y,coef,n,xmax,ymax,0)
      cmax=coef(imx)+(coef(imx+1)-coef(imx))/(x(imx+1)-x(imx))*(xmax-x(imx))
      result=result+
     &  (xmax-x(imx))*0.5d0
     &  *(ymax+y(imx))
     &  -(xmax-x(imx))**3/24.d0
     &  *(cmax+coef(imx))

      do i=i1,i2

        result=result+
     &    (x(i+1)-x(i))*0.5d0
     &    *(y(i)+y(i+1))
     &    -(x(i+1)-x(i))**3/24.d0
     &    *(coef(i)+coef(i+1))

      enddo

      istat=0

      return
      end
