*CMZ :  2.68/00 25/05/2012  11.03.55  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  10.48.03  by  Michael Scheer
*CMZ :  2.63/05 14/08/2009  12.21.34  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_spline_integral_buff(x,y,n,result,nbuff,margin
     &  ,coef,work1,work2,work3,work4,istat)
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

c---  calculates integral of y(x) via splines

      implicit none

      integer ibuff,n,istat,nbuff,margin,mbuff,i1,i2,nbfull

      double precision x(n),y(n),result,
     &  coef(n),work1(n),work2(n),work3(n),work4(n),sum

      istat=-1
      result=0.0d0

      mbuff=nbuff+2*margin

      if (n.le.mbuff.or.3*margin.ge.n) then
        call util_spline_integral(x,y,n,result
     &    ,coef,work1,work2,work3,work4)
        return
      endif

      nbfull=n/nbuff

      call util_spline_integral_window(
     &  x,y,nbuff+margin,x(1),x(nbuff),result
     &  ,coef,work1,work2,work3,work4,-1,istat)

      if (istat.ne.0) then
        result=0.0d0
        return
      endif

      i2=nbuff

      do ibuff=1,nbfull-2
        i1=i2
        i2=i1+nbuff
        call util_spline_integral_window(
     &    x(i1-margin),y(i1-margin),mbuff,x(i1),x(i2),sum
     &    ,coef,work1,work2,work3,work4,-1,istat)
        if (istat.ne.0) then
          result=0.0d0
          return
        endif
        result=result+sum
      enddo !nbfull

      call util_spline_integral_window(
     &  x(i2),y(i2),n-i2+1,x(i2),x(n),sum
     &  ,coef,work1,work2,work3,work4,-1,istat)

      if (istat.ne.0) then
        result=0.0d0
        return
      endif

      result=result+sum

      return
      end
