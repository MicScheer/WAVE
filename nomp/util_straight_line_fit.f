*CMZ :  2.69/00 26/10/2012  09.35.15  by  Michael Scheer
*-- Author :    Michael Scheer   25/10/2012
      subroutine util_straight_line_fit(n,x,y,e,a,b,chi2,erra,errb,istat)
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

c Fit a and b, such that f = a*x+b and Sum(i,((a*xi+b-yi-fi)/ei)**2)=min
c See also grafit.kumac

c The uncertainties erra and errb of the fitted parameters a (slope) and b
c (offset) agree to the fits of MINUIT, expept when the error e(i) are zero.
c Then the errors are such, that chi2/ndf is scaled to unity

      integer istat,n,i

      double precision x(n),y(n),e(n),chi2,wsum,wxy,wx,wy,wx2,
     &  d,ep,erra,errb,a,b,w,esuma

      wx=0.0d0
      wy=0.0d0
      wxy=0.0d0
      wx2=0.0d0
      wsum=0.0d0
      chi2=0.0d0
      erra=0.0d0
      errb=0.0d0
      istat=0

      if (n.lt.2) then
        istat=-1
        return
      else if (n.eq.2) then
        if (x(2).eq.x(1)) then
          istat=-2
          return
        endif
        a=(y(2)-y(1))/(x(2)-x(1))
        b=y(2)-a*x(2)
        return
      endif

      ep=1.0d0
      esuma=0.0d0
      do i=1,n
        if (e(i).lt.0.0d0) then
          istat=1
          return
        else if (e(i).eq.0.0d0) then
          ep=0.0d0
        endif
        esuma=esuma+abs(e(i))
      enddo

      if(esuma.ne.0.0d0.and.ep.eq.0.0d0) then
        istat=2
        return
      endif

      do i=1,n
        if (esuma.ne.0.0d0) then
          w=1.0d0/e(i)**2
        else
          w=1.0d0
        endif
        wx=wx+w*x(i)
        wy=wy+w*y(i)
        wx2=wx2+w*x(i)*x(i)
        wxy=wxy+w*x(i)*y(i)
        wsum=wsum+w
      enddo

      d=wsum*wx2-wx**2

      if (d.ne.0.0d0) then

        a=(wsum*wxy-wx*wy)/d
        b=(wx2*wy-wx*wxy)/d

        if (esuma.ne.0.0d0) then
          do i=1,n
            chi2=chi2+((a*x(i)+b-y(i))/e(i))**2
          enddo
          erra=sqrt(wsum/d)
          errb=sqrt(wx2/d)
        else
          do i=1,n
            chi2=chi2+(a*x(i)+b-y(i))**2
          enddo
          erra=sqrt(wsum/d*chi2/(n-2))
          errb=sqrt(wx2/d*chi2/(n-2))
        endif

      else
        istat=-1
        return
      endif

      return
      end
