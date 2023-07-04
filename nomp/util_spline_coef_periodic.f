*CMZ :  4.01/02 17/03/2023  21.52.12  by  Michael Scheer
*CMZ :  2.66/18 10/12/2010  14.04.32  by  Michael Scheer
*CMZ : 00.00/07 10/12/2010  13.21.53  by  Michael Scheer
*-- Author :    Michael Scheer   10/12/2010
      subroutine util_spline_coef_periodic(x,y,n,ypp,aa,bb,cc,c,cn,ifail)
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

c--- calculates spline coefficients for periodic function
c--- the interval must be closed, i.e. x(n)-x(1)=periodelength and y(n)=y(1)

c--   input:

c-       n: number of x,y-values, must be at least five
c-       x: array of x-values
c-       y: array of y-values

c--   ouput:

c-       ypp:   spline-coefficients
c-     ifail:   error status

c--   workingspace: aa(n),bb(n),cc(n),c(n),cn(n)

      implicit none

      integer n,j,n1,n2,ifail
      real*8  x(n),y(n),ypp(n),aa(n),bb(n),cc(n),c(n),cn(n)

      if(n.lt.5) then
        ifail=-1
        return
      endif

      if(y(n).ne.y(1)) then
        ifail=-2
        return
      endif

      cn=0.0d0 !letzte Spalte der Matrix

c      y(n)=y(1)
c      ypp(n)=ypp(1)

c n=5
c      aa(1)*ypp(4)+bb(1)*ypp(1)+cc(1)*ypp(2)=c(1)
c      aa(2)*ypp(1)+bb(2)*ypp(2)+cc(2)*ypp(3)=c(2)
c      aa(3)*ypp(2)+bb(3)*ypp(3)+cc(3)*ypp(4)=c(3)
c      aa(4)*ypp(3)+bb(4)*ypp(4)+cc(4)*ypp(1)=c(4)

      n1=n-1
      n2=n-2

      aa(1)=(x(n)-x(n1))/6.d0
      bb(1)=(x(2)-x(1)+(x(n)-x(n1)))/3.d0
      cc(1)=(x(2)-x(1))/6.d0
      c(1)=(y(2)-y(1))/(x(2)-x(1))
     &    -(y(n)-y(n1))/(x(n)-x(n1))

      do j=2,n1
          aa(j)=(x(j  )-x(j-1))/6.d0
          bb(j)=(x(j+1)-x(j-1))/3.d0
          cc(j)=(x(j+1)-x(j  ))/6.d0
          c(j)=(y(j+1)-y(j  ))/(x(j+1)-x(j  ))
     &          -(y(j  )-y(j-1))/(x(j  )-x(j-1))
      enddo !j

c Auf Dreiecksmatrix bringen

c Oberste Zeile

      cc(1)=cc(1)/bb(1)
      aa(1)=aa(1)/bb(1)
      c(1) = c(1)/bb(1)
      bb(1)=1.0d0

c 2. Zeile, d.h. die erste regulaere, cn(j) ist die letzte Spalte der Matrix

      bb(2)=bb(2)/aa(2)
      cc(2)=cc(2)/aa(2)
      c(2) = c(2)/aa(2)
      aa(2)=1.0d0

      aa(2)=0.0d0
      bb(2)=bb(2)-cc(1)
      cn(2)=-aa(1)
      c(2) = c(2)-c(1)

      cc(2)=cc(2)/bb(2)
      cn(2)=cn(2)/bb(2)
      c(2) = c(2)/bb(2)

      bb(2)=1.0d0

c Nun die hoeheren bis n-3

      do j=3,n-3

        bb(j)=bb(j)/aa(j)
        cc(j)=cc(j)/aa(j)
        c(j) = c(j)/aa(j)

        aa(j)=0.0d0
        bb(j)=bb(j)-cc(j-1)
        cn(j)=-cn(j-1)
        c(j) = c(j)-c(j-1)

        cc(j)=cc(j)/bb(j)
        cn(j)=cn(j)/bb(j)
        c(j) = c(j)/bb(j)

        bb(j)=1.0d0

      enddo

c vorletzte Zeile

      bb(n2)=bb(n2)/aa(n2)
      cc(n2)=cc(n2)/aa(n2)
      c(n2) = c(n2)/aa(n2)
      aa(n2)=1.0d0

      aa(n2)=0.0d0
      bb(n2)=bb(n2)-cc(n2-1)
      cc(n2)=cc(n2)-cn(n2-1)
      c(n2) = c(n2)-c(n2-1)

      cc(n2)=cc(n2)/bb(n2)
      c(n2) = c(n2)/bb(n2)

      bb(n2)=1.0d0

c Letzte Zeile, benutze ypp als Puffer

      ypp=0.0d0
      ypp(n2)=aa(n1)/cc(n1)
      ypp(n1)=bb(n1)/cc(n1)
      c(n1)=c(n1)/cc(n1)
      cc(n1)=1.0d0
      ypp(1)=cc(n1)

c Oberste Zeile abziehen

      ypp(1)=0.0d0
      ypp(2)=-cc(1)
      ypp(n1)=ypp(n1)-aa(1)
      c(n1)=c(n1)-c(1)

      do j=2,n2

        c(n1)=c(n1)/ypp(j)
        ypp=ypp/ypp(j)

        ypp(j)=ypp(j)-bb(j)
        ypp(j+1)=ypp(j+1)-cc(j)
        ypp(n1)=ypp(n1)-cn(j)
        c(n1)=c(n1)-c(j)

      enddo

      c(n1)=c(n1)/ypp(n1)

c Ernten

      ypp(n1)=c(n1)

c Vorletzte Zeile

      bb(n2)=bb(n2)/cc(n2)
      c(n2)=c(n2)/cc(n2)
      cc(n2)=1.0d0

      c(n2)=c(n2)-c(n1)
      cc(n2)=0.0d0

      c(n2)=c(n2)/bb(n2)
      bb(n2)=1.0d0

      ypp(n2)=c(n2)

c Letzte Spalte nullen und normieren, letzte und vorletzte sind bereits fertig

      bb(1)=bb(1)/aa(1)
      cc(1)=cc(1)/aa(1)
      c(1)=c(1)/aa(1)
      aa(1)=1.0d0

      c(1)=c(1)-c(n1)
      aa(1)=0.0d0

      cc(1)=cc(1)/bb(1)
      c(1)=c(1)/bb(1)
      bb(1)=1.0d0

c Regulaere Zeilen
      do j=2,n-3

        bb(j)=bb(j)/cn(j)
        cc(j)=cc(j)/cn(j)
        c(j)=c(j)/cn(j)
        cn(j)=1.0d0

        c(j)=c(j)-c(n1)
        cn(j)=0.0d0

        cc(j)=cc(j)/bb(j)
        c(j)=c(j)/bb(j)
        bb(j)=1.0d0

      enddo

      do j=n-3,2,-1
         ypp(j)=c(j)-cc(j)*ypp(j+1)
      enddo

      ypp(1)=(c(1)-cc(1)*ypp(2)-aa(1)*ypp(n1))/bb(1)
      ypp(n)=ypp(1)

      ifail=0

      return
      end
