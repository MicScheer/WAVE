*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  0.99/00 16/02/2004  17.22.39  by  Michael Scheer
*CMZ :  0.00/08 21/01/2004  17.45.18  by  Michael Scheer
*CMZ :  0.00/06 14/01/2004  16.32.20  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  16.08.27  by  Michael Scheer
*CMZ :  0.00/04 19/12/2003  18.32.08  by  Michael Scheer
*-- Author :    Michael Scheer   19/12/2003
      subroutine areim(x,a,b,z,are,aim)

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

      implicit none

      double precision x,a,b,z,are,aim,arg1,x2,z2,a2,b2,y,az2,y2,b3

      x2=x*x
      z2=z*z
      a2=a*a
      b2=b*b
      b3=b2*b
      az2=a2*z2
      y=a*x+b
      y2=y*y

      arg1=sqrt(y2+x2+z2)

      are=2.0d0*z2*(
     &  arg1*((a*x-b)*(az2+b2))
     &  +(((a2+1.d0)*x2-b2+z2)*a2-b2)*z2
     &  +(( a2-1.d0)*x2-b2)*b2)

      aim=2.0d0*z*(
     &  arg1*((az2+b*y)*a*z2+b3*x)
     &  + ((y*a2+2.0d0*b)*z2 + b*(y2+2.d0*x2))*a*z2
     &  +y*b3*x)

      return
      end
