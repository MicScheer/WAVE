*CMZ :  4.01/02 17/03/2023  21.47.36  by  Michael Scheer
*CMZ :  2.57/04 01/02/2006  15.08.35  by  Michael Scheer
*CMZ :  2.57/00 03/11/2005  16.29.31  by  Michael Scheer
*CMZ :  2.56/00 11/10/2005  15.29.20  by  Michael Scheer
*-- Author :    Michael Scheer   10/10/2005
      subroutine thickapp(icirc,al,ar,alpha0,corr)
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

      double precision pi
      parameter (PI=3.14159265359d0)

      double precision al,ar,corr,ys,xs,f,f0,phi,alpha0

      integer icirc

      if (icirc.ne.0) then

C Berechnet effektive Oeffnung durch einen Kollimator
C Ansatz: Schnittflaeche zweier verschobener Kreise

        xs=alpha0*al/2.0d0
        ys=sqrt(ar**2-xs**2)
        phi=2.0d0*acos(xs/ar)

        f0=pi*ar**2
        f=2.0d0*(f0*phi/(2.0d0*pi)-ys*xs)

        corr=f/f0

      else  !(icirc.ne.0) then

c Naeherung: Der cos wird vernachlaessigt, d.h. z.B. bei al=0 haengt die Apertur nicht
c         vom Winkel ab, was falsch ist, da sie mit cos(alpha) kleiner wird.

        corr=max(1.0d0-tan(alpha0)*al/(2.0d0*ar),0.0d0)

      endif  !(icirc.ne.0) then

      return
      end
