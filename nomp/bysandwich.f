*CMZ :  4.01/07 18/01/2025  16.41.01  by  Michael Scheer
*CMZ :  4.00/11 26/07/2021  09.08.58  by  Michael Scheer
*CMZ :  3.06/00 11/02/2019  12.49.34  by  Michael Scheer
*CMZ :  3.04/00 19/01/2018  16.33.13  by  Michael Scheer
*CMZ :  3.03/02 16/02/2016  12.18.47  by  Michael Scheer
*CMZ :  3.01/00 15/07/2013  08.04.32  by  Michael Scheer
*CMZ :  2.68/03 07/08/2012  13.09.30  by  Michael Scheer
*CMZ :  2.66/07 04/12/2009  16.11.19  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.02/00 19/12/97  16.15.53  by  Michael Scheer
*CMZ :  1.01/00 28/10/97  12.14.09  by  Michael Scheer
*CMZ : 00.01/08 01/04/95  16.54.24  by  Michael Scheer
*CMZ : 00.01/07 10/03/95  11.22.55  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.21.20  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.03  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.42  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE bysandwich(nmag,im,XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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

      INTEGER nmag,IM

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION XIN,YIN,ZIN,BX,BY,BZ,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION dist,fint,gap,fringe,fa,fb,fc,xrs(3),
     &  strength,pin(3),center(3),pout(3),vnin(3),vnout(3),b(3),edge(2),x,y,z,x2,x3,x4,x5,y2,y3

      integer istatus,modus

      x=xin
      y=yin
      z=zin

      BXOUT=0.0d0
      BYOUT=0.0d0
      BZOUT=0.0d0

      AXOUT=0.0d0
      AYOUT=0.0d0
      AZOUT=0.0d0

      BX=0.0d0
      BY=0.0d0
      BZ=0.0d0

      dist=dot_product([pmag(1,im),pmag(2,im), pmag(3,im)]-[x,y,z],[pmag(4,im),pmag(5,im), pmag(6,im)])

      if (dist.ge.0.0d0) return

      modus=int(pmag(13,im))
      strength=pmag(14,im)
      fringe=pmag(15,im)
      fa=pmag(16,im)
      fb=pmag(17,im)
      fc=pmag(18,im)

      if (dist.ge.-fringe) then
        x=-dist
        if (modus.eq.0) then
          byout=1.0d0
        else if (modus.eq.1) then
          bxout=y*fa
          byout=x*fa
        else if (modus.eq.3) then
          x2=x*x
          x3=x2*x
          bxout=(2.0d0*fa*x+3.0d0*fb*x2)*y
          byout=fa*x2+fb*x3+y**2*(-fa-3.0d0*fb*x)
        else if (modus.eq.5) then
          x2=x*x
          x3=x2*x
          x4=x3*x
          x5=x4*x
          y2=y*y
          y3=y2*y
          bxout=y*(3.0d0*Fa*x2+4.0d0*fb*x3+5.0d0*fc*x4)
     &      +y3*(-fa-4.0d0*fb*x-10.0d0*fc*x2) !This term is not Maxwell conform
          !byout=(fa*x3+fb*x4+fc*x5)+y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2) ! The sign seems to be wrong in the manual
          byout=(fa*x3+fb*x4+fc*x5)-y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2)
        endif
        byout=byout*strength
        bxout=-bxout*strength
        return
      endif

      dist=dot_product([pmag(7,im),pmag(8,im), pmag(9,im)]-[x,y,z],[pmag(10,im),pmag(11,im), pmag(12,im)])
      if (dist.lt.0.0d0) return

      if (dist.le.fringe) then
        x=dist
        if (modus.eq.0) then
          byout=1.0d0
        else if (modus.eq.1) then
          bxout=y*fa
          byout=x*fa
        else if (modus.eq.3) then
          x2=x*x
          x3=x2*x
          bxout=(2.0d0*fa*x+3.0d0*fb*x2)*y
          byout=fa*x2+fb*x3+y**2*(-fa-3.0d0*fb*x)
        else if (modus.eq.5) then
          x2=x*x
          x3=x2*x
          x4=x3*x
          x5=x4*x
          y2=y*y
          y3=y2*y
          bxout=y*(3.0d0*Fa*x2+4.0d0*fb*x3+5.0d0*fc*x4)
     &      +y3*(-fa-4.0d0*fb*x-10.0d0*fc*x2) !This term is not Maxwell conform
          !byout=(fa*x3+fb*x4+fc*x5)+y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2) ! The sign seems to be wrong in the manual
          byout=(fa*x3+fb*x4+fc*x5)-y2*x*(3.0d0*fa+6.0d0*fb*x+10.0d0*fc*x2)
        endif
        byout=byout*strength
        bxout=bxout*strength
      else
        byout=strength
      endif

      return
      end
