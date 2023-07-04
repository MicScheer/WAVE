*CMZ :  3.06/00 15/02/2019  14.44.39  by  Michael Scheer
*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  3.02/00 09/10/2014  14.11.09  by  Michael Scheer
*CMZ :  2.52/05 16/08/2004  15.08.37  by  Michael Scheer
*-- Author :    Michael Scheer   16/08/2004
      subroutine bpolypl2(xpl,ypl,col,ixyz)
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

*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      implicit none

*KEEP,bforce.
      include 'bforce.cmn'
*KEND.

      real xpl(2), ypl(2), x(2), y(2), col,
     &  rmtyp20,rmtyp24,rmtyp31

      integer ixyz

      data rmtyp20/20./
      data rmtyp24/-9999./
      data rmtyp31/31./

      if (nbforcx*nbforcy*nbforcz.eq.0) return

      call muwk(0,0)

      call mgset('PLCI',col)
      call mgset('PMCI',col)
      call mgset('MTYP',rmtyp24)
      call mgset('MSCF',0.1)

      if (ixyz.eq.12) then
        x(1)=torqcenxmm
        y(1)=torqcenymm
      else if (ixyz.eq.13) then
        x(1)=torqcenxmm
        y(1)=torqcenzmm
      else if (ixyz.eq.23) then
        x(1)=torqcenzmm
        y(1)=torqcenymm
      endif
      call mpm(1,x,y)
      call mgset('PMCI',1.)

      x(1)=xpl(1)
      x(2)=xpl(2)
      y(1)=ypl(1)
      y(2)=ypl(1)
      call mpl(2,x,y)

      x(1)=xpl(1)
      x(2)=xpl(1)
      y(1)=ypl(1)
      y(2)=ypl(2)
      call mpl(2,x,y)

      x(1)=xpl(1)
      x(2)=xpl(2)
      y(1)=ypl(2)
      y(2)=ypl(2)
      call mpl(2,x,y)

      x(1)=xpl(2)
      x(2)=xpl(2)
      y(1)=ypl(1)
      y(2)=ypl(2)
      call mpl(2,x,y)

      call muwk(0,0)

      return
      end
