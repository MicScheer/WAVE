*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 19/12/2021  10.13.31  by  Michael Scheer
*CMZ :  3.02/04 12/12/2014  15.28.39  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  2.69/00 24/10/2012  16.55.08  by  Michael Scheer
*CMZ :  2.67/04 14/05/2012  13.05.02  by  Michael Scheer
*CMZ :  2.67/01 21/02/2012  10.22.54  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  16.23.54  by  Michael Scheer
*-- Author :    Michael Scheer   19/01/2012
      subroutine hbook2m(id,ctitle,nxbins,xmin,xmax,nybins,ymin,ymax,vmx)
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

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

c creates 2D-histo for root.

      integer id,nxbins,nybins,k,i
      real xmin,xmax,ymin,ymax,vmx
      character(*) ctitle
      character c1
      character(8) cname,cname2




        call mh_book2(id,ctitle,nxbins,dble(xmin),dble(xmax),
     &    nybins,dble(ymin),dble(ymax))



      return
      end
