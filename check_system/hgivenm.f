*CMZ :  3.02/00 11/10/2014  21.21.59  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  11.16.25  by  Michael Scheer
*CMZ :  0.01/03 11/06/2014  11.35.22  by  Michael Scheer
*CMZ :  0.01/02 06/06/2014  16.22.45  by  Michael Scheer
*CMZ :  0.01/01 01/06/2014  13.01.56  by  Michael Scheer
*CMZ :  0.01/00 29/04/2014  12.02.22  by  Michael Scheer
*CMZ :  0.00/02 29/04/2014  09.10.41  by  Michael Scheer
*CMZ :  0.00/01 28/04/2014  20.38.18  by  Michael Scheer
*-- Author :    Michael Scheer   26/04/2014
      subroutine hgivenm(id,title,nvar,chtag,rlow,rhigh)
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

*KEEP,ntupinfo.
      include 'ntupinfo.cmn'
*KEND.

      real rlow(nvar),rhigh(nvar)
      real*8 rlow8(nvar),rhigh8(nvar)
      character(*) title
      character(*) chtag(1)
      character(1024) cline
      integer nvar,id,i
      integer*8 nvar8

      return
      end
