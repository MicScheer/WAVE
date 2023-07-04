*CMZ :  3.02/03 03/11/2014  10.42.03  by  Michael Scheer
*CMZ :  2.68/05 02/10/2012  09.01.01  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  11.58.18  by  Michael Scheer
*CMZ :  2.66/12 24/06/2010  12.50.52  by  Michael Scheer
*CMZ :  2.66/08 12/03/2010  15.02.37  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  14.08.46  by  Michael Scheer
*-- Author :    Michael Scheer   24/02/2010
      subroutine ubunch(xbunch,ybunch,zbunch,ypbunch,zpbunch,gambunch,
     &  dtphase)
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

      use bunchmod

      implicit none

      double precision dtphase,xbunch,ybunch,zbunch,
     &  ypbunch,zpbunch,gambunch

      integer ical

      data ical/0/

      if (ical.eq.0) then

        xbunch=xbunch
        ybunch=ybunch
        zbunch=zbunch
        ypbunch=ypbunch
        zpbunch=zpbunch
        gambunch=gambunch
        dtphase=dtphase

        ical=1
      endif !ical

      return
      end

