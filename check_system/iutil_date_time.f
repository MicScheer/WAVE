*CMZ :  4.00/15 02/06/2022  09.35.26  by  Michael Scheer
*CMZ :  3.01/02 25/02/2014  14.52.02  by  Michael Scheer
*CMZ :  2.70/06 07/01/2013  13.42.44  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  13.41.59  by  Michael Scheer
*CMZ :  2.66/20 18/10/2011  08.06.47  by  Michael Scheer
*CMZ :  2.66/18 01/12/2010  10.11.23  by  Michael Scheer
*CMZ :  2.66/13 23/11/2010  09.59.43  by  Michael Scheer
*CMZ :  2.63/00 10/01/2008  12.40.35  by  Michael Scheer
*CMZ :  2.58/01 17/01/2007  10.59.07  by  Michael Scheer
*CMZ :  2.57/04 13/01/2006  11.10.13  by  Michael Scheer
*CMZ :  2.57/03 09/12/2005  11.19.05  by  Michael Scheer
*CMZ :  2.52/15 03/01/2005  15.51.05  by  Michael Scheer
*CMZ :  2.47/22 03/12/2003  10.29.33  by  Michael Scheer
*CMZ :  2.46/02 21/01/2003  16.21.35  by  Michael Scheer
*CMZ :  2.41/09 14/08/2002  17.26.05  by  Michael Scheer
*CMZ :  2.39/00 03/01/2002  12.30.43  by  Michael Scheer
*CMZ :  2.20/01 02/01/2001  11.38.55  by  Michael Scheer
*CMZ :  2.15/00 01/05/2000  11.48.08  by  Michael Scheer
*CMZ :  2.13/04 24/01/2000  17.57.30  by  Michael Scheer
*CMZ :  2.13/00 06/10/99  16.32.34  by  Michael Scheer
*CMZ :  1.03/06 01/07/98  10.29.18  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.57.32  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.41  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine iutil_date_time

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
      ! calculates date and time as numbers and corresponding string
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

      ! real: util_time
      ! integer i_util_time, i_util_date,idatetime(8)
      ! integer*8 i_util_date_time
      ! string: ch_i_util_date_time, dtday,dttime,dtzone

      integer i,k
      character(6) c6
      call date_and_time(dtday,dttime,dtzone,idatetime)

      read(dttime,*) util_time

      i_util_time=nint(util_time)
      write(c6,'(i6)') i_util_time

      ch_i_util_date_time=trim(dtday) // c6
      read(ch_i_util_date_time,'(i6)') i_util_date
      read(ch_i_util_date_time,*) i_util_date_time
      k=0
      do i=1,len_trim(ch_i_util_date_time)
        if (ch_i_util_date_time(i:i).ne.' ') then
          k=k+1
          ch_i_util_date_time(k:k)=ch_i_util_date_time(i:i)
        endif
      enddo

      return
      end
