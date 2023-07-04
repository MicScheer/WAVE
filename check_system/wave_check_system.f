*CMZ :  4.01/00 14/03/2023  10.42.13  by  Michael Scheer
*CMZ :  4.00/17 15/11/2022  10.46.52  by  Michael Scheer
*CMZ :  4.00/07 28/05/2020  14.09.00  by  Michael Scheer
*CMZ :  4.00/04 05/08/2019  11.45.37  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  11.52.16  by  Michael Scheer
*CMZ :  3.07/01 21/03/2019  15.16.59  by  Michael Scheer
*CMZ :  3.05/01 04/05/2018  14.52.22  by  Michael Scheer
*CMZ :  3.05/00 26/04/2018  13.14.44  by  Michael Scheer
*CMZ :  3.03/02 29/02/2016  16.23.51  by  Michael Scheer
*CMZ :  3.03/01 05/11/2015  13.09.06  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  11.46.09  by  Michael Scheer
*CMZ :  3.01/00 17/07/2013  12.16.14  by  Michael Scheer
*CMZ :  3.00/01 20/03/2013  10.21.47  by  Michael Scheer
*CMZ :  3.00/00 14/03/2013  10.40.48  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.41.36  by  Michael Scheer
*CMZ :  2.70/00 29/11/2012  16.23.27  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  16.44.43  by  Michael Scheer
*CMZ :  2.63/02 13/03/2008  15.29.00  by  Michael Scheer
*CMZ :  2.57/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.47/23 17/02/2004  13.10.49  by  Michael Scheer
*CMZ :  2.47/09 20/05/2003  14.43.57  by  Michael Scheer
*CMZ :  2.47/08 15/05/2003  16.39.38  by  Michael Scheer
*CMZ :  2.47/03 12/03/2003  15.45.33  by  Michael Scheer
*CMZ :  2.38/00 12/12/2001  16.46.28  by  Michael Scheer
*CMZ :  2.17/00 02/11/2000  16.39.03  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  11.58.00  by  Michael Scheer
*CMZ :  2.13/04 24/01/2000  15.33.01  by  Michael Scheer
*CMZ : 00.01/10 27/08/96  16.19.56  by  Michael Scheer
*CMZ : 00.01/09 09/10/95  17.58.12  by  Michael Scheer
*CMZ : 00.00/07 19/05/94  12.06.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  12.01.00  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  17.00.43  by  Michael Scheer

*-- Author : Michael Scheer

*KEEP,GPLHINT.
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

!+seq,waveenv.

      PROGRAM wave_check_system

      !use waveenv

      implicit none

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.
c+seq,platform.

      character(1024) chrunwave

      print*,""
      print*,"--- Calling wavesystem"
      call wavesystemcheck
      print*,"--- Returned from wavesystem"

      chrunwave=trim(chwavehome) // chpathsep  // "bin" //  chpathsep // "wave.exe"

      print*,""
      print*,"--- Command to run WAVE is:"
      print*,trim(chrunwave)
      print*,""

      end program wave_check_system
