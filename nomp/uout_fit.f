*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.13/03 17/12/99  11.45.46  by  Michael Scheer
*CMZ :  1.00/00 10/07/97  13.57.51  by  Michael Scheer
*CMZ : 00.01/10 02/06/96  12.04.24  by  Michael Scheer
*CMZ : 00.01/06 20/02/95  16.06.00  by  Michael Scheer
*CMZ : 00.01/05 01/02/95  15.24.06  by  Michael Scheer
*CMZ : 00.01/04 26/01/95  15.55.24  by  Michael Scheer
*CMZ : 00.00/07 24/05/94  09.48.34  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UOUT_FIT
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.
      OPEN(UNIT=99,FILE='UOUT.FIT',STATUS='OLD')
         WRITE(99,*)SPEC(1)
         WRITE(99,*)WTRA2I
         WRITE(99,*)WTRA(2,1,NCO)
         WRITE(99,*)WTRA(3,1,NCO)
         WRITE(99,*)WTRA(2,2,NCO)/WTRA(1,2,NCO)
         WRITE(99,*)WTRA(3,2,NCO)/WTRA(1,2,NCO)
      CLOSE(99)

      RETURN
      END
