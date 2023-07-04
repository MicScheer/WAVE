*CMZ :  4.00/04 17/05/2019  11.51.09  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.48/04 16/03/2004  10.43.05  by  Michael Scheer
*CMZ :  2.48/03 03/03/2004  12.49.39  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  18.04.15  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.08.41  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  10.26.40  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.28.07  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_WARTE(IBATCH)
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

      INTEGER IBATCH
      CHARACTER(5) ANS

      WRITE(6,*)
      WRITE(6,*)'PAUSE CAUSED BY ROUTINE UTIL_WARTE'
      WRITE(6,*)
      READ(5,'(1A5)',ERR=99) ANS

      IF (IBATCH.NE.0) RETURN


      IF (ANS.EQ.'BLITZ'.OR.ANS.EQ.'blitz') THEN
         call sleep(1)
      ELSE IF (ANS.EQ.'SHORT'.OR.ANS.EQ.'short') THEN
         call sleep(15)
      ELSE IF (ANS.EQ.'LONG'.OR.ANS.EQ.'long') THEN
        call sleep(300)
      ELSE IF (ANS.EQ.'PAUSE'.OR.ANS.EQ.'pause') THEN
        call sleep(3600)
      ENDIF


99    RETURN
      END
