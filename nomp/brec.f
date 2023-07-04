*CMZ :  4.00/07 04/06/2020  17.17.29  by  Michael Scheer
*CMZ :  3.06/00 18/02/2019  17.12.04  by  Michael Scheer
*CMZ :  3.03/02 07/03/2016  09.53.53  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/08 01/04/95  16.52.27  by  Michael Scheer
*CMZ : 00.01/07 03/03/95  17.01.49  by  Michael Scheer
*-- Author :    Michael Scheer   03/03/95
      SUBROUTINE BREC(XIN,YIN,ZIN,BX,BY,BZ,AX,AY,AZ)
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

C     TO CALCULATED MAGNETIC FIELD OF REC-STRUCTURE WITH REC_FIELD ROUTINES

      use ompmod

      IMPLICIT NONE

      DOUBLE PRECISION XIN,YIN,ZIN,X,Y,Z,BX,BY,BZ,AX,AY,AZ
      INTEGER ICAL

*KEEP,klotz.
      include 'klotz.cmn'
*KEND.

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN
        CALL REC_INIT
        ICAL=1
      ENDIF

      X=XIN*1000.0D0
      Y=YIN*1000.0D0
      Z=ZIN*1000.0D0

      AX=0.0
      AY=0.0
      AZ=0.0

c      if (irecsolve.ne.0) call rec_solve

      if (iomp.eq.0) then
        CALL REC_BFELD(X,Y,Z,BX,BY,BZ)
      else
        CALL REC_BFELD_OMP(X,Y,Z,BX,BY,BZ)
      endif

      RETURN
      END
