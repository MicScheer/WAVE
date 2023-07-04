*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.56/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.41/10 29/04/2004  15.29.30  by  Michael Scheer
*CMZ :  2.36/00 07/11/2001  14.23.12  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.12/03 07/07/99  16.16.04  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.36.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.17  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.49  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ABSCOEF(FREQ,ABSMU,ABSDEN,ABSCOM,IFILTER,ICAL)
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

C--- READS ABSORPTION COEFFICIENT FROM FILE CALCULATES COEFFICIENT FOR
C    GIVEN FREQUENCY

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      INTEGER I,ICAL,ICALO,NCOEF,NDAFRQP,IERR,IFILTER

      PARAMETER(NDAFRQP=10000)

      CHARACTER(65) ABSCOM
      DOUBLE PRECISION ABSDEN,AFREQ(NDAFRQP),ABSCO(NDAFRQP),FREQ,ABSMU

      DATA ICALO/-1/

C--- READ DATA FILE

      IF (ICAL.NE.ICALO) THEN

        OPEN(UNIT=LUNABS,FILE=FILEABS,STATUS='OLD',FORM='FORMATTED')

        READ(LUNABS,'(A65)') ABSCOM
        READ(LUNABS,*) ABSDEN
        READ(LUNABS,*) NCOEF

        IF (NCOEF.GT.NDAFRQP) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN ABSCOEF ***'
          WRITE(LUNGFO,*)'DIMENSION EXCEEDED'
          WRITE(LUNGFO,*)'INCREASE PARAMETER NDAFRQP IN THIS ROUTINE'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN ABSCOEF ***'
          WRITE(6,*)'DIMENSION EXCEEDED'
          WRITE(6,*)'INCREASE PARAMETER NDAFRQP IN THIS ROUTINE'
          WRITE(6,*)
          STOP

        ENDIF !NCOEF

        DO I=1,NCOEF
          READ(LUNABS,*) AFREQ(I),ABSCO(I)
        ENDDO !I

        CLOSE(LUNABS)

        ICALO=ICAL

      ENDIF !ICAL

      IF (IFILTER.EQ.1) THEN
        CALL ABSNOSPLI(AFREQ,ABSCO,NCOEF,FREQ,ABSMU,IERR,1)
      ELSE
        CALL ABSNOSPLI(AFREQ,ABSCO,NCOEF,FREQ,ABSMU,IERR,-1)
      ENDIF

      IF (IERR.NE.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN ABSCOEF ***'
        WRITE(LUNGFO,*)'CALL TO SR ABSNOSPLI FAILED'
        WRITE(LUNGFO,*)'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE'
        WRITE(LUNGFO,*)FILEABS
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)

        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN ABSCOEF ***'
        WRITE(6,*)'CALL TO SR ABSNOSPLI FAILED'
        WRITE(6,*)'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE'
        WRITE(6,*)FILEABS
        WRITE(6,*)
        WRITE(6,*)

        STOP

      ENDIF   !IERR

      RETURN
      END
