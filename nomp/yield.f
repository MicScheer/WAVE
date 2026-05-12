*CMZ :          12/05/2026  08.10.39  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/03 02/05/2008  14.41.00  by  Michael Scheer
*CMZ :  2.52/13 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.41/10 29/04/2004  15.29.30  by  Michael Scheer
*CMZ :  2.36/00 07/11/2001  14.28.25  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.16/04 24/06/2000  21.14.53  by  Michael Scheer
*CMZ :  2.16/01 15/06/2000  17.07.13  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.12/03 07/07/99  16.16.04  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.36.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.17  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.49  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE YIELD(IEFFI,FREQ,ABSMU,ABSCOM)
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

C--- READS YIELD COEFFICIENTS FROM FILE CALCULATES COEFFICIENT FOR
C    GIVEN FREQUENCY

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      DOUBLE PRECISION, dimension(:), allocatable ::  AFREQ,ABSCO,Y2,AA,BB,CC,C

      INTEGER :: I,NCOEF,NDAFRQP,IERR,MODE,IEFFI,ical=0,icomment,istat

      CHARACTER(65) ABSCOM, cline

      DOUBLE PRECISION FREQ,ABSMU,x

      SAVE

C--- READ DATA FILE

      IERR=0

      IF (ICAL.EQ.0) THEN

        OPEN(UNIT=LUNEFF,FILE=FILEFF,STATUS='OLD',FORM='FORMATTED')

        READ(LUNABS,'(A65)') ABSCOM
        READ(LUNABS,'(A65)') cline

        icomment=1

        do i=1,65
          if (cline(i:i).eq.'.') then
            icomment=0
            exit
          endif
        enddo

        if (icomment.eq.1) then
          backspace(lunabs)
          READ(LUNABS,*) NCOEF
        else
          abscom='No Commment on file ' // trim(fileff)
          rewind(lunabs)
          ncoef=0
          do while(.true.)
            read(lunabs,*,iostat=istat) x
            if (istat.ne.0) exit
            ncoef=ncoef+1
          enddo
          rewind(lunabs)
        endif

        allocate(AFREQ(ncoef),ABSCO(ncoef),
     &  Y2(ncoef),AA(ncoef),BB(ncoef),CC(ncoef),C(ncoef))

        DO I=1,NCOEF
          READ(LUNABS,*) AFREQ(I),ABSCO(I)
        ENDDO !I

        CLOSE(LUNABS)

        IF (IEFFI.GT.0) THEN
          CALL UTIL_SPLINE_COEF(AFREQ,ABSCO,NCOEF,-9999.0d0,-9999.0d0,
     &      Y2,AA,BB,CC,C)
        ENDIF   !IEFFI

      ENDIF !ICAL

      IF (IEFFI.GT.0) THEN

        IF (ICAL.EQ.0) THEN
          MODE=-1
        ELSE
          MODE=0
        ENDIF

        CALL UTIL_SPLINE_INTER(AFREQ,ABSCO,Y2,NCOEF,FREQ,ABSMU,MODE)

        IF (IERR.NE.0) THEN

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN YIELD ***'
          WRITE(LUNGFO,*)'CALL TO SR UTIL_SPLINE_INTER FAILED'
          WRITE(LUNGFO,*)'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE'
          WRITE(LUNGFO,*)FILEFF
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)

          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN YIELD ***'
          WRITE(6,*)'CALL TO SR UTIL_SPLINE_INTER FAILED'
          WRITE(6,*)'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE'
          WRITE(6,*)FILEFF
          WRITE(6,*)
          WRITE(6,*)

          STOP

        ENDIF   !IERR

      ELSE    !IEFFI

        IF (IEFFI.EQ.-1) THEN  !IEFFI
          CALL ABSNOSPLI(AFREQ,ABSCO,NCOEF,FREQ,ABSMU,IERR,1)
        ELSE   !IEFFI
          CALL ABSNOSPLI(AFREQ,ABSCO,NCOEF,FREQ,ABSMU,IERR,-1)
        ENDIF      !IEFFI

        IF (IERR.NE.0) THEN

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN YIELD ***'
          WRITE(LUNGFO,*)'CALL TO ABSNOSPLI FAILED'
          WRITE(LUNGFO,*)'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE'
          WRITE(LUNGFO,*)FILEFF
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)

          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN YIELD ***'
          WRITE(6,*)'CALL TO ABSNOSPLI FAILED'
          WRITE(6,*)'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE'
          WRITE(6,*)FILEFF
          WRITE(6,*)
          WRITE(6,*)

          STOP

        ENDIF   !IERR

      ENDIF   !IEFFI

      ICAL=1

      RETURN
      END
