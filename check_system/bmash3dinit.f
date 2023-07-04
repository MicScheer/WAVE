*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.21.18  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  18.22.03  by  Michael Scheer
*CMZ :  1.00/00 28/07/97  15.09.21  by  Michael Scheer
*CMZ : 00.02/00 14/11/96  17.55.17  by  Michael Scheer
*CMZ : 00.01/09 05/10/95  16.34.35  by  Michael Scheer
*CMZ :  1.00/01 25/09/95  14.26.06  by  Michael Scheer
*CMZ :  1.00/00 21/09/95  17.14.08  by  Michael Scheer
*-- Author :    Michael Scheer   21/09/95
      SUBROUTINE BMASH3DINIT

      IMPLICIT NONE
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

*KEEP,bmash.
      include 'bmash.cmn'
*KEND.

      CHARACTER(80) FILECOEFF

      INTEGER NCOEF,IORD

      INTEGER IX,IY,IZ,IS,NS,IR

      DATA FILECOEFF/'WI:BMASH.COEF'/

      MORD=MORDP
      MFIT=MFITP
      MTOT=MTOTP
      NPOI=NPOIP

      NCOEF=0

        OPEN(UNIT=99,FILE=FILECOEFF,STATUS='OLD')


100       READ(99,*,END=900) IX,IY,IZ,NS

          DO IORD=1,MORDP
         IF(IX+IY+IZ.EQ.IORD) GOTO 50
          ENDDO
          GOTO 100

50        NCOEF=NCOEF+1

          IF (NCOEF.GT.MTOTP) THEN
             WRITE(6,*)
             WRITE(6,*)
     &   '*** ERROR IN BMASH3DINIT: DIMENSION MTOTP EXCEEDED ***'
             WRITE(6,*)
             STOP
          ENDIF

          IF (NS.GT.NSTAKP) THEN
             WRITE(6,*)
             WRITE(6,*)
     &   '*** ERROR IN BMASH3DINIT: DIMENSION NSTAKP EXCEEDED ***'
             WRITE(6,*)
             STOP
          ENDIF

          INDEX(1,NCOEF)=IX
          INDEX(2,NCOEF)=IY
          INDEX(3,NCOEF)=IZ

          IF (NS.GT.0) THEN
         INDEX(4,NCOEF)=NS
         DO IS=1,NS
             READ(99,*)(FSTAK(IR,IS,NCOEF),IR=1,4)
         ENDDO
          ELSE
         INDEX(4,NCOEF)=1
         FSTAK(1,1,NCOEF)=IX
         FSTAK(2,1,NCOEF)=IY
         FSTAK(3,1,NCOEF)=IZ
         FSTAK(4,1,NCOEF)=1.D0
          ENDIF   !IS

      GOTO 100

900   CLOSE(99)

      IF (NCOEF.NE.MFITP) THEN
         WRITE(6,*)
     &'*** ERROR IN BMASH3DINIT: WRONG NUMBER OF COEFFICIENTS ***'
         WRITE(6,*)'CHECK ',FILECOEFF
         WRITE(6,*)
         STOP
      ENDIF

      RETURN
      END
