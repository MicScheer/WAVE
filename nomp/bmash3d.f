*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  22.21.26  by  Michael Scheer
*CMZ :  1.00/00 29/07/97  10.15.30  by  Michael Scheer
*CMZ : 00.02/00 21/11/96  14.40.10  by  Michael Scheer
*CMZ :  1.00/03 27/09/95  17.15.24  by  Michael Scheer
*CMZ :  1.00/01 25/09/95  18.20.26  by  Michael Scheer
*CMZ :  1.00/00 21/09/95  17.13.57  by  Michael Scheer
*-- Author :    Michael Scheer   21/09/95
      SUBROUTINE BMASH3D(X,Y,Z,BX,BY,BZ,IFAIL)
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

C---  TO FIT 3D POTENTIAL V=SUM( CG(I,J,K) * X**(I-1) * Y**(J-1) * Z**(K-1))
C     OF A MAGNETIC FIELD B=(BX,BY,BZ)=-GRAD(V)
C

C--- OUTPUT:

C     IFAIL : FAILURE FLAG

      IMPLICIT NONE


*KEEP,bmash.
      include 'bmash.cmn'
*KEND.

      INTEGER IFAIL,I,NCINDG

      DOUBLE PRECISION X(NPOIP),Y(NPOIP),Z(NPOIP)
     &                  ,BX(NPOIP),BY(NPOIP),BZ(NPOIP)

      INTEGER ICAL,ICOEF,IX,IY,IZ,IS,NS
     &         ,JX,JY,JZ
      INTEGER ICG(MORDP+1,MORDP+1,MORDP+1)

      DOUBLE PRECISION VC(MFITP)

      DATA ICAL/0/

      DO IZ=1,MORDP+1
      DO IY=1,MORDP+1
      DO IX=1,MORDP+1
          CG(IX,IY,IZ)=0.D0
          ICG(IX,IY,IZ)=0
      ENDDO
      ENDDO
      ENDDO

      IF (ICAL.EQ.0) THEN
          CALL BMASH3DINIT
      ENDIF !ICAL

      CALL BMASH3DFIT(X,Y,Z,BX,BY,BZ,VC,IFAIL)

C--- GET FITTED COEFFICIENTS

      DO ICOEF=1,MFITP
          IX=INDEX(1,ICOEF)
          IY=INDEX(2,ICOEF)
          IZ=INDEX(3,ICOEF)
C         IF (ICAL.EQ.0) READ(5,*)VC(ICOEF)
          CG(IX+1,IY+1,IZ+1)=VC(ICOEF)
          ICG(IX+1,IY+1,IZ+1)=1
      ENDDO

C--- CALCULATE OTHER COEFFICIENTS

      DO ICOEF=1,MFITP
            NS=INDEX(4,ICOEF)
          IF (NS.GT.1) THEN
             IX=INDEX(1,ICOEF)
             IY=INDEX(2,ICOEF)
             IZ=INDEX(3,ICOEF)
             DO IS=1,NS
            JX=NINT(FSTAK(1,IS,ICOEF))
            JY=NINT(FSTAK(2,IS,ICOEF))
            JZ=NINT(FSTAK(3,IS,ICOEF))
            IF (IX.NE.JX.OR.IY.NE.JY.OR.IZ.NE.JZ) THEN
            ICG(JX+1,JY+1,JZ+1)=2
            CG(JX+1,JY+1,JZ+1)=CG(JX+1,JY+1,JZ+1)
     &         +FSTAK(4,IS,ICOEF)*CG(IX+1,IY+1,IZ+1)
                 ENDIF
             ENDDO   !IS
         ENDIF
      ENDDO

      IF (ICAL.EQ.0) THEN

         NCINDG=0
              DO IZ=1,MORDP+1
              DO IY=1,MORDP+1
              DO IX=1,MORDP+1
                  IF (ICG(IX,IY,IZ).NE.0) THEN
            NCINDG=NCINDG+1
            ICINDG(1,NCINDG)=IX
            ICINDG(2,NCINDG)=IY
            ICINDG(3,NCINDG)=IZ
             ENDIF
              ENDDO
              ENDDO
              ENDDO

         IF (NCINDG.NE.MTOTP) THEN
             WRITE(6,*)
             WRITE(6,*)
     & '*** ERROR IN BMASH3D: WRONG NUMBER OF COEFFICIENTS'
             WRITE(6,*)
             STOP
         ENDIF

         ICAL=1

      ENDIF

      DO I=1,NCINDG
          CINDG(I)=CG(ICINDG(1,I),ICINDG(2,I),ICINDG(3,I))
      ENDDO

      RETURN
      END
