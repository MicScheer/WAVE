*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  13.48.25  by  Michael Scheer
*CMZ : 00.01/03 18/11/94  11.39.59  by  Michael Scheer
*CMZ :  0.00/03 09/11/94  18.29.44  by  Michael Scheer
*CMZ :  0.00/02 02/11/94  12.17.57  by  Michael Scheer
*-- Author :    Michael Scheer   31/10/94
      SUBROUTINE BRCOEF(LUN,FILE
     &                 ,XLX,YLY,ZLZ
     &                 ,XKX,YKY,ZKZ
     &                 ,NFIRSTX,NORDX,NSTEPX
     &                 ,NFIRSTY,NORDY,NSTEPY
     &                 ,NDIMQ,Q,QA0,QA,IFHALBA,GAP2PI,WIDTH)

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

      CHARACTER(50) FILE

      INTEGER LUN,IORD,JORD,IORD1,JORD1,NDIMQ
     &                 ,NFIRSTX,NORDX,NSTEPX
     &                 ,NFIRSTY,NORDY,NSTEPY
     &                 ,I,J
     &                 ,IFHALBA

      INTEGER IORDX,IORDY,IDUMX,IDUMY

      DIMENSION Q(NDIMQ,NDIMQ)
      DOUBLE PRECISION XLX,YLY,ZLZ,XKX,YKY,ZKZ,Q
      DOUBLE PRECISION GAP2PI,WIDTH
      DOUBLE PRECISION QA0(NDIMQ),QA(NDIMQ,NDIMQ)

      DO I=1,NDIMQ
      DO J=2,NDIMQ
             Q(I,J)=0.
      ENDDO
      ENDDO


      OPEN(UNIT=LUN,FILE=FILE,STATUS='OLD')


         READ(LUN,*)IFHALBA
         READ(LUN,*)XLX
         READ(LUN,*)YLY
         READ(LUN,*)ZLZ
         READ(LUN,*)XKX
         READ(LUN,*)YKY
         READ(LUN,*)ZKZ
         READ(LUN,*)GAP2PI
         READ(LUN,*)WIDTH
         READ(LUN,*)NFIRSTX,NORDX,NSTEPX
         READ(LUN,*)NFIRSTY,NORDY,NSTEPY

         IF (IFHALBA.EQ.0) THEN
           DO IORD=NFIRSTX,NORDX,NSTEPX
           DO JORD=NFIRSTY,NORDY,NSTEPY
               READ(LUN,*)IORD1,JORD1,Q(IORD1+1,JORD1+1)
           ENDDO
           ENDDO
           NFIRSTX=NFIRSTX+1    !FORTRAN-INDICES
           NORDX=NORDX+1
           NFIRSTY=NFIRSTY+1
           NORDY=NORDY+1
         ELSE !IFHALBA
           IORD=0
           DO IORDY=NFIRSTY,NORDY,NSTEPY
               iord=iord+1
               READ(LUN,*)IDUMY,QA0(IDUMY)
           DO IORDX=NFIRSTX,NORDX,NSTEPX
               iord=iord+1
               READ(LUN,*)IDUMY,IDUMX,QA(IDUMY,IDUMX)
           ENDDO
           ENDDO
         ENDIF !IFHALBA
      CLOSE(LUN)

      RETURN
      END
