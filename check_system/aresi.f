*CMZ :  3.00/00 11/03/2013  10.37.07  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.08.30  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.39  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.54  by  Michael Scheer
*-- Author :
C****************************************************************
      SUBROUTINE ARESI(AIN,MAXTRAP,NORDNGP,
     &                   RESXM,RESXAVE,RESX,
     &                   RESYM,RESYAVE,RESY,
     &                   RESPXM,RESPXAVE,RESPX,
     &                   RESPYM,RESPYAVE,RESPY,RESALL)
C****************************************************************

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

C--- RECALCULATES TRACK USING THE GENERATING FUNCTION AND COMPARES WITH
C    TRACKING RESULTS

      IMPLICIT NONE

      INTEGER MAXTRAP,NORDNGP,NKOEF,I,J,K,L,M

*KEEP,genfun.
      include 'genfun.cmn'
*KEND.

      DOUBLE PRECISION AIN(NORDNG,NORDNG,NORDNG,NORDNG)
      DOUBLE PRECISION RESX2,RESY2,RESPX2,RESPY2,RESALL2

      DOUBLE PRECISION
     &                   RESXM,RESXAVE,RESX,
     &                   RESYM,RESYAVE,RESY,
     &                   RESPXM,RESPXAVE,RESPX,
     &                   RESPYM,RESPYAVE,RESPY,RESALL

*KEEP,ttracks.
      include 'ttracks.cmn'
*KEND.

      IF(MAXTRA.NE.MAXTRAP) STOP '*** SR ARESI: MAXTRA FALSCH ***'
      IF(NORDNGP.NE.NORDNG) STOP '*** SR ARESI: NORDNG FALSCH ***'

      RESXAVE=0.D0
      RESX2=0.D0
      RESXM=0.D0
      RESPXAVE=0.D0
      RESPX2=0.D0
      RESPXM=0.D0
      RESYAVE=0.D0
      RESY2=0.D0
      RESYM=0.D0
      RESPYAVE=0.D0
      RESPY2=0.D0
      RESPYM=0.D0
      RESALL2=0.D0

      DO L=1,NORDNG
        DO K=1,NORDNG
          DO J=1,NORDNG
            DO I=1,NORDNG

              ADUM(I,J,K,L)=AIN(I,J,K,L)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO M=1,MTRAJ

          CALL ERZTRA    (XIC(M),PXF(M),YIC(M),PYF(M),
     &           XFN(M),PXIN(M),YFN(M),PYIN(M))

          RESXF(M)=XFN(M)-XFC(M)
          RESYF(M)=YFN(M)-YFC(M)
          RESPXI(M)=PXIN(M)-PXI(M)
          RESPYI(M)=PYIN(M)-PYI(M)

          RESALL2=RESALL2+RESXF(M)**2+RESYF(M)**2+RESPXI(M)**2+RESPYI(M)**2

          RESXAVE=RESXAVE+RESXF(M)
          RESX2=RESX2+RESXF(M)*RESXF(M)
          IF(DABS(RESXF(M)).GT.DABS(RESXM)) RESXM=RESXF(M)

          RESPXAVE=RESPXAVE+RESPXI(M)
          RESPX2=RESPX2+RESPXI(M)*RESPXI(M)
          IF(DABS(RESPXI(M)).GT.DABS(RESPXM)) RESPXM=RESPXI(M)

          RESYAVE=RESYAVE+RESYF(M)
          RESY2=RESY2+RESYF(M)*RESYF(M)
          IF(DABS(RESYF(M)).GT.RESYM) RESYM=RESYF(M)

          RESPYAVE=RESPYAVE+RESPYI(M)
          RESPY2=RESPY2+RESPYI(M)*RESPYI(M)
          IF(DABS(RESPYI(M)).GT.RESPYM) RESPYM=RESPYI(M)

      ENDDO

      RESXAVE=RESXAVE/MTRAJ
      RESX=DSQRT(RESX2/MTRAJ)
      RESPXAVE=RESPXAVE/MTRAJ
      RESPX=DSQRT(RESPX2/MTRAJ)

      RESYAVE=RESYAVE/MTRAJ
      RESY=DSQRT(RESY2/MTRAJ)
      RESPYAVE=RESPYAVE/MTRAJ
      RESPY=DSQRT(RESPY2/MTRAJ)

      RESALL=DSQRT(RESALL2/MTRAJ)

      RETURN
      END
