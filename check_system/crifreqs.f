*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.16/08 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/03 16/12/99  10.51.00  by  Michael Scheer
*CMZ :  2.12/04 20/08/99  11.22.54  by  Michael Scheer
*CMZ :  2.12/03 29/07/99  14.01.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.24  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.50  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CRIFREQS(IOBSV)
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C     INTEGRATES SPECTRAL SINGLE OBSERVATION POINT
C     OVER ALL PHOTON ENERGIES TO CALCULATE CRITICAL ENERGY

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      INTEGER IFREQ,ISOUR,IOBSV
      DOUBLE PRECISION FBUFF(NDFREQP),DFBUFF(NDFREQP),SBUFF(NDFREQP),DSUM,SUM,SUM2

      IF(NFREQ.LT.2) RETURN

      DO IFREQ=1,NFREQ
          FBUFF(IFREQ)=0.0
      DO ISOUR=1,NSOURCE
          FBUFF(IFREQ)=FBUFF(IFREQ)
     &                +SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
     &                /BANWID*ECHARGE1
      ENDDO !ISOUR
      ENDDO !IFREQ

      DO IFREQ=2,NFREQ
          DFBUFF(IFREQ)=FREQ(IFREQ)-FREQ(IFREQ-1)
      ENDDO !IFREQ

      SUM=0.0
      SBUFF(1)=0.0
      DO IFREQ=2,NFREQ
          DSUM=(FBUFF(IFREQ)+FBUFF(IFREQ-1))/2.*DFBUFF(IFREQ)
          SUM=SUM+DSUM
          SBUFF(IFREQ)=SUM
      ENDDO !IFREQ

      SUM2=SUM/2.
      DO IFREQ=2,NFREQ
          IF (SBUFF(IFREQ).GT.SUM2) GOTO 100
      ENDDO !IFREQ

100   CONTINUE

      IF (SBUFF(IFREQ)-SBUFF(IFREQ-1).NE.0.) THEN
           FREQC=FREQ(IFREQ-1)
     &       +(SUM2-SBUFF(IFREQ-1))/(SBUFF(IFREQ)-SBUFF(IFREQ-1))
     &       *DFBUFF(IFREQ)
      ELSE
           FREQC=0.0
      ENDIF

      IF (FREQC.EQ.0.0) FREQC=FREQ(NFREQ/2+1)

      IF (FREQC.NE.0.0) WELLENC=WTOE1/FREQC

      RETURN
      END
