*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.56.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.20  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.34  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CONVUN
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

C--- CONVERTS FREQUENCES ON INPUT FILE FROM EV TO NM OR VICE VERSA

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IFREQ
      DOUBLE PRECISION FBUFF(NDFREQP)

      DO IFREQ=1,NFREQ
          FBUFF(IFREQ)=WTOE1/FREQ(IFREQ)
      ENDDO

      DO IFREQ=1,NFREQ
             FREQ(IFREQ)=FBUFF(NFREQ+1-IFREQ)
      ENDDO


      IF (FREQC.NE.0.0) FREQC=WTOE1/FREQC
      IF (FREQCF.NE.0.0) FREQCF=WTOE1/FREQCF

      RETURN
      END
