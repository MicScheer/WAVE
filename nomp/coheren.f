*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.16/08 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.46.45  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.11  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.01  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE COHEREN
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IFREQ,ISOUR
      DOUBLE PRECISION OSZILL,OSZMAX,OSZMIN
      DOUBLE PRECISION DTIME,DTIMEP,DTIM

      OSZMAX=-1.D30
      OSZMIN= 1.D30

      DO IFREQ=1,NFREQ
          DO ISOUR=1,NSOURCE

          DTIME=SOURCET(2,ISOUR)-SOURCET(1,ISOUR)
          DTIMEP=(SOURCEE(1,1,ISOUR)-SOURCEA(1,1,ISOUR))/CLIGHT1

          DTIM=DTIME-DTIMEP

         OSZILL=FREQ(IFREQ)/HBAREV1*DTIM/2.D0/PI1
         IF (OSZILL.LT.OSZMIN) OSZMIN=OSZILL
         IF (OSZILL.GT.OSZMAX) OSZMAX=OSZILL

          ENDDO   !ISOUR
      ENDDO !IFREQ

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR COHEREN:'
      WRITE(LUNGFO,*)
     &'     Minimum and maximum number of phases  between'
      WRITE(LUNGFO,*)
     &'     photons and electrons:'
      WRITE(LUNGFO,*)'     ',SNGL(OSZMIN),SNGL(OSZMAX)
      WRITE(LUNGFO,*)

      RETURN
      END
