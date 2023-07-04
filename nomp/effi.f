*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.36/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.46  by  Michael Scheer
*CMZ :  2.16/01 15/06/2000  16.16.24  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/03 17/12/99  10.44.11  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  20.12.40  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.46  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE EFFI
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

C--- READS PHOTO YIELD FROM FILE AND APPLIES TO SPECTRUM

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER ISOUR,IOBSV,IFREQ

      DOUBLE PRECISION AEFF

C- ABSORPTION COEFFICIENT (SPLINE INTERPOLATION)

      DO IFREQ=1,NFREQ
          CALL YIELD(IEFFI,FREQ(IFREQ),AEFF,EFFCOM)
          EFF(IFREQ)=AEFF
      ENDDO !IFREQ

      DO ISOUR=1,NSOURCE
      DO IOBSV=1,NOBSV
      DO IFREQ=1,NFREQ

          ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
            SPEC(ILIOBFR)=
     &      SPEC(ILIOBFR)*EFF(IFREQ)

      ENDDO !IFREQ
      ENDDO !IOBSV
      ENDDO !ISOUR

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &   '     *** SR EFFI CALLED, I.E. SPECTRUM IS GIVEN AS SEEN BY DETECTOR ***'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     FILE WITH PHOTO YIELD:'
      WRITE(LUNGFO,*)'     ',FILEFF
      WRITE(LUNGFO,*)'     COMMENT ON DATA FILE:'
      WRITE(LUNGFO,*)'     ',EFFCOM
      WRITE(LUNGFO,*)

      RETURN
      END
