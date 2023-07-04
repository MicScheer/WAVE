*CMZ :  4.00/10 07/09/2020  15.46.05  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.56/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.36/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.44  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/03 17/12/99  10.44.11  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  20.12.40  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.46  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE FILTER
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

C--- READS ABSORPTION COEFFICIENT FROM FILE AND APPLIES FILTER TO SPECTRUM

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

      DOUBLE PRECISION AMU

      INTEGER ISOUR,IOBSV,IFREQ,ICAL

      DATA ICAL/0/

C- ABSORPTION COEFFICIENT (SPLINE INTERPOLATION)

      ICAL=ICAL+1

      DO IFREQ=1,NFREQ
        CALL ABSCOEF(FREQ(IFREQ),AMU,ABSDEN,ABSCOM,IFILTER,ICAL)
        ABSMU(IFREQ)=AMU
        EXPMU(IFREQ)=DEXP(-ABSMU(IFREQ)*ABSTHI*ABSDEN)
        IF (EXPMU(IFREQ).GT.1.D0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN FILTER ***'
          WRITE(LUNGFO,*)'REDUCTION FACTOR GREATER THAN ONE'
          WRITE(LUNGFO,*)'PHOTON ENERGY:',FREQ(IFREQ)
          WRITE(LUNGFO,*)'REDUCTION FACTOR:',ABSMU(IFREQ)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN FILTER ***'
          WRITE(6,*)'REDUCTION FACTOR GREATER THAN ONE'
          WRITE(6,*)'PHOTON ENERGY:',FREQ(IFREQ)
          WRITE(6,*)'REDUCTION FACTOR:',EXPMU(IFREQ)
          WRITE(6,*)'MASS ABSORPTION COEFFICIENT [KG/M**2]:',ABSMU(IFREQ)
          WRITE(6,*)
          WRITE(6,*)
          STOP
        ENDIF
      ENDDO !IFREQ

      DO ISOUR=1,NSOURCE
        DO IOBSV=1,NOBSV
          DO IFREQ=1,NFREQ

            ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
            SPEC(ILIOBFR)=SPEC(ILIOBFR)*EXPMU(IFREQ)

          ENDDO !IFREQ
        ENDDO !IOBSV
      ENDDO !ISOUR

      if (istokes.ne.0) then
        DO IOBSV=1,NOBSV
          DO IFREQ=1,NFREQ
            IOBFR=IOBSV+NOBSV*(IFREQ-1)
            stokes(1:4,iobfr)=stokes(1:4,iobfr)*expmu(ifreq)
          ENDDO !IFREQ
        ENDDO !IOBSV
      endif

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     *** SR FILTER CALLED, I.E. SPECTRUM IS GIVEN BEHIND ABSORBER ***'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     FILE WITH ABSORPTION COEFFICIENTS:'
      WRITE(LUNGFO,*)'     ',FILEABS
      WRITE(LUNGFO,*)'     COMMENT ON DATA FILE:'
      WRITE(LUNGFO,*)'     ',ABSCOM
      WRITE(LUNGFO,*)'     THICKNESS OF FILTER:',ABSTHI
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      RETURN
      END
