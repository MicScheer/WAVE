*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.53/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.52/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.28  by  Michael Scheer
*CMZ :  1.01/02 12/12/97  10.25.25  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.02  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECINTF
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
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C     INTEGRATES POWER SPECTRA OVER ALL PHOTON ENERGIES FOR FOLDED INTENSITY

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
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEND.

      INTEGER IY,IZ,ISOUR,IOBSV,IFREQ,IWARNW,IWARNS
      DOUBLE PRECISION S2(NDFREQP),RESULT,SIMPLE

      IWARNS=0
      IWARNW=0

      IF (NFREQ0.LE.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING SR SPECINTF ***'
        WRITE(LUNGFO,*)'INTEGRATION NOT POSSIBLE, ONLY ONE FREQUENCY'
        WRITE(LUNGFO,*)'CHANGE PARAMETER NFREQ2P IN NAMELIST CONTRL'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING SR SPECINTF ***'
        WRITE(6,*)'INTEGRATION NOT POSSIBLE, ONLY ONE FREQUENCY'
        WRITE(6,*)'CHANGE PARAMETER NFREQ2P IN NAMELIST CONTRL'
        WRITE(6,*)
        RETURN
      ENDIF

C--- LOOP OVER ALL SOURCES

      WFLUXTIF=0.0
      DO IOBSV=1,NOBSV
        SPECTOTIF(IOBSV)=0.0
      ENDDO !IOBSV

      DO ISOUR=1,NSOURCE

C--- LOOP OVER ALL OBSERVATION POINTS

        IF (IPIN.NE.0) THEN

          DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY
            DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ

              IOBSV=(IY-1)*NOBSVZ+IZ

C--- FILL INTEGRATION BUFFER

              DO IFREQ=NFREQ0M,NFREQ0P
                S2(IFREQ-NFREQ0M+1)=
     &            SPECF(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
     &            /BANWID*ECHARGE1
              ENDDO  !IFREQ

C--- DO INTEGRATION OF BUFFER

              CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)

              IF (IWARNW.EQ.0.AND.RESULT.NE.0.0
     &        .AND.DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
     &'** WARNING SPECINTF: PROBLEM WITH INTEGRATION OF SPECTUM (SPEC)'
                WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
                WRITE(6,*)
                WRITE(6,*)
     &'** WARNING SPECINTF: PROBLEM WITH INTEGRATION OF SPECTRUM (SPEC)'
                WRITE(LUNGFO,*)
     &            '(FLUX THROUGH PINHOLE)'
                WRITE(6,*)'CHECK RESULTS CAREFULLY'
                IWARNW=1
              ENDIF

              ILIOB=ISOUR+NSOURCE*(IOBSV-1)
              SPECIF(ILIOB)=RESULT
              SPECTOTIF(IOBSV)=SPECTOTIF(IOBSV)+SPECIF(ILIOB)

            ENDDO !IZ
          ENDDO   !IY


C--- INTEGRATION OF FLUX THROUGH PINHOLE

          DO IFREQ=NFREQ0M,NFREQ0P
            S2(IFREQ-NFREQ0M+1)=
     &        WFLUXF(ISOUR+NSOURCE*(IFREQ-1))/BANWID*ECHARGE1
          ENDDO   !IFREQ

          CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)

          CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)

          IF (IWARNW.EQ.0.AND.RESULT.NE.0.0
     &    .AND.DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &'*** WARNING SPECINTF: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(LUNGFO,*)
     &        '(FLUX THROUGH PINHOLE)'
            WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
            WRITE(6,*)
            WRITE(6,*)
     &'*** WARNING SPECINTF: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(6,*)
     &        '(FLUX THROUGH PINHOLE)'
            WRITE(6,*)'CHECK RESULTS CAREFULLY'
            IWARNW=1
          ENDIF

          WFLUXIF(ISOUR)=RESULT
          WFLUXTIF=WFLUXTIF+WFLUXIF(ISOUR)

        ELSE    !IPIN

          DO IOBSV=1,NOBSV

C--- FILL INTEGRATION BUFFER

            DO IFREQ=NFREQ0M,NFREQ0P
              S2(IFREQ-NFREQ0M+1)=
     &          SPECF(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
     &          /BANWID*ECHARGE1
            ENDDO !IFREQ

C--- DO INTEGRATION OF BUFFER

            CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)

            IF (IWARNS.EQ.0.AND.RESULT.NE.0.0
     &      .AND.DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &'*** WARNING SPECINTF: PROBLEM WITH INTEGRATION OF SPECTRUM'
              WRITE(LUNGFO,*)
     &          '(FLUX THROUGH PINHOLE)'
              WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
              WRITE(6,*)
              WRITE(6,*)
     &'*** WARNING SPECINTF: PROBLEM WITH INTEGRATION OF SPECTRUM'
              WRITE(6,*)'CHECK RESULTS CAREFULLY'
              IWARNS=1
            ENDIF

            ILIOB=ISOUR+NSOURCE*(IOBSV-1)
            SPECIF(ILIOB)=RESULT
            SPECTOTIF(IOBSV)=SPECTOTIF(IOBSV)+SPECIF(ILIOB)

          ENDDO   !IOBSV

        ENDIF    !IPIN

      ENDDO !ISOUR

      RETURN
      END
