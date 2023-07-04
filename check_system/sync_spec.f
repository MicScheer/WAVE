*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.41/10 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.46  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.13/03 17/12/99  11.43.31  by  Michael Scheer
*CMZ :  2.13/00 15/11/99  17.43.14  by  Michael Scheer
*CMZ :  1.03/03 23/02/98  14.54.49  by  Michael Scheer
*CMZ :  1.03/02 23/02/98  14.30.43  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  18.47.22  by  Michael Scheer
*CMZ : 00.00/07 25/05/94  17.53.18  by  Michael Scheer
*-- Author :    Michael Scheer   18/05/94
      SUBROUTINE SYNC_SPEC
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

C     ROUTINE TO READ FILE OF SPECTRUM OF SYNC AND TO CONVERT TO DOSE
C     FILE: DES:SYNC_SP.DAT

C     MODIFICATION: IDESYNC_N.EQ.-9999. SIMPLE VERSION FOR EGS4-OUTPUT

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      CHARACTER(65) COMLAT,COMMAG
      INTEGER NSPEC,NBIN,IBIN,JSPEC,ISEG,ILAY,ISEL
      REAL*4 CORR,EDUM,CDUM

      IF (IDESYNC.NE.-9999) THEN

      IPIN=0
      IHPIN=0
      NOBSV=1
      NOBSVZ=1
      NOBSVY=1
      OBS1X=9999.
      OBS1Y=9999.
      OBS1Z=9999.
      OBSV(1,1)=OBS1X
      OBSV(2,1)=OBS1Y
      OBSV(3,1)=OBS1Z
      OBSVZ(1)=OBS1Z
      OBSVY(1)=OBS1Y

      WRITE(6,*)' '
      WRITE(6,*)'      *** SR SYNC_SPEC:'
      WRITE(6,*)' '

      OPEN(UNIT=90,FILE='DES:SYNC_SP.DAT',STATUS='OLD')
          READ(90,'(1A65)')COMLAT
          READ(90,'(1A65)')COMMAG
          READ(90,*)EDUM,CDUM
          IF (DABS(DMYENERGY-EDUM)/EDUM.GT.1.D-3
     &        .OR.
     &        DABS(DMYCUR-CDUM)/CDUM.GT.1.D-3) THEN
              WRITE(6,*)' '
              WRITE(6,*)
     &'*** WARNING SR SYNC_SPEC: ENERGY OR CURRENT ON DATA FILE NOT CONSISTENT WITH ORIGINAL VALUE ON FILE WAVE.IN'
              WRITE(6,*)' OLD VALUES OVERWRITTEN'
              WRITE(6,*)' '
              WRITE(LUNGFO,*)' '
              WRITE(LUNGFO,*)
     &'*** WARNING SR SYNC_SPEC: ENERGY OR CURRENT ON DATA FILE NOT CONSISTENT WITH ORIGINAL VALUE ON FILE WAVE.IN'
              WRITE(LUNGFO,*)' OLD VALUES OVERWRITTEN'
              WRITE(LUNGFO,*)' '
              DMYENERGY=EDUM
              DMYCUR=CDUM
          ENDIF
          DMYGAMMA=DMYENERGY/EMASSG1
          READ(90,*)CORR,BANWID
          READ(90,*)NSPEC,NBIN
          WRITE(6,*)' '
          WRITE(6,*)' Comments on file DES:SYNC_SP.DAT:'
          WRITE(6,*)COMLAT
          WRITE(6,*)COMMAG
          WRITE(6,*)' '
          IF (IDESYNC.LT.0) THEN
100         WRITE(6,*)' Number of spectra:',NSPEC
            WRITE(6,*)' Which one do you want?'
            READ(5,*)ISEL
            IF (ISEL.GT.NSPEC.OR.ISEL.LE.0) GOTO 100
200         WRITE(6,*)' Enter area for dose calculation [mm**2]:'
            READ(5,*)AREAM2
            IF (AREAM2.LE.0.) GOTO 200
            AREAM2=AREAM2*1.E-6
          ELSE
            ISEL=IDESYNC
            IF (ISEL.GT.NSPEC) THEN
               WRITE(6,*)' '
               WRITE(6,*)
     &' *** ERROR IN SYNC_SP: IDESYNC OUT OF RANGE, CHECK INPUT FILE *** '
               WRITE(6,*)' '
               WRITE(LUNGFO,*)' '
               WRITE(LUNGFO,*)
     &         ' *** ERROR IN SYNC_SP: IDESYNC OR ISEL OUT OF RANGE ***'
               WRITE(LUNGFO,*)' '
               WRITE(LUNGFO,*)' '
               STOP
            ENDIF
          ENDIF
          DO JSPEC=1,ISEL
          READ(90,*)ISEG,ILAY
          WRITE(6,*)
          WRITE(6,*)' SEGMENT, LAYER:',ISEG,ILAY
          WRITE(6,*)
          DO IBIN=1,NBIN
             IF (IBIN.GT.NDFREQP) THEN
                WRITE(6,*)
     &          '*** SR SYNC_SP: DIMENSION NDFREQP EXCEEDED ***'
                STOP
             ENDIF
            ILIOBFR=1+NSOURCE*NOBSV*(IBIN-1)
             READ(90,*)FREQ(IBIN),SPEC(ILIOBFR)
             FREQ(IBIN)=FREQ(IBIN)*1000.
             SPEC(ILIOBFR)=SPEC(ILIOBFR)/ECHARGE1/AREAM2
          ENDDO !IBIN
          ENDDO   !JSPEC
      CLOSE(90)
      NFREQ=NBIN

      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)'      *** SR SYNC_SPEC:'
      WRITE(LUNGFO,*)' '

      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)'      Comments on file DES:SYNC_SP.DAT:'
      WRITE(LUNGFO,*)'      ',COMLAT
      WRITE(LUNGFO,*)'      ',COMMAG
      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)'      BEAM ENERGY, CURRENT:',EDUM,CDUM
      WRITE(LUNGFO,*)'      SEGMENT, LAYER:      ',ISEG,ILAY
      WRITE(LUNGFO,*)'      AREA [mm**2]:        ',SNGL(AREAM2)*1.E6
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      ELSE    !IDESYNC.EQ.-9999

      IPIN=0
      IHPIN=0
      NOBSV=1
      NOBSVZ=1
      NOBSVY=1
      OBS1X=9999.
      OBS1Y=9999.
      OBS1Z=9999.
      OBSV(1,1)=OBS1X
      OBSV(2,1)=OBS1Y
      OBSV(3,1)=OBS1Z
      OBSVZ(1)=OBS1Z
      OBSVY(1)=OBS1Y

      WRITE(6,*)' '
      WRITE(6,*)'      *** SR SYNC_SPEC:'
      WRITE(6,*)' '

      OPEN(UNIT=90,FILE='KEK:KEK_EGS4_DOSE.WAV',STATUS='OLD')

          READ(90,'(1A65)')COMLAT
          READ(90,*)EDUM,CDUM

          IF (DABS(DMYENERGY-EDUM)/EDUM.GT.1.D-3
     &        .OR.
     &        DABS(DMYCUR-CDUM)/CDUM.GT.1.D-3) THEN
              WRITE(6,*)' '
              WRITE(6,*)
     &'*** WARNING SR SYNC_SPEC: ENERGY OR CURRENT ON DATA FILE NOT CONSISTENT WITH ORIGINAL VALUE ON FILE WAVE.IN'
              WRITE(6,*)' OLD VALUES OVERWRITTEN'
              WRITE(6,*)' '
              WRITE(LUNGFO,*)' '
              WRITE(LUNGFO,*)
     &'*** WARNING SR SYNC_SPEC: ENERGY OR CURRENT ON DATA FILE NOT CONSISTENT WITH ORIGINAL VALUE ON FILE WAVE.IN'
              WRITE(LUNGFO,*)' OLD VALUES OVERWRITTEN'
              WRITE(LUNGFO,*)' '
              DMYENERGY=EDUM
              DMYCUR=CDUM
          ENDIF

          DMYGAMMA=DMYENERGY/EMASSG1
          CORR=1.D0
        BANWID=1.D0
          NSPEC=1
          ISEL=1

          READ(90,*)NBIN,AREAM2
        AREAM2=AREAM2*1.E-4
           WRITE(6,*)' '
           WRITE(6,*)' Comments on file KEK:KEK_EGS4_DOSE.WAV:'
           WRITE(6,*)COMLAT
           WRITE(6,*)' '

             IF (IBIN.GT.NDFREQP) THEN
                WRITE(6,*)
     &          '*** SR SYNC_SP: DIMENSION NDFREQP EXCEEDED ***'
                STOP
             ENDIF

          DO IBIN=1,NBIN
            ILIOBFR=1+NSOURCE*NOBSV*(IBIN-1)
             READ(90,*)FREQ(IBIN),SPEC(ILIOBFR)
             FREQ(IBIN)=FREQ(IBIN)*1000.
             SPEC(ILIOBFR)=SPEC(ILIOBFR)/AREAM2
          ENDDO !IBIN

      CLOSE(90)
      NFREQ=NBIN

      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)'      *** SR SYNC_SPEC:'
      WRITE(LUNGFO,*)' '

      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)'      Comments on file KEK:KEK_EGS4_DOSE.WAV:'
      WRITE(LUNGFO,*)'      ',COMLAT
      WRITE(LUNGFO,*)' '
      WRITE(LUNGFO,*)'      BEAM ENERGY, CURRENT:',EDUM,CDUM
      WRITE(LUNGFO,*)'      AREA [mm**2]:        ',SNGL(AREAM2)*1.E6
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      ENDIF   !IDESYNC.EQ.-9999

      RETURN
      END
