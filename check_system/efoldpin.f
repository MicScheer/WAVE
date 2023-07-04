*CMZ :  4.00/07 06/04/2020  10.44.37  by  Michael Scheer
*CMZ :  4.00/04 10/05/2019  14.43.12  by  Michael Scheer
*CMZ :  3.03/02 03/12/2015  13.35.35  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.61/02 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.52/16 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/00 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.51/00 13/05/2004  12.01.48  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ : 00.01/06 01/02/95  16.22.16  by  Michael Scheer
*CMZ : 00.01/04 25/01/95  16.36.51  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.19.11  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.43  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE EFOLDPIN
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
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- THE ROUTINE FOLDS THE STOKES VECTORS INSIDE THE PINHOLE
C      WITH A GAUSSIAN TO TAKE THE BEAM ENERGY SPREAD INTO ACCOUNT

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEND.

      INTEGER IFREQ
      INTEGER IOBSV
      REAL*4 F,S0,S1,S2,S3

      INTEGER NFOLDP,kmode
      PARAMETER (NFOLDP=1000)

      DOUBLE PRECISION DF,F3SIG
      DOUBLE PRECISION S0E(NDFREQP),S2E(NDFREQP),S3E(NDFREQP),S1E(NDFREQP)
      DOUBLE PRECISION S0EF(NDFREQP),S2EF(NDFREQP),S3EF(NDFREQP),S1EF(NDFREQP)

      IF (IPIN.NE.0.AND.IPIN.NE.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN EFOLDPIN: IPIN is not 1  ***'
        WRITE(LUNGFO,*)'leaving IEFOLDPIN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN EFOLDPIN: IPIN is not 1  ***'
        WRITE(6,*)'leaving IEFOLDPIN'
        WRITE(6,*)
        RETURN
      ENDIF

      IF (iabs(IEFOLD).NE.1.and.iefold.ne.3) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN EFOLDPIN: IEFOLD not in [-1,1,3] ***'
        WRITE(LUNGFO,*)'leaving IEFOLDPIN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN EFOLDPIN: IEFOLD not in [-1,1,3] ***'
        WRITE(6,*)'leaving IEFOLDPIN'
        WRITE(6,*)
        RETURN
      ENDIF

      IF (IFOLD.NE.0.AND.IFOLD.NE.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING IN EFOLDPIN: IFOLD is neither 0 or nor 1  ***'
        WRITE(LUNGFO,*)'leaving IEFOLDPIN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
     &    '*** WARNING IN EFOLDPIN: IFOLD is neither 0 or nor 1  ***'
        WRITE(6,*)'leaving IEFOLDPIN'
        WRITE(6,*)
        RETURN
      ENDIF

      DO IOBSV=1,NOBSV

        DO IFREQ=1,NFREQ
          F=FREQ(IFREQ)
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          S0=STOKES(1,IOBFR)
          S1=STOKES(2,IOBFR)
          S2=STOKES(3,IOBFR)
          S3=STOKES(4,IOBFR)
          S0E(IFREQ)=S0
          S1E(IFREQ)=S1
          S2E(IFREQ)=S2
          S3E(IFREQ)=S3
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S1E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S1EF(IFREQ))
            kmode=1
          ELSE
            S1EF(IFREQ)=9.999D-11
          ENDIF
          STOKESE(2,IOBFR)=S1EF(IFREQ)
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S2E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S2EF(IFREQ))
            kmode=1
          ELSE
            S2EF(IFREQ)=9.999D-11
          ENDIF
          STOKESE(3,IOBFR)=S2EF(IFREQ)
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S3E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S3EF(IFREQ))
            kmode=1
          ELSE
            S3EF(IFREQ)=9.999D-11
          ENDIF
          STOKESE(4,IOBFR)=S3EF(IFREQ)
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S0E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S0EF(IFREQ))
            kmode=1
          ELSE
            S0EF(IFREQ)=9.999D-11
          ENDIF
          STOKESE(1,IOBFR)=S0EF(IFREQ)
        ENDDO !IFREQ

      ENDDO !NOBSV

      IF (IFOLD.NE.0) THEN

C STOKES DISTRIBUTION IN PINHOLE (FOLDED)

      DO IOBSV=1,NOBSV

        DO IFREQ=1,NFREQ
          F=FREQ(IFREQ)
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          S0=STOKESF(1,IOBFR)
          S1=STOKESF(2,IOBFR)
          S2=STOKESF(3,IOBFR)
          S3=STOKESF(4,IOBFR)
          S0E(IFREQ)=S0
          S1E(IFREQ)=S1
          S2E(IFREQ)=S2
          S3E(IFREQ)=S3
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S0E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S0EF(IFREQ))
            kmode=1
          ELSE
            S0EF(IFREQ)=9.999D-11
          ENDIF
          STOKESEF(1,IOBFR)=S0EF(IFREQ)
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S1E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S1EF(IFREQ))
            kmode=1
          ELSE
            S1EF(IFREQ)=9.999D-11
          ENDIF
          STOKESEF(2,IOBFR)=S1EF(IFREQ)
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S2E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S2EF(IFREQ))
            kmode=1
          ELSE
            S2EF(IFREQ)=9.999D-11
          ENDIF
          STOKESEF(3,IOBFR)=S2EF(IFREQ)
        ENDDO !IFREQ

        kmode=-1
        DO IFREQ=NFREQEM+1,NFREQ-NFREQEP
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          DF=FREQ(IFREQ)*ESPREAD*2.D0*NSIGE
          F3SIG=FREQ(IFREQ)*ESPREAD*2.D0
          IF(FREQ(IFREQ)-DF.GE.FREQ(1)
     &        .AND.
     &        FREQ(IFREQ)+DF.LE.FREQ(NFREQ)) THEN
            CALL EFOLD_GAUSS
     &        (NFREQ,FREQ,S3E,F3SIG,NSIGE,FREQ(IFREQ),kmode,S3EF(IFREQ))
            kmode=1
          ELSE
            S3EF(IFREQ)=9.999D-11
          ENDIF
          STOKESEF(4,IOBFR)=S3EF(IFREQ)
        ENDDO !IFREQ

      ENDDO !NOBSV

      ENDIF !(IFOLD.NE.0)

      RETURN
      END
