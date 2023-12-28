*CMZ :  4.01/04 27/11/2023  13.42.16  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.39.06  by  Michael Scheer
*CMZ :  3.05/09 10/02/2005  12.58.04  by  Michael Scheer
*CMZ :  2.52/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.36/00 07/11/2001  14.17.58  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  17.23.02  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.16/05 04/08/2000  11.46.38  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.28  by  Michael Scheer
*CMZ :  2.13/02 08/12/99  18.50.05  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  13.48.27  by  Michael Scheer
*CMZ :  1.01/01 10/12/97  13.26.03  by  Michael Scheer
*CMZ :  1.00/00 30/06/97  11.25.46  by  Michael Scheer
*CMZ : 00.02/05 18/03/97  15.39.20  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.24.12  by  Michael Scheer
*CMZ : 00.00/07 03/06/94  10.14.21  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.24  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.00  by  Michael Scheer
*-- Author :
      SUBROUTINE SPECINT

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C     INTEGRATES POWER SPECTRA OVER ALL PHOTON ENERGIES

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,source.
      include 'source.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      CHARACTER(50) CDUM
      INTEGER IY,IZ,ISOUR,IOBSV,IFREQ,IERR
      INTEGER ICAL,NMU,IMU,IWARNW,IWARNS
      DOUBLE PRECISION S2(NDFREQP),RESULT,SIMPLE
      DOUBLE PRECISION EMUDUM(1000),AMUDUM(1000),DENDUM

      DATA ICAL/0/

      IF (NFREQ0.LE.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING SR SPECINT ***'
        WRITE(LUNGFO,*)'INTEGRATION NOT POSSIBLE, ONLY ONE FREQUENCY'
        WRITE(LUNGFO,*)'CHANGE PARAMETER NFREQ2P IN NAMELIST CONTRL'
        WRITE(LUNGFO,*)
c        WRITE(6,*)
c        WRITE(6,*)'*** WARNING SR SPECINT ***'
c        WRITE(6,*)'INTEGRATION NOT POSSIBLE, ONLY ONE FREQUENCY'
c        WRITE(6,*)'CHANGE PARAMETER NFREQ2P IN NAMELIST CONTRL'
c        WRITE(6,*)
        ISPECINT=0
        RETURN
      ENDIF

C--- LOOP OVER ALL SOURCES

      WFLUXTI=0.0
      DO IOBSV=1,NOBSV
        SPECTOTI(IOBSV)=0.0
        ENEDOS(IOBSV)=0.0
      ENDDO !IOBSV

      IWARNS=0
      IWARNW=0
      DO ISOUR=1,NSOURCE

C--- LOOP OVER ALL OBSERVATION POINTS

        IF (IPIN.NE.0) THEN

          DO IY=1,NOBSVY
            DO IZ=1,NOBSVZ

              IOBSV=(IY-1)*NOBSVZ+IZ

C--- FILL INTEGRATION BUFFER

              DO IFREQ=NFREQ0M,NFREQ0P
                S2(IFREQ-NFREQ0M+1)=
     &            SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
     &            /BANWID*ECHARGE1
              ENDDO !IFREQ

C--- DO INTEGRATION OF BUFFER

              CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)

              IF (IWARNS.EQ.0.AND.RESULT.NE.0.0.AND.
     &              DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTUM (SPEC)'
                WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
                WRITE(6,*)
                WRITE(6,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM (SPEC)'
                WRITE(6,*)'CHECK RESULTS CAREFULLY'
                IWARNS=1
              ENDIF

              ILIOB=ISOUR+NSOURCE*(IOBSV-1)
              SPECI(ILIOB)=RESULT
              SPECTOTI(IOBSV)=SPECTOTI(IOBSV)+SPECI(ILIOB)

            ENDDO !IZ
          ENDDO !IY


C--- INTEGRATION OF FLUX THROUGH PINHOLE

          DO IFREQ=NFREQ0M,NFREQ0P
            S2(IFREQ-NFREQ0M+1)=
     &        WFLUX(ISOUR+NSOURCE*(IFREQ-1))/BANWID*ECHARGE1
          ENDDO !IFREQ

          CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)

          IF (IWARNW.EQ.0.AND.RESULT.NE.0.0.AND.
     &        DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &'*** WARNING SR SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(LUNGFO,*)
     &        '(FLUX THROUGH PINHOLE)'
            WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
            WRITE(6,*)
            WRITE(6,*)
     &'*** WARNING SR SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(6,*)
     &        '(FLUX THROUGH PINHOLE)'
            WRITE(6,*)'CHECK RESULTS CAREFULLY'
            IWARNW=1
          ENDIF

          WFLUXI(ISOUR)=RESULT
          WFLUXTI=WFLUXTI+WFLUXI(ISOUR)

        ELSE   !IPIN


          DO IOBSV=1,NOBSV

C--- FILL INTEGRATION BUFFER

            DO IFREQ=NFREQ0M,NFREQ0P
              S2(IFREQ-NFREQ0M+1)=
     &          SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
     &          /BANWID*ECHARGE1
            ENDDO !IFREQ

C--- DO INTEGRATION OF BUFFER

            CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)
            IF (IWARNS.EQ.0.AND.RESULT.NE.0.0.AND.
     &          DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM (SPEC)'
              WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
              WRITE(6,*)
              WRITE(6,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM (SPEC)'
              WRITE(6,*)'CHECK RESULTS CAREFULLY'
              IWARNS=1
            ENDIF

            ILIOB=ISOUR+NSOURCE*(IOBSV-1)
            SPECI(ILIOB)=RESULT
            SPECTOTI(IOBSV)=SPECTOTI(IOBSV)+SPECI(ILIOB)

          ENDDO !IOBSV

        ENDIF    !IPIN

      ENDDO !ISOUR

C--- DOSE CALCULATIONS

      IF (IDOSE.NE.0) THEN

        IWARNW=0
        IWARNS=0

        IF (ICAL.EQ.0) THEN
          OPEN(UNIT=99,FILE='ABSORPDOSE.RP',STATUS='OLD')

          READ(99,'(A50)') CDUM
          READ(99,*) DENDUM
          READ(99,*) NMU
          IF (NMU.GT.1000) STOP '*** SR SPECINT: NMU EXCEEDED 1000 ***'
          DO IMU=1,NMU
            READ(99,*) EMUDUM(IMU),AMUDUM(IMU)
          ENDDO !NMU
          DO IFREQ=NFREQ0M,NFREQ0P
            CALL ABSNOSPLI
     &        (EMUDUM,AMUDUM,NMU,FREQ(IFREQ-NFREQ0M+1),
     &        ABSMUEN(IFREQ-NFREQ0M+1),IERR,IDOSE)
            IF (IERR.NE.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** ERROR IN SPECINT ***'
              WRITE(LUNGFO,*)'CALL TO SR ABSNOSPLI FAILED'
              WRITE(LUNGFO,*)
     &'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE ABSORPDOSE.RP'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(6,*)
              WRITE(6,*)'*** ERROR IN SPECINT ***'
              WRITE(6,*)'CALL TO SR ABSNOSPLI FAILED'
              WRITE(6,*)
     &'CHECK PHOTONENERGIES IN NAMELIST FREQN AND FILE ABSORPDOSE.RP'
              WRITE(6,*)
              WRITE(6,*)
              STOP
            ENDIF   !IERR
          ENDDO
          CLOSE(99)
          ICAL=1
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &'     SPECINT: MASS ENERGY ABSORPTION COEFFS READ FROM FILE'
          WRITE(LUNGFO,*)'     COMMENT ON FILE IS:'
          WRITE(LUNGFO,*)'     ',CDUM
          WRITE(LUNGFO,*)
        ENDIF !ICAL

        DO IOBSV=1,NOBSV

C--- FILL INTEGRATION BUFFER FOR ABSORBED DOSE

          DO IFREQ=NFREQ0M,NFREQ0P
            S2(IFREQ-NFREQ0M+1)=
     &        SPECTOT(IOBSV+NOBSV*(IFREQ-1))/BANWID*ECHARGE1
            S2(IFREQ-NFREQ0M+1)=
     &        S2(IFREQ-NFREQ0M+1)*ABSMUEN(IFREQ-NFREQ0M+1)
          ENDDO !IFREQ

C--- DO INTEGRATION OF BUFFER

          CALL SPBUFINT(FREQ,S2,NFREQ0,RESULT,SIMPLE)
          IF (IWARNS.EQ.0.AND.RESULT.NE.0.0.AND.
     &        DABS((RESULT-SIMPLE)/RESULT).GT.1.D-1) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(LUNGFO,*)
     &        '(ABSORBED DOSE)'
            WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
            WRITE(6,*)
            WRITE(6,*)
     &'*** WARNING SR SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(6,*)
     &        '(ABSORBED DOSE)'
            WRITE(6,*)'CHECK RESULTS CAREFULLY'
            IWARNS=1
          ENDIF

          ENEDOS(IOBSV)=ENEDOS(IOBSV)+RESULT

        ENDDO !IOBSV

        IF (IPIN.NE.0) THEN

          IWARNW=0
          DO IFREQ=NFREQ0M,NFREQ0P
            S2(IFREQ-NFREQ0M+1)=
     &        WFLUXT(IFREQ-NFREQ0M+1)/BANWID*ECHARGE1/PINW/PINH
            S2(IFREQ-NFREQ0M+1)=
     &        S2(IFREQ-NFREQ0M+1)*ABSMUEN(IFREQ-NFREQ0M+1)
          ENDDO !IFREQ

C--- DO INTEGRATION OF BUFFER

          CALL SPBUFINT(FREQ,S2,NFREQ0,PINDOS,SIMPLE)

          IF (IWARNW.EQ.0.AND.PINDOS.NE.0.0.AND.
     &        DABS((PINDOS-SIMPLE)/PINDOS).GT.1.D-1) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(LUNGFO,*)
     &        '(ABSORBED ENERGY DOSE THROUGH PINHOLE)'
            WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
            WRITE(6,*)
            WRITE(6,*)
     &'*** WARNING SPECINT: PROBLEM WITH INTEGRATION OF SPECTRUM'
            WRITE(6,*)
     &        '(ABSORBED ENERGY DOSE THROUGH PINHOLE)'
            WRITE(6,*)'CHECK RESULTS CAREFULLY'
            IWARNW=1
          ENDIF

        ENDIF !IPIN

      ENDIF !IDOSE

      IF (IPINCIRC.EQ.0) THEN
        DO ISOUR=1,NSOURCE
          DO IZ=1,NOBSVZ
            CALL BLENDSPECIV(ISOUR,IZ)
          ENDDO
        ENDDO
      ENDIF !(IPINCIRC.EQ.0)

      RETURN
      END
