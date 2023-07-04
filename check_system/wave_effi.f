*CMZ :  2.57/05 10/03/2006  10.57.34  by  Michael Scheer
*CMZ :  2.52/13 14/12/2004  09.21.02  by  Michael Scheer
*-- Author :    Michael Scheer   09/12/2004
      PROGRAM WAVEEFFI

        IMPLICIT NONE

        DOUBLE PRECISION FREQLOW,FREQHIG,DFREQ,FREQL,FREQH,
     &    X1,Y1,Z1,POW1,POWS1,POWSF1,
     &    X2,Y2,Z2,POW2,POWS2,POWSF2

        EXTERNAL LNBLNK
        INTEGER NINTFREQ,NDFREQ,IFREQ,NFREQ,ILOOP,NLOOP,LAST,LNBLNK

        CHARACTER(8) DAY
        CHARACTER(10) TIM

        CHARACTER(128) CVAL(3),CHSUB(3)
        CHARACTER(128) CLINE

        PRINT*
        PRINT*,'Program WAVEEFFI'
        PRINT*,'Runs WAVE in a loop to sum up power yield'
        PRINT*
        PRINT*
        PRINT*,'*** Set in wave.in.perl:'
        PRINT*
        PRINT*,'   IEFFI=[0,1,-2,-1]'
        PRINT*,'   FILEPOW=wave.pow'
        PRINT*
        PRINT*,'   FREQLOW=9999.'
        PRINT*,'   FREQHIG=9999.'
        PRINT*,'   NINTFREQ=9999'
        PRINT*
        PRINT*,'Enter lower and upper photon energy:'
        READ(5,*)FREQL,FREQH
17      PRINT*
        PRINT*,'Total number of photon energies and number per WAVE run:'
        READ(5,*)NFREQ,NDFREQ

        IF (NDFREQ.LT.2) THEN
          PRINT*,'Error: Number per run must be greater than 1'
          GOTO 17
        ENDIF

        NLOOP=1+(NFREQ-NDFREQ)/(NDFREQ-1)

        IF (NLOOP.LT.2) THEN
          STOP '*** Error in WAVE_EFFI: Total number lower than 2'
        ENDIF

        IF (NDFREQ+(NLOOP-1)*(NDFREQ-1).LT.NFREQ) NLOOP=NLOOP+1

        PRINT*
        PRINT*,'   Number of loops to be done:',NLOOP
        PRINT*

        DFREQ=(FREQH-FREQL)/(NFREQ-1)
        FREQLOW=FREQL

        IFREQ=1

        DO ILOOP=1,NLOOP

          IF (IFREQ+(NDFREQ-1).GT.NFREQ) THEN
            NDFREQ=NFREQ-IFREQ+1
          ENDIF

          IFREQ=IFREQ+(NDFREQ-1)

          NINTFREQ=NDFREQ
          FREQHIG=FREQLOW+(NINTFREQ-1)*DFREQ

          IF (NINTFREQ.LE.1) GOTO 99

          CALL DATE_AND_TIME(DAY,TIM)
          PRINT*, ILOOP,'of',NLOOP,SNGL(FREQLOW),SNGL(FREQHIG),'  ',
     &      TIM(1:2),':',TIM(3:4),':',TIM(5:6),' '

          WRITE(CVAL(1),*)FREQLOW
          CHSUB(1)=' FREQLOW='//CVAL(1)
          WRITE(CVAL(2),*)FREQHIG
          CHSUB(2)=' FREQHIG='//CVAL(2)
          WRITE(CVAL(3),*)NINTFREQ
          CHSUB(3)=' NINTFREQ='//CVAL(3)

          OPEN(UNIT=20,FILE='wave.in.perl',status='old',readonly)
          CALL SYSTEM('rm -f wave.in')

          OPEN(UNIT=21,FILE='wave.in',status='new')

10        CONTINUE

            READ(20,'(A)',END=90)CLINE

            LAST = LNBLNK(CLINE)

            CALL CCOSUB(CLINE,LAST,CLINE,1,128,'FREQLOW=9999.',CHSUB(1))
            CALL CCOSUB(CLINE,LAST,CLINE,1,128,'FREQHIG=9999.',CHSUB(2))
            CALL CCOSUB(CLINE,LAST,CLINE,1,128,'NINTFREQ=9999',CHSUB(3))

            LAST = LNBLNK(CLINE)
            WRITE(21,'(A)')CLINE(1:LAST)

          GOTO 10

90        CLOSE(21)
          CLOSE(20)

          FREQLOW=FREQHIG

C------------

          CALL SYSTEM('rm -f wave.pow')
          CALL SYSTEM('nice -n 2 bin/wave.exe > wave.log')

          IF (ILOOP.EQ.1) THEN

            CALL SYSTEM('mv -f wave.pow wave.pow.sum')

          ELSE

            CALL SYSTEM('rm -f wave.effi')
            CALL SYSTEM('rm -f wave.effi_f')
            CALL SYSTEM('mv -f wave.pow.sum wave.pow.dum')
            CALL SYSTEM('rm -f wave.pow.sum')

            OPEN(UNIT=20,FILE='wave.pow.sum',status='new')
            OPEN(UNIT=21,FILE='wave.pow.dum',status='old',readonly)
            OPEN(UNIT=22,FILE='wave.pow',status='old',readonly)
            OPEN(UNIT=23,FILE='wave.effi_f',status='new')
            OPEN(UNIT=24,FILE='wave.effi',status='new')

2           CONTINUE

              READ(21,*,END=92)X1,Y1,Z1,POW1,POWS1,POWSF1
              READ(22,*)X2,Y2,Z2,POW2,POWS2,POWSF2

              IF (X1.NE.X2.OR.Y1.NE.Y2.OR.Z1.NE.Z2)
     &          STOP '*** Error: wave.pow does not match'

              POW1=POW1+POW2
              POWS1=POWS1+POWS2
              POWSF1=POWSF1+POWSF2

              WRITE(20,'(6(1PE13.5))')X1,Y1,Z1,POW1,POWS1,POWSF1
              WRITE(23,'(3(1PE13.5))')Z1,Y1,POWSF1
              WRITE(24,'(3(1PE13.5))')Z1,Y1,POWS1

            GOTO 2

92          CLOSE(20)
            CLOSE(21)
            CLOSE(22)
            CLOSE(23)
            CLOSE(24)

          ENDIF !ILOOP.EQ.1

        ENDDO !LOOP

99    CONTINUE

      PRINT*
      PRINT*,
     &  '--- Results written to wave.pow.sum, wave.effi, and wave.effi_f ---'
      PRINT*

      STOP
      END
