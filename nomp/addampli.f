*CMZ :  4.00/15 27/04/2022  08.29.38  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  10.39.59  by  Michael Scheer
*CMZ :  3.01/00 17/06/2013  08.45.46  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.69/02 08/11/2012  09.57.15  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  11.58.46  by  Michael Scheer
*CMZ :  2.66/12 24/06/2010  12.50.52  by  Michael Scheer
*CMZ :  2.66/08 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  09.23.32  by  Michael Scheer
*CMZ :  2.54/07 27/11/2009  15.57.31  by  Michael Scheer
*CMZ :  2.41/10 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.36/00 08/11/2001  14.44.28  by  Michael Scheer
*CMZ :  2.20/01 19/02/2001  12.17.18  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.40.08  by  Michael Scheer
*CMZ :  2.16/04 19/06/2000  12.30.42  by  Michael Scheer
*CMZ :  2.15/00 18/05/2000  11.40.29  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  11.45.20  by  Michael Scheer
*CMZ :  2.13/03 17/01/2000  16.17.57  by  Michael Scheer
*CMZ :  2.10/01 22/03/99  10.04.12  by  Michael Scheer
*CMZ :  2.00/02 12/01/99  16.57.12  by  Michael Scheer
*CMZ :  2.00/00 11/01/99  13.31.50  by  Michael Scheer
*-- Author :    Michael Scheer   05/01/99
      SUBROUTINE ADDAMPLI

*KEEP,GPLHINT.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      use bunchmod

C--- SUBROUTINE ADDAMPLI TO READ, WRITE AND TREAT ARRAY REAIMA OF
C--- COMPLEXE FIELD AMPLITUDES.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,cmzsel.
      include 'cmzsel.cmn'
*KEND.

      CHARACTER(80) CODEAMP,C80
      CHARACTER(134) FILESUPER
      CHARACTER CLAST

      INTEGER IO,IFR,IXYZ,IERROR,IREP,NTOTIN,NTOT2IN
      INTEGER IREPMXP,JCMZNOCMPLX,nelec
      PARAMETER (IREPMXP=1000)

      COMPLEX*16 DPHASE,DMODU,AX,AY,AZ,DDMODU(IREPMXP),DMODU0
     &  ,AX0,AY0,AZ0
      COMPLEX*8 APOLH,APOLR,APOLL,APOL45

      DOUBLE PRECISION AXR,AXI,AYR,AYI,AZR,AZI,PHI0,PHI,TPHASE,OMEGA,RANRMS,
     &  R0,R02,H2,H2R2,GAMMA21,TWOPI,OMEGAR,TMODULATOR,DTMOD
     &  ,DTPHASE,FREQR

      double precision xub,yub,zub,ypub,zpub,gammaub

      REAL STOK1,STOK2,STOK3,STOK4
      REAL XRAN(IREPMXP),xranmar(1),wbuff(10000),rr(2)

      INTEGER IOBSV

      INTEGER IGETLASTCHAR,ILAST
      EXTERNAL IGETLASTCHAR

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      SUBROUTINE ADDAMPLI:'
      WRITE(LUNGFO,*)

      iampobsv=icbrill
      codeamp=code

      DTPHASE=(WTRA2I+(1.D0/(DMYGAMMA*DMYGAMMA))
     &  *(XSTOP-XSTART)/2.D0)/CLIGHT1
      FREQR=2.D0*PI1/DTPHASE*HBAREV1

      IF (IMAMPLI.GT.0.AND.AMPFREQ.EQ.-9999.D0) THEN
        AMPFREQ=FREQR
      ENDIF !(IMAMPLI.GT.0.AND.AMPFREQ.EQ.-9999.D0)

      if (iamprep.ge.0) then

        IF (ABS(FREQR-AMPFREQ)/FREQR.GT.1.E-6) THEN
          WRITE(LUNGFO,*)'*** WARNING IN ADDAMPLI:'
          WRITE(LUNGFO,*)'1. harmonical and AMPFREQ differ.'
          WRITE(LUNGFO,*)SNGL(FREQR),SNGL(AMPFREQ)
          WRITE(6,*)'*** WARNING IN ADDAMPLI:'
          WRITE(6,*)'1. harmonical and AMPFREQ differ.'
          WRITE(6,*)SNGL(FREQR),SNGL(AMPFREQ)
        ENDIF  !(ABS(FREQR-AMPFREQ)/FREQR.GT.1.E-6)

C CURRENT SETTING MUST AGREE TO SETTING ON FILE WITH AMPLITUDES
C TO HAVE CORRECTLY DEFINED ALLOCATABLE ARRAYS

        nelec=0

      else

        nelec=-iamprep

      endif !iamprep.ge.0

      CALL AMPCHECK(
     &  NSOURCE,NOBSV,NFREQ,IFREQ2P,
     &  NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,
     &  MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY,
     &  PINW,PINH,PINR,IPIN,IF1DIM,IPINCIRC,AMPFREQ,iamprep,
     &  ibunch,iubunch,bunchlen,
     &  IERROR)

      IF (iamprep.ge.0.and.IAMPREP.LT.2
     &    .AND.(AMPPHI(1).NE.0.D0.OR.AMPSHIFT(1).NE.0.D0)) THEN
        C80=
     &    '*** WARNING IN ADDAMPLI: AMPPHI(1) or AMPHSHIFT(1) not zero'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)C80
        WRITE(6,*)
        WRITE(6,*)C80
      ENDIF !(AMPPHI(1).NE.0.D0.OR.AMPHSHIFT(1).NE.0.D0)

      GAMMA21=1.D0/DMYGAMMA/DMYGAMMA
      TWOPI=2.D0*PI1
      OMEGAR=AMPFREQ/HBAREV1

      IAMPREAD=1

      IF (NSOURCE.NE.1) THEN

        WRITE(LUNGFO,*)
     &    ' *** ERROR IN ADDAMPLI: NUMBER OF SOURCES MUST BE ONE!! ***'

        WRITE(6,*)
     &    ' *** ERROR IN ADDAMPLI: NUMBER OF SOURCES MUST BE ONE!! ***'

        STOP '--- PROGRAM ABORTED ---'

      ENDIF !(NSOURCE.NE.1)

      IF (IMAMPLI.GT.0.AND.IAMPSKIP.NE.0) THEN

        WRITE(LUNGFO,*)
     &    ' *** WARNING IN ADDAMPLI: IAMPSKIP.NE.0 .AND. IMAMPLI.GT.0 ***'
        WRITE(LUNGFO,*)
     &    ' *** CHECK NAMELIST $AMPLIN IN INPUT FILE WAVE.IN ***'

        WRITE(6,*)
     &    ' *** WARNING IN ADDAMPLI: IAMPSKIP.NE.0 .AND. IMAMPLI.GT.0 ***'
        WRITE(6,*)
     &    ' *** CHECK NAMELIST $AMPLIN IN INPUT FILE WAVE.IN ***'

      ENDIF !(IAMPSKIP.NE.0)

      IF (IMAMPLI.LT.0.AND.IAMPSKIP.EQ.0) THEN

        WRITE(LUNGFO,*)
     &    ' *** WARNING IN ADDAMPLI: IAMPSKIP.EQ.0 .AND. IMAMPLI.LT.0 ***'
        WRITE(LUNGFO,*)
     &    ' *** CHECK NAMELIST $AMPLIN IN INPUT FILE WAVE.IN ***'

        WRITE(6,*)
     &    ' *** WARNING IN ADDAMPLI: IAMPSKIP.EQ.0 .AND. IMAMPLI.LT.0 ***'
        WRITE(6,*)
     &    ' *** CHECK NAMELIST $AMPLIN IN INPUT FILE WAVE.IN ***'

      ENDIF !(IMAMPLI.LT.0.AND.IAMPSKIP.EQ.0)

      IF (IMAMPLI.EQ.1.OR.IMAMPLI.EQ.3) THEN

C--- CREATE NEW FILE AND WRITE TO IT

        OPEN(UNIT=LUNAMPLI,FILE=FILEAMPLI,STATUS='unknown'
     &    ,FORM='FORMATTED')

        WRITE(LUNAMPLI,*)ICODE,' ',CODE
        WRITE(LUNAMPLI,*)ICMZNOCMPLX
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)NSOURCE,NOBSV,NFREQ,IFREQ2P
        WRITE(LUNAMPLI,*)NOBSVZ,NOBSVY,MOBSVZ,MOBSVY
        WRITE(LUNAMPLI,*)MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY
        WRITE(LUNAMPLI,*)IPIN,IF1DIM,IPINCIRC
        WRITE(LUNAMPLI,*)IBUNCH,IUBUNCH,BUNCHLEN
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)PINCEN
        WRITE(LUNAMPLI,*)PINW,PINH,PINR
        WRITE(LUNAMPLI,*)OBSVDZ,OBSVDY
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)SPECNOR
        WRITE(LUNAMPLI,*)VPOLA(1)
        WRITE(LUNAMPLI,*)VPOLA(2)
        WRITE(LUNAMPLI,*)VPOLA(3)
        WRITE(LUNAMPLI,*)VSTOKES(1,1)
        WRITE(LUNAMPLI,*)VSTOKES(1,2)
        WRITE(LUNAMPLI,*)VSTOKES(1,3)
        WRITE(LUNAMPLI,*)VSTOKES(2,1)
        WRITE(LUNAMPLI,*)VSTOKES(2,2)
        WRITE(LUNAMPLI,*)VSTOKES(2,3)
        WRITE(LUNAMPLI,*)VSTOKES(3,1)
        WRITE(LUNAMPLI,*)VSTOKES(3,2)
        WRITE(LUNAMPLI,*)VSTOKES(3,3)
        WRITE(LUNAMPLI,*)VSTOKES(4,1)
        WRITE(LUNAMPLI,*)VSTOKES(4,2)
        WRITE(LUNAMPLI,*)VSTOKES(4,3)
        WRITE(LUNAMPLI,*)

        WRITE(LUNAMPLI,*)(OBSVZ(IO),IO=1,NOBSVZ)
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)(OBSVY(IO),IO=1,NOBSVY)
        WRITE(LUNAMPLI,*)

        DO IO=1,NOBSV
          WRITE(LUNAMPLI,*)(OBSV(IXYZ,IO),IXYZ=1,3)
        ENDDO

        WRITE(LUNAMPLI,*)AMPFREQ
        WRITE(LUNAMPLI,*)

        DO IFR=1,NFREQ
          WRITE(LUNAMPLI,*)
          WRITE(LUNAMPLI,*)FREQ(IFR)
          DO IO=1,NOBSV
            IOBFR=IO+NOBSV*(IFR-1)
            WRITE(LUNAMPLI,*)
     &        REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR)
            WRITE(LUNAMPLI,*)
     &        REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR)
            WRITE(LUNAMPLI,*)
     &        REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR)
          ENDDO
        ENDDO
        WRITE(LUNAMPLI,*)

        WRITE(LUNGFO,*)'      array REAIMA written to file'
        WRITE(LUNGFO,*)'      ',FILEAMPLI
        WRITE(LUNGFO,*)

        CLOSE(LUNAMPLI)

        IF (IMAMPLI.EQ.3) THEN
C--- WRITE FILES FOR PROGRAM PHASE OF JOHANNES BAHRDT
          CALL PHASE_BAHRDT
        ENDIF

      ELSE IF (IMAMPLI.EQ.2) THEN

C--- WRITE TO END OF EXISTING FILE

        OPEN(UNIT=LUNAMPLI,FILE=FILEAMPLI,STATUS='OLD'
     &    ,FORM='FORMATTED',ACCESS='APPEND')

        WRITE(LUNAMPLI,*)ICODE,' ',CODE
        WRITE(LUNAMPLI,*)ICMZNOCMPLX
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)NSOURCE,NOBSV,NFREQ,IFREQ2P
        WRITE(LUNAMPLI,*)NOBSVZ,NOBSVY,MOBSVZ,MOBSVY
        WRITE(LUNAMPLI,*)MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY
        WRITE(LUNAMPLI,*)IPIN,IF1DIM,IPINCIRC
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)PINCEN
        WRITE(LUNAMPLI,*)PINW,PINH,PINR
        WRITE(LUNAMPLI,*)OBSVDZ,OBSVDY
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)SPECNOR
        WRITE(LUNAMPLI,*)VPOLA(1)
        WRITE(LUNAMPLI,*)VPOLA(2)
        WRITE(LUNAMPLI,*)VPOLA(3)
        WRITE(LUNAMPLI,*)VSTOKES(1,1)
        WRITE(LUNAMPLI,*)VSTOKES(1,2)
        WRITE(LUNAMPLI,*)VSTOKES(1,3)
        WRITE(LUNAMPLI,*)VSTOKES(2,1)
        WRITE(LUNAMPLI,*)VSTOKES(2,2)
        WRITE(LUNAMPLI,*)VSTOKES(2,3)
        WRITE(LUNAMPLI,*)VSTOKES(3,1)
        WRITE(LUNAMPLI,*)VSTOKES(3,2)
        WRITE(LUNAMPLI,*)VSTOKES(3,3)
        WRITE(LUNAMPLI,*)VSTOKES(4,1)
        WRITE(LUNAMPLI,*)VSTOKES(4,2)
        WRITE(LUNAMPLI,*)VSTOKES(4,3)
        WRITE(LUNAMPLI,*)

        WRITE(LUNAMPLI,*)(OBSVZ(IO),IO=1,NOBSVZ)
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)(OBSVY(IO),IO=1,NOBSVY)
        WRITE(LUNAMPLI,*)

        DO IO=1,NOBSV
          WRITE(LUNAMPLI,*)(OBSV(IXYZ,IO),IXYZ=1,3)
        ENDDO

        WRITE(LUNAMPLI,*)AMPFREQ
        WRITE(LUNAMPLI,*)

        DO IFR=1,NFREQ
          WRITE(LUNAMPLI,*)
          WRITE(LUNAMPLI,*)FREQ(IFR)
          DO IO=1,NOBSV
            IOBFR=IO+NOBSV*(IFR-1)
            WRITE(LUNAMPLI,*)
     &        REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR)
            WRITE(LUNAMPLI,*)
     &        REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR)
            WRITE(LUNAMPLI,*)
     &        REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR)
          ENDDO
        ENDDO

        WRITE(LUNAMPLI,*)

        WRITE(LUNGFO,*)'      array REAIMA written to file'
        WRITE(LUNGFO,*)'      ',FILEAMPLI
        WRITE(LUNGFO,*)

        CLOSE(LUNAMPLI)

      ELSE IF (IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3
     &    .or.
     &    imampli.eq.0.and.iampskip.eq.0.and.iamprep.lt.0) THEN

C--- READ REAIMA FROM FILE, OVERWRITTING ARRAY

        IF (IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3) then

          OPEN(UNIT=LUNAMPLI,FILE=FILEAMPLI,STATUS='OLD'
     &      ,FORM='FORMATTED')

          WRITE(LUNGFO,*)'      reading spectra from file'
          WRITE(LUNGFO,*)'      ',FILEAMPLI
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '      length for phaseshifts [m], scaling factors, '//
     &      'jobnumbers and comments:'
          WRITE(LUNGFO,*)

          IAMPREAD=0

          AXR=0.0D0
          AYR=0.0D0
          AZR=0.0D0

          AXI=0.0D0
          AYI=0.0D0
          AZI=0.0D0

          DO IFR=1,NFREQ
            DO IO=1,NOBSV

              IOBFR=IO+NOBSV*(IFR-1)
              AFREQ(1,IOBFR)=(0.D0,0.D0)
              AFREQ(2,IOBFR)=(0.D0,0.D0)
              AFREQ(3,IOBFR)=(0.D0,0.D0)

            ENDDO !NOBSV
          ENDDO   !NFREQ

        end IF !(IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3) then

10      CONTINUE

        IF (IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3) then

          READ(LUNAMPLI,'(A80)',END=90)CODEAMP
          READ(LUNAMPLI,*)JCMZNOCMPLX
          IF (JCMZNOCMPLX.NE.ICMZNOCMPLX) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*) '*** ERROR IN ADDAMPLI ***'
            WRITE(LUNGFO,*)
     &        'DIFFERENT WAVE VERSION CONCERNING COMPLEX NUMBERS'
            WRITE(LUNGFO,*) 'USED FOR WRITING AND READING FILEAMPLI'
            WRITE(LUNGFO,*) '*** PROGRAM WAVE ABORTED ***'
            STOP
          ENDIF
          READ(LUNAMPLI,*)
          READ(LUNAMPLI,*)NSOURCE,NOBSV,NFREQ,IFREQ2P
          READ(LUNAMPLI,*)NOBSVZ,NOBSVY,MOBSVZ,MOBSVY
          READ(LUNAMPLI,*)MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY
          READ(LUNAMPLI,*)IPIN,IF1DIM,IPINCIRC
          read(LUNAMPLI,*)IBUNCH,IUBUNCH,BUNCHLEN
          READ(LUNAMPLI,*)
          READ(LUNAMPLI,*)PINCEN
          READ(LUNAMPLI,*)PINW,PINH,PINR

          IAMPOBSV=0
          DO IOBSV=1,NOBSV
            IF (ABS(OBSV(2,IOBSV)).LT.1.D-15
     &          .AND.ABS(OBSV(3,IOBSV)).LT.1.D-15) THEN
              IF (IAMPOBSV.NE.0) THEN
                WRITE(LUNGFO,*)
     &            ' *** ERROR IN ADDAMPLI: on-axis observation '//
     &            'point not unique ***'
                WRITE(6,*)
     &            ' *** ERROR IN ADDAMPLI: on-axis observation '//
     &            'point not unique ***'
                STOP '--- PROGRAM ABORTED ---'
              ENDIF  !IAMPOBSV.NE.0
              IAMPOBSV=IOBSV
            ENDIF !ABS(OBSV(3,IOBSV)).LT.1.D-15
          ENDDO   !NOBSV

          IF (IAMPOBSV.EQ.0) THEN
            WRITE(LUNGFO,*)
     &        ' *** ERROR IN ADDAMPLI: no on-axis observation point found ***'
            WRITE(6,*)
     &        ' *** ERROR IN ADDAMPLI: no on-axis observation point found ***'
            STOP '--- PROGRAM ABORTED ---'
          ENDIF   !IAMPOBSV.NE.0

          READ(LUNAMPLI,*)OBSVDZ,OBSVDY
          READ(LUNAMPLI,*)
          READ(LUNAMPLI,*)SPECNOR
          READ(LUNAMPLI,*)VPOLA(1)
          READ(LUNAMPLI,*)VPOLA(2)
          READ(LUNAMPLI,*)VPOLA(3)
          READ(LUNAMPLI,*)VSTOKES(1,1)
          READ(LUNAMPLI,*)VSTOKES(1,2)
          READ(LUNAMPLI,*)VSTOKES(1,3)
          READ(LUNAMPLI,*)VSTOKES(2,1)
          READ(LUNAMPLI,*)VSTOKES(2,2)
          READ(LUNAMPLI,*)VSTOKES(2,3)
          READ(LUNAMPLI,*)VSTOKES(3,1)
          READ(LUNAMPLI,*)VSTOKES(3,2)
          READ(LUNAMPLI,*)VSTOKES(3,3)
          READ(LUNAMPLI,*)VSTOKES(4,1)
          READ(LUNAMPLI,*)VSTOKES(4,2)
          READ(LUNAMPLI,*)VSTOKES(4,3)
          READ(LUNAMPLI,*)
          READ(LUNAMPLI,*)(OBSVZ(IO),IO=1,NOBSVZ)
          READ(LUNAMPLI,*)
          READ(LUNAMPLI,*)(OBSVY(IO),IO=1,NOBSVY)
          READ(LUNAMPLI,*)

          DO IO=1,NOBSV
            READ(LUNAMPLI,*)(OBSV(IXYZ,IO),IXYZ=1,3)
          ENDDO

          READ(LUNAMPLI,*)AMPFREQ
          READ(LUNAMPLI,*)

          DO IFR=1,NFREQ
            READ(LUNAMPLI,*)
            READ(LUNAMPLI,*)FREQ(IFR)
            DO IO=1,NOBSV
              IOBFR=IO+NOBSV*(IFR-1)
              READ(LUNAMPLI,*)
     &          REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR)
              READ(LUNAMPLI,*)
     &          REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR)
              READ(LUNAMPLI,*)
     &          REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR)
            ENDDO !NOBSV
          ENDDO   !NFREQ

          READ(LUNAMPLI,*)

          CALL AMPCHECK(
     &      NSOURCE,NOBSV,NFREQ,IFREQ2P,
     &      NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,
     &      MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY,
     &      PINW,PINH,PINR,IPIN,IF1DIM,IPINCIRC,AMPFREQ,iamprep,
     &      ibunch,iubunch,bunchlen,
     &      IERROR)

          IF (IERROR.NE.0) THEN
            WRITE(LUNGFO,*)
     &        ' *** ERROR IN ADDAMPLI: data on file incompatible ***'
            WRITE(6,*)
     &        ' *** ERROR IN ADDAMPLI: data on file incompatible ***'
            STOP '--- PROGRAM ABORTED ---'
          ENDIF   !IERROR

          IAMPREAD=IAMPREAD+1

        end if !(IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3

        IF (IAMPREP.GT.1.AND.IAMPREAD.GT.1) THEN
          WRITE(LUNGFO,*)
     &      ' *** ERROR IN ADDAMPLI: IAMPREP > 1, but more than ***'
          WRITE(LUNGFO,*)
     &      ' *** one device on amplitude file ***'
          WRITE(6,*)
     &      ' *** ERROR IN ADDAMPLI: IAMPREP > 1, but more than ***'
          WRITE(6,*)
     &      ' *** one device on amplitude file ***'
          STOP '--- PROGRAM ABORTED ---'
        ENDIF  !(IAMPREP.GT.0.AND.IAMPREAD.GT.1)

        IF ((IAMPREAD.GT.1.OR.IAMPREP.GT.1).AND.
     &      (AMPPHI(IAMPREAD).EQ.0.D0
     &      .OR.AMPSHIFT(IAMPREAD).EQ.0.D0)) THEN
          C80=
     &      '*** WARNING IN ADDAMPLI: AMPPHI(N) or AMPHSHIFT(N) zero'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)C80
          WRITE(LUNGFO,*)'N=',IAMPREAD
          WRITE(6,*)
          WRITE(6,*)C80
          WRITE(6,*)'N=',IAMPREAD
        ENDIF  !(AMPPHI(N).EQ.0.D0

        IF (IAMPREAD.GT.IAMPDIMP) THEN

          WRITE(LUNGFO,*)
     &      ' *** ERROR IN ADDAMPLI: DIMENSION IAMPDIMP EXCEEDED ***'
          WRITE(6,*)
     &      ' *** ERROR IN ADDAMPLI: INCREASE IAMPDIMP IN AMPLI.CMN ***'

          STOP '--- PROGRAM ABORTED ---'

        ENDIF  !(IAMPREAD.GT.IAMPDIMP)

        WRITE(LUNGFO,*)
     &    '      ',IAMPREAD,AMPSHIFT(IAMPREAD),AMPSCALE(IAMPREAD)
        WRITE(LUNGFO,*)'      ',CODEAMP
        WRITE(LUNGFO,*)

C PHASE ADVANCE OF MODULATOR IN TIME

        dtmod=0.0d0
        IF (AMPPHI(IAMPREAD).LT.0.D0) THEN
          DTMOD=-AMPPHI(IAMPREAD)/CLIGHT1
        ELSE if (omegar.ne.0.0d0) then !(AMPPHI(IAMPREAD).LT.0.D0)
          DTMOD=TWOPI*AMPPHI(IAMPREAD)/OMEGAR
        ENDIF  !(AMPPHI(IAMPREAD).LT.0.D0)

        TMODULATOR=DTMOD-AMPSHIFT(IAMPREAD)/2.D0/CLIGHT1*GAMMA21

        IF (TMODULATOR.LT.0.D0) THEN
          WRITE(LUNGFO,*)
     &      '*** ERROR IN ADDAMPLI: phaseshift lower than zero for device ',IAMPREAD
          WRITE(6,*)
     &      '*** ERROR IN ADDAMPLI: phaseshift lower than zero for device ',IAMPREAD
          STOP '--- PROGRAM ABORTED ---'
        ENDIF  !(TMODULATOR.LT.0.D0)

        R0=OBSV(1,IAMPOBSV)
        R02=R0*R0

        DO IFR=1,NFREQ

C--- ADJUST PHASE TO GIVEN OBSERVATION POINT AND COMPONENT

          IOBFR=IAMPOBSV+NOBSV*(IFR-1)

          if (iamprep.gt.0) then

            IF (REAIMA(IAMPCOMP,1,IOBFR).EQ.0.D0
     &          .AND.REAIMA(IAMPCOMP,2,IOBFR).EQ.0.D0) THEN
              PHI0=-9999.
            ELSE
              PHI0=ATAN2
     &          (REAIMA(IAMPCOMP,2,IOBFR)
     &          ,REAIMA(IAMPCOMP,1,IOBFR))
            ENDIF

            IF (PHI0.EQ.-9999.) THEN

              WRITE(LUNGFO,*)
     &          ' *** ERROR IN ADDAMPLI: PHASE ADJUSTMENT FAILED ***'
              WRITE(LUNGFO,*)
     &          ' TRY DIFFERENT IAMPCOMP IN NAMELIST AMPLIN'
              WRITE(LUNGFO,*)'IAMPCOMP NOW:',IAMPCOMP
              WRITE(LUNGFO,*)'IFREQ, IOBSV:',IFR,IAMPOBSV
              WRITE(LUNGFO,*)'REAIMA(IAMPCOMP,1,IOBFR):',
     &          REAIMA(IAMPCOMP,1,IOBFR)
              WRITE(LUNGFO,*)'REAIMA(IAMPCOMP,2,IOBFR):',
     &          REAIMA(IAMPCOMP,2,IOBFR)

              WRITE(6,*)
     &          ' *** ERROR IN ADDAMPLI: PHASE ADJUSTMENT FAILED ***'
              WRITE(6,*)
     &          ' TRY DIFFERENT IAMPCOMP IN NAMELIST AMPLIN'
              WRITE(6,*)'IAMPCOMP NOW:',IAMPCOMP
              WRITE(6,*)'IFREQ, IOBSV:',IFR,IAMPOBSV
              WRITE(6,*)'REAIMA(IAMPCOMP,1,IOBFR):',
     &          REAIMA(IAMPCOMP,1,IOBFR)
              WRITE(6,*)'REAIMA(IAMPCOMP,2,IOBFR):',
     &          REAIMA(IAMPCOMP,2,IOBFR)

              WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
              WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
              STOP

            ENDIF !(PHI0.EQ.0.D0)

          else ! iamprep>0
            phi0=0.0d0
          ENDIF   ! iamprep>0

          DO IO=1,NOBSV

            if (iamprep.gt.0) then
              H2= OBSV(2,IO)**2+OBSV(3,IO)**2
              H2R2=H2/R02
              TPHASE=TMODULATOR+AMPSHIFT(IAMPREAD)*(H2R2+GAMMA21)
     &          /2.D0/CLIGHT1
            else
              tphase=0.0d0
            endif

            OMEGA=FREQ(IFR)/HBAREV1
            PHI=TPHASE*OMEGA

            DPHASE=CDEXP (DCMPLX(0.D0,PHI-PHI0))
            DMODU=DCMPLX(AMPSCALE(IAMPREAD),0.D0)*DPHASE

            IOBFR=IO+NOBSV*(IFR-1)

            AXR=REAIMA(1,1,IOBFR)
            AYR=REAIMA(2,1,IOBFR)
            AZR=REAIMA(3,1,IOBFR)
            AXI=REAIMA(1,2,IOBFR)
            AYI=REAIMA(2,2,IOBFR)
            AZI=REAIMA(3,2,IOBFR)

            AX=DCMPLX(AXR,AXI)*DMODU
            AY=DCMPLX(AYR,AYI)*DMODU
            AZ=DCMPLX(AZR,AZI)*DMODU

            IOBFR=IO+NOBSV*(IFR-1)
            AFREQ(1,IOBFR)=AFREQ(1,IOBFR)+AX
            AFREQ(2,IOBFR)=AFREQ(2,IOBFR)+AY
            AFREQ(3,IOBFR)=AFREQ(3,IOBFR)+AZ

          ENDDO   !NOBSV
        ENDDO  !NFREQ

        IF (IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3) then
          GOTO 10
        ENDIF

90      continue
        IF (IMAMPLI.EQ.-1.OR.IMAMPLI.EQ.-3) then
          CLOSE(LUNAMPLI)
        endif

C--- SPECTRUM ARRAYS{

        IF (IAMPREP.NE.0) THEN

          IF (IAMPREP.GT.IREPMXP) THEN
            WRITE(LUNGFO,*)
     &        ' *** ERROR IN ADDAMPLI: DIMENSION IAMPMXP EXCEEDED ***'
            WRITE(6,*)
     &        ' *** ERROR IN ADDAMPLI: DIMENSION IAMPMXP EXCEEDED ***'
            STOP '--- PROGRAM ABORTED ---'
          ENDIF

          if (iampseed.gt.0) then

            IF (IAMPSEED.NE.0) THEN
              CALL RMARIN(IAMPSEED,NTOTIN,NTOT2IN) !CERN V113
            ENDIF

            CALL RNORML(XRAN,IAMPREP,rr)

            DO IREP=1,IAMPREP
              XRAN(IREP)=AMPRAN*XRAN(IREP)
              RANRMS=RANRMS+XRAN(IREP)**2
            ENDDO
            RANRMS=SQRT(RANRMS/IAMPREP)

          end if !(iampseed.gt.0) then

          DO IFR=1,NFREQ
            DO IO=1,NOBSV

              IOBFR=IO+NOBSV*(IFR-1)
              REAIMA(1,1,IOBFR)=0.D0
              REAIMA(2,1,IOBFR)=0.D0
              REAIMA(3,1,IOBFR)=0.D0
              REAIMA(1,2,IOBFR)=0.D0
              REAIMA(2,2,IOBFR)=0.D0
              REAIMA(3,2,IOBFR)=0.D0

            ENDDO !NOBSV
          ENDDO   !NFREQ

          R0=OBSV(1,IAMPOBSV)
          R02=R0*R0

          if (iamprep.gt.0) then

            DO IFR=1,NFREQ
              DO IO=1,NOBSV

cerror? 10.3.2010              if (iampseed.gt.0) then
                if (iampread.gt.0) then

                  H2= OBSV(2,IO)**2+OBSV(3,IO)**2
                  H2R2=H2/R02

                  TPHASE=TMODULATOR+AMPSHIFT(IAMPREAD)*(H2R2+GAMMA21)
     &              /2.D0/CLIGHT1

                  OMEGA=FREQ(IFR)/HBAREV1
                  PHI=TPHASE*OMEGA

                  DMODU=CDEXP(DCMPLX(0.D0,PHI))
                  DMODU0=DMODU

                endif !iampread.gt.0

                IOBFR=IO+NOBSV*(IFR-1)

                AX=AFREQ(1,IOBFR)
                AY=AFREQ(2,IOBFR)
                AZ=AFREQ(3,IOBFR)
                AX0=AX
                AY0=AY
                AZ0=AZ

                DO IREP=1,IAMPREP

                  IF (AMPRAN.NE.0.D0) THEN
                    PHI=XRAN(IREP)*TWOPI/OMEGAR*OMEGA
                    DDMODU(IREP)=CDEXP(DCMPLX(0.D0,PHI))
                  ENDIF   !AMPRAN

                  IOBFR=IO+NOBSV*(IFR-1)
                  REAIMA(1,1,IOBFR)=REAIMA(1,1,IOBFR)+DREAL(AX)
                  REAIMA(2,1,IOBFR)=REAIMA(2,1,IOBFR)+DREAL(AY)
                  REAIMA(3,1,IOBFR)=REAIMA(3,1,IOBFR)+DREAL(AZ)
                  REAIMA(1,2,IOBFR)=REAIMA(1,2,IOBFR)+DIMAG(AX)
                  REAIMA(2,2,IOBFR)=REAIMA(2,2,IOBFR)+DIMAG(AY)
                  REAIMA(3,2,IOBFR)=REAIMA(3,2,IOBFR)+DIMAG(AZ)

                  IF (AMPRAN.NE.0.D0) THEN

                    DMODU=DMODU0*DDMODU(IREP)

                    AX0=AX0*DMODU0
                    AY0=AY0*DMODU0
                    AZ0=AZ0*DMODU0
                    AX=AX0*DDMODU(IREP)
                    AY=AY0*DDMODU(IREP)
                    AZ=AZ0*DDMODU(IREP)

                  ELSE    !(AMPRAN.NE.0.D0)

                    AX=AX*DMODU
                    AY=AY*DMODU
                    AZ=AZ*DMODU

                  ENDIF   !(AMPRAN.NE.0.D0)

                ENDDO   !IAMPREP

              ENDDO !IFR=1,NFREQ
            ENDDO !IO=1,NOBSV

          else !if (iamprep.gt.0) then

            DO IREP=1,nelec

              if (iampcoh.eq.0) then
                CALL RNORML(xranmar,1,rr)
c                call ranmar(xranmar,1)
                tphase=ampcohsig*xranmar(1)/clight1
              else if (iampcoh.eq.2.or.iampcoh.eq.3.or.iampcoh.eq.4) then
                call bunch(tphase)
              else if (iampcoh.eq.-1) then
                call ubunch(xub,yub,zub,ypub,zpub,gammaub,tphase)
              else
                write(lungfo,*)' *** Error in ADDAMPLI: Bad value for IAMPCOH!'
                write(6,*)' *** Error in ADDAMPLI: Bad value for IAMPCOH!'
                stop '*** Program WAVE aborted ***'
              endif

              if (irep.le.10000) then
                wbuff(irep)=tphase*clight1*1.0e9
              endif

              DO IFR=1,NFREQ
                DO IO=1,NOBSV

                  IOBFR=IO+NOBSV*(IFR-1)

                  AX0=AFREQ(1,IOBFR)
                  AY0=AFREQ(2,IOBFR)
                  AZ0=AFREQ(3,IOBFR)

                  OMEGA=FREQ(IFR)/HBAREV1
                  PHI=TPHASE*OMEGA

                  DMODU=CDEXP(DCMPLX(0.D0,PHI))

                  !(AX0,AY0,AZ0)*DMODU since we deal with absolute phases here
                  AX=AX0*DMODU
                  AY=AY0*DMODU
                  AZ=AZ0*DMODU

                  IOBFR=IO+NOBSV*(IFR-1)

                  REAIMA(1,1,IOBFR)=REAIMA(1,1,IOBFR)+DREAL(AX)
                  REAIMA(2,1,IOBFR)=REAIMA(2,1,IOBFR)+DREAL(AY)
                  REAIMA(3,1,IOBFR)=REAIMA(3,1,IOBFR)+DREAL(AZ)
                  REAIMA(1,2,IOBFR)=REAIMA(1,2,IOBFR)+DIMAG(AX)
                  REAIMA(2,2,IOBFR)=REAIMA(2,2,IOBFR)+DIMAG(AY)
                  REAIMA(3,2,IOBFR)=REAIMA(3,2,IOBFR)+DIMAG(AZ)

                ENDDO !IFR=1,NFREQ
              ENDDO !IO=1,NOBSV

            ENDDO   !IAMPREP

            if (ampbunchcharge.ne.0.0d0) then
              reaima=reaima*sqrt(abs(AMPBUNCHCHARGE)/echarge1/nelec**2)
            else
              reaima=reaima/sqrt(dble(nelec))
            endif

            open(unit=99,file='wave_bunch.dat')
            do irep=1,min(nelec,10000)
              write(99,*)wbuff(irep)
            enddo
            close(99)

          end if !(iamprep.gt.0) then

          DO IFR=1,NFREQ
            DO IO=1,NOBSV

              IOBFR=IO+NOBSV*(IFR-1)
              AXR=REAIMA(1,1,IOBFR)
              AYR=REAIMA(2,1,IOBFR)
              AZR=REAIMA(3,1,IOBFR)
              AXI=REAIMA(1,2,IOBFR)
              AYI=REAIMA(2,2,IOBFR)
              AZI=REAIMA(3,2,IOBFR)


              AX=DCMPLX(AXR,AXI)
              AY=DCMPLX(AYR,AYI)
              AZ=DCMPLX(AZR,AZI)

              IOBFR=IO+NOBSV*(IFR-1)
              AFREQ(1,IOBFR)=AX
              AFREQ(2,IOBFR)=AY
              AFREQ(3,IOBFR)=AZ

            ENDDO   !NOBSV
          ENDDO   !NFREQ

        ENDIF   !(IAMPREP.NE.0)

        DO IFR=1,NFREQ
          DO IO=1,NOBSV

            IOBFR=IO+NOBSV*(IFR-1)

            AXR=DREAL(AFREQ(1,IOBFR))
            AYR=DREAL(AFREQ(2,IOBFR))
            AZR=DREAL(AFREQ(3,IOBFR))

            AXI=DIMAG(AFREQ(1,IOBFR))
            AYI=DIMAG(AFREQ(2,IOBFR))
            AZI=DIMAG(AFREQ(3,IOBFR))

            REAIMA(1,1,IOBFR)=AXR
            REAIMA(2,1,IOBFR)=AYR
            REAIMA(3,1,IOBFR)=AZR
            REAIMA(1,2,IOBFR)=AXI
            REAIMA(2,2,IOBFR)=AYI
            REAIMA(3,2,IOBFR)=AZI


            SPEC(1+NSOURCE*(IO-1+NOBSV*(IFR-1)))=
     &        (AXR*AXR+AXI*AXI
     &        +AYR*AYR+AYI*AYI
     &        +AZR*AZR+AZI*AZI
     &        )*SPECNOR

            IF (ISTOKES.NE.0) THEN

              AX=DCMPLX(AXR,AXI)
              AY=DCMPLX(AYR,AYI)
              AZ=DCMPLX(AZR,AZI)

              APOLH=
     &          AX*CONJG(VSTOKES(1,1))
     &          +AY*CONJG(VSTOKES(1,2))
     &          +AZ*CONJG(VSTOKES(1,3))

              APOLR=
     &          AX*CONJG(VSTOKES(2,1))
     &          +AY*CONJG(VSTOKES(2,2))
     &          +AZ*CONJG(VSTOKES(2,3))

              APOLL=
     &          AX*CONJG(VSTOKES(3,1))
     &          +AY*CONJG(VSTOKES(3,2))
     &          +AZ*CONJG(VSTOKES(3,3))

              APOL45=
     &          AX*CONJG(VSTOKES(4,1))
     &          +AY*CONJG(VSTOKES(4,2))
     &          +AZ*CONJG(VSTOKES(4,3))

              STOK1=
     &          REAL(APOLR*CONJG(APOLR))+
     &          REAL(APOLL*CONJG(APOLL))

              STOK2=-STOK1+
     &          2.*REAL(APOLH*CONJG(APOLH))

              STOK3=
     &          2.*REAL(APOL45*CONJG(APOL45))-
     &          STOK1

              STOK4=
     &          REAL(APOLR*CONJG(APOLR))-
     &          REAL(APOLL*CONJG(APOLL))


              IOBFR=IO+NOBSV*(IFR-1)
              STOKES(1,IOBFR)=STOK1*SPECNOR
              STOKES(2,IOBFR)=STOK2*SPECNOR
              STOKES(3,IOBFR)=STOK3*SPECNOR
              STOKES(4,IOBFR)=STOK4*SPECNOR

            ENDIF !ISTOKES

          ENDDO   !NOBSV
        ENDDO  !NFREQ

        IF (IMAMPLI.EQ.-3) THEN
C--- WRITE FILES FOR PROGRAM PHASE OF JOHANNES BAHRDT
          CALL PHASE_BAHRDT
        ENDIF

C--- SPECTRUM ARRAYS}

      ELSE  !IMAMPLI

        WRITE(LUNGFO,*)
     &    ' *** ERROR IN ADDAMPLI: UNDEFINED MODE REQUESTED ***'
        WRITE(LUNGFO,*)
     &    ' CHECK IMAMPLI IN NAMELIST $AMPLIN IN INPUT FILE'
        WRITE(6,*)
     &    ' *** ERROR IN ADDAMPLI: UNDEFINED MODE REQUESTED ***'
        WRITE(6,*)
     &    ' CHECK IMAMPLI IN NAMELIST $AMPLIN IN INPUT FILE'

        STOP '--- PROGRAM ABORTED ---'

      ENDIF !IMAMPLI

      IF (NSOURCE.NE.1) THEN

        WRITE(LUNGFO,*)
     &    ' *** ERROR IN ADDAMPLI: NUMBER OF SOURCES MUST BE ONE!! ***'

        WRITE(6,*)
     &    ' *** ERROR IN ADDAMPLI: NUMBER OF SOURCES MUST BE ONE!! ***'

        STOP '--- PROGRAM ABORTED ---'

      ENDIF !(NSOURCE.NE.1)

      WRITE(LUNGFO,*)
     &  '      MODE (IMAMPLI, IAMPTERM): ',IMAMPLI,IAMPTERM
      WRITE(LUNGFO,*)
     &  '             IAMPCOMP: ',IAMPCOMP
      WRITE(LUNGFO,*)
     &  '             IAMPREP,IAMPSUP:   ',IAMPREP,IAMPSUP
      WRITE(LUNGFO,*)
     &  '             AMPRAN, IAMPSEED:  '
     &  ,SNGL(AMPRAN),IAMPSEED
      WRITE(LUNGFO,*)
     &  '             AMPFREQ,AMPPHI:    ',SNGL(AMPFREQ),AMPPHI(IAMPREAD)
      WRITE(LUNGFO,*)
     &  '             DETOUR OF e- [m]:  ',DTMOD*CLIGHT1
      WRITE(LUNGFO,*)

      IF (IAMPSUP.EQ.1) THEN

C--- CREATE NEW FILE AND WRITE RESULTING REAIMA TO IT

        ILAST=IGETLASTCHAR(1,128,FILEAMPLI,CLAST)
        FILESUPER=FILEAMPLI(1:ILAST)//'_SUPER'

        OPEN(UNIT=LUNAMPLI,FILE=FILESUPER,STATUS='unknown'
     &    ,FORM='FORMATTED')

        WRITE(LUNAMPLI,*)ICODE,' ',CODE
        WRITE(LUNAMPLI,*)ICMZNOCMPLX
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)NSOURCE,NOBSV,NFREQ,IFREQ2P
        WRITE(LUNAMPLI,*)NOBSVZ,NOBSVY,MOBSVZ,MOBSVY
        WRITE(LUNAMPLI,*)MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY
        WRITE(LUNAMPLI,*)IPIN,IF1DIM,IPINCIRC
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)PINCEN
        WRITE(LUNAMPLI,*)PINW,PINH,PINR
        WRITE(LUNAMPLI,*)OBSVDZ,OBSVDY
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)SPECNOR
        WRITE(LUNAMPLI,*)VPOLA(1)
        WRITE(LUNAMPLI,*)VPOLA(2)
        WRITE(LUNAMPLI,*)VPOLA(3)
        WRITE(LUNAMPLI,*)VSTOKES(1,1)
        WRITE(LUNAMPLI,*)VSTOKES(1,2)
        WRITE(LUNAMPLI,*)VSTOKES(1,3)
        WRITE(LUNAMPLI,*)VSTOKES(2,1)
        WRITE(LUNAMPLI,*)VSTOKES(2,2)
        WRITE(LUNAMPLI,*)VSTOKES(2,3)
        WRITE(LUNAMPLI,*)VSTOKES(3,1)
        WRITE(LUNAMPLI,*)VSTOKES(3,2)
        WRITE(LUNAMPLI,*)VSTOKES(3,3)
        WRITE(LUNAMPLI,*)VSTOKES(4,1)
        WRITE(LUNAMPLI,*)VSTOKES(4,2)
        WRITE(LUNAMPLI,*)VSTOKES(4,3)
        WRITE(LUNAMPLI,*)

        WRITE(LUNAMPLI,*)(OBSVZ(IO),IO=1,NOBSVZ)
        WRITE(LUNAMPLI,*)
        WRITE(LUNAMPLI,*)(OBSVY(IO),IO=1,NOBSVY)
        WRITE(LUNAMPLI,*)

        DO IO=1,NOBSV
          WRITE(LUNAMPLI,*)(OBSV(IXYZ,IO),IXYZ=1,3)
        ENDDO

        WRITE(LUNAMPLI,*)AMPFREQ
        WRITE(LUNAMPLI,*)

        DO IFR=1,NFREQ
          WRITE(LUNAMPLI,*)
          WRITE(LUNAMPLI,*)FREQ(IFR)
          DO IO=1,NOBSV
            IOBFR=IO+NOBSV*(IFR-1)
            WRITE(LUNAMPLI,*)
     &        REAIMA(1,1,IOBFR),REAIMA(1,2,IOBFR)
            WRITE(LUNAMPLI,*)
     &        REAIMA(2,1,IOBFR),REAIMA(2,2,IOBFR)
            WRITE(LUNAMPLI,*)
     &        REAIMA(3,1,IOBFR),REAIMA(3,2,IOBFR)
          ENDDO
        ENDDO
        WRITE(LUNAMPLI,*)

        WRITE(LUNGFO,*)'      array REAIMA written to file'
        WRITE(LUNGFO,*)'      ',FILEAMPLI//'_SUPER'
        WRITE(LUNGFO,*)

        CLOSE(LUNAMPLI)

      ENDIF !IAMPSUP

      IF (IAMPTERM.NE.0) THEN
        IABEND=6
        iroottrees=0
        WRITE(LUNGFO,*)
     &    '*** SR ADDAMPLI: GOING TO TERMINATE WAVE DUE TO FLAG IAMPTERM ***'
        WRITE(6,*)
     &    '*** SR ADDAMPLI: GOING TO TERMINATE WAVE DUE TO FLAG IAMPTERM ***'
      ENDIF

      IF (AMPRAN.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      Rel. rms phase error (from generated errors): '
     &    ,RANRMS
      ENDIF

      RETURN
      END
