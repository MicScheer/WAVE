*CMZ :  3.07/00 16/03/2019  15.22.09  by  Michael Scheer
*CMZ :  3.03/02 27/02/2017  13.51.45  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  16.22.59  by  Michael Scheer
*CMZ :  2.52/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.35/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.34/09 24/09/2001  12.08.47  by  Michael Scheer
*CMZ :  2.34/01 25/06/2001  14.39.42  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.14.41  by  Michael Scheer
*CMZ :  2.31/00 23/04/2001  18.27.11  by  Michael Scheer
*CMZ :  2.20/09 23/03/2001  11.01.07  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.08.07  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.27  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.41  by  Michael Scheer
*CMZ :  1.00/00 31/07/97  17.36.12  by  Michael Scheer
*CMZ : 00.02/05 18/03/97  15.48.44  by  Michael Scheer
*CMZ : 00.02/04 26/02/97  12.07.43  by  Michael Scheer
*CMZ : 00.02/00 10/12/96  18.07.52  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  12.58.16  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.49.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.04  by  Michael Scheer
*-- Author :
      SUBROUTINE BLENDE_circ_omp(ISOUR,kfreq,ith)

*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

      use circpinmod

C--- INTEGRATES THE SPLINES THAT INTERPOLATE THE INTENSITY INSIDE THE PINHOLE

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


      INTEGER kfreq,ISOUR,IY,IZ,IOBSV,IIY,ICAL,IR,IP,NR,MR,MP,KDUM,IERR
      INTEGER JCAL,IWBLEN,IDUM,ith
      INTEGER IWRPHIS,IWRPHIF,IWSOUR,IWFREQ,ILIOBFR1

      DOUBLE PRECISION DSUM,RPHI,DIA
      DOUBLE PRECISION SUMZ(NDOBSVYP),S2(NDOBSVYP),SUM,OBSVYF(NDOBSVYP)
      DOUBLE PRECISION SUMY(NDOBSVZP)
      DOUBLE PRECISION SUMZP(NDOBSVYP),S2P(NDOBSVYP),SUMP,DSUMP
      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &        ,FPHI(NDOBSVYP)
     &        ,SZ(NDOBSVZP),SY(NDOBSVYP),s(ndobsvzp)
     &        ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y


      DATA ICAL/0/,JCAL/0/,iwblen/0/

      save ical,jcal,iwfreq,iwsour,iwblen

      ierr=0

      if (ipin.eq.3) then
        ILIOBFR=ISOUR+NSOURCE*(kfreq-1)
        WFLUX(ILIOBFR)=SPEC(ILIOBFR)*pinr*twopi1
        print*,"***** IPIN=3!!"
        return
      endif

c      IWBLEN=0

      IF (IRPHI.NE.0) THEN !INTEGRATION WITH RESPECT TO POLAR COORDINATES
        print*,"*** Obsolete option IRPHI not implemented for OMP version of WAVE ***"
        write(lungfo,*)"*** Obsolete option IRPHI not implemented for OMP version of WAVE ***"
        stop
      ENDIF !IRPHI

      WFLUX(ISOUR+NSOURCE*(kfreq-1))=SUM

      IF (irphi.eq.0.and.(IWBLEN.NE.0.OR.IERR.NE.0)) THEN
        IF (JCAL.EQ.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      *** SUBROUTINE blende_circ_omp:'
          WRITE(LUNGFO,*)
     &      '      LINES INDICATED BY * SHOW A RAW ESTIMATE OF ERRORS DUE TO'
          WRITE(LUNGFO,*)
     &      '      SPLINE FAILURE IF REL. ERROR .GT. 1E-5 (FIRST NUMBER IS SOURCE)'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'      *** SUBROUTINE blende_circ_omp:'
          WRITE(6,*)
     &      '      LINES INDICATED BY * SHOW A RAW ESTIMATE OF ERRORS DUE TO'
          WRITE(6,*)
     &      '      SPLINE FAILURE IF REL. ERROR .GT. 1E-5 (FIRST NUMBER IS SOURCE)'
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)
     &      '      source, energy, flux, flux+error, ratio:'
          JCAL=1
        ENDIF   !JCAL

        IF (SUMP.NE.0.D0) THEN
          DSUM=SUM/SUMP
        ELSE
          DSUM=-9999.
        ENDIF

        IF (DABS(DSUM-1.D0).GT.1.D-5) THEN
          WRITE(LUNGFO,*)'*',ISOUR,
     &      SNGL(FREQ(kfreq)),SNGL(SUM),SNGL(SUMP),SNGL(DSUM)
          WRITE(6,*)'*',ISOUR,
     &      SNGL(FREQ(kfreq)),SNGL(SUM),SNGL(SUMP),SNGL(DSUM)
        ENDIF
      ENDIF !IWBLEN

      DO IY=1,NOBSVY
        DO IZ=1,NOBSVZ
          IOBSV=(IY-1)*NOBSVZ+IZ
          IOBFR=IOBSV+NOBSV*(kfreq-1)
          SPECTOT(IOBFR)=SPECTOT(IOBFR)+
     &      SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1)))
        ENDDO   !IZ
      ENDDO   !IY

      WFLUXT(kfreq)=WFLUXT(kfreq)+WFLUX(ISOUR+NSOURCE*(kfreq-1))

      RETURN
      END
