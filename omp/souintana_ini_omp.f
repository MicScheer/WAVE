*CMZ :  3.06/00 26/02/2019  11.09.34  by  Michael Scheer
*CMZ :  3.05/04 27/06/2018  13.53.00  by  Michael Scheer
*CMZ :  3.05/02 15/05/2018  09.37.03  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE souintana_ini_omp(ISOUR,IOBSV,INSIDE)

*KEEP,gplhint.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.
      use bunchmod

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,trackf90.
      include 'trackf90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,souintanac.
      include 'souintanac.cmn'
*KEND.

      double precision x1,droix
      integer isour,iobsv,inside,iroi

      CHARACTER(5) CTUP(NTUPP)

      data ctup /'t','x','y','z','rx','ry','rz','rt','p','expr','expi','roi'
     &  ,'iob','ie','yob','zob','bet1n','om','dt','by2','isou'
     &  ,'spec','reax','imax','reay','imay','reaz','imaz','dom1',
     &  'betx','bety','betz','betxp','betyp','betzp','nx','ny','nz'/

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       SUBROUTINE souintana_ini_omp:'
      WRITE(LUNGFO,*)

      IF (NFREQ.GT.NDFREQ) THEN
        WRITE(LUNGFO,*)
     &    '*** ERROR IN souintana_ini_omp: NUMBER OF MAXIMUM PHOTON ENERGIES EXCEEDED'
        WRITE(LUNGFO,*)
     &    'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
     &    '*** ERROR IN souintana_ini_omp: NUMBER OF MAXIMUM PHOTON ENERGIES EXCEEDED'
        WRITE(6,*)
     &    'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF    !(NFREQ.GT.NDFREQP)

      IF (ISPECMODE.EQ.1) THEN
        DTIM00=DTMCO
      ELSE
        DTIM00=DTIM0
      ENDIF

      DTIM01=1.D0/DTIM00

      IWARNbet1n=0
      IWARNROI=0

      C=CLIGHT1
      C1=1.0D0/CLIGHT1

      DOM=(FREQ(2)-FREQ(1))/HBAREV1
      OM=FREQ(1)/HBAREV1
      ZIDOM=ZI*DOM
      ZIOM=ZI*OM
      ZIC=ZI*CLIGHT1

      IF (IWFILINT.LT.0) THEN
        CALL hbookm(NIDSOURCE,'RADIATION INTEGRALS$',NTUPP
     &    , '//WAVE',1024,CTUP)
      ENDIF !(IWFILINT.LT.0)

      jobunch=icbrill
      if (iobunch.ne.-9999) then
        jobunch=iobunch
      endif

      nphsp=nbunch
      nelec=neinbunch*nphsp

c Flux density is normalized to number of electrons per bunch or bunch charge
c and dmycurr. The field is normalized such, that flux dens = ABS(field)**2
      if (ibunch.ne.0.and.bunchcharge.ne.0.0d0) then
        sqnbunch=nbunch
        sqnphsp=sqrt(bunchcharge/echarge1)
     &    *neinbunch
     &    /(bunchcharge/echarge1)
        bunnor=1.0d0/nbunch
      else
        sqnbunch=nbunch
        sqnphsp=sqrt(dble(neinbunch))
        bunnor=1.0d0/nbunch
      endif

      if (iobsv.ne.jobunch) then
        print*,'****************************************************'
        print*,
     &    '*** SEVERE ERROR IN souintana_ini_omp: IOBSV.NE.iobunch FOR FIRST CALL'
        print*,'*** Programm WAVE aborted ***'
        print*,'****************************************************'
        stop
      endif

      X1=xelec

      IF (ISPECMODE.EQ.1) THEN
        XENDSOU=DWX(MCO)    !FINAL X
      ELSE
        XENDSOU=SOURCEEO(1,1,ISOUR)    !FINAL X
      ENDIF

      IF (NROI.LT.0) THEN
        DROIX=(XENDSOU-X1)/(NROIA-1)
        DO IROI=1,NROIA
          ROIX(IROI)=X1+(IROI-1)*DROIX
          ROIP(IROI)=1.0D0
        ENDDO
      ENDIF   !(NROI.LT.0)

      ROIX(1)=ROIX(1)-1.0D-6
      ROIX(NROIA)=ROIX(NROIA)+1.0D-6

      ipoiroi=0
      if (nroia.eq.0) then
        roi(1,1)=-1.0d30
        roi(1,2)=+1.0d30
        roi(2,1)=1.0d0
      endif

      DO IROI=1,NROIA
        roi(1,iroi)=roix(iroi)
        roi(2,iroi)=roip(iroi)
      ENDDO

      RETURN
      END subroutine souintana_ini_omp
