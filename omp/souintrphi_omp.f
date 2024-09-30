*CMZ :          30/09/2024  14.51.49  by  Michael Scheer
*CMZ :  4.01/04 14/11/2023  11.33.48  by  Michael Scheer
*CMZ :  4.01/03 02/06/2023  08.38.01  by  Michael Scheer
*CMZ :  4.00/17 28/11/2022  17.49.12  by  Michael Scheer
*CMZ :  4.00/13 07/11/2021  17.22.27  by  Michael Scheer
*CMZ :  3.07/01 21/03/2019  12.26.41  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  12.15.26  by  Michael Scheer
*CMZ :  3.06/00 28/02/2019  17.18.16  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.12.29  by  Michael Scheer
*CMZ :  3.05/04 27/06/2018  14.00.13  by  Michael Scheer
*CMZ :  3.05/03 22/05/2018  07.27.47  by  Michael Scheer
*-- Author :    Michael Scheer   04/09/2009
      SUBROUTINE SOUINTRPHI_OMP(ISOUR,INSIDE)
*KEEP,GPLHINT.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,workf90u.
      include 'workf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,amplif90u.
      include 'amplif90u.cmn'
*KEND.

      use bunchmod
      use souintmod
      use omp_lib

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE
C---- RESULTS ARE STORE IN AFREQRPHI AND SPECPOWRPHI

         IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,primkin.
      include 'primkin.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

      double precision, dimension (:,:,:), allocatable :: dphsou
      double precision, dimension (:), allocatable :: frq
      integer, dimension (:), allocatable :: iinside,kinside,linside

      COMPLEX*16 ZIOM,ZI,ZIDOM,ZONE,ZICR1,ZIC,baff(3),daff(3)
      COMPLEX*16 EXPOM,DEXPOMPH1,DEXPOMPH,DDEXPOMPH,DEXPOM,EXPOMV2
      COMPLEX*16 DMODU,DMODU0,DDMODU,AX,AY,AZ,AX0,AY0,AZ0,
     &  bx0,by0,bz0
      COMPLEX*16 APOL,APOLH,APOLR,APOLL,APOL45

      DOUBLE PRECISION T0,T1,T2,TENDSOU,X0,X1,X2,X10,Y1,Y2,Z1,Z2,XENDSOU,R0
     &  ,T,DT,DT2,DT0,DTIM00,DTIM01,VXP,VYP,VZP,TENDSOU1
     &  ,R02
c     &  ,H2,H2R2
     &  ,PHI,FREQR,CORRR0,R00,R2,POW
     &  ,X2B,Y2B,Z2B
     &  ,DGAMMA,DGAMSUM,BETA,GAMGAM,GAMGAM0,AMPDT,sqnphsp,sqnbunch,
     &  sqbunnor
     &  ,are(6),aim(6),
     &  xn1,slopein,slope,drn1,drn2,zn1,yn1,wi

      DOUBLE PRECISION STOK1,STOK2,STOK3,STOK4

      DOUBLE PRECISION VX1,VY1,VZ1,BX1,BY1,BZ1
      DOUBLE PRECISION VX2,VY2,VZ2,BX2,BY2,BZ2,AX2D,AY2D,AZ2D
      DOUBLE PRECISION ECDUM,BS,BSQ,ECMAXS,bs1
      DOUBLE PRECISION TS,DPHASE,DPHSOUR(2,2),phase
      DOUBLE PRECISION C1,OM,DOM,GAMMA

      DOUBLE PRECISION BX,BY,BZ,RX,RY,RZ,PX,PY,PZ,RNBX,RNBY,RNBZ
      DOUBLE PRECISION R1,RNX,RNY,RNZ,DOM1,DOM2,BET1N,DUM11,R,BPX,BPY,BPZ
      DOUBLE PRECISION WGANG,OPANG,rn
      double precision br2,rnr2,br4,rnr4,b3,yp2zp2i,
c     &  yp2zp2ia,
     &  f(3),yp(3),ypp,a(3),fdt(3),filo,fihi,dfdt

      DOUBLE PRECISION RARG(5),C

      DOUBLE PRECISION DROIX,DTPHASE,DXEXI,CENXEXI,roi(nroip)
      DOUBLE PRECISION BET1NO,XRPHI,YRPHI,ZRPHI,speck

      double precision fillb(29)

      INTEGER INSIDE
      INTEGER ISOUR,IOBSV,kfreq,JFREQ,IZAEHL,NZAEHL,IX10,I,ICAL,ICOMP
      INTEGER*8
     &  NZAEHL10,MZAEHL,kzaehl,iizaehl,ir1,ir2
      INTEGER ICSPL,IROI,II,IZTOTS,IWARNBET1N,LSTEP,
     &  MCOUNT,NCOUNT,NCOUNT10,N10,ICOUNT,
     &  jobsv,nelec,norad,iwarnwi,jliob,jliobfr,jobfr,
     &  jvelofield,nfrq,ith,kungfo

      INTEGER NTUPP,IC
      PARAMETER (NTUPP=38)
      REAL*8 FILLT(NTUPP)
      CHARACTER(5) CTUP(NTUPP)

      data ctup /'t','x','y','z','rx','ry','rz','rt','p','expr','expi','roi'
     &  ,'iob','ie','yob','zob','bet1n','om','dt','by2','isou'
     &  ,'spec','reax','imax','reay','imay','reaz','imaz','dom1',
     &  'betx','bety','betz','betxp','betyp','betzp','nx','ny','nz'/
     &
      DATA ICAL/0/
      DATA ZI/(0.0D0,1.0D0)/
      DATA ZONE/(1.0D0,0.0D0)/
      DATA IWARNBET1N/0/

      save ical

      allocate(frq(nfreq))
      allocate(kinside(nobsvrphi),linside(nobsvrphi),iinside(nobsvrphi),
     &  dphsou(2,2,max(1,mthreads)))

      linside=0
      iinside=0
      kinside=0

      IF (ICAL.EQ.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       SUBROUTINE SOUINTRPHI_OMP:'
        WRITE(LUNGFO,*)'       (triggered by MPINR)'
        WRITE(LUNGFO,*)

        if(ibunch.ne.0) then
          write(lungfo,*)
          write(lungfo,*)
     &      '*** Warning in SOUINTRPHI_OMP: Use of IBUNCH may result in problems,'
          write(lungfo,*)
     &      '*** since radiation cone might be not centered and lacking cylindrical symmetry '
          write(lungfo,*)
          write(6,*)
          write(6,*)
     &      '*** Warning in SOUINTRPHI_OMP: Use of IBUNCH may result in problems,'
          write(6,*)
     &      '*** since radiation cone might be not centered and lacking cylindrical symmetry '
          write(6,*)
        endif

        IF (NFREQ.GT.NDFREQ) THEN
          WRITE(LUNGFO,*)
     &      '*** ERROR IN SOUINTRPHI_OMP: NUMBER OF MAXIMUM PHOTON ENERGIES EXCEEDED'
          WRITE(LUNGFO,*)
     &      'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
          WRITE(6,*)
     &      '*** ERROR IN SOUINTRPHI_OMP: NUMBER OF MAXIMUM PHOTON ENERGIES EXCEEDED'
          WRITE(6,*)
     &      'INCREASE PARAMETER NDFREQP IN CMPARA.CMN'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
          STOP
        ENDIF    !(NFREQ.GT.NDFREQP)

        DO II=1,NSOURCE
          DO I=1,NROIA
            IWARNROI(I,II)=0
          ENDDO
        ENDDO

        IF (ISPECMODE.EQ.1) THEN
          DTIM00=DTMCO
        ELSE
          DTIM00=DTIM0
        ENDIF

        DTIM01=1.D0/DTIM00

        C=CLIGHT1
        C1=1.D0/CLIGHT1

        DOM=(FREQ(2)-FREQ(1))/HBAREV1
        OM=FREQ(1)/HBAREV1
        ZIDOM=ZI*DOM
        ZIOM=ZI*OM
        ZIC=ZI*CLIGHT1

        IF (IWFILINT.LT.0) THEN
          CALL hbookm(NIDSOURCE,'RADIATION INTEGRALS$',NTUPP
     &      , '//WAVE',1024,CTUP)
        ENDIF !(IWFILINT.LT.0)

        nphsp=nbunch
        nelec=neinbunch*nphsp

c Flux density is normalized to number of electrons per bunch or bunch charge
c and dmycurr. The field is normalized such, that flux dens = ABS(field)**2
        if (ibunch.ne.0.and.bunchcharge.ne.0.0d0) then
          sqnbunch=nbunch
          sqnphsp=sqrt(bunchcharge/echarge1)
     &      *neinbunch
     &      /(bunchcharge/echarge1)
          bunnor=1.0d0/nbunch
        else
          sqnbunch=nbunch
          sqnphsp=sqrt(dble(neinbunch))
          bunnor=1.0d0/nbunch
        endif

        if (ibunphase.eq.0) then
          sqbunnor=sqrt(dfloat(neinbunch))/sqrt(dfloat(nbunch))
        else
          sqbunnor=1.0d0/sqrt(dfloat(nbunch))
        endif

      ENDIF !ICAL

      IF (ielec.eq.1) THEN

        WRITE(LUNGFO,*)'            SOURCE NUMBER',ISOUR,':'
        WRITE(LUNGFO,*)

        do kfreq=1,nfreq
          ampzmax(kfreq)=0.0d0
          azcos(kfreq)=1.0d0
          azsin(kfreq)=0.0d0
        enddo
        kobs=iobunch

        X1=xelec

        IF (NROI.LT.0) THEN
          DROIX=(XENDSOU-X1)/(NROIA-1)
          DO IROI=1,NROIA
            ROIX(IROI)=X1+(IROI-1)*DROIX
            ROIP(IROI)=1.0D0
          ENDDO
        ENDIF   !(NROI.LT.0)

        ROIX(1)=ROIX(1)-1.0D-6
        ROIX(NROIA)=ROIX(NROIA)+1.0D-6

        DO IROI=1,NROIA
          IPOIROI(IROI)=0
          if (ical.eq.0) then
            roi(iroi)=roix(iroi)
          endif
        ENDDO
      ENDIF !ielec.eq.1

      phaserphi=0.0d0
      LSTEP=0
      DGAMSUM=0.0D0

      gamma=egamma
      beta=dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))

c error? 16feb07        WGANG2=(WGWINFC/GAMMA)**2+2.0D0/(GAMMA**2*(1.0D0+DMYBETA))
C ERROR 11jan08        WGANG2=(WGWINFC/GAMMA)**2
      WGANG=WGWINFC/GAMMA

      ICSPL=0

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS     T=-R0*C1
      T=0.0D0 !WICHTIG HIER WEGEN TENDSOU-T WEITER UNTEN

      R0=OBSVRPHI(1,1)-SOURCEAO(1,1,ISOUR)

      IF (ISPECMODE.EQ.1) THEN
        T0=DWT(1)
        T1=T0
        T2=DWT(MCO)
        XENDSOU=DWX(MCO)    !FINAL X
      ELSE
        T0=SOURCET(1,ISOUR)
        T1=T0
        T2=SOURCET(2,ISOUR)
        XENDSOU=SOURCEEO(1,1,ISOUR)    !FINAL X
      ENDIF

      TENDSOU=T2-T1

      X1=xelec
      Y1=yelec
      Z1=zelec

      VX1=vxelec
      VY1=vyelec
      VZ1=vzelec

      BX1=SOURCEAO(1,4,ISOUR)
      BY1=SOURCEAO(2,4,ISOUR)
      BZ1=SOURCEAO(3,4,ISOUR)
      BS1=SQRT(BX1**2+BY1**2+BZ1**2)

      IZTOTS=0

      X0=X1
      X2=X1
      X10=(XENDSOU-X0)/10.1D0

      NZAEHL=NLPOIO
      DT0=TENDSOU/NZAEHL

      DT=DT0

      if (ielec.eq.1) then

        KZAEHL=0

        IR1=-1
        DO IROI=1,NROIA
          IF (ROIX(IROI).GT.X1.AND.ROIX(IROI).LT.XENDSOU.AND.IR1.EQ.-1) THEN
            IR1=IROI
            GOTO 11
          ENDIF
        ENDDO

11      DO IROI=1,NROIA
          IR2=IROI
          IF (roi(IROI).GT.XENDSOU) THEN
            roi(IROI)=XENDSOU
            IR2=IR2-1
            IF (roi(IR2).LT.X1) THEN
              roi(IR2)=X1
            ENDIF
            GOTO 12
          ENDIF
        ENDDO

12      CONTINUE

        KZAEHL=KZAEHL+NZAEHL*ROIP(IR2)*(XENDSOU-roi(IR2))/(XENDSOU-X1)

        IF (IR1.NE.-1) THEN

          KZAEHL=KZAEHL+NZAEHL*ROIP(IR1-1)*(roi(IR1)-X1)/(XENDSOU-X1)

          DO IROI=IR1,IR2-1
            IF (roi(IROI).GT.X1.OR.roi(IROI)+1.LT.XENDSOU) THEN
              KZAEHL=KZAEHL+NZAEHL*ROIP(IROI)*(roi(IROI+1)-roi(IROI))/(XENDSOU-X1)
            ELSE IF (roi(IROI).GT.X1.OR.roi(IROI)+1.LT.XENDSOU) THEN
              KZAEHL=KZAEHL+NZAEHL*ROIP(IROI)*(roi(IROI+1)-roi(IROI))/(XENDSOU-X1)
            ENDIF
          ENDDO

        ENDIF

      endif !ielec.eq.1

      IF (X1.LT.ROIX(1).OR.XENDSOU.GT.ROIX(NROIA)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN SOUINTRPHI_OMP: X OUTSIDE ROIS ***'
        WRITE(LUNGFO,*)'CHECK NAMELIST $ROIN'
        WRITE(LUNGFO,*)' *** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN SOUINTRPHI_OMP: X OUTSIDE ROIS ***'
        WRITE(6,*)'CHECK NAMELIST $ROIN'
        WRITE(6,*)' *** PROGRAM WAVE ABORTED ***'
        STOP
      ENDIF   !IROI

      X2=X1
      Y2=Y1
      Z2=Z1

      VX2=VX1
      VY2=VY1
      VZ2=VZ1

      BX2=BX1
      BY2=BY1
      BZ2=BZ1
      BS=BS1

C--- LOOP OVER STEPS

      DO IROI=1,NROIA
        IPOIROI(IROI)=0
      ENDDO

      IROI=1
      DO I=1,NROIA
        IF (X1.GE.ROIX(I)) THEN
          IROI=I
        ENDIF !(X1.GE.ROIX(I))
      ENDDO   !IROI

      DT=DT0/ROIP(IROI)

      NZAEHL=MAX(5,NINT((TENDSOU-T)/DT))
      DT=(TENDSOU-T)/NZAEHL

      TENDSOU1=TENDSOU-DT
      DT2=DT/2.D0

C- CHECK STEPS SIZE

      IF (IWARNROI(IROI,ISOUR).EQ.0) THEN
        IF (DT.GT.DTIM00) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** WARNING IN SOUINTRPHI_OMP, SOURCE, ROI:',ISOUR,IROI
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
          WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
          WRITE(6,*)
          WRITE(6,*)
     &      '*** WARNING IN SOUINTRPHI_OMP, SOURCE, ROI:',ISOUR,IROI
          WRITE(6,*)
          WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
          WRITE(6,*)'SIZE FOR TRAJECTORY!'
          WRITE(6,*)
          WRITE(6,*)
     &      'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
          WRITE(6,*)
          IWARNROI(IROI,ISOUR)=1
        ENDIF !DT
      ENDIF !IWARNROI

      IROI=IROI+1

      IZAEHL=0 !LOOP COUNTER

      nutrack=ielec
      nustep=izaehl

      if (ielec.eq.1) then
        iizaehl=0 !total number of steps in souintana
        NZAEHL10=KZAEHL*nelec*nobsvrphi/10
        MZAEHL=NZAEHL10
        IX10=1
      endif

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS     T=-R0*C1

      T=-DT
      TS=-DT

      expom1rphi=ZONE
      DEXPOMPH1=ZONE

      IF (ifreq2p.EQ.0) THEN
        DO JFREQ=1,NFREQ
          EXPOM2P0(1,JFREQ)=ZONE
        ENDDO
      ENDIF

      afferphi=(0.0D0,0.0D0)
      yp2zp2i=0.0d0
c      yp2zp2ia=0.0d0
      f=0.0d0

1000  IZAEHL=IZAEHL+1

      nustep=izaehl

      IF (ISOUR.eq.1.and.IIZAEHL.GE.MZAEHL) THEN
        CALL date_and_time(dtday,dttime,dtzone,idatetime)
        WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
        IX10=IX10+1
        if (ix10.eq.10) then
          mzaehl=NZAEHL10*9.9
        else
          MZAEHL=MZAEHL+NZAEHL10
        endif
      ENDIF

      IF (IROI.LE.NROIA) THEN

        IF (X2.GE.ROIX(IROI)) THEN

          DT=DT0/ROIP(IROI)
          NZAEHL=NINT((TENDSOU-T)/DT)

          IF (ISPECMODE.EQ.1) THEN
            DT=(TENDSOU-T)/(NZAEHL-1)
          ELSE
            DT=(TENDSOU-T)/NZAEHL
          ENDIF

          TENDSOU1=TENDSOU-DT

          DT2=DT/2.D0

          IF (IWARNROI(IROI,ISOUR).EQ.0) THEN

            IF (DT.GT.DTIM00) THEN

              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &          '*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &          'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
              WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &          'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
              WRITE(6,*)
              WRITE(6,*)
     &          '*** WARNING IN SOUINTWIGS, SOURCE, ROI:',ISOUR,IROI
              WRITE(6,*)
              WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
              WRITE(6,*)'SIZE FOR TRAJECTORY!'
              WRITE(6,*)
              WRITE(6,*)
     &          'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
              WRITE(6,*)

              IWARNROI(IROI,ISOUR)=1

            ENDIF !DT

          ENDIF !IWARNROI

          IROI=IROI+1

        ENDIF   !X2

      ENDIF   !IROI

      IPOIROI(IROI)=IPOIROI(IROI)+1

      T=T+DT

      IF (LSTEP.EQ.1) THEN

        IF (X2.LE.XENDSOU) THEN

          DT=(MIN(XENDSOU,XIEND)-X2)/VX2
          DT2=DT/2.0D0

        ELSE

          TS=TS-DT
          T=T-DT

          DT=(MIN(XENDSOU,XIEND)-X1)/VX2
          DT2=DT/2.0D0

          X2=X1
          Y2=Y1
          Z2=Z1

          VX2=VX1
          VY2=VY1
          VZ2=VZ1

          BX2=BX1
          BY2=BY1
          BZ2=BZ1

        ENDIF

      ENDIF

      X1=X2
      Y1=Y2
      Z1=Z2

      VX1=VX2
      VY1=VY2
      VZ1=VZ2

      BX1=BX2
      BY1=BY2
      BZ1=BZ2
      BS1=BS

      IF (ISPECMODE.NE.1) THEN

C GET MAGNETIC FIELD {

        X2B=X1+VX1*DT2
        Y2B=Y1+VY1*DT2
        Z2B=Z1+VZ1*DT2
        norad=0
        if (ibmasksp.ne.0) then
          ibmasksp=-abs(ibmasksp)
          call mybfeld(x2b,y2b,z2b,bx2,by2,bz2,ax2d,ay2d,az2d)
          if ((bx2**2+by2**2+bz2**2).ne.0.0d0) then
            norad=1
          endif
          ibmasksp=-ibmasksp
        endif

        X2=WSOU(1,1,IZAEHL)
        Y2=WSOU(2,1,IZAEHL)
        Z2=WSOU(3,1,IZAEHL)

        VX2=WSOU(1,2,IZAEHL)
        VY2=WSOU(2,2,IZAEHL)
        VZ2=WSOU(3,2,IZAEHL)

        VXP=WSOU(1,3,IZAEHL)
        VYP=WSOU(2,3,IZAEHL)
        VZP=WSOU(3,3,IZAEHL)

        DT=   wsou(1,4,IZAEHL)
        BETA= wsou(2,4,IZAEHL)
        GAMMA=wsou(3,4,IZAEHL)

        bX2=WSOU(1,5,IZAEHL)
        bY2=WSOU(2,5,IZAEHL)
        bZ2=WSOU(3,5,IZAEHL)

        BX=VX2*C1
        BY=VY2*C1
        BZ=VZ2*C1

        BPX=VXP*C1
        BPY=VYP*C1
        BPZ=VZP*C1

C MOVE ONE STEP }

      ELSE  !ISPECMODE

        CALL WAVE_TRACK_INTER(TS,X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,BS,ICSPL,
     &    GAMMA)

        norad=0
        if (ibmasksp.ne.0) then
          ibmasksp=-abs(ibmasksp)
          call mybfeld(x2b,y2b,z2b,bx2,by2,bz2,ax2d,ay2d,az2d)
          if ((bx2**2+by2**2+bz2**2).ne.0.0d0) then
            norad=1
          endif
          ibmasksp=-ibmasksp
        endif

        IF (IENELOSS.NE.0) THEN
          BETA=DSQRT((1.0D0-1.0D0/GAMMA)*(1.D0+1.0D0/GAMMA))
        ENDIF

        BSQ=BS*BS
        BY2=BSQ

        BX=VX2*C1
        BY=VY2*C1
        BZ=VZ2*C1

        BPX=VXP*C1
        BPY=VYP*C1
        BPZ=VZP*C1

      ENDIF !ISPECMODE

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

      jvelofield=ivelofield
      nfrq=nfreq
      frq=freq(1:nfreq)
      kungfo=lungfo

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& firstprivate(r0,x1,y1,z1,c,c1,om,dom,x2,y2,z2,t,zic,beta,gamma,pi1)
!$OMP& firstprivate(bx,by,bz,bpx,bpy,bpz)
!$OMP& firstprivate(WGANG,nelec,zone,zi,HBAREV1,ZIDOM)
!$OMP& firstprivate(specnor,bunnor,sqnphsp,dexpbunch)
!$OMP& firstprivate(dt,iroi,nidsource,jvelofield,nfrq,kungfo)
!$OMP& SHARED(nobsvrphi,iobunch,isour,ielec,izaehl,nsource)
!$OMP& SHARED(iinside,linside,kinside)
!$OMP& SHARED(phaserphi,expom1rphi,SPECPOWRPHI,afferphi,dphsou)
!$OMP& SHARED(norad,frq,ifreq2p,JWFILINT,IWFILINT,xiend,xianf)
!$OMP& SHARED(IWARNBET1N,BET1NO,obsvrphi)

      ith=OMP_GET_THREAD_NUM()+1

!$OMP DO

      DO JOBSV=1,NOBSVRPHI

        if (jobsv.eq.1.and.iobunch.ne.-9999) then
          iobsv=iobunch
        else if (jobsv.eq.iobunch) then
          iobsv=1
        else
          iobsv=jobsv
        endif

        jliob=ISOUR+NSOURCE*(IOBSV-1)

        XRPHI=OBSVRPHI(1,IOBSV)
        YRPHI=OBSVRPHI(2,IOBSV)*SIN(OBSVRPHI(3,IOBSV))
        ZRPHI=OBSVRPHI(2,IOBSV)*COS(OBSVRPHI(3,IOBSV))

        r=sqrt((xrphi-x1)**2+((yrphi-y1)**2+(zrphi-z1)**2))
        PHASE=(r-r0)*c1 ! needed for phase of field amplitude

        if (izaehl.eq.1) then
          phaserphi(iobsv)=phase
          expom1rphi(iobsv)=cdexp(dcmplx(0.0d0,phaserphi(iobsv)*om))
        endif

        RX=XRPHI-X2
        RY=YRPHI-Y2
        RZ=ZRPHI-Z2

        R=SQRT(RX*RX+RY*RY+RZ*RZ)
        R1=1.D0/R
        ZICR1=ZIC*R1

        RNX=RX*R1
        RNY=RY*R1
        RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

        BET1N=(1.0D0-BX*RNX)-BY*RNY-BZ*RNZ
c 20090928{
        br2=by**2+bz**2
        rnr2=rny**2+rnz**2
        b3=beta**3
        br4=br2**2
        rnr4=rnr2**2

        if(br2.lt.1.0d-4.and.rnr2.lt.1.0d-4) then
          bet1n=
     &      1.0d0/(1.0+beta)/gamma**2
     &      +beta*(rnr2/2.0d0
     &      +rnr4/8.0d0)
     &      +(br2/2.0d0
     &      -br2*rnr2/4.0d0
     &      -br2*rnr4/16.0d0)/beta
     &      +b3*br4*(1.0d0/8.0d0
     &      -rnr2/16.0d0
     &      -rnr4/64.0d0)
     &      -by*rny
     &      -bz*rnz
        endif
c }20090928

        OPANG=BX/BETA*RNX+BY/BETA*RNY+BZ/BETA*RNZ

        IF (ABS(OPANG).LE.1.0D0) THEN
          OPANG=ACOS(OPANG)
        ELSE IF (OPANG.GT.1.0D0) THEN
          OPANG=0.0D0
        ELSE
          OPANG=-PI1
        ENDIF

        DUM11=1.0D0/BET1N
        DOM1=1.0D0/(R*BET1N*BET1N)

        IF (IOBSV.EQ.1) THEN
          IF (IZAEHL.EQ.1) THEN
            BET1NO=BET1N
          ELSE IF (iundulator.eq.0.and.(BET1N-BET1NO)/BET1N.GT.0.05.AND.IWARNBET1N.EQ.0) THEN
            WRITE(kungfo,*)
            WRITE(kungfo,*)'*** WARNING IN SOUINTRPHI_OMP  ***'
            WRITE(kungfo,*)'DISCONTINUITY IN INTEGRAND'
            WRITE(kungfo,*)
     &        'Check results carefully, change BMOVECUT, MYINUM, NLPOI etc.'
            WRITE(kungfo,*)
            WRITE(kungfo,*)'ISTEP,X,BET1N,BET1NO:',IZAEHL,SNGL(X1),SNGL(BET1N),SNGL(BET1NO)
            WRITE(kungfo,*)'FURTHER WARNINGS ARE SUPPRESSED!'
            WRITE(kungfo,*)
            WRITE(6,*)
            WRITE(6,*)'*** WARNING IN SOUINTRPHI_OMP  ***'
            WRITE(6,*)'DISCONTINUITY IN INTEGRAND'
            WRITE(6,*)
     &        'Check results carefully, change BMOVECUT, MYINUM, NLPOI etc.'
            WRITE(6,*)
            WRITE(6,*)'ISTEP,X,BET1N,BET1NO:',IZAEHL,SNGL(X1),SNGL(BET1N),SNGL(BET1NO)
            WRITE(6,*)
            WRITE(6,*)'FURTHER WARNINGS ARE SUPPRESSED!'
            WRITE(6,*)
            IWARNBET1N=1
          ENDIF
          BET1NO=BET1N
        ENDIF

        RNBX=RNX-BX
        RNBY=RNY-BY
        RNBZ=RNZ-BZ

        PX=(RNBY*BPZ-RNBZ*BPY)
        PY=(RNBZ*BPX-RNBX*BPZ)
        PZ=(RNBX*BPY-RNBY*BPX)

        IF (jvelofield.EQ.0) THEN !2 WEGEN POWER
          DOM2=C*DOM1*R1/GAMMA**2
          RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
          RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
          RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
        ELSE IF (jvelofield.EQ.1) THEN
          RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
          RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
          RARG(3)=(RNX*PY-RNY*PX)*DOM1
        ELSE IF (jvelofield.LT.0) THEN
          DOM2=C*DOM1*R1/GAMMA**2
          RARG(1)=(RNX-BX)*DOM2
          RARG(2)=(RNY-BY)*DOM2
          RARG(3)=(RNZ-BZ)*DOM2
        ELSE   !jvelofield
          WRITE(6,*)
     &      '*** ERROR IN SOUINTRPHI_OMP: BAD VALUE OF jvelofield  ***'
          WRITE(6,*) '*** PROGRAM WAVE ABORTED  ***'
          STOP
        ENDIF  !jvelofield

        IF (iinside(iobsv).EQ.0.AND.OPANG.LE.WGANG) THEN

          if (iobsv.eq.1) then
            DPHSOU(1,1,ith)=BET1N*DT*frq(1)/HBAREV1
            DPHSOU(1,2,ith)=BET1N*DT*frq(nfrq)/HBAREV1
          endif

          iinside(iobsv)=1
          kinside(iobsv)=1
          linside(iobsv)=linside(iobsv)+1

          IF (linside(iobsv).GT.1) THEN
            WRITE(kungfo,*)
            WRITE(kungfo,*)'*** WARNING IN SOUINTRPHI_OMP  ***'
            WRITE(kungfo,*)'*** SOURCE:',ISOUR
            WRITE(kungfo,*)'STRANGE SOURCE, CONTAINS SEVERAL SOURCES'
            WRITE(kungfo,*)'SOURCE AND OBSERVATION POINT:'
            WRITE(kungfo,*)
     &        ISOUR,OBSVRPHI(1,IOBSV),OBSVRPHI(2,IOBSV),OBSVRPHI(3,IOBSV)
            WRITE(kungfo,*)
     &        'RESULTS OF SPECTRUM CALCULATIONS MAY BE UNRELIABLE'
            WRITE(kungfo,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
            WRITE(6,*)
            WRITE(6,*)'*** WARNING IN SOUINTRPHI_OMP  ***'
            WRITE(6,*)'*** SOURCE:',ISOUR
            WRITE(6,*)'*** STRANGE SOURCE, CONTAINS SEVERAL SOURCES'
            WRITE(6,*)'SOURCE AND OBSERVATION POINT:',
     &        ISOUR,OBSVRPHI(1,IOBSV),OBSVRPHI(2,IOBSV),OBSVRPHI(3,IOBSV)
            WRITE(6,*)
            WRITE(6,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
            WRITE(6,*)'WARNING OF SPECTRUM CALCULATIONS ARE UNRELIABLE'
            linside(iobsv)=linside(iobsv)-1   !SUPRESS LOTS OF WARNINGS
          ENDIF   !linside(ith)
        ELSE IF (iinside(iobsv).EQ.1.AND.OPANG.GT.WGANG) THEN
          iinside(iobsv)=0
        ENDIF   !iinside(iobsv)

        IF (iinside(iobsv).NE.0) THEN

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS      RARG(4)=T+R*C1

          DPHASE=BET1N*DT

          RARG(4)=phaserphi(iobsv)
          RARG(5)=(RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11

          if (norad.ne.0) rarg=0.0d0

C REAL PART OF INTEGRAND }

C COMPLEX PART OF INTEGRAND {

C    ASSUMES FREQ(I+1)=2*FREQ(I)   FOR ifreq2p=2
C    OR FREQ(I+1)=FREQ(I)+DELTA    FOR ifreq2p>2

C--- LOOP OVER ALL FREQUENCES

          kfreq=1
          if (nelec.gt.1) then
            dexpbunch=phexp(kfreq)
          else
            dexpbunch=(1.0d0,0.0d0)
          endif

          IFROB=kfreq+nfrq*(IOBSV-1)

          OM=frq(kfreq)/HBAREV1
          ZIOM=ZI*OM

          EXPOM=expom1rphi(iobsv)
          DEXPOMPH1=EXP(ZIOM*DPHASE)
          DEXPOMPH=DEXPOMPH1

          IF(ifreq2p.GT.2) THEN
            DEXPOM=EXP(ZIDOM*phaserphi(iobsv))
            DDEXPOMPH=EXP(ZIDOM*DPHASE)
          ELSE IF(ifreq2p.EQ.0) THEN
            EXPOM2P0(2,kfreq)=EXP(ZIOM*DPHASE)
            EXPOM=EXPOM2P0(1,kfreq)
          ENDIF   !ifreq2p

          IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN

            SPECPOWRPHI(jliob)=SPECPOWRPHI(jliob)+RARG(5)*DT

            DO ICOMP=1,3
              daff(icomp)=RARG(ICOMP)/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)*DEXPbunch/sqnphsp
              afferphi(icomp,ifrob)=afferphi(icomp,ifrob)+daff(icomp)
            ENDDO   !ICOMP

            baff(1)=(rny*daff(3)-rnz*daff(2))
            baff(2)=(rnz*daff(1)-rnx*daff(3))
            baff(3)=(rnx*daff(2)-rny*daff(1))

            afferphi(4:6,ifrob)=afferphi(4:6,ifrob)+baff(1:3)/clight1

          ENDIF   !XIANF

          IF (IWFILINT.NE.0) THEN
            IF (MOD(IZAEHL,JWFILINT).EQ.0) THEN
              IF (IWFILINT.LT.0) THEN
                FILLT(1)=T
                FILLT(2)=X2
                FILLT(3)=Y2
                FILLT(4)=Z2
                FILLT(5)=RARG(1)
                FILLT(6)=RARG(2)
                FILLT(7)=RARG(3)
                FILLT(8)=RARG(4)
                FILLT(9)=RARG(5)
                FILLT(10)=dREAL(EXPOM)
                FILLT(11)=dIMAG(EXPOM)
                FILLT(12)=IROI-1
                FILLT(13)=IOBSV
                FILLT(14)=kfreq
                FILLT(17)=BET1N
                FILLT(18)=OM
                FILLT(19)=DT
                FILLT(20)=BY2
                FILLT(21)=ISOUR
                FILLT(15)=YRPHI
                FILLT(16)=ZRPHI
                FILLT(22)=
     &            (
     &            REAL(afferphi(1,ifrob))*REAL(afferphi(1,ifrob))
     &            +IMAG(afferphi(1,ifrob))*IMAG(afferphi(1,ifrob))
     &            +REAL(afferphi(2,ifrob))*REAL(afferphi(2,ifrob))
     &            +IMAG(afferphi(2,ifrob))*IMAG(afferphi(2,ifrob))
     &            +REAL(afferphi(3,ifrob))*REAL(afferphi(3,ifrob))
     &            +IMAG(afferphi(3,ifrob))*IMAG(afferphi(3,ifrob))
     &            )*SPECNOR*bunnor
                FILLT(23)=dREAL(afferphi(1,ifrob))*SPECNOR*bunnor
                FILLT(24)=dIMAG(afferphi(1,ifrob))*SPECNOR*bunnor
                FILLT(25)=dREAL(afferphi(2,ifrob))*SPECNOR*bunnor
                FILLT(26)=dIMAG(afferphi(2,ifrob))*SPECNOR*bunnor
                FILLT(27)=dREAL(afferphi(3,ifrob))*SPECNOR*bunnor
                FILLT(28)=dIMAG(afferphi(3,ifrob))*SPECNOR*bunnor
                FILLT(29)=DOM1

                FILLT(30)=bx
                FILLT(31)=by
                FILLT(32)=bz
                FILLT(33)=bpx
                FILLT(34)=bpy
                FILLT(35)=bpz
c                ef(1:3)=real(afferphi(1:3,ifrob))
c                bf(1:3)=real(afferphi(4:6,ifrob))
c                rnx=ef(2)*bf(3)-ef(3)*bf(2)
c                rny=ef(3)*bf(1)-ef(1)*bf(3)
c                rnz=ef(1)*bf(2)-ef(2)*bf(1)
                rnx=real(
     &            afferphi(2,ifrob)*conjg(afferphi(6,ifrob))-
     &            afferphi(3,ifrob)*conjg(afferphi(5,ifrob)))
                rny=real(
     &            afferphi(3,ifrob)*conjg(afferphi(4,ifrob))-
     &            afferphi(1,ifrob)*conjg(afferphi(6,ifrob)))
                rnz=real(
     &            afferphi(1,ifrob)*conjg(afferphi(5,ifrob))-
     &            afferphi(2,ifrob)*conjg(afferphi(4,ifrob)))
                rn=sqrt(rnx**2+rny**2+rnz**2)
                FILLT(36)=rnx/rn
                FILLT(37)=rny/rn
                FILLT(38)=rnz/rn

                CALL hfm(NIDSOURCE,FILLT)

              ELSE IF (ISOUR.EQ.IWFILINT.AND.IOBSV.EQ.1) THEN

                WRITE(LUNINT,*) IZAEHL,kfreq,X2
                WRITE(LUNINT,*) (RARG(1),IC=1,3)
                WRITE(LUNINT,*) RARG(4)*OM,RARG(5)
                WRITE(LUNINT,*)REAL(EXPOM),IMAG(EXPOM)
                WRITE(LUNINT,*)RARG(1)*REAL(EXPOM),RARG(1)*IMAG(EXPOM)
                WRITE(LUNINT,*)RARG(2)*REAL(EXPOM),RARG(2)*IMAG(EXPOM)
                WRITE(LUNINT,*)RARG(3)*REAL(EXPOM),RARG(3)*IMAG(EXPOM)

              ENDIF !IWFILINT.LT.0
            ENDIF !JFILINT
          ENDIF !IWFILINT.NE.0

          DO kfreq=2,nfrq

            IFROB=kfreq+nfrq*(IOBSV-1)

            IF (ifreq2p.GT.2) THEN
              OM=OM+DOM
              EXPOM=EXPOM*DEXPOM
              DEXPOMPH=DEXPOMPH*DDEXPOMPH
            ELSE IF(ifreq2p.EQ.2) THEN
              OM=OM*2.0D0
              EXPOM=EXPOM*EXPOM
              DEXPOMPH=DEXPOMPH*DEXPOMPH
            ELSE IF(ifreq2p.EQ.0) THEN
              OM=frq(kfreq)/HBAREV1
              ZIOM=ZI*OM
              EXPOM2P0(2,kfreq)=EXP(ZIOM*DPHASE)
              EXPOM=EXPOM2P0(1,kfreq)
              DEXPOMPH=EXPOM2P0(2,kfreq)
            ELSE
              OM=frq(kfreq)/HBAREV1
              ZIOM=ZI*OM
              DEXPOMPH=EXP(ZIOM*DPHASE)
            ENDIF

            if (nelec.gt.1) then
              dexpbunch=phexp(kfreq)
            endif

            IF (X2.GE.XIANF.AND.X2.LE.XIEND.and.norad.eq.0) THEN
              EXPOMV2=1.0D0/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)

              DO ICOMP=1,3
                daff(icomp)=RARG(ICOMP)*EXPOMV2*DEXPbunch/sqnphsp
c                print*,izaehl,ifrob,icomp,RARG(ICOMP),EXPOMV2,DEXPbunch,sqnphsp
c                print*,daff
                afferphi(ICOMP,ifrob)=afferphi(ICOMP,ifrob)+daff(icomp)
              ENDDO

              baff(1)=(rny*daff(3)-rnz*daff(2))
              baff(2)=(rnz*daff(1)-rnx*daff(3))
              baff(3)=(rnx*daff(2)-rny*daff(1))

              afferphi(4:6,ifrob)=afferphi(4:6,ifrob)+baff(1:3)/clight1

            ENDIF !XIEND

            IF (IWFILINT.NE.0) THEN
              IF (MOD(IZAEHL,JWFILINT).EQ.0) THEN
                IF (IWFILINT.LT.0) THEN
                  FILLT(1)=T
                  FILLT(2)=X2
                  FILLT(3)=Y2
                  FILLT(4)=Z2
                  FILLT(5)=RARG(1)
                  FILLT(6)=RARG(2)
                  FILLT(7)=RARG(3)
                  FILLT(8)=RARG(4)
                  FILLT(9)=RARG(5)
                  FILLT(10)=dREAL(EXPOM)
                  FILLT(11)=dIMAG(EXPOM)
                  FILLT(12)=IROI-1
                  FILLT(13)=IOBSV
                  FILLT(14)=kfreq
                  FILLT(17)=BET1N
                  FILLT(18)=OM
                  FILLT(19)=DT
                  FILLT(20)=BY2
                  FILLT(21)=ISOUR
                  FILLT(15)=YRPHI
                  FILLT(16)=ZRPHI
                  FILLT(22)=
     &              (
     &              REAL(afferphi(1,ifrob))*REAL(afferphi(1,ifrob))
     &              +IMAG(afferphi(1,ifrob))*IMAG(afferphi(1,ifrob))
     &              +REAL(afferphi(2,ifrob))*REAL(afferphi(2,ifrob))
     &              +IMAG(afferphi(2,ifrob))*IMAG(afferphi(2,ifrob))
     &              +REAL(afferphi(3,ifrob))*REAL(afferphi(3,ifrob))
     &              +IMAG(afferphi(3,ifrob))*IMAG(afferphi(3,ifrob))
     &              )*SPECNOR*bunnor
                  FILLT(23)=dREAL(afferphi(1,ifrob))*SPECNOR*bunnor
                  FILLT(24)=dIMAG(afferphi(1,ifrob))*SPECNOR*bunnor
                  FILLT(25)=dREAL(afferphi(2,ifrob))*SPECNOR*bunnor
                  FILLT(26)=dIMAG(afferphi(2,ifrob))*SPECNOR*bunnor
                  FILLT(27)=dREAL(afferphi(3,ifrob))*SPECNOR*bunnor
                  FILLT(28)=dIMAG(afferphi(3,ifrob))*SPECNOR*bunnor
                  FILLT(29)=DOM1
                  FILLT(30)=bx
                  FILLT(31)=by
                  FILLT(32)=bz
                  FILLT(33)=bpx
                  FILLT(34)=bpy
                  FILLT(35)=bpz
c                ef(1:3)=real(afferphi(1:3,ifrob))
c                bf(1:3)=real(afferphi(4:6,ifrob))
c                rnx=ef(2)*bf(3)-ef(3)*bf(2)
c                rny=ef(3)*bf(1)-ef(1)*bf(3)
c                rnz=ef(1)*bf(2)-ef(2)*bf(1)
                  rnx=real(
     &              afferphi(2,ifrob)*conjg(afferphi(6,ifrob))-
     &              afferphi(3,ifrob)*conjg(afferphi(5,ifrob)))
                  rny=real(
     &              afferphi(3,ifrob)*conjg(afferphi(4,ifrob))-
     &              afferphi(1,ifrob)*conjg(afferphi(6,ifrob)))
                  rnz=real(
     &              afferphi(1,ifrob)*conjg(afferphi(5,ifrob))-
     &              afferphi(2,ifrob)*conjg(afferphi(4,ifrob)))
                  rn=sqrt(rnx**2+rny**2+rnz**2)
                  FILLT(36)=rnx/rn
                  FILLT(37)=rny/rn
                  FILLT(38)=rnz/rn

                  CALL hfm(NIDSOURCE,FILLT)

                ELSE IF (ISOUR.EQ.IWFILINT.AND.IOBSV.EQ.1) THEN

                  WRITE(LUNINT,*) IZAEHL,kfreq,X2
                  WRITE(LUNINT,*) (RARG(1),IC=1,3)
                  WRITE(LUNINT,*) RARG(4)*OM,RARG(5)
                  WRITE(LUNINT,*)REAL(EXPOM),IMAG(EXPOM)
                  WRITE(LUNINT,*)RARG(1)*REAL(EXPOM),RARG(1)*IMAG(EXPOM)
                  WRITE(LUNINT,*)RARG(2)*REAL(EXPOM),RARG(2)*IMAG(EXPOM)
                  WRITE(LUNINT,*)RARG(3)*REAL(EXPOM),RARG(3)*IMAG(EXPOM)

                ENDIF !IWFILINT.LT.0
              ENDIF !JWFILINT
            ENDIF !IWFILINT.NE.0

          ENDDO   !LOOP OVER ALL FREQUENCES
        ENDIF   !iinside(iobsv)

C COMPLEX PART OF INTEGRAND }

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

        phase=phaserphi(iobsv)
        phaserphi(iobsv)=phaserphi(iobsv)+DPHASE
        expom1rphi(iobsv)=expom1rphi(iobsv)*DEXPOMPH1

        IF(ifreq2p.EQ.0) THEN

          DO JFREQ=1,nfrq
            OM=frq(JFREQ)/HBAREV1
            ZIOM=ZI*OM
            EXPOM2P0(1,JFREQ)=EXPOM2P0(1,JFREQ)*EXPOM2P0(2,JFREQ)
          ENDDO
        ENDIF

        IF (iinside(iobsv).NE.0.and.iobsv.eq.1) THEN
          DPHSOU(2,1,ith)=BET1N*DT*frq(1)/HBAREV1
          DPHSOU(2,2,ith)=BET1N*DT*frq(nfrq)/HBAREV1
        ENDIF

      ENDDO !IOBSV=1,NOBSVRPHI

!$OMP END DO
!$OMP END PARALLEL

      IIZAEHL=IIZAEHL+nobsvrphi !total step counter

      TS=TS+DT

C--- END OF LOOP OVER TIME STEPS

c      yp2zp2ia=yp2zp2ia
c     &  +((vy1/vx1)**2+(vy2/vx2)**2+(vz1/vx1)**2+(vz2/vx2)**2)*beta*clight1*dt2

      f(3)=((vy2/vx2)**2+(vz2/vx2)**2)
      fdt(3)=dt

      if (lstep.eq.1) then
        yp(1)=(f(2)-f(1))/fdt(2)
        yp(3)=(f(3)-f(2))/fdt(3)
        yp(2)=(yp(3)+yp(1))/2.0d0
        ypp=(yp(3)-yp(1))/(fdt(2)+fdt(3))*2.0d0
        a(3)=ypp/2.0d0
        a(2)=yp(2)-2.0d0*a(3)*fdt(2)
        a(1)=f(2)-a(2)*fdt(2)-a(3)*fdt(2)**2
        dfdt=fdt(2)+fdt(3)
        fihi=a(1)*dfdt+a(2)/2.0d0*dfdt**2+a(3)/3.0d0*dfdt**3
        dfdt=0.0d0
        filo=a(1)*dfdt+a(2)/2.0d0*dfdt**2+a(3)/3.0d0*dfdt**3
        yp2zp2i=yp2zp2i+fihi-filo
        yp2zp2i=yp2zp2i*beta*clight1
      else if (izaehl.ge.3) then
        yp(1)=(f(2)-f(1))/fdt(2)
        yp(3)=(f(3)-f(2))/fdt(3)
        yp(2)=(yp(3)+yp(1))/2.0d0
        ypp=(yp(3)-yp(1))/(fdt(2)+fdt(3))*2.0d0
        a(3)=ypp/2.0d0
        a(2)=yp(2)-2.0d0*a(3)*fdt(2)
        a(1)=f(2)-a(2)*fdt(2)-a(3)*fdt(2)**2
        dfdt=fdt(2)+fdt(3)
        fihi=a(1)*dfdt+a(2)/2.0d0*dfdt**2+a(3)/3.0d0*dfdt**3
        dfdt=fdt(2)
        filo=a(1)*dfdt+a(2)/2.0d0*dfdt**2+a(3)/3.0d0*dfdt**3
        yp2zp2i=yp2zp2i+fihi-filo
      endif

      f(1)=f(2)
      fdt(1)=fdt(2)
      f(2)=f(3)
      fdt(2)=fdt(3)

c      stop "Ende"
      if (ispecmode.eq.2) then
        if (izaehl.lt.ipoisou(isour)) goto 1000
      else
        IF (X2.LT.XENDSOU-VX2*DT.AND.X2.LT.(XIEND-VX2*DT).AND.LSTEP.EQ.0)
     &    GOTO 1000
        IF (LSTEP.EQ.0) THEN
          LSTEP=1
          GOTO 1000
        ENDIF
      endif

      IF (IAMPLI.LT.0) THEN

        if (nelec.eq.1) then
          print*,' '
          print*,'Starting phase repetition, counting from one to 10 to show progress:'
          print*,' '
        endif

        DXEXI=MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &    -MAX(SOURCEAO(1,1,ISOUR),XIANF)
        if (ampr2corr.eq.-9999.0d0) ampr2corr=dxexi
        CENXEXI=(MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &    +MAX(SOURCEAO(1,1,ISOUR),XIANF))/2.D0
        GAMGAM0=(SOURCEG(1,1,ISOUR)*(egamma/dmygamma))**2
        GAMGAM=(
     &    (SOURCEG(1,1,ISOUR)+SOURCEG(2,2,ISOUR))*(egamma/dmygamma)
     &    )**2

        XRPHI=OBSVRPHI(1,1)

        slopein=sqrt(vyin**2+vzin**2)/vxin
        slope=sqrt(vyelec**2+vzelec**2)/vxelec

        if (myinum.gt.nlpoi/dxexi) then
          WI=(WTRA2IS(ISOUR)
     &      -DXEXI/2.0D0*slopein**2) !wi is detour for on-axis particle
     &      *(dmygamma/egamma)**2
        else
          if (iwarnwi.eq.0) then
            write(lungfo,*)
            write(lungfo,*)'*** Warning in SOUINTANA:'
            write(lungfo,*)'*** MYINUM is rather small with respect to NLPOI'
            write(lungfo,*)'*** Length of trajectories are now calculated by simple'
            write(lungfo,*)'*** integration with SOUINTANA, which might be poor'
            write(lungfo,*)
            write(lungfo,*)
            write(6,*)'*** Warning in SOUINTANA:'
            write(6,*)'*** MYINUM is rather small with respect to NLPOI'
            write(6,*)'*** Length of trajectories are now calculated by simple'
            write(6,*)'*** integration with SOUINTANA, which might be poor'
            write(6,*)
            iwarnwi=1
          endif
          wi=(yp2zp2i/2.0d0
     &      -DXEXI/2.0D0*slopein**2) !wi is detour for on-axis particle
     &      *(dmygamma/egamma)**2
        endif

        xn1=cenxexi
        yn1=(xn1-cenxexi)*vyelec/vxelec
        zn1=(xn1-cenxexi)*vyelec/vxelec

        drn2=(
     &    (yn1+dxexi*vyelec/vxelec)**2+
     &    (zn1+dxexi*vzelec/vxelec)**2
     &    )/
     &    (2.0d0*(xrphi-xn1-dxexi))

        drn1=(
     &    yn1**2+
     &    zn1**2
     &    )/
     &    (2.0d0*(xrphi-xn1))

        DTPHASE=(
     &    WI+DXEXI*(slope**2/2.0d0+1.0d0/(2.0D0*GAMGAM0))
     &    +drn2-drn1)
     &    /CLIGHT1*GAMGAM0/GAMGAM

        AMPDT=AMPSHIFT(1)/CLIGHT1/2.0D0/GAMGAM0
        FREQR=2.0D0*PI1/DTPHASE*HBAREV1

        ICOUNT=0
        NCOUNT=NFREQ*NOBSVRPHI*ABS(IAMPLI)
        NCOUNT10=NCOUNT/10
        MCOUNT=NCOUNT10
        N10=0

      ENDIF !IAMPLI

      DO kfreq=1,NFREQ

        DO IOBSV=1,NOBSVRPHI

          jliobfr=ISOUR+NSOURCE*(IOBSV-1+NOBSVRPHI*(kfreq-1))
          IFROB=kfreq+NFREQ*(IOBSV-1)
          jobfr=IOBSV+NOBSVRPHI*(kfreq-1)

          IF (IAMPLI.LT.0) THEN

            YRPHI=OBSVRPHI(2,IOBSV)*SIN(OBSVRPHI(3,IOBSV))
            ZRPHI=OBSVRPHI(2,IOBSV)*COS(OBSVRPHI(3,IOBSV))

            OM=frq(kfreq)/HBAREV1

            AX0=afferphi(1,ifrob)
            AY0=afferphi(2,ifrob)
            AZ0=afferphi(3,ifrob)

            AX=AX0
            AY=AY0
            AZ=AZ0

            BX0=afferphi(4,ifrob)
            BY0=afferphi(5,ifrob)
            BZ0=afferphi(6,ifrob)

            BX=BX0
            BY=BY0
            BZ=BZ0

            afferphi(1:6,ifrob)=(0.0D0,0.0D0)

            R0=XRPHI-CENXEXI
            R02=R0*R0
            R00=R0

            xn1=cenxexi
            yn1=(xn1-cenxexi)*vyelec/vxelec
            zn1=(xn1-cenxexi)*vzelec/vxelec

            drn2=(
     &        (yn1+dxexi*vyelec/vxelec-yrphi)**2+
     &        (zn1+dxexi*vzelec/vxelec-zrphi)**2
     &        )/
     &        (2.0d0*(xrphi-xn1-dxexi))

            drn1=(
     &        (yn1-yrphi)**2+
     &        (zn1-zrphi)**2
     &        )/
     &        (2.0d0*(xrphi-xn1))

          DTPHASE=(
     &      WI+DXEXI*(slope**2/2.0d0+1.0d0/(2.0D0*GAMGAM0))
     &      +drn2-drn1)
     &      /CLIGHT1*GAMGAM0/GAMGAM
     &      +AMPDT

            PHI=2.D0*PI1*frq(kfreq)*ECHARGE1/HPLANCK1*DTPHASE

            DMODU=EXP(ZI*PHI)
            DMODU0=DMODU
            DDMODU=ZONE

            DO I=1,-IAMPLI

              R0=xrphi+DXEXI/2.D0*(-IAMPLI-2*(I-1)-1)-CENXEXI
              CORRR0=R00/R0
            !corrects for mistake of averaging over 1/r2, if e.g.
            !the repeated device is a long undulator
     &        *(R0/(R0-ampr2corr/2.0d0))**2
              R02=R0*R0

              xn1=cenxexi-dxexi/2.d0*(-iampli-2*(i-1)-1)
     &          *((R0-ampr2corr/2.0d0)/R0)**2 !empirically, due to depth of field
              yn1=(xn1-cenxexi)*vyelec/vxelec
              zn1=(xn1-cenxexi)*vzelec/vxelec

              drn2=(
     &          (yn1+dxexi*vyelec/vxelec-yrphi)**2+
     &          (zn1+dxexi*vzelec/vxelec-zrphi)**2
     &          )/
     &          (2.0d0*(xrphi-xn1-dxexi))

              drn1=(
     &          (yn1-yrphi)**2+
     &          (zn1-zrphi)**2
     &          )/
     &          (2.0d0*(xrphi-xn1))

              DTPHASE=(
     &          WI+DXEXI*(slope**2/2.0d0+1.0d0/(2.0D0*GAMGAM0))
     &          +drn2-drn1)
     &          /CLIGHT1*GAMGAM0/GAMGAM
     &          +AMPDT

              PHI=2.D0*PI1*frq(kfreq)*ECHARGE1/HPLANCK1*DTPHASE

              DMODU=EXP(ZI*PHI)
              DMODU0=DMODU
              DDMODU=ZONE

              afferphi(1,ifrob)=afferphi(1,ifrob)+AX
              afferphi(2,ifrob)=afferphi(2,ifrob)+AY
              afferphi(3,ifrob)=afferphi(3,ifrob)+AZ

              afferphi(4,ifrob)=afferphi(1,ifrob)+BX
              afferphi(5,ifrob)=afferphi(2,ifrob)+BY
              afferphi(6,ifrob)=afferphi(3,ifrob)+BZ

              IF (AMPRAN.NE.0.D0) THEN
                PHI=2.D0*PI1*XRANA(I)/FREQR*frq(kfreq)
                DDMODU=EXP(ZI*PHI)
              ENDIF   !(AMPRAN.NE.0.D0)

              AX0=AX0*DMODU0
              AY0=AY0*DMODU0
              AZ0=AZ0*DMODU0

              AX=AX0*CORRR0
              AY=AY0*CORRR0
              AZ=AZ0*CORRR0

              DMODU=DMODU0*DDMODU
              AX=AX*DMODU
              AY=AY*DMODU
              AZ=AZ*DMODU

              BX0=BX0*DMODU0
              BY0=BY0*DMODU0
              BZ0=BZ0*DMODU0

              BX=BX0*CORRR0
              BY=BY0*CORRR0
              BZ=BZ0*CORRR0

              BX=BX*DMODU
              BY=BY*DMODU
              BZ=BZ*DMODU

              IF (kfreq.EQ.1) THEN
                jliob=ISOUR+NSOURCE*(IOBSV-1)
                IF (I.EQ.1) THEN
                  POW=SPECPOWRPHI(jliob)
                  SPECPOWRPHI(jliob)=0.0D0
                ENDIF !(I.EQ.1) THEN
                R02=(OBSVRPHI(1,IOBSV)-CENXEXI)**2
     &            +OBSVRPHI(2,IOBSV)**2+OBSVRPHI(3,IOBSV)**2
                R2=(OBSVRPHI(1,IOBSV)-CENXEXI-DXEXI*(I-ABS(IAMPLI)/2+1))**2
     &            +OBSVRPHI(2,IOBSV)**2+OBSVRPHI(3,IOBSV)**2
                SPECPOWRPHI(jliob)=SPECPOWRPHI(jliob)+POW*R02/R2
     &            *R2/(sqrt(R2)-ampr2corr/2.0d0)**2/nelec
              ENDIF !kfreq.EQ.1

              ICOUNT=ICOUNT+1
              IF (nelec.eq.1.and.ICOUNT.EQ.MCOUNT) THEN
                N10=N10+1
                CALL date_and_time(dtday,dttime,dtzone,idatetime)
                WRITE(6,*)' ',N10,ICOUNT/(NFREQ*NOBSVRPHI),' ',
     &            dttime(1:2),':',dttime(3:4),':',dttime(5:6)
                MCOUNT=MCOUNT+NCOUNT10
                IF (N10.EQ.9) MCOUNT=NCOUNT
              ENDIF

            ENDDO !IAMPLI

          ENDIF   !(IAMPLI.LT.0)

          if (ihbunch.ne.0.and.iobsv.eq.1) then
            fillb(1)=ibun
            fillb(2)=isub
            fillb(3)=ielec
            fillb(4)=bunchx
            fillb(5)=xelec
            fillb(6)=yelec
            fillb(7)=zelec
            fillb(8)=vyelec/vxelec
            fillb(9)=vzelec/vxelec
            fillb(10)=x2
            fillb(11)=y2
            fillb(12)=z2
            fillb(13)=vy2/vx2
            fillb(14)=vz2/vx2
            fillb(15)=egamma*emassg1
            fillb(16)=gamma*emassg1
            fillb(17)=OBSVRPHI(1,IOBSV)
            fillb(18)=OBSVRPHI(2,IOBSV)*SIN(OBSVRPHI(3,IOBSV))
            fillb(19)=OBSVRPHI(2,IOBSV)*COS(OBSVRPHI(3,IOBSV))
            fillb(20)=kfreq
            fillb(21)=frq(kfreq)
            speck=
     &        DREAL(
     &        afferphi(1,IFROB)*CONJG(afferphi(1,IFROB))
     &        +afferphi(2,IFROB)*CONJG(afferphi(2,IFROB))
     &        +afferphi(3,IFROB)*CONJG(afferphi(3,IFROB))
     &        )*specnor*bunnor
            fillb(22)=speck*nelec

            if (istokes.ne.0) then

              APOLH=
     &          afferphi(1,IFROB)*CONJG(VSTOKES(1,1))
     &          +afferphi(2,IFROB)*CONJG(VSTOKES(1,2))
     &          +afferphi(3,IFROB)*CONJG(VSTOKES(1,3))

              APOLR=
     &          afferphi(1,IFROB)*CONJG(VSTOKES(2,1))
     &          +afferphi(2,IFROB)*CONJG(VSTOKES(2,2))
     &          +afferphi(3,IFROB)*CONJG(VSTOKES(2,3))

              APOLL=
     &          afferphi(1,IFROB)*CONJG(VSTOKES(3,1))
     &          +afferphi(2,IFROB)*CONJG(VSTOKES(3,2))
     &          +afferphi(3,IFROB)*CONJG(VSTOKES(3,3))

              APOL45=
     &          afferphi(1,IFROB)*CONJG(VSTOKES(4,1))
     &          +afferphi(2,IFROB)*CONJG(VSTOKES(4,2))
     &          +afferphi(3,IFROB)*CONJG(VSTOKES(4,3))

              STOK1=
     &          APOLR*CONJG(APOLR)+
     &          APOLL*CONJG(APOLL)

              STOK2=-STOK1+
     &          2.0d0*APOLH*CONJG(APOLH)

              STOK3=
     &          2.0d0*APOL45*CONJG(APOL45)-
     &          STOK1

              STOK4=
     &          APOLR*CONJG(APOLR)-
     &          APOLL*CONJG(APOLL)

              fillb(23)=stok1*specnor*bunnor*nelec
              fillb(24)=stok2*specnor*bunnor*nelec
              fillb(25)=stok3*specnor*bunnor*nelec
              fillb(26)=stok4*specnor*bunnor*nelec

            else
              fillb(23)=fillb(22)
              fillb(24:26)=0.0d0
            endif !istokes

            call hfm(nidbunch,fillb)
            fillb(27)=specpowrphi(isour)*pownor*bunnor*nelec
            fillb(28)=isour
            fillb(29)=t2
          endif ! iobsv=1

          if (
     &        (ibunphase.eq.1.and.ielec.eq.1.and.iobunch.eq.-9999)
     &        .or.
     &        isub.eq.neinbunch
     &        ) then
            are(1:6)=dreal(afferphi(1:6,ifrob))
            aim(1:6)=dimag(afferphi(1:6,ifrob))
            ampz(kfreq)=sqrt(are(3)**2+aim(3)**2)
            if (ampz(kfreq).gt.ampzmax(kfreq)) then
              ampzmax(kfreq)=ampz(kfreq)
              kobs(kfreq)=iobsv
            endif
          endif

        ENDDO !NOBSVRPHI

        if (ibunphase.eq.1
     &    .or.isub.eq.neinbunch
     &    ) then
          ifrob=kfreq+nfreq*(kobs(kfreq)-1)
          are(1:6)=dreal(afferphi(1:6,ifrob))
          aim(1:6)=dimag(afferphi(1:6,ifrob))
          ampz(kfreq)=sqrt(are(3)**2+aim(3)**2)
          if (ampz(kfreq).ne.0.0d0) then
            azcos(kfreq)=are(3)/ampz(kfreq)
            azsin(kfreq)=aim(3)/ampz(kfreq)
          else
            azcos(kfreq)=1.0d0
            azsin(kfreq)=0.0d0
          endif
        endif

        DO IOBSV=1,NOBSVRPHI

          ifrob=kfreq+nfreq*(iobsv-1)

          if (ibunphase.eq.1) then
            are(1:6)=dreal(afferphi(1:6,ifrob))
            aim(1:6)=dimag(afferphi(1:6,ifrob))

            afferphi(1:6,ifrob)=dcmplx(
     &        azcos(kfreq)*are+azsin(kfreq)*aim,
     &        -azsin(kfreq)*are+azcos(kfreq)*aim
     &        )
          endif

          unphrphi(1:6,ifrob)=unphrphi(1:6,ifrob)+afferphi(1:6,ifrob)

          if (isub.eq.neinbunch) then

cold *** This does not work, we really need to call cyltocart for each bunch
cold            if (ibunphase.eq.0) then
cold
coldc destroy phase information between bunches by synchronization
cold
cold
ccold if ibunphase.ne.0, we have done it already
cold              are(1:3)=dreal(unphrphi(1:3,ifrob))
cold              aim(1:3)=dimag(unphrphi(1:3,ifrob))
cold
cold              unphrphi(1:3,ifrob)=dcmplx(
cold     &          azcos(kfreq)*are+azsin(kfreq)*aim,
cold     &          -azsin(kfreq)*are+azcos(kfreq)*aim
cold     &          )
cold            endif

            afreqrphi(1:6,ifrob)=afreqrphi(1:6,ifrob)
     &        +unphrphi(1:6,ifrob)
cold     &        *sqbunnor

            unphrphi(1,IFROB)=(0.0d0,0.0d0)
            unphrphi(2,IFROB)=(0.0d0,0.0d0)
            unphrphi(3,IFROB)=(0.0d0,0.0d0)

          endif !isub.eq.neinbunch

        ENDDO !NOBSVRPHI

      ENDDO !kfreq

cold      if (ielec.ne.nelec) return
      if (isub.ne.neinbunch) return

      call cyltocart(isour)

      do iobsv=1,nobsv

        do kfreq=1,nfreq

          jliobfr=isour+nsource*(iobsv-1+nobsv*(kfreq-1))
          ifrob=kfreq+nfreq*(iobsv-1)
          jobfr=iobsv+nobsv*(kfreq-1)

          om=frq(kfreq)/hbarev1

          if(speccut.gt.0.0d0) then
            if (ispecmode.eq.1) ecmaxs=ecmax(isour)
            if(frq(kfreq).gt.speccut*ecdipev1*dmyenergy**2*ecmaxs) then
              afreq(1:6,ifrob)=(0.0d0,0.0d0)
            endif
          endif

          afreq(1:3,ifrob)=afreq(1:3,ifrob)*reflec(1:3)
          afreq(4:6,ifrob)=afreq(4:6,ifrob)*reflec(1:3)

          if (ipola.eq.0) then

            spec(jliobfr)=spec(jliobfr)+
     &        dreal(
     &        afreq(1,ifrob)*conjg(afreq(1,ifrob))
     &        +afreq(2,ifrob)*conjg(afreq(2,ifrob))
     &        +afreq(3,ifrob)*conjg(afreq(3,ifrob))
     &        )*specnor*bunnor

            reaima(1:3,1,jobfr)=reaima(1:3,1,jobfr)+
     &        dreal(afreq(1:3,ifrob))/sqnbunch
            reaima(1:3,2,jobfr)=reaima(1:3,2,jobfr)+
     &        dimag(afreq(1:3,ifrob))/sqnbunch

            reaima(6:8,1,jobfr)=reaima(6:8,1,jobfr)+
     &        dreal(afreq(4:6,ifrob))/sqnbunch
            reaima(6:8,2,jobfr)=reaima(6:8,2,jobfr)+
     &        dimag(afreq(4:6,ifrob))/sqnbunch

          else    !ipola

            apol=
     &        afreq(1,ifrob)*conjg(vpola(1))
     &        +afreq(2,ifrob)*conjg(vpola(2))
     &        +afreq(3,ifrob)*conjg(vpola(3))

            spec(jliobfr)=spec(jliobfr)+
     &        dreal(apol*conjg(apol))*specnor*bunnor

            reaima(1:3,1,jobfr)=reaima(1:3,1,jobfr)+
     &        dreal(afreq(1:3,ifrob))/sqnbunch
            reaima(1:3,2,jobfr)=reaima(1:3,2,jobfr)+
     &        dimag(afreq(1:3,ifrob))/sqnbunch

            reaima(6:8,1,jobfr)=reaima(6:8,1,jobfr)+
     &        dreal(afreq(4:6,ifrob))/sqnbunch
            reaima(6:8,2,jobfr)=reaima(6:8,2,jobfr)+
     &        dimag(afreq(4:6,ifrob))/sqnbunch

          endif   !ipola

          if (istokes.ne.0) then

            apolh=
     &        afreq(1,ifrob)*conjg(vstokes(1,1))
     &        +afreq(2,ifrob)*conjg(vstokes(1,2))
     &        +afreq(3,ifrob)*conjg(vstokes(1,3))

            apolr=
     &        afreq(1,ifrob)*conjg(vstokes(2,1))
     &        +afreq(2,ifrob)*conjg(vstokes(2,2))
     &        +afreq(3,ifrob)*conjg(vstokes(2,3))

            apoll=
     &        afreq(1,ifrob)*conjg(vstokes(3,1))
     &        +afreq(2,ifrob)*conjg(vstokes(3,2))
     &        +afreq(3,ifrob)*conjg(vstokes(3,3))

            apol45=
     &        afreq(1,ifrob)*conjg(vstokes(4,1))
     &        +afreq(2,ifrob)*conjg(vstokes(4,2))
     &        +afreq(3,ifrob)*conjg(vstokes(4,3))

            stok1=
     &        apolr*conjg(apolr)+
     &        apoll*conjg(apoll)

            stok2=-stok1+
     &        2.*apolh*conjg(apolh)

            stok3=
     &        2.*apol45*conjg(apol45)-
     &        stok1

            stok4=
     &        apolr*conjg(apolr)-
     &        apoll*conjg(apoll)

            stokes(1,jobfr)=stokes(1,jobfr)+
     &        stok1*specnor*bunnor

            stokes(2,jobfr)=stokes(2,jobfr)+
     &        stok2*specnor*bunnor

            stokes(3,jobfr)=stokes(3,jobfr)+
     &        stok3*specnor*bunnor

            stokes(4,jobfr)=stokes(4,jobfr)+
     &        stok4*specnor*bunnor

          endif   !istokes

          afreq(1:6,ifrob)=(0.0d0,0.0d0)

        enddo !kfreq

        jliob=isour+nsource*(iobsv-1)
        specpow(jliob)=specpow(jliob)*pownor

      enddo !nobsv


      if (mpinr.ne.0) then

c only used for Ntuple 5700 so far, 28.3.2012

        do iobsv=1,nobsvrphi

          do kfreq=1,nfreq

            jliobfr=isour+nsource*(iobsv-1+nobsvrphi*(kfreq-1))
            ifrob=kfreq+nfreq*(iobsv-1)
            jobfr=iobsv+nobsvrphi*(kfreq-1)

            om=frq(kfreq)/hbarev1

            if(speccut.gt.0.0d0) then
              if(frq(kfreq).gt.speccut*ecdipev1*dmyenergy**2*ecmaxs) then
                afreqrphi(1:6,ifrob)=(0.0d0,0.0d0)
              endif
            endif

            afreqrphi(1:3,ifrob)=afreqrphi(1,ifrob)*reflec(1:3)
            afreqrphi(4:6,ifrob)=afreqrphi(4:6,ifrob)*reflec(1:3)

            if (ipola.eq.0) then

              specrphi(jliobfr)=specrphi(jliobfr)+
     &          dreal(
     &          afreqrphi(1,ifrob)*conjg(afreqrphi(1,ifrob))
     &          +afreqrphi(2,ifrob)*conjg(afreqrphi(2,ifrob))
     &          +afreqrphi(3,ifrob)*conjg(afreqrphi(3,ifrob))
     &          )*specnor*bunnor

              reaimarphi(1:3,1,jobfr)=reaimarphi(1:3,1,jobfr)+
     &          dreal(afreq(1:3,ifrob))/sqnbunch
              reaimarphi(1:3,2,jobfr)=reaimarphi(1:3,2,jobfr)+
     &          dimag(afreq(1:3,ifrob))/sqnbunch

              reaimarphi(6:8,1,jobfr)=reaimarphi(6:8,1,jobfr)+
     &          dreal(afreq(4:6,ifrob))/sqnbunch
              reaimarphi(6:8,2,jobfr)=reaimarphi(6:8,2,jobfr)+
     &          dimag(afreq(4:6,ifrob))/sqnbunch

            else    !ipola

              apol=
     &          afreqrphi(1,ifrob)*conjg(vpola(1))
     &          +afreqrphi(2,ifrob)*conjg(vpola(2))
     &          +afreqrphi(3,ifrob)*conjg(vpola(3))

              specrphi(jliobfr)=specrphi(jliobfr)+
     &          dreal(apol*conjg(apol))*specnor*bunnor

              reaimarphi(1:3,1,jobfr)=reaimarphi(1:3,1,jobfr)+
     &          dreal(afreq(1:3,ifrob))/sqnbunch
              reaimarphi(1:3,2,jobfr)=reaimarphi(1:3,2,jobfr)+
     &          dimag(afreq(1:3,ifrob))/sqnbunch

              reaimarphi(6:8,1,jobfr)=reaimarphi(6:8,1,jobfr)+
     &          dreal(afreq(4:6,ifrob))/sqnbunch
              reaimarphi(6:8,2,jobfr)=reaimarphi(6:8,2,jobfr)+
     &          dimag(afreq(4:6,ifrob))/sqnbunch

            endif   !ipola

            afreqrphi(1:6,ifrob)=(0.0d0,0.0d0)

          enddo !kfreq

      enddo !nobsvrphi

      endif !(mpinr.ne.0) then

      dphsour=0.0d0
      do ith=1,max(1,mthreads)
        dphsour(:,:)=dphsour(:,:)+dphsou(:,:,ith)
      enddo

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '       phase advance per step at beginning and end of source for'
      WRITE(LUNGFO,*)
     &  '       lowest and highest photon energy at selected observation point:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       beginning:',SNGL(DPHSOUR(1,1)),SNGL(DPHSOUR(1,2))
      WRITE(LUNGFO,*)'       end:      ',SNGL(DPHSOUR(2,1)),SNGL(DPHSOUR(2,2))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ROIs (boundary, precision, points):'
      WRITE(LUNGFO,*)

      DO IROI=1,NROIA-1
        WRITE(LUNGFO,*)
     &    IROI,SNGL(ROIX(IROI)),SNGL(ROIP(IROI)),IPOIROI(IROI+1)
      ENDDO
      WRITE(LUNGFO,*)
     &  NROI,SNGL(ROIX(NROIA))

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       SOURCE, TOTAL NUMBER OF STEPS:',ISOUR,IZAEHL
      WRITE(LUNGFO,*)'       (controlled by NLPOI and namelist $ROIN)'
      WRITE(LUNGFO,*)

      inside=0
      do iobsv=1,nobsvrphi
        if (iinside(iobsv).ne.0.0d0) then
          inside=1
          exit
        endif
      enddo

      deallocate(frq,kinside,dphsou,linside)

      ICAL=1

      RETURN
      END
