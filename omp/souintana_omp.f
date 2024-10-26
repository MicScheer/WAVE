*CMZ :          11/05/2024  10.21.49  by  Michael Scheer
*CMZ :  4.01/05 19/04/2024  15.05.35  by  Michael Scheer
*CMZ :  4.01/04 26/11/2023  16.52.38  by  Michael Scheer
*CMZ :  4.01/03 11/06/2023  11.04.36  by  Michael Scheer
*CMZ :  4.00/17 28/11/2022  17.48.39  by  Michael Scheer
*CMZ :  4.00/15 13/03/2022  19.00.20  by  Michael Scheer
*CMZ :  4.00/13 06/12/2021  13.18.24  by  Michael Scheer
*CMZ :  3.08/01 03/04/2019  15.40.06  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  12.15.25  by  Michael Scheer
*CMZ :  3.06/00 26/02/2019  10.53.54  by  Michael Scheer
*CMZ :  3.05/13 19/09/2018  13.46.31  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.12.29  by  Michael Scheer
*CMZ :  3.05/04 27/06/2018  13.51.56  by  Michael Scheer
*CMZ :  3.05/03 22/05/2018  07.12.30  by  Michael Scheer
*CMZ :  3.05/02 15/05/2018  15.27.15  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE souintana_omp(ISOUR,IOBSV,INSIDE)

*KEEP,gplhint.
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
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

      use bunchmod
      use wbetaf90m
      use souintmod

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE
C---- RESULTS ARE STORE IN AFREQ AND SPECPOW

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
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,tralin.
      include 'tralin.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,specdip.
      include 'specdip.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEEP,souintanac.
      include 'souintanac.cmn'
*KEND.

c      double precision ampzmax(ndfreqp)
c      integer kobs(ndfreqp)
c      common/kobsc/ampzmax,kobs

      REAL*8 FSPEC(31)

      double precision rn_cross_beta(3), rn_cross_rn_cross_beta(3)
      double precision h2,ddist,dist0,dist02
     & ,vn,dgamma

c      COMPLEX*16 ZIOM,ZI,ZIDOM,ZONE,ZICR1,ZIC
      COMPLEX*16 EXPOM1,EXPOM,DEXPOMPH1,DEXPOMPH,DDEXPOMPH,DEXPOM,EXPOMV2,
     &  zicr1,baff(3),daff(3)
      COMPLEX*16 APOL,APOLH,APOLR,APOLL,APOL45
     &  ,DMODU,DMODU0,DDMODU,AX,AY,AZ,AX0,AY0,AZ0,bx0,by0,bz0,bxc,byc,bzc

      double precision apz(nfreq),azc(nfreq),azs(nfreq)

      DOUBLE PRECISION T0,T1,T2,TENDSOU,X0,X1,X2,X10,Y1,Y2,Z1,Z2,R0
     &  ,T,DT,DT2,DT0,VXP,VYP,VZP,TENDSOU1
     &  ,R02
c     &  ,H2,H2R2
     &  ,PHI,FREQR,CORRR0,R00,R2,POW,powpow
     &  ,X2B,Y2B,Z2B
     &  ,DGAMSUM,BETA,GAMGAM,GAMGAM0,AMPDT
     &  ,xn1,slopein,slope,drn1,drn2,zn1,yn1,wi,
     &  zz,yy,zzp,yyp,zzi,yyi,yypi,zzpi,
     &  yeleco,zeleco,zpeleco,ypeleco
c     &  ,ef(3),bf(3)

      DOUBLE PRECISION VX1,VY1,VZ1,BX1,BY1,BZ1
      DOUBLE PRECISION VX2,VY2,VZ2,BX2,BY2,BZ2,AX2D,AY2D,AZ2D
      DOUBLE PRECISION ECDUM,BS,BSQ,BS1
      DOUBLE PRECISION TS,DPHASE,DPHSOUR(2,2)
      DOUBLE PRECISION GAMMA

      DOUBLE PRECISION BX,BY,BZ,RX,RY,RZ,PX,PY,PZ,rn,RNBX,RNBY,RNBZ
      DOUBLE PRECISION R1,RNX,RNY,RNZ,DOM1,DOM2,BET1N,DUM11,R,BPX,BPY,BPZ
      DOUBLE PRECISION WGANG,OPANG

      DOUBLE PRECISION RARG(5),PHASE
      DOUBLE PRECISION DROIX,DTPHASE,DXEXI,CENXEXI
      DOUBLE PRECISION STOK1,STOK2,STOK3,STOK4,BET1NO
      double precision br2,rnr2,br4,rnr4,b3
     &  ,are(6),aim(6),yp2zp2i,robsv,phiobsv,
     &  f(3),yp(3),ypp,a(3),fdt(3),filo,fihi,dfdt,xobsv,yobsv,zobsv,speck,
     &  tfmh(2,2),tfmv(2,2),tfm1(2,2),
     &  w22(2,2),dum22(2,2),
     &  tfmhi(2,2),tfmvi(2,2),
     &  tfmdehi(2,2),tfmdevi(2,2),
     &  tfmdeh(2,2),tfmdev(2,2),
     &  tfmhc(2,2),tfmvc(2,2),
     &  tfmdehc(2,2),tfmdevc(2,2),
     &  tfmhtoti(2,2),tfmvtoti(2,2),
     &  tfmdehtoti(2,2),tfmdevtoti(2,2)
     &  ,rq,cpsi,alpha,spsi,rm,betafun,psi
     &  ,alpha0(2),beta0(2)

      real*8 fillb(41)
      real rnrn(2)
      INTEGER IINSIDE,JINSIDE,INSIDE,iw2(2),ifail
      INTEGER ISOUR,isourold,IOBSV,kfreq,JFREQ,IZAEHL,NZAEHL,I,ICAL,ICOMP
      INTEGER ICSPL,IROI,II,IZTOTS,LSTEP,IR1,IR2
      integer job,jfrob,norad,iwarnwi,kfrob,jliobfr,jliob,jobfr

      INTEGER IC
      INTEGER jpin
      common /souintc/ jpin

      REAL*8 FILLT(NTUPP)
      CHARACTER(5) CTUP(NTUPP)

      DOUBLE PRECISION wth,wta,
     &                H6,H26,A2,A21H6,A3AH26,B,B2,B21H6,B3BH26,DT10
      integer icount,mode,klo,khi,k
      data ctup /'t','x','y','z','rx','ry','rz','rt','p','expr','expi','roi'
     &  ,'iob','ie','yob','zob','bet1n','om','dt','by2','isou'
     &  ,'spec','reax','imax','reay','imay','reaz','imaz','dom1',
     &  'betx','bety','betz','betxp','betyp','betzp','nx','ny','nz'/

      DATA isourold/0/
      DATA ICAL/0/
      DATA tfm1(1,1),tfm1(1,2),tfm1(2,1),tfm1(2,2)/1.0d0,0.0d0,0.0d0,1.0d0/
c      DATA ZI/(0.0D0,1.0D0)/
c      DATA ZONE/(1.0D0,0.0D0)/
c      DATA IWARNBET1N/0/
      ypeleco=vyelec/vxelec
      zpeleco=vzelec/vxelec
      zeleco=zelec
      yeleco=yelec

c Transfermatrices

      if (ibunch.ne.0.and.iampli.lt.0.and.isour.ne.isourold) then

        tfmh=tfm1
        tfmv=tfm1
        tfmhc=tfm1
        tfmvc=tfm1
        tfmhtoti=tfm1
        tfmvtoti=tfm1

        tfmdeh=tfm1
        tfmdev=tfm1
        tfmdehc=tfm1
        tfmdevc=tfm1
        tfmdehtoti=tfm1
        tfmdevtoti=tfm1

        x1=sourceao(1,1,isour)
        y1=sourceao(2,1,isour)
        z1=sourceao(3,1,isour)

        x2=sourceeo(1,1,isour)

        alpha0(1)=-wbetasub(3,1)/2.d0
        alpha0(2)=-wbetasub(5,1)/2.d0
        beta0(1)=wbetasub(2,1)
        beta0(2)=wbetasub(4,1)

        if (alpha0(1).gt.0.001) then
          write(6,*)' '
          write(6,*)'*** Warning in souintana_omp: Derivative of hori. beta function'
          write(6,*)'beginnning of the source is greater than 0.001!'
          write(6,*)' '
        endif

        if (alpha0(2).gt.0.001) then
          write(lungfo,*)' '
          write(lungfo,*)'*** Warning in souintana_omp: Derivative of vert. beta function'
          write(lungfo,*)'in source center is greater than 0.001!'
          write(lungfo,*)'source center:',(x1+x2)/2.0d0
          write(lungfo,*)'beta, alpha:',beta0(2),alpha0(2)
          write(lungfo,*)'Maybe it is useful, to set IBL0CUT'
          write(lungfo,*)' '
          write(6,*)' '
          write(6,*)'*** Warning in souintana_omp: Derivative of vert. beta function'
          write(6,*)'in source center is greater than 0.001!'
          write(6,*)'source center:',(x1+x2)/2.0d0
          write(6,*)'beta, alpha:',beta0(2),alpha0(2)
          write(6,*)'Maybe it is useful, to set IBL0CUT'
          write(6,*)' '
        endif

        alpha=-wbetasub(3,3)/2.0d0
        betafun=  wbetasub(2,3)
        psi=   wbetasub(8,3)

        cpsi=cos(psi)
        spsi=sin(psi)
        rq=sqrt(betafun/beta0(1))
        rm=sqrt(betafun*beta0(1))

        tfmhc(1,1) = rq * (cpsi+alpha0(1)*spsi)
        tfmhc(1,2) = rm * spsi
        tfmhc(2,1)=
     &    ((alpha0(1)-alpha)*cpsi - (1.0d0+alpha0(1)*alpha)*spsi) / rm
        tfmhc(2,2)=
     &    (cpsi-alpha*spsi) / rq

        alpha=-wbetasub(5,3)/2.d0
        betafun=  wbetasub(4,3)
        psi=   wbetasub(9,3)

        cpsi=cos(psi)
        spsi=sin(psi)
        rq=sqrt(betafun/beta0(2))
        rm=sqrt(betafun*beta0(2))

        tfmvc(1,1) = rq * (cpsi+alpha0(2)*spsi)
        tfmvc(1,2) = rm * spsi
        tfmvc(2,1)=
     &    ((alpha0(2)-alpha)*cpsi - (1.0d0+alpha0(2)*alpha)*spsi) / rm
        tfmvc(2,2)=
     &    (cpsi-alpha*spsi) / rq

        alpha=-wbetasub(3,3)/2.d0
        betafun=  wbetasub(2,3)
        psi=   wbetasub(8,3)

        cpsi=cos(psi)
        spsi=sin(psi)
        rq=sqrt(betafun/beta0(1))
        rm=sqrt(betafun*beta0(1))

        tfmh(1,1) = rq * (cpsi+alpha0(1)*spsi)
        tfmh(1,2) = rm * spsi
        tfmh(2,1)=
     &    ((alpha0(1)-alpha)*cpsi - (1.0d0+alpha0(1)*alpha)*spsi) / rm
        tfmh(2,2)=
     &    (cpsi-alpha*spsi) / rq

        alpha=-wbetasub(5,3)/2.d0
        betafun=  wbetasub(4,3)
        psi=   wbetasub(9,3)

        cpsi=cos(psi)
        spsi=sin(psi)
        rq=sqrt(betafun/beta0(2))
        rm=sqrt(betafun*beta0(2))

        tfmv(1,1) = rq * (cpsi+alpha0(2)*spsi)
        tfmv(1,2) = rm * spsi
        tfmv(2,1)=
     &    ((alpha0(2)-alpha)*cpsi - (1.0d0+alpha0(2)*alpha)*spsi) / rm
        tfmv(2,2)=
     &    (cpsi-alpha*spsi) / rq

        w22=tfmh

        dum22(1,1)=1.0d0
        dum22(1,2)=0.0d0
        dum22(2,1)=0.0d0
        dum22(2,2)=1.0d0

        call deqinv(2,w22,2,iw2,ifail,2,dum22)

        if (ifail.ne.0) then
          write(6,*)'*** Error in souintana_omp: Matrix invertation failed'
          write(6,*)'Please check horizontal beta functions.'
          write(lungfo,*)'*** Error in souintana_omp: Matrix invertation failed'
          write(lungfo,*)'Please check horizontal beta functions.'
          stop '*** Program WAVE aborted ***'
        endif

        tfmhi=w22

        w22=tfmv

        dum22(1,1)=1.0d0
        dum22(1,2)=0.0d0
        dum22(2,1)=0.0d0
        dum22(2,2)=1.0d0

        call deqinv(2,w22,2,iw2,ifail,2,dum22)

        if (ifail.ne.0) then
          write(6,*)'*** Error in souintana_omp: Matrix invertation failed'
          write(6,*)'Please check vertical beta functions.'
          write(lungfo,*)'*** Error in souintana_omp: Matrix invertation failed'
          write(lungfo,*)'Please check vertical beta functions.'
          stop '*** Program WAVE aborted ***'
        endif

        tfmvi=w22

        tfmdehi=tfmhi
        tfmdevi=tfmvi

        tfmhtoti=tfmhi
        tfmvtoti=tfmvi

        if (iampli.lt.0) then
          do i=1,-iampli/2-1
            call util_matrix_multiplication(2,2,2,tfmhtoti,tfmhi,tfmhtoti,w22)
            call util_matrix_multiplication(2,2,2,tfmvtoti,tfmvi,tfmvtoti,w22)
          enddo
          tfmdehtoti=tfmhtoti
          tfmdevtoti=tfmhtoti
        endif

      endif !isour

      if (inside.ne.-3.or.ielec.eq.1) then
        xobsv=obsv(1,iobsv)
        yobsv=obsv(2,iobsv)
        zobsv=obsv(3,iobsv)
      else
        call grndmm(rnrn,2)  !s. 39
        xobsv=obsv(1,iobsv)
        yobsv=pincen(2)-pinh/2.0d0+rnrn(1)*pinh
        zobsv=pincen(3)-pinw/2.0d0+rnrn(2)*pinw
        if (ipincirc.ne.0) then
          yobsv=1.0d30
          zobsv=1.0d30
          do while (sqrt(yobsv**2+zobsv**2).gt.pinr)
            call grndmm(rnrn,2)  !s. 39
            yobsv=(rnrn(1)-0.5)*2.0d0*pinr
            zobsv=(rnrn(2)-0.5)*2.0d0*pinr
          enddo
          yobsv=pincen(2)+yobsv
          zobsv=pincen(3)+zobsv
        endif
      endif

      IF (jpin.ne.3.and.ielec.eq.1.and.IOBSV.EQ.jobunch
     &    .or.jpin.eq.3.and.ielec.eq.1) THEN

        WRITE(LUNGFO,*)'            SOURCE NUMBER',ISOUR,':'
        WRITE(LUNGFO,*)

        ampzmax(1:nfreq)=0.0d0
        azcos(1:nfreq)=1.0d0
        azsin(1:nfreq)=0.0d0

        X1=xelec

      ENDIF !IF (ielec.eq.1.and.IOBSV.EQ.jobunch) THEN

      jliob=ISOUR+NSOURCE*(IOBSV-1)
c?6.11.      if (jpin.ne.3.or.jpin.eq.3.and.ielec.eq.1) SPECPOW(jliob)=0.0D0

      LSTEP=0
      DGAMSUM=0.0D0

      gamma=egamma
      beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))

      WGANG=WGWINFC/GAMMA

      ICSPL=0

      if (inside.ne.-3) then
        INSIDE=1
        iinside=-1
      endif
      IINSIDE=0
      JINSIDE=0

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS     T=-R0*C1
      T=0.0D0 !WICHTIG HIER WEGEN TENDSOU-T WEITER UNTEN

      IF (ISPECMODE.EQ.1) THEN
        T0=DWT(1)
        T1=T0
        T2=DWT(MCO)
c        XENDSOU=DWX(MCO)    !FINAL X
      ELSE
        T0=SOURCET(1,ISOUR)
        T1=T0
        T2=SOURCET(2,ISOUR)
c        XENDSOU=SOURCEEO(1,1,ISOUR)    !FINAL X
      ENDIF

      TENDSOU=T2-T1

      IF (X1.LT.roi(1,1).OR.XENDSOU.GT.roi(1,NROIA)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN souintana_omp: X OUTSIDE ROIS ***'
        WRITE(LUNGFO,*)'CHECK NAMELIST $ROIN'
        WRITE(LUNGFO,*)' *** PROGRAM WAVE ABORTED ***'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN souintana_omp: X OUTSIDE ROIS ***'
        WRITE(6,*)'CHECK NAMELIST $ROIN'
        WRITE(6,*)' *** PROGRAM WAVE ABORTED ***'
        STOP
      ENDIF   !IROI

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
c      DT0=TENDSOU/NZAEHL
      DT0=TENDSOU/dble(NZAEHL-1)

      DT=DT0
      ical=ical+1
      !if (ical.gt.1.and.iobsv.gt.nobsvz/2+1) call util_break

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

      IROI=1
      DO I=1,NROIA
        IF (X1.GE.roi(1,I)) THEN
          IROI=I
        ENDIF !(X1.GE.roi(1,I))
      ENDDO   !IROI

      DT=DT0/roi(2,IROI)

      NZAEHL=MAX(5,NINT((TENDSOU-T)/DT))
      DT=(TENDSOU-T)/NZAEHL

      TENDSOU1=TENDSOU-DT
      DT2=DT/2.D0

C- CHECK STEPS SIZE

      IF (IWARNROI(IROI,ISOUR).EQ.0) THEN
        IF (DT-DTIM00.ge.dtim00*0.001) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** WARNING IN souintana_omp, SOURCE, ROI:',ISOUR,IROI
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
          WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
          WRITE(LUNGFO,*)
          write(lungfo,*)'Step size for source point:',dt*clight1
          write(lungfo,*)'Step size for trajectory:',dtim00*clight1
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
          WRITE(6,*)
          WRITE(6,*)
     &      '*** WARNING IN souintana_omp, SOURCE, ROI:',ISOUR,IROI
          WRITE(6,*)
          WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
          WRITE(6,*)'SIZE FOR TRAJECTORY!'
          WRITE(6,*)
          write(6,*)'Step size for source point:',dt*clight1
          write(6,*)'Step size for trajectory:',dtim00*clight1
          WRITE(6,*)
          WRITE(6,*)
     &      'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
          WRITE(6,*)
          IWARNROI(IROI,ISOUR)=1
        ENDIF !DT
      ENDIF !IWARNROI

      IROI=IROI+1

      IZAEHL=0 !LOOP COUNTER for each track

      nutrack=ielec
      nustep=izaehl

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS     T=-R0*C1

      T=-DT
      TS=-DT

c20.11.2023
c20.11.2023      R0=XOBSV-SOURCEAO(1,1,ISOUR)
      R0=obsv(1,icbrill)-SOURCEAO(1,1,ISOUR)
      h2=((yobsv-y1)**2+(zobsv-z1)**2)/(xobsv-x1)**2
      if (h2.lt.0.01) then
        r=abs(xobsv-x1)*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &    -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
      else
        r=sqrt((xobsv-x1)**2+((yobsv-y1)**2+(zobsv-z1)**2))
      endif
c      phase=(obsv(1,icbrill)-x1)*c1
c20.11.2023
      PHASE=(r-r0)*c1 ! needed for phase of field amplitude
      EXPOM1=ZONE
      DEXPOMPH1=ZONE

      IF (ifreq2p.EQ.0) THEN
        DO JFREQ=1,NFREQ
          EXPOM2P0(1,JFREQ)=ZONE
        ENDDO
      ENDIF

      yp2zp2i=0.0D0
      f=0.0d0

      powpow=0.0d0

1000  IZAEHL=IZAEHL+1
      !all util_break
      nustep=izaehl
c      IIZAEHL=IIZAEHL+1 !total step counter

      IF (IROI.LE.NROIA) THEN

        IF (X2.GE.roi(1,IROI)) THEN

          DT=DT0/roi(2,IROI)
          NZAEHL=NINT((TENDSOU-T)/DT)

          IF (ISPECMODE.EQ.1) THEN
            DT=(TENDSOU-T)/(NZAEHL-1)
          ELSE
            DT=(TENDSOU-T)/NZAEHL
          ENDIF

          TENDSOU1=TENDSOU-DT

          DT2=DT/2.D0

          IF (IWARNROI(IROI,ISOUR).EQ.0) THEN

            IF (DT-DTIM00.ge.dtim00*0.001) THEN

              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &          '*** WARNING IN souintana_omp, SOURCE, ROI:',ISOUR,IROI
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &          'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
              WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
              WRITE(LUNGFO,*)
              write(lungfo,*)'Step size for source point:',dt*clight1
              write(lungfo,*)'Step size for trajectory:',dtim00*clight1
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
     &          'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
              WRITE(6,*)
              WRITE(6,*)
     &          '*** WARNING IN souintana_omp, SOURCE, ROI:',ISOUR,IROI
              WRITE(6,*)
              WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
              WRITE(6,*)'SIZE FOR TRAJECTORY!'
              WRITE(6,*)
              write(6,*)'Step size for source point:',dt*clight1
              write(6,*)'Step size for trajectory:',dtim00*clight1
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

      if (ibun.eq.1.and.isub.eq.1) IPOIROI(IROI)=IPOIROI(IROI)+1

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
          BS=BS1

        ENDIF !X2

      ENDIF !LSTEP

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

c{wave_track_inter, inline

c        CALL WAVE_TRACK_INTER(TS,X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,BS,ICSPL,
c     &    GAMMA)

        IF (ICOUNT.EQ.0) THEN
          MODE=0
          DT=(DWT(2)-DWT(1))
          DT10=DT*1.D-10
          DO I=2,MCO
         IF (ABS(DWT(I)-DWT(I-1)-DT).GT.DT10) THEN
           MODE=1
           GOTO 19
                ENDIF
          ENDDO
19        KLO=1
          KHI=MCO
          ICOUNT=1
        ENDIF

      IF (MODE.EQ.1) THEN

        IF (KLO.GE.MCO.OR.KLO.LT.1.OR.KHI.GT.MCO.OR.KHI.LT.2) THEN
          KLO=1
          KHI=MCO
        ENDIF

          IF (T.GE.DWT(KLO).AND.T.LT.DWT(KLO+1)) THEN
            KHI=KLO+1
            GOTO 2
          ELSE IF (T.LT.DWT(KLO).OR.T.GE.DWT(KHI)) THEN
            KLO=1
            KHI=MCO
          ENDIF

          K=1
111       K=K*2
          KHI=KLO+K
          IF (KHI.GE.MCO) GOTO 122
          IF (T.GT.DWT(KHI)) THEN
            KLO=KHI
            GOTO 111
          ELSE
            GOTO 1
          ENDIF

122       KHI=MCO

1         IF (KHI-KLO.GT.1) THEN
            K=(KHI+KLO)/2
            IF(DWT(K).GT.T)THEN
              KHI=K
            ELSE
              KLO=K
            ENDIF
            GOTO 1
          ENDIF

        ELSE !MODE

          IF (T.GE.DWT(1).AND.T.LT.DWT(MCO)) THEN
            KLO=T/DT+1
            KHI=KLO+1
            IF (KHI.GT.MCO) THEN
              KHI=MCO
              KLO=KHI-1
            ENDIF
          ELSE IF (T.LT.DWT(1)) THEN
            KLO=1
            KHI=2
          ELSE IF (T.GE.DWT(MCO)) THEN
            KLO=MCO-1
            KHI=MCO
          ENDIF

        ENDIF !MODE

2       wtH=DWT(KHI)-DWT(KLO)

        IF (wtH.EQ.0.) THEN
          WRITE(6,*) '*** ERROR IN WAVE_TRACK_INTER: BAD INPUT ***'
          STOP
        ENDIF

        H6=wtH/6.D0
        H26=H6*wtH
        wta=(DWT(KHI)-T)/wtH
        A2=wta*wta
        A3AH26=(A2-1.D0)*wta*H26
        A21H6=(-3.D0*A2+1.D0)*H6
        B=(T-DWT(KLO))/wtH
        B2=B*B
        B21H6=(3.D0*B2-1.D0)*H6
        B3BH26=(B2-1.D0)*B*H26

        X2=wta*DWX(KLO)+B*DWX(KHI)+A3AH26*DWX2P(KLO)+B3BH26*DWX2P(KHI)
        VX2=(-DWX(KLO)+DWX(KHI))/wtH+A21H6*DWX2P(KLO)+B21H6*DWX2P(KHI)
        VXP=wta*DWX2P(KLO)+B*DWX2P(KHI)

        Y2=wta*DWY(KLO)+B*DWY(KHI)+A3AH26*DWY2P(KLO)+B3BH26*DWY2P(KHI)
        VY2=(-DWY(KLO)+DWY(KHI))/wtH+A21H6*DWY2P(KLO)+B21H6*DWY2P(KHI)
        VYP=wta*DWY2P(KLO)+B*DWY2P(KHI)

        Z2=wta*DWZ(KLO)+B*DWZ(KHI)+A3AH26*DWZ2P(KLO)+B3BH26*DWZ2P(KHI)
        VZ2=(-DWZ(KLO)+DWZ(KHI))/wtH+A21H6*DWZ2P(KLO)+B21H6*DWZ2P(KHI)
        VZP=wta*DWZ2P(KLO)+B*DWZ2P(KHI)

        BS=wta*DWB(KLO)+B*DWB(KHI)+A3AH26*DWB2P(KLO)+B3BH26*DWB2P(KHI)

        GAMMA=(TRAGAM(KLO)+TRAGAM(KHI))/2.0D0

c}wave_track_inter, inline
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
          BETA=DSQRT((1.0D0-1.0D0/GAMMA)*(1.0D0+1.0D0/GAMMA))
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

      RX=XOBSV-X2
      RY=YOBSV-Y2
      RZ=ZOBSV-Z2

      R=SQRT(RX*RX+RY*RY+RZ*RZ)

      R1=1.0D0/R
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
     &    1.0d0/(1.0+beta)/gamma**2
     &    +beta*(rnr2/2.0d0
     &    +rnr4/8.0d0)
     &    +(br2/2.0d0
     &    -br2*rnr2/4.0d0
     &    -br2*rnr4/16.0d0)/beta
     &    +b3*br4*(1.0d0/8.0d0
     &    -rnr2/16.0d0
     &    -rnr4/64.0d0)
     &    -by*rny
     &    -bz*rnz
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

      DUM11=1.D0/BET1N
      DOM1=1.D0/(R*BET1N*BET1N)

      IF (IZAEHL.EQ.1) THEN
        BET1NO=BET1N
      ELSE IF (iundulator.eq.0.and.(BET1N-BET1NO)/BET1N.GT.0.05.AND.IWARNBET1N.EQ.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN souintana_omp  ***'
        WRITE(LUNGFO,*)'DISCONTINUITY IN INTEGRAND'
        WRITE(LUNGFO,*)
     &    'Check results carefully, change BMOVECUT, MYINUM, NLPOI etc.'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'IELEC,ISTEP,X,BET1N,BET1NO:',IELEC,IZAEHL,SNGL(X1),SNGL(BET1N),SNGL(BET1NO)
        WRITE(LUNGFO,*)'FURTHER WARNINGS ARE SUPPRESSED!'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN souintana_omp  ***'
        WRITE(6,*)'DISCONTINUITY IN INTEGRAND'
        WRITE(6,*)
     &    'Check results carefully, change BMOVECUT, MYINUM, NLPOI etc.'
        WRITE(6,*)
        WRITE(6,*)'IELEC,ISTEP,X,BET1N,BET1NO:',IELEC,IZAEHL,SNGL(X1),SNGL(BET1N),SNGL(BET1NO)
        WRITE(6,*)
        WRITE(6,*)'FURTHER WARNINGS ARE SUPPRESSED!'
        WRITE(6,*)
        IWARNBET1N=1
      ENDIF

      BET1NO=BET1N

      RNBX=RNX-BX
      RNBY=RNY-BY
      RNBZ=RNZ-BZ

      PX=(RNBY*BPZ-RNBZ*BPY)
      PY=(RNBZ*BPX-RNBX*BPZ)
      PZ=(RNBX*BPY-RNBY*BPX)

      IF (IVELOFIELD.EQ.0) THEN !2 WEGEN POWER
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
        RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
      ELSE IF (IVELOFIELD.EQ.1) THEN
        RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
        RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
        RARG(3)=(RNX*PY-RNY*PX)*DOM1
      ELSE IF (IVELOFIELD.LT.0) THEN
        DOM2=C*DOM1*R1/GAMMA**2
        RARG(1)=(RNX-BX)*DOM2
        RARG(2)=(RNY-BY)*DOM2
        RARG(3)=(RNZ-BZ)*DOM2
      ELSE  !IVELOFIELD
        WRITE(6,*)
     &    '*** ERROR IN souintana_omp: BAD VALUE OF IVELOFIELD  ***'
        WRITE(6,*) '*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF !IVELOFIELD

      IF (IINSIDE.EQ.0.AND.OPANG.LE.WGANG) THEN
        if (ibun.eq.1.and.isub.eq.neinbunch) then
          DPHSOUR(1,1)=BET1N*DT*FREQ(1)/HBAREV1
          DPHSOUR(1,2)=BET1N*DT*FREQ(NFREQ)/HBAREV1
        endif
        IINSIDE=1
        INSIDE=1
        JINSIDE=JINSIDE+1
        IF (JINSIDE.GT.1.and.ielec.eq.1) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN souintana_omp  ***'
          WRITE(LUNGFO,*)'*** SOURCE:',ISOUR
          WRITE(LUNGFO,*)'STRANGE SOURCE, CONTAINS SEVERAL SOURCES'
          WRITE(LUNGFO,*)'SOURCE AND OBSERVATION POINT:'
          WRITE(LUNGFO,*)ISOUR,XOBSV,YOBSV,ZOBSV
          WRITE(LUNGFO,*)
     &      'RESULTS OF SPECTRUM CALCULATIONS MAY BE UNRELIABLE'
          WRITE(LUNGFO,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN souintana_omp  ***'
          WRITE(6,*)'*** SOURCE:',ISOUR
          WRITE(6,*)'*** STRANGE SOURCE, CONTAINS SEVERAL SOURCES'
          WRITE(6,*)'SOURCE AND OBSERVATION POINT:'
          WRITE(6,*)ISOUR,XOBSV,YOBSV,ZOBSV
          WRITE(6,*)'*** CHECK COLLIMATOR, PINHOLE, WGWINFC ... ***'
          WRITE(6,*)'WARNING OF SPECTRUM CALCULATIONS ARE UNRELIABLE'
          JINSIDE=JINSIDE-1   !SUPRESS LOTS OF WARNINGS
        ENDIF  !JINSIDE
      ELSE IF (IINSIDE.EQ.1.AND.OPANG.GT.WGANG) THEN
        IINSIDE=0
      ENDIF   !IINSIDE

      IF (IINSIDE.NE.0) THEN

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS      RARG(4)=T+R*C1

        DPHASE=BET1N*DT

        RARG(4)=PHASE
        RARG(5)=(RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11/
     &    (nphsp*neinbunch)

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

        kfrob=kfreq+NFREQ*(IOBSV-1)

        OM=FREQ(kfreq)/HBAREV1
        ZIOM=ZI*OM

        if (izaehl.eq.1) then
          EXPOM1=CDEXP(DCMPLX(0.D0,phase*OM))
        endif

        EXPOM=EXPOM1
        DEXPOMPH1=EXP(ZIOM*DPHASE)
        DEXPOMPH=DEXPOMPH1

        IF(ifreq2p.GT.2) THEN
          DEXPOM=EXP(ZIDOM*PHASE)
          DDEXPOMPH=EXP(ZIDOM*DPHASE)
        ELSE IF(ifreq2p.EQ.0) THEN
          EXPOM2P0(2,kfreq)=EXP(ZIOM*DPHASE)
          EXPOM=EXPOM2P0(1,kfreq)
        ENDIF  !ifreq2p

        !all util_break
        !expom=exp(zidom*)

        IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN

          IF (ECMAXS.LT.BS) ECMAXS=BS

          pow=rarg(5)*dt
          powpow=powpow+pow
          SPECPOW(jliob)=SPECPOW(jliob)+pow

          DO ICOMP=1,3
            daff(icomp)=
     &        RARG(ICOMP)/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)*DEXPbunch/sqnphsp
            affe(icomp,kfrob)=affe(icomp,kfrob)+daff(icomp)
c            if (iobsv.eq.icbrill.and.kmode.eq.1.and.icomp.eq.3) then
c              print*,izaehl,t,affe(3,kfrob)*conjg(affe(3,kfrob))
c            endif
          ENDDO   !ICOMP

c          baff(1)=conjg(rny*daff(3)-rnz*daff(2))
c          baff(2)=conjg(rnz*daff(1)-rnx*daff(3))
c          baff(3)=conjg(rnx*daff(2)-rny*daff(1))

          baff(1)=(rny*daff(3)-rnz*daff(2))
          baff(2)=(rnz*daff(1)-rnx*daff(3))
          baff(3)=(rnx*daff(2)-rny*daff(1))

          !baff(1)=dcmplx(rny*dreal(daff(3))-rnz*dreal(daff(2)),0.0d0)
          !baff(2)=dcmplx(rnz*dreal(daff(1))-rny*dreal(daff(3)),0.0d0)
          !baff(3)=dcmplx(rnx*dreal(daff(2))-rnx*dreal(daff(1)),0.0d0)

          affe(4:6,kfrob)=affe(4:6,kfrob)+baff(1:3)/clight1

        ENDIF  !XIANF

        IF (ibun.eq.1.and.isub.eq.neinbunch) THEN
          IF (IWFILINT.NE.0) THEN
            IF (MOD(IZAEHL,JWFILINT).EQ.0) THEN
              IF (IWFILINT.LT.0) THEN
                FILLT(1)=T+DT
                FILLT(2)=X2
                FILLT(3)=Y2
                FILLT(4)=Z2
                FILLT(5)=RARG(1)
                FILLT(6)=RARG(2)
                FILLT(7)=RARG(3)
                FILLT(8)=RARG(4)
                FILLT(9)=RARG(5)
                FILLT(10)=DREAL(EXPOM*DEXPBUNCH)
                FILLT(11)=DIMAG(EXPOM*DEXPBUNCH)
                FILLT(12)=IROI-1
                FILLT(13)=IOBSV
                FILLT(14)=kfreq
                FILLT(15)=YOBSV
                FILLT(16)=ZOBSV
                if (abs(fillt(15)).lt.1.0d-15) fillt(15)=0.0d0
                if (abs(fillt(16)).lt.1.0d-15) fillt(16)=0.0d0
                FILLT(17)=BET1N
                FILLT(18)=OM
                FILLT(19)=DT
                FILLT(20)=BY2
                FILLT(21)=ISOUR
                FILLT(22)=
     &            (
     &            DREAL(affe(1,kfrob))*DREAL(affe(1,kfrob))
     &            +DIMAG(affe(1,kfrob))*DIMAG(affe(1,kfrob))
     &            +DREAL(affe(2,kfrob))*DREAL(affe(2,kfrob))
     &            +DIMAG(affe(2,kfrob))*DIMAG(affe(2,kfrob))
     &            +DREAL(affe(3,kfrob))*DREAL(affe(3,kfrob))
     &            +DIMAG(affe(3,kfrob))*DIMAG(affe(3,kfrob))
     &            )*specnor*bunnor
                FILLT(23)=DREAL(affe(1,kfrob))*specnor*bunnor
                FILLT(24)=DIMAG(affe(1,kfrob))*specnor*bunnor
                FILLT(25)=DREAL(affe(2,kfrob))*specnor*bunnor
                FILLT(26)=DIMAG(affe(2,kfrob))*specnor*bunnor
                FILLT(27)=DREAL(affe(3,kfrob))*specnor*bunnor
                FILLT(28)=DIMAG(affe(3,kfrob))*specnor*bunnor
                FILLT(29)=DOM1
                FILLT(30)=bx
                FILLT(31)=by
                FILLT(32)=bz
                FILLT(33)=bpx
                FILLT(34)=bpy
                FILLT(35)=bpz
c                ef(1:3)=real(affe(1:3,kfrob))
c                bf(1:3)=real(affe(4:6,kfrob))
c                rnx=ef(2)*bf(3)-ef(3)*bf(2)
c                rny=ef(3)*bf(1)-ef(1)*bf(3)
c                rnz=ef(1)*bf(2)-ef(2)*bf(1)
                rnx=real(
     &            affe(2,kfrob)*conjg(affe(6,kfrob))-
     &            affe(3,kfrob)*conjg(affe(5,kfrob)))
                rny=real(
     &            affe(3,kfrob)*conjg(affe(4,kfrob))-
     &            affe(1,kfrob)*conjg(affe(6,kfrob)))
                rnz=real(
     &            affe(1,kfrob)*conjg(affe(5,kfrob))-
     &            affe(2,kfrob)*conjg(affe(4,kfrob)))
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
                WRITE(LUNINT,*)RARG(1)*REAL(EXPOM*DEXPBUNCH),RARG(1)*IMAG(EXPOM*DEXPBUNCH)
                WRITE(LUNINT,*)RARG(2)*REAL(EXPOM*DEXPBUNCH),RARG(2)*IMAG(EXPOM*DEXPBUNCH)
                WRITE(LUNINT,*)RARG(3)*REAL(EXPOM*DEXPBUNCH),RARG(3)*IMAG(EXPOM*DEXPBUNCH)

              ENDIF !IWFILINT.LT.0
            ENDIF !JFILINT
          ENDIF !IWFILINT.NE.0
        ENDIF !ibun

        DO kfreq=2,NFREQ

          kfrob=kfreq+NFREQ*(IOBSV-1)

          IF    (ifreq2p.GT.2) THEN
            OM=OM+DOM
            EXPOM=EXPOM*DEXPOM
            DEXPOMPH=DEXPOMPH*DDEXPOMPH
          ELSE IF(ifreq2p.EQ.2) THEN
            OM=OM*2.0D0
            EXPOM=EXPOM*EXPOM
            DEXPOMPH=DEXPOMPH*DEXPOMPH
          ELSE IF(ifreq2p.EQ.0) THEN
            OM=FREQ(kfreq)/HBAREV1
            ZIOM=ZI*OM
            EXPOM2P0(2,kfreq)=EXP(ZIOM*DPHASE)
            EXPOM=EXPOM2P0(1,kfreq)
            DEXPOMPH=EXPOM2P0(2,kfreq)
          ELSE
            OM=FREQ(kfreq)/HBAREV1
            ZIOM=ZI*OM
            DEXPOMPH=EXP(ZIOM*DPHASE)
          ENDIF

          if (nelec.gt.1) then
            dexpbunch=phexp(kfreq)
          endif

          IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN

            EXPOMV2=1.0D0/BET1N/OM*EXPOM*(ZONE-DEXPOMPH)

            DO ICOMP=1,3
              daff(icomp)=RARG(ICOMP)*EXPOMV2*DEXPbunch/sqnphsp
              affe(icomp,kfrob)=affe(icomp,kfrob)+daff(icomp)
            ENDDO

c            baff(1)=conjg(rny*daff(3)-rnz*daff(2))
c            baff(2)=conjg(rnz*daff(1)-rnx*daff(3))
c            baff(3)=conjg(rnx*daff(2)-rny*daff(1))

            baff(1)=(rny*daff(3)-rnz*daff(2))
            baff(2)=(rnz*daff(1)-rnx*daff(3))
            baff(3)=(rnx*daff(2)-rny*daff(1))

          !baff(1)=dcmplx(rny*dreal(daff(3))-rnz*dreal(daff(2)),0.0d0)
          !baff(2)=dcmplx(rnz*dreal(daff(1))-rny*dreal(daff(3)),0.0d0)
          !baff(3)=dcmplx(rnx*dreal(daff(2))-rnx*dreal(daff(1)),0.0d0)

            affe(4:6,kfrob)=affe(4:6,kfrob)+baff(1:3)/clight1

          ENDIF !XIEND

          IF (ibun.eq.1.and.isub.eq.neinbunch) then
            IF (IWFILINT.NE.0) THEN
              IF (MOD(IZAEHL,JWFILINT).EQ.0) THEN
                IF (IWFILINT.LT.0) THEN
                  FILLT(1)=T+DT
                  FILLT(2)=X2
                  FILLT(3)=Y2
                  FILLT(4)=Z2
                  FILLT(5)=RARG(1)
                  FILLT(6)=RARG(2)
                  FILLT(7)=RARG(3)
                  FILLT(8)=RARG(4)
                  FILLT(9)=RARG(5)
                  FILLT(10)=DREAL(EXPOM*DEXPBUNCH)
                  FILLT(11)=DIMAG(EXPOM*DEXPBUNCH)
                  FILLT(12)=IROI-1
                  FILLT(13)=IOBSV
                  FILLT(14)=kfreq
                  FILLT(15)=YOBSV
                  FILLT(16)=ZOBSV
                  if (abs(fillt(15)).lt.1.0d-15) fillt(15)=0.0d0
                  if (abs(fillt(16)).lt.1.0d-15) fillt(16)=0.0d0
                  FILLT(17)=BET1N
                  FILLT(18)=OM
                  FILLT(19)=DT
                  FILLT(20)=BY2
                  FILLT(21)=ISOUR
                  FILLT(22)=
     &              (
     &              DREAL(affe(1,kfrob))*DREAL(affe(1,kfrob))
     &              +DIMAG(affe(1,kfrob))*DIMAG(affe(1,kfrob))
     &              +DREAL(affe(2,kfrob))*DREAL(affe(2,kfrob))
     &              +DIMAG(affe(2,kfrob))*DIMAG(affe(2,kfrob))
     &              +DREAL(affe(3,kfrob))*DREAL(affe(3,kfrob))
     &              +DIMAG(affe(3,kfrob))*DIMAG(affe(3,kfrob))
     &              )*specnor*bunnor
                  FILLT(23)=DREAL(affe(1,kfrob))*specnor*bunnor
                  FILLT(24)=DIMAG(affe(1,kfrob))*specnor*bunnor
                  FILLT(25)=DREAL(affe(2,kfrob))*specnor*bunnor
                  FILLT(26)=DIMAG(affe(2,kfrob))*specnor*bunnor
                  FILLT(27)=DREAL(affe(3,kfrob))*specnor*bunnor
                  FILLT(28)=DIMAG(affe(3,kfrob))*specnor*bunnor
                  FILLT(29)=DOM1
                  FILLT(30)=bx
                  FILLT(31)=by
                  FILLT(32)=bz
                  FILLT(33)=bpx
                  FILLT(34)=bpy
                  FILLT(35)=bpz
c                ef(1:3)=real(affe(1:3,kfrob))
c                bf(1:3)=real(affe(4:6,kfrob))
c                rnx=ef(2)*bf(3)-ef(3)*bf(2)
c                rny=ef(3)*bf(1)-ef(1)*bf(3)
c                rnz=ef(1)*bf(2)-ef(2)*bf(1)
                  rnx=real(
     &              affe(2,kfrob)*conjg(affe(6,kfrob))-
     &              affe(3,kfrob)*conjg(affe(5,kfrob)))
                  rny=real(
     &              affe(3,kfrob)*conjg(affe(4,kfrob))-
     &              affe(1,kfrob)*conjg(affe(6,kfrob)))
                  rnz=real(
     &              affe(1,kfrob)*conjg(affe(5,kfrob))-
     &              affe(2,kfrob)*conjg(affe(4,kfrob)))
                  rn=sqrt(rnx**2+rny**2+rnz**2)
                  FILLT(36)=rnx/rn
                  FILLT(37)=rny/rn
                  FILLT(38)=rnz/rn

                  CALL hfm(NIDSOURCE,FILLT)

                ELSE IF (ISOUR.EQ.IWFILINT.AND.IOBSV.EQ.1) THEN

                  WRITE(LUNINT,*) IZAEHL,kfreq,X2
                  WRITE(LUNINT,*) (RARG(1),IC=1,3)
                  WRITE(LUNINT,*) RARG(4)*OM,RARG(5)
                  WRITE(LUNINT,*)REAL(EXPOM*DEXPBUNCH),IMAG(EXPOM*DEXPBUNCH)
                  WRITE(LUNINT,*)RARG(1)*REAL(EXPOM*DEXPBUNCH),RARG(1)*IMAG(EXPOM*DEXPBUNCH)
                  WRITE(LUNINT,*)RARG(2)*REAL(EXPOM*DEXPBUNCH),RARG(2)*IMAG(EXPOM*DEXPBUNCH)
                  WRITE(LUNINT,*)RARG(3)*REAL(EXPOM*DEXPBUNCH),RARG(3)*IMAG(EXPOM*DEXPBUNCH)

                ENDIF !IWFILINT.LT.0
              ENDIF !JWFILINT
            ENDIF !IWFILINT.NE.0
          ENDIF !ibun

        ENDDO   !LOOP OVER ALL FREQUENCES

      ENDIF   !IINSIDE

C COMPLEX PART OF INTEGRAND }

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

      PHASE=PHASE+DPHASE
      EXPOM1=EXPOM1*DEXPOMPH1

      IF(ifreq2p.EQ.0) THEN

        DO JFREQ=1,NFREQ
          OM=FREQ(JFREQ)/HBAREV1
          ZIOM=ZI*OM
          EXPOM2P0(1,JFREQ)=EXPOM2P0(1,JFREQ)*EXPOM2P0(2,JFREQ)
        ENDDO
      ENDIF

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

      if (ibun.eq.1.and.isub.eq.neinbunch) then

        IF (IINSIDE.NE.0) THEN
          DPHSOUR(2,1)=BET1N*DT*FREQ(1)/HBAREV1
          DPHSOUR(2,2)=BET1N*DT*FREQ(NFREQ)/HBAREV1
        ENDIF

C- STORE NUMBER OF POINTS FOR INTEGRATION

        IF (IOBSV.EQ.ICBRILL) IZTOT(ISOUR)=IZAEHL

      endif !ibun.eq.1

      IF (IAMPLI.LT.0) THEN

        DXEXI=MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &    -MAX(SOURCEAO(1,1,ISOUR),XIANF)
        if (ampr2corr.eq.-9999.0d0) ampr2corr=dxexi
        CENXEXI=(MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &    +MAX(SOURCEAO(1,1,ISOUR),XIANF))/2.D0
c        GAMGAM0=(SOURCEG(1,1,ISOUR))**2
c        GAMGAM=((SOURCEG(1,1,ISOUR)+SOURCEG(2,2,ISOUR)))**2
        GAMGAM0=(SOURCEG(1,1,ISOUR)*(egamma/dmygamma))**2
        GAMGAM=(
     &    (SOURCEG(1,1,ISOUR)+SOURCEG(2,2,ISOUR))*(egamma/dmygamma)
     &    )**2

c        DTPHASE=(WTRA2IS(ISOUR)+(1.0D0/GAMGAM0)*DXEXI/2.D0)/CLIGHT1
c     &    *GAMGAM0/GAMGAM

        slopein=sqrt(vyin**2+vzin**2)/vxin
        slope=sqrt(vyelec**2+vzelec**2)/vxelec

        if (myinum.gt.nlpoi/dxexi) then
          WI=(WTRA2IS(ISOUR)
     &      -DXEXI/2.0D0*slopein**2) !wi is detour for on-axis particle
     &      *(dmygamma/egamma)**2
        else
          if (iwarnwi.eq.0) then
            write(lungfo,*)
            write(lungfo,*)'*** Warning in souintana_omp:'
            write(lungfo,*)'*** MYINUM is rather small with respect to NLPOI'
            write(lungfo,*)'*** Length of trajectories are now calculated by simple'
            write(lungfo,*)'*** integration with souintana_omp, which might be poor'
            write(lungfo,*)
            write(lungfo,*)
            write(6,*)'*** Warning in souintana_omp:'
            write(6,*)'*** MYINUM is rather small with respect to NLPOI'
            write(6,*)'*** Length of trajectories are now calculated by simple'
            write(6,*)'*** integration with souintana_omp, which might be poor'
            write(6,*)
            iwarnwi=1
          endif
          wi=(yp2zp2i/2.0d0
     &      -DXEXI/2.0D0*slopein**2) !wi is detour for on-axis particle
     &      *(dmygamma/egamma)**2
        endif

        xn1=cenxexi
        yn1=(xn1-cenxexi)*vyelec/vxelec
        zn1=(xn1-cenxexi)*vzelec/vxelec

        drn2=(
     &    (yn1+dxexi*vyelec/vxelec)**2+
     &    (zn1+dxexi*vzelec/vxelec)**2
     &    )/
     &    (2.0d0*(xobsv-xn1-dxexi))

        drn1=(
     &    yn1**2+
     &    zn1**2
     &    )/
     &    (2.0d0*(xobsv-xn1))

        DTPHASE=(
     &    WI+DXEXI*(slope**2/2.0d0+1.0d0/(2.0D0*GAMGAM0))
     &    +drn2-drn1)
     &    /CLIGHT1*GAMGAM0/GAMGAM

        AMPDT=AMPSHIFT(1)/CLIGHT1/2.0D0/GAMGAM0
        FREQR=2.D0*PI1/DTPHASE*HBAREV1
        POW=SPECPOW(jliob)

        if (jpin.ne.3.or.jpin.eq.3.and.ielec.eq.1) SPECPOW(jliob)=0.0D0

        DO I=1,-IAMPLI
          R02=(XOBSV-CENXEXI)**2+YOBSV**2+ZOBSV**2
          R2=(XOBSV-CENXEXI-DXEXI*(I-ABS(IAMPLI)/2+1))**2
     &      +YOBSV**2+ZOBSV**2
          SPECPOW(jliob)=SPECPOW(jliob)+POW*R02/R2
     &      *R2/(sqrt(R2)-ampr2corr/2.0d0)**2/nelec
        ENDDO

      ENDIF  !endif iampli.lt.0

      DO kfreq=1,NFREQ

        jliobfr=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(kfreq-1))
        kfrob=kfreq+NFREQ*(IOBSV-1)
        jobfr=IOBSV+NOBSV*(kfreq-1)

        OM=FREQ(kfreq)/HBAREV1

        IF (IAMPLI.LT.0) THEN

          AX0=affe(1,kfrob)
          AY0=affe(2,kfrob)
          AZ0=affe(3,kfrob)

          AX=AX0
          AY=AY0
          AZ=AZ0

          BX0=affe(4,kfrob)
          BY0=affe(5,kfrob)
          BZ0=affe(6,kfrob)

          BXc=BX0
          BYc=BY0
          BZc=BZ0

          affe(1:6,kfrob)=(0.0D0,0.0D0)

          R0=OBSV(1,NOBSV/2+1)-CENXEXI
          R02=R0*R0
          R00=R0

c          H2=(YOBSV-vyelec)**2+(ZOBSV-vzelec)**2
c          H2R2=H2/R02
c
c          DTPHASE=(WTRA2IS(ISOUR)+(H2R2+1.0D0/GAMGAM0)*DXEXI/2.D0)/CLIGHT1
c     &      *GAMGAM0/GAMGAM
c     &      +AMPDT

          xn1=cenxexi
          yn1=(xn1-cenxexi)*vyelec/vxelec
          zn1=(xn1-cenxexi)*vzelec/vxelec

          drn2=(
     &      (yn1+dxexi*vyelec/vxelec-yobsv)**2+
     &      (zn1+dxexi*vzelec/vxelec-zobsv)**2
     &      )/
     &      (2.0d0*(xobsv-xn1-dxexi))

          drn1=(
     &      (yn1-yobsv)**2+
     &      (zn1-zobsv)**2
     &      )/
     &      (2.0d0*(xobsv-xn1))

          DTPHASE=(
     &      WI+DXEXI*(slope**2/2.0d0+1.0d0/(2.0D0*GAMGAM0))
     &      +drn2-drn1)
     &      /CLIGHT1*GAMGAM0/GAMGAM
     &      +AMPDT

          PHI=2.D0*PI1*FREQ(kfreq)*ECHARGE1/HPLANCK1*DTPHASE

          DMODU=EXP(ZI*PHI)
          DMODU0=DMODU
          DDMODU=ZONE

          if (ibunch.ne.0) then

            zzi=zeleco
            yyi=yeleco
            zzpi=zpeleco
            yypi=ypeleco

c  transform to beginning of first section
            zz=tfmhtoti(1,1)*zzi+tfmhtoti(1,2)*zzpi
            zzp=tfmhtoti(2,1)*zzi+tfmhtoti(2,2)*zzpi
            yy=tfmvtoti(1,1)*yyi+tfmvtoti(1,2)*yypi
            yyp=tfmvtoti(2,1)*yyi+tfmvtoti(2,2)*yypi

          endif !ibunch

          DO I=1,-IAMPLI

            if (ibunch.ne.0) then

c  transform to center of section, no closed orbit!
                zelec=tfmhc(1,1)*zz+tfmhc(1,2)*zzp
                zpelec=tfmhc(2,1)*zz+tfmhc(2,2)*zzp
                yelec=tfmvc(1,1)*yy+tfmvc(1,2)*yyp
                ypelec=tfmvc(2,1)*yy+tfmvc(2,2)*yyp

c  transform to beginning of next section

              zzi=zz
              yyi=yy
              zzpi=zzp
              yypi=yyp

              zz=tfmh(1,1)*zzi+tfmh(1,2)*zzpi
              zzp=tfmh(2,1)*zzi+tfmh(2,2)*zzpi
              yy=tfmv(1,1)*yyi+tfmv(1,2)*yypi
              yyp=tfmv(2,1)*yyi+tfmv(2,2)*yypi

            endif !ibunch

            R0=OBSV(1,NOBSV/2+1)+DXEXI/2.D0*(-IAMPLI-2*(I-1)-1)-CENXEXI
            CORRR0=R00/R0
            !corrects for mistake of averaging over 1/r2, if e.g.
            !the repeated device is a long undulator
     &        *(R0/(R0-ampr2corr/2.0d0))**2

            R02=R0*R0
c            H2=(YOBSV)**2+(ZOBSV)**2
c            H2R2=H2/R02

c            DTPHASE=(WTRA2IS(ISOUR)+(H2R2+1.0D0/GAMGAM0)*DXEXI/2.D0)/CLIGHT1
c     &        *GAMGAM0/GAMGAM
c     &        +AMPDT

c            slope=sqrt(vyelec**2+vzelec**2)/vxelec

            slope=sqrt(ypelec**2+zpelec**2)*gamgam0/gamgam

            slope=sqrt(ypelec**2+zpelec**2)

            xn1=cenxexi-dxexi/2.d0*(-iampli-2*(i-1)-1)
     &          /(R0/(R0-ampr2corr/2.0d0))**2 !empirically
            yn1=(xn1-cenxexi)*ypelec
            zn1=(xn1-cenxexi)*zpelec

            yn1=yelec
            zn1=zelec

            drn2=(
     &        (yn1+dxexi*ypelec-yobsv)**2+
     &        (zn1+dxexi*zpelec-zobsv)**2
     &        )/
     &        (2.0d0*(xobsv-xn1-dxexi))

            drn1=(
     &        (yn1-yobsv)**2+
     &        (zn1-zobsv)**2
     &        )/
     &        (2.0d0*(xobsv-xn1))

            DTPHASE=(
     &        WI+DXEXI*(slope**2/2.0d0+1.0d0/(2.0D0*GAMGAM0))
     &        +drn2-drn1)
     &        /CLIGHT1*GAMGAM0/GAMGAM
     &        +AMPDT

            PHI=2.D0*PI1*FREQ(kfreq)*ECHARGE1/HPLANCK1*DTPHASE

            DMODU=EXP(ZI*PHI)
            DMODU0=DMODU
            DDMODU=ZONE

              affe(1,kfrob)=affe(1,kfrob)+AX
              affe(2,kfrob)=affe(2,kfrob)+AY
              affe(3,kfrob)=affe(3,kfrob)+AZ

              affe(4,kfrob)=affe(4,kfrob)+bXc
              affe(5,kfrob)=affe(5,kfrob)+bYc
              affe(6,kfrob)=affe(6,kfrob)+bZc

            IF (AMPRAN.NE.0.0D0) THEN
              PHI=2.D0*PI1*XRANA(I)/FREQR*FREQ(kfreq)
              DDMODU=EXP(ZI*PHI)
            ENDIF   !(AMPRAN.NE.0.0D0)

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

            BXc=BX0*CORRR0
            BYc=BY0*CORRR0
            BZc=BZ0*CORRR0

            BXc=BXc*DMODU
            BYc=BYc*DMODU
            BZc=BZc*DMODU

          ENDDO !IAMPLI

          zelec=zeleco
          yelec=yeleco

        ENDIF  !(IAMPLI.LT.0)

        if (jpin.eq.3) then

          FSPEC(1)=ISOUR
          FSPEC(2)=IOBSV
          FSPEC(3)=xobsv
          FSPEC(4)=yobsv
          FSPEC(5)=zobsv
          if (abs(fspec(4)).lt.1.0d-15) fspec(4)=0.0d0
          if (abs(fspec(5)).lt.1.0d-15) fspec(5)=0.0d0
          FSPEC(6)=FREQ(kfreq)
          FSPEC(7)=
     &      (
     &      DREAL(affe(1,kfrob))*DREAL(affe(1,kfrob))
     &      +DIMAG(affe(1,kfrob))*DIMAG(affe(1,kfrob))
     &      +DREAL(affe(2,kfrob))*DREAL(affe(2,kfrob))
     &      +DIMAG(affe(2,kfrob))*DIMAG(affe(2,kfrob))
     &      +DREAL(affe(3,kfrob))*DREAL(affe(3,kfrob))
     &      +DIMAG(affe(3,kfrob))*DIMAG(affe(3,kfrob))
     &      )*specnor*bunnor
          FSPEC(8)=1
          FSPEC(9)=1
          FSPEC(10)=kfreq
          FSPEC(11)=dreal(affe(1,kfrob))*sqrt(specnor*bunnor)
          FSPEC(12)=dimag(affe(1,kfrob))*sqrt(specnor*bunnor)
          FSPEC(13)=dreal(affe(2,kfrob))*sqrt(specnor*bunnor)
          FSPEC(14)=dimag(affe(2,kfrob))*sqrt(specnor*bunnor)
          FSPEC(15)=dreal(affe(3,kfrob))*sqrt(specnor*bunnor)
          FSPEC(16)=dimag(affe(3,kfrob))*sqrt(specnor*bunnor)
          FSPEC(17)=0.0d0
          FSPEC(18)=0.0d0
          FSPEC(19)=0.0d0
          FSPEC(20)=0.0d0

          FSPEC(21)=dreal(affe(4,kfrob))*sqrt(specnor*bunnor)
          FSPEC(22)=dimag(affe(4,kfrob))*sqrt(specnor*bunnor)
          FSPEC(23)=dreal(affe(5,kfrob))*sqrt(specnor*bunnor)
          FSPEC(24)=dimag(affe(5,kfrob))*sqrt(specnor*bunnor)
          FSPEC(25)=dreal(affe(6,kfrob))*sqrt(specnor*bunnor)
          FSPEC(26)=dimag(affe(6,kfrob))*sqrt(specnor*bunnor)
          FSPEC(27)=0.0d0
          FSPEC(28)=0.0d0
          FSPEC(29)=0.0d0
          FSPEC(30)=0.0d0

          if (ispecdip.le.0) then
            cenxexi=(min(sourceeo(1,1,isour),xiend)
     &        +max(sourceao(1,1,isour),xianf))/2.d0
          else
            cenxexi=x0dip(isour)
          endif

          dist0=pincen(1)-cenxexi
          dist02=dist0**2

          h2=(yobsv**2+zobsv**2)/dist02
          if (h2.lt.0.01) then
            ddist=dist0*(h2/2.0d0-h2**2/8.0d0)
          else
            ddist=dist0*(sqrt(1.0d0+h2)-1.0d0)
          endif

          dphase=ddist/freq(kfreq)*wtoe1*1.0d9*twopi1

          FSPEC(31)=dphase

          call hfm(nidspec,fspec)

          FSPEC(1)=xobsv
          FSPEC(2)=yobsv
          FSPEC(3)=zobsv
          if (abs(fspec(2)).lt.1.0d-15) fspec(2)=0.0d0
          if (abs(fspec(3)).lt.1.0d-15) fspec(3)=0.0d0
          FSPEC(4)=powpow*pownor
          FSPEC(5)=0.0d0
          FSPEC(6)=1.0d0
          FSPEC(7)=1.0d0
          FSPEC(8)=iobsv
          FSPEC(9)=ISOUR

          call hfm(nidpow,fspec)

          if (istokes.ne.0) then

            APOLH=
     &        affe(1,kfrob)*CONJG(VSTOKES(1,1))
     &        +affe(2,kfrob)*CONJG(VSTOKES(1,2))
     &        +affe(3,kfrob)*CONJG(VSTOKES(1,3))

            APOLR=
     &        affe(1,kfrob)*CONJG(VSTOKES(2,1))
     &        +affe(2,kfrob)*CONJG(VSTOKES(2,2))
     &        +affe(3,kfrob)*CONJG(VSTOKES(2,3))

            APOLL=
     &        affe(1,kfrob)*CONJG(VSTOKES(3,1))
     &        +affe(2,kfrob)*CONJG(VSTOKES(3,2))
     &        +affe(3,kfrob)*CONJG(VSTOKES(3,3))

            APOL45=
     &        affe(1,kfrob)*CONJG(VSTOKES(4,1))
     &        +affe(2,kfrob)*CONJG(VSTOKES(4,2))
     &        +affe(3,kfrob)*CONJG(VSTOKES(4,3))

            STOK1=
     &        APOLR*CONJG(APOLR)+
     &        APOLL*CONJG(APOLL)

            STOK2=-STOK1+
     &        2.0d0*APOLH*CONJG(APOLH)

            STOK3=
     &        2.0d0*APOL45*CONJG(APOL45)-
     &        STOK1

            STOK4=
     &        APOLR*CONJG(APOLR)-
     &        APOLL*CONJG(APOLL)

            FSPEC(1)=IOBSV
            FSPEC(2)=xobsv
            FSPEC(3)=yobsv
            FSPEC(4)=zobsv
            if (abs(fspec(3)).lt.1.0d-15) fspec(3)=0.0d0
            if (abs(fspec(4)).lt.1.0d-15) fspec(4)=0.0d0
            FSPEC(5)=FREQ(kfreq)
            FSPEC(6)=STOK1*specnor*bunnor
            FSPEC(7)=STOK2*specnor*bunnor
            FSPEC(8)=STOK3*specnor*bunnor
            FSPEC(9)=STOK4*specnor*bunnor
            FSPEC(10)=1.0d0
            FSPEC(11)=1.0d0
            FSPEC(12)=kfreq

            CALL hfm(NIDSTOK,FSPEC)

          ENDIF   !ISTOKES

        endif !ipin.eq.3

        if (ihbunch.ne.0.and.iobsv.eq.icbrill) then
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
            fillb(17)=xobsv
            fillb(18)=yobsv
            fillb(19)=zobsv
            if (abs(fillb(18)).lt.1.0d-15) fillb(18)=0.0d0
            if (abs(fillb(19)).lt.1.0d-15) fillb(19)=0.0d0
            fillb(20)=kfreq
            fillb(21)=freq(kfreq)
            speck=
     &        DREAL(
     &        affe(1,kfrob)*CONJG(affe(1,kfrob))
     &        +affe(2,kfrob)*CONJG(affe(2,kfrob))
     &        +affe(3,kfrob)*CONJG(affe(3,kfrob))
     &        )*specnor*bunnor
            fillb(22)=speck*nelec

            if (istokes.ne.0) then

              APOLH=
     &          affe(1,kfrob)*CONJG(VSTOKES(1,1))
     &          +affe(2,kfrob)*CONJG(VSTOKES(1,2))
     &          +affe(3,kfrob)*CONJG(VSTOKES(1,3))

              APOLR=
     &          affe(1,kfrob)*CONJG(VSTOKES(2,1))
     &          +affe(2,kfrob)*CONJG(VSTOKES(2,2))
     &          +affe(3,kfrob)*CONJG(VSTOKES(2,3))

              APOLL=
     &          affe(1,kfrob)*CONJG(VSTOKES(3,1))
     &          +affe(2,kfrob)*CONJG(VSTOKES(3,2))
     &          +affe(3,kfrob)*CONJG(VSTOKES(3,3))

              APOL45=
     &          affe(1,kfrob)*CONJG(VSTOKES(4,1))
     &          +affe(2,kfrob)*CONJG(VSTOKES(4,2))
     &          +affe(3,kfrob)*CONJG(VSTOKES(4,3))

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

          fillb(27)=powpow*pownor*nelec
          fillb(28)=isour
          fillb(29)=t2

          fillb(30)=dreal(affe(1,kfrob))
          fillb(31)=dimag(affe(1,kfrob))
          fillb(32)=dreal(affe(2,kfrob))
          fillb(33)=dimag(affe(2,kfrob))
          fillb(34)=dreal(affe(3,kfrob))
          fillb(35)=dimag(affe(3,kfrob))
          fillb(36)=dreal(affe(4,kfrob))
          fillb(37)=dimag(affe(4,kfrob))
          fillb(38)=dreal(affe(5,kfrob))
          fillb(39)=dimag(affe(5,kfrob))
          fillb(40)=dreal(affe(6,kfrob))
          fillb(41)=dimag(affe(6,kfrob))

          call hfm(nidbunch,fillb)
        endif

        if (ibunphase.eq.1) then

          are(1:6)=dreal(affe(1:6,kfrob))
          aim(1:6)=dimag(affe(1:6,kfrob))

          if (iobunch.eq.-9999.and.ielec.eq.1) then
            ampz(kfreq)=sqrt(are(3)**2+aim(3)**2)
            if (ampz(kfreq).gt.ampzmax(kfreq)) then
              ampzmax(kfreq)=ampz(kfreq)
              kobs(kfreq)=iobsv
            endif
          else if (iobsv.eq.iobunch) then
            ampz(kfreq)=sqrt(are(3)**2+aim(3)**2)
            if (ampz(kfreq).ne.0.0d0) then
              azcos(kfreq)=are(3)/ampz(kfreq)
              azsin(kfreq)=aim(3)/ampz(kfreq)
            else
              azcos(kfreq)=1.0d0
              azsin(kfreq)=0.0d0
            endif
          endif

          if (iobsv.eq.nobsv.and.kfreq.eq.nfreq) then

            do jfreq=1,nfreq

              jfrob=jfreq+nfreq*(kobs(jfreq)-1)
              are(1:6)=dreal(affe(1:6,jfrob))
              aim(1:6)=dimag(affe(1:6,jfrob))

              ampz(jfreq)=sqrt(are(3)**2+aim(3)**2)
              if (ampz(jfreq).ne.0.0d0) then
                azcos(jfreq)=are(3)/ampz(jfreq)
                azsin(jfreq)=aim(3)/ampz(jfreq)
              else
                azcos(jfreq)=1.0d0
                azsin(jfreq)=0.0d0
              endif

            enddo !jfreq

            do job=1,nobsv
              do jfreq=1,nfreq

                jfrob=jfreq+nfreq*(job-1)

                are(1:6)=dreal(affe(1:6,jfrob))
                aim(1:6)=dimag(affe(1:6,jfrob))

                affe(1:6,jfrob)=dcmplx(
     &            azcos(jfreq)*are+azsin(jfreq)*aim,
     &            -azsin(jfreq)*are+azcos(jfreq)*aim
     &            )

                AFREQ(1:6,jfrob)=AFREQ(1:6,jfrob)+affe(1:6,jfrob)

                affe(1:6,jfrob)=(0.0D0,0.0D0)

              enddo
            enddo

          endif !iobsv.eq.nobsv

        else !ibunphase

          if (iobsv.eq.nobsv.and.kfreq.eq.nfreq) then
            AFREQ=AFREQ+affe
            affe=(0.0D0,0.0D0)
          endif

        endif !ibunphase

      ENDDO !kfreq

      if (isub.eq.neinbunch.and.iobsv.eq.nobsv) then

        do job=1,nobsv
          do kfreq=1,nfreq

            kfrob=kfreq+nfreq*(job-1)
            jliobfr=isour+nsource*(job-1+nobsv*(kfreq-1))
            jobfr=job+nobsv*(kfreq-1)

            IF(SPECCUT.GT.0.0D0) THEN
              ECMAXS=ECMAX(ISOUR)
              IF(FREQ(kfreq).GT.SPECCUT*ecdipev1*DMYENERGY**2*ECMAXS) THEN
                AFREQ(1:6,kfrob)=(0.0D0,0.0D0)
              ENDIF
            ENDIF

            AFREQ(1:3,kfrob)=AFREQ(1:3,kfrob)*REFLEC(1:3)
            AFREQ(4:6,kfrob)=AFREQ(4:6,kfrob)*REFLEC(1:3)

            IF (IPOLA.EQ.0) THEN

              speck=
     &          DREAL(
     &          AFREQ(1,kfrob)*CONJG(AFREQ(1,kfrob))
     &          +AFREQ(2,kfrob)*CONJG(AFREQ(2,kfrob))
     &          +AFREQ(3,kfrob)*CONJG(AFREQ(3,kfrob))
     &          )*specnor*bunnor

              SPEC(jliobfr)=SPEC(jliobfr)+speck

              REAIMA(1:3,1,jobfr)=REAIMA(1:3,1,jobfr)+
     &          DREAL(AFREQ(1:3,kfrob))/sqnbunch
              REAIMA(1:3,2,jobfr)=REAIMA(1:3,2,jobfr)+
     &          DIMAG(AFREQ(1:3,kfrob))/sqnbunch

              REAIMA(6:8,1,jobfr)=REAIMA(6:8,1,jobfr)+
     &          DREAL(AFREQ(4:6,kfrob))/sqnbunch
              REAIMA(6:8,2,jobfr)=REAIMA(6:8,2,jobfr)+
     &          DIMAG(AFREQ(4:6,kfrob))/sqnbunch

            ELSE    !IPOLA

              APOL=
     &          AFREQ(1,kfrob)*CONJG(VPOLA(1))
     &          +AFREQ(2,kfrob)*CONJG(VPOLA(2))
     &          +AFREQ(3,kfrob)*CONJG(VPOLA(3))

              SPEC(jliobfr)=SPEC(jliobfr)+
     &          DREAL(APOL*CONJG(APOL))*specnor*bunnor

              REAIMA(1:3,1,jobfr)=REAIMA(1:3,1,jobfr)+
     &          DREAL(AFREQ(1:3,kfrob))/sqnbunch
              REAIMA(1:3,2,jobfr)=REAIMA(1:3,2,jobfr)+
     &          DIMAG(AFREQ(1:3,kfrob))/sqnbunch

              REAIMA(6:8,1,jobfr)=REAIMA(6:8,1,jobfr)+
     &          DREAL(AFREQ(4:6,kfrob))/sqnbunch
              REAIMA(6:8,2,jobfr)=REAIMA(6:8,2,jobfr)+
     &          DIMAG(AFREQ(4:6,kfrob))/sqnbunch

            ENDIF   !IPOLA

            IF (ISTOKES.NE.0) THEN

              APOLH=
     &          AFREQ(1,kfrob)*CONJG(VSTOKES(1,1))
     &          +AFREQ(2,kfrob)*CONJG(VSTOKES(1,2))
     &          +AFREQ(3,kfrob)*CONJG(VSTOKES(1,3))

              APOLR=
     &          AFREQ(1,kfrob)*CONJG(VSTOKES(2,1))
     &          +AFREQ(2,kfrob)*CONJG(VSTOKES(2,2))
     &          +AFREQ(3,kfrob)*CONJG(VSTOKES(2,3))

              APOLL=
     &          AFREQ(1,kfrob)*CONJG(VSTOKES(3,1))
     &          +AFREQ(2,kfrob)*CONJG(VSTOKES(3,2))
     &          +AFREQ(3,kfrob)*CONJG(VSTOKES(3,3))

              APOL45=
     &          AFREQ(1,kfrob)*CONJG(VSTOKES(4,1))
     &          +AFREQ(2,kfrob)*CONJG(VSTOKES(4,2))
     &          +AFREQ(3,kfrob)*CONJG(VSTOKES(4,3))

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

              STOKES(1,jobfr)=STOKES(1,jobfr)+
     &          STOK1*specnor*bunnor

              STOKES(2,jobfr)=STOKES(2,jobfr)+
     &          STOK2*specnor*bunnor

              STOKES(3,jobfr)=STOKES(3,jobfr)+
     &          STOK3*specnor*bunnor

              STOKES(4,jobfr)=STOKES(4,jobfr)+
     &          STOK4*specnor*bunnor

            ENDIF !ISTOKES

            AFREQ(1,kfrob)=(0.0d0,0.0d0)
            AFREQ(2,kfrob)=(0.0d0,0.0d0)
            AFREQ(3,kfrob)=(0.0d0,0.0d0)

          enddo !kfreq
        enddo !job

      endif !isub.eq.neinbunch

      if (ibun.eq.nbunch.and.isub.eq.neinbunch) then
        jliob=ISOUR+NSOURCE*(IOBSV-1)
        SPECPOW(jliob)=SPECPOW(jliob)*POWNOR
      endif

      IF (
     &    jpin.ne.0.and.jpin.ne.3.and.IOBSV.EQ.jobunch
     &    .or.
     &    (jpin.eq.3.or.jpin.eq.0).and.ielec.eq.1
     &    ) THEN


        WRITE(LUNGFO,*)
     &    '       phase advance per step at beginning and end of source for'
        if (jpin.ne.3) then
          WRITE(LUNGFO,*)
     &      '       lowest and highest photon energy at selected observation point:'
        else
          WRITE(LUNGFO,*)
     &      '       lowest and highest photon energy for first electron:'
        endif
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       beginning:',SNGL(DPHSOUR(1,1)),SNGL(DPHSOUR(1,2))
        WRITE(LUNGFO,*)'       end:      ',SNGL(DPHSOUR(2,1)),SNGL(DPHSOUR(2,2))
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       ROIs (boundary, precision, points):'
        WRITE(LUNGFO,*)

        DO IROI=1,NROIA-1
          WRITE(LUNGFO,*)
     &      IROI,SNGL(roi(1,IROI)),SNGL(roi(2,IROI)),IPOIROI(IROI+1)
        ENDDO
        WRITE(LUNGFO,*)
     &    NROI,SNGL(roi(1,NROIA))

      ENDIF !IOBSV

      IF (
     &    jpin.ne.0.and.jpin.ne.3.and.IOBSV.EQ.NOBSV
     &    .or.
     &    (jpin.eq.0.or.jpin.eq.3).and.ielec.eq.1) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'       SOURCE, TOTAL NUMBER OF STEPS:',ISOUR,IZAEHL
        WRITE(LUNGFO,*)'       (controlled by NLPOI and namelist $ROIN)'
        WRITE(LUNGFO,*)
      ENDIF

      isourold=isour

      RETURN
      END
