*CMZ :  4.01/03 01/06/2023  06.50.18  by  Michael Scheer
*CMZ :  4.00/16 09/09/2022  17.17.34  by  Michael Scheer
*CMZ :  4.00/15 26/03/2022  11.38.18  by  Michael Scheer
*CMZ :  4.00/13 06/12/2021  12.45.17  by  Michael Scheer
*CMZ :  4.00/07 28/04/2020  21.30.02  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  12.15.26  by  Michael Scheer
*CMZ :  3.06/00 26/02/2019  17.28.29  by  Michael Scheer
*CMZ :  3.05/03 17/05/2018  17.16.20  by  Michael Scheer
*CMZ :  3.05/02 15/05/2018  15.28.23  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE souintall_omp(ISOUR)
*KEEP,gplhint.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,amplif90u.
      include 'amplif90u.cmn'
*KEND.

      use bunchmod
      use ompmod
      use wbetaf90m
      use souintmod

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE

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
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,tralin.
      include 'tralin.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,strings.
      include 'strings.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEEP,souintanac.
      include 'souintanac.cmn'
*KEND.

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RESRAN
      double precision, dimension(:), allocatable :: apz,azs,azc
      DOUBLE PRECISION :: RMS,RESPOW,vx1,vy1,vz1,vn,schwingungen,
     &  yy,zz,yyp,zzp,
     &  tfmh(2,2),tfmv(2,2),tfmhi(2,2),tfmvi(2,2),
     &  tfmdeh(2,2),tfmdev(2,2),tfmdehi(2,2),tfmdevi(2,2),tfm1(2,2),
     &  w22(2,2),dum22(2,2),
     &  x1,y1,z1,x2,wx1,dxb
     &  ,alpha0(2),beta0(2),alphah,alphav,
     &  xub,yub,zub,ypub,zpub,gammaub,tphase,
     &  beta0h,s0h,beta0v,s0v,gammah,gammav,
     &  xco,yco,zco,ypco,zpco,egammaco,dpp=0.0d0,
     &  xfirst,yfirst,zfirst,ypfirst,zpfirst,egammafirst
     &  ,rm,rq,alpha,beta
     &  ,psi,spsi,cpsi

      real xran(5)

      INTEGER IOBSV,ISOUR,JX10,JDX10,INSIDE,nutracko
      INTEGER NTOTIN,NTOT2IN,IREP,KINSIDE
      INTEGER JEFOLD,JFOLD,ICAL,kfreq,ICOUNT,kobsv,
     &  iw2(2),ifail,i1,icen,ierr,i2,lzaehl

      INTEGER jpin
      common/souintc/jpin
c      common/kobsc/ampzmax,kobs

      DATA tfm1(1,1),tfm1(1,2),tfm1(2,1),tfm1(2,2)/1.0d0,0.0d0,0.0d0,1.0d0/
      DATA tfmhi(1,1),tfmhi(1,2),tfmhi(2,1),tfmhi(2,2)/1.0d0,0.0d0,0.0d0,1.0d0/
      DATA tfmvi(1,1),tfmvi(1,2),tfmvi(2,1),tfmvi(2,2)/1.0d0,0.0d0,0.0d0,1.0d0/
      DATA ICAL/0/

      save ical

      allocate(ampzmax(nfreq),kobs(nfreq))

      IF (ICAL.EQ.0) THEN

c        print*,"******************* OMP ****************"
c        print*," "
c        print*,"--- SOUINTALL_OMP CALLED ---"
c        print*," "

        if (bunchcharge.eq.0.0d0) then
          bunchcharge=neinbunch*echarge1
        endif

        if (ibunch.eq.0) then
          ibunphase=0
          iobunch=0
        else if (abs(ibunch).ne.1) then
                write(LUNGFO,*)
     &            '*** Error in SOUINTALL_OMP: Bad value of IBUNCH (not -1,0,1)!  ***'
                write(6,*)
     &            '*** Error in SOUINTALL_OMP: Bad value of IBUNCH (not -1,0,1)!  ***'
                stop '*** Program WAVE aborted ***'
        endif

        if (iwbunch.gt.0) then
          open(unit=22,file='wave_phasespace_bunch.dat',status='old',
     &      iostat=ierr)
          if (ierr.eq.0) close(22,status='delete')
          open(unit=22,file='wave_phasespace_bunch.dat',status='new')
        endif

        if (mpinr.eq.0) then
          if (ibunphase.eq.0) iobunch=icbrill
          jobunch=icbrill
          if (iobunch.eq.0) then
            iobunch=icbrill
          endif
        else
          if (ibunphase.eq.0) iobunch=1
          jobunch=1
          if (iobunch.eq.0) then
            iobunch=1
          endif
        endif

        if (ibunch.eq.0) then
          neinbunch=1
          nbunch=1
          if (iobunch.eq.-9999) then
            if (mpinr.ne.0) then
              jobunch=icbrill
            else
              jobunch=1
            endif
          else
            jobunch=iobunch
          endif
        else if (iobunch.ne.-9999) then
          jobunch=iobunch
        endif

        nelec=nbunch*neinbunch

        IF (IAMPSEED.NE.0) CALL RMARIN(IAMPSEED,NTOTIN,NTOT2IN)

        if (jobunch.gt.nobsv) then
          write(6,*)' '
          write(6,*)' *** WARNING IN SOUINTALL_OMP: IOBUNCH.GT.NOBSV'
          write(6,*)' SET ACCORDING TO SELECTED POINT'
          write(6,*)' '
          write(lungfo,*)' '
          write(lungfo,*)' *** WARNING IN SOUINTALL_OMP: IOBUNCH.GT.NOBSV'
          write(lungfo,*)' SET ACCORDING TO SELECTED POINT'
          write(lungfo,*)' '
        endif

        allocate(wbetasub(16,3))

        ICAL=1
      ENDIF !ICAL

      xub=sourceao(1,1,isour)
      yub=sourceao(2,1,isour)
      zub=sourceao(3,1,isour)
      ypub=sourceao(2,2,isour)/sourceao(1,2,isour)
      zpub=sourceao(3,2,isour)/sourceao(1,2,isour)
      gammaub=sourceg(1,1,isour)

      xco=xub
      yco=yub
      zco=zub
      ypco=ypub
      zpco=zpub
      egammaco=gammaub

      xfirst=xco
      yfirst=yco
      zfirst=zco
      ypfirst=ypco
      zpfirst=zpco
      egammafirst=egammaco

      xelec=xco
      yelec=yco
      zelec=zco
      ypelec=ypco
      zpelec=zpco
      egamma=egammaco

      if (ibunch.ne.0) then

        if (iemit.eq.0.and.ilintra.lt.0) then
          nutracko=nutrack
          allocate(wbeta(16,nco))
          allocate(wbetak(3,nco))
c calculate transfer matrix from entrance to center
          call wbetfn
          nutrack=nutracko
        endif

        if (ilintra.lt.0) then

          tfmh=tfm1
          tfmv=tfm1
          tfmdeh=tfm1
          tfmdev=tfm1

          x1=sourceao(1,1,isour)
          y1=sourceao(2,1,isour)
          z1=sourceao(3,1,isour)
          x2=sourceeo(1,1,isour)

          dxb=(wbeta(1,nco)-wbeta(1,1))/(nco-1)

          wx1=wbeta(1,1)

          i1=nint((x1-wx1)/dxb)
          i2=nint((x2-wx1)/dxb)

          if (i1.lt.1) then
            if (i1.lt.0) then
              write(6,*)
     &          '*** Error in SOUTINALL: No horizontal beta function for first point of source found ***'
              write(6,*)
     &          'Taking value from x = ',wbeta(1,1)
            endif
            i1=1
          endif

          if (i2.gt.nco) then
              write(6,*)
     &          '*** Error in SOUTINALL: No horizontal beta function for last point of source found ***'
              write(6,*)
     &          'Taking value from x = ',wbeta(1,nco)
            i2=nco
          endif

          if (xlintra.eq.-9999.0d9) then
            icen=nint(((x1+x2)/2.0d0-wx1)/dxb)
          else
            icen=nint(xlintra/dxb)+1
            if (icen.le.i1.or.icen.gt.i2) then
              write(6,*)
     &          '*** Error in SOUINTALL_OMP: Bad XLINTRA, please check input *** '
              write(6,*)
     &          '*** XLINTRA MUST BE INSIDE SOURCE! *** '
              stop '*** WAVE aborted'
            endif
          endif

          wbetasub(1:16,1)=wbeta(1:16,i1)
          wbetasub(1:16,2)=wbeta(1:16,icen)
          wbetasub(1:16,3)=wbeta(1:16,i2)

          alpha0(1)=-wbeta(3,i1)/2.d0
          alpha0(2)=-wbeta(5,i1)/2.d0
          beta0(1)=wbeta(2,i1)
          beta0(2)=wbeta(4,i1)

          alpha=-wbeta(3,icen)/2.d0
          beta=  wbeta(2,icen)
          psi=   wbeta(8,icen)

          cpsi=cos(psi)
          spsi=sin(psi)
          rq=sqrt(beta/beta0(1))
          rm=sqrt(beta*beta0(1))
          tfmh(1,1) = rq * (cpsi+alpha0(1)*spsi)
          tfmh(1,2) = rm * spsi
          tfmh(2,1)=
     &      ((alpha0(1)-alpha)*cpsi - (1.0d0+alpha0(1)*alpha)*spsi) / rm
          tfmh(2,2)=
     &      (cpsi-alpha*spsi) / rq

          alpha=-wbeta(5,icen)/2.d0
          beta=  wbeta(4,icen)
          psi=   wbeta(9,icen)

          cpsi=cos(psi)
          spsi=sin(psi)
          rq=sqrt(beta/beta0(2))
          rm=sqrt(beta*beta0(2))
          tfmv(1,1) = rq * (cpsi+alpha0(2)*spsi)
          tfmv(1,2) = rm * spsi
          tfmv(2,1)=
     &      ((alpha0(2)-alpha)*cpsi - (1.0d0+alpha0(2)*alpha)*spsi) / rm
          tfmv(2,2)=
     &      (cpsi-alpha*spsi) / rq

          tfmdeh=tfmh
          tfmdev=tfmv

        endif !(ilintra.eq.0) then

        if(ilintra.gt.0) then

          open(unit=99,file='wave_lintra.dat',status='old')
          call util_skip_comment(99)
          read(99,*)tfmh(1,1),tfmh(1,2)
          read(99,*)tfmh(2,1),tfmh(2,2)
          read(99,*)tfmv(1,1),tfmv(1,2)
          read(99,*)tfmv(2,1),tfmv(2,2)
          read(99,*)wbetasub(1:16,1)
          read(99,*)wbetasub(1:16,2)
          read(99,*)wbetasub(1:16,3)
          close(99)

          tfmdeh=tfmh
          tfmdev=tfmv

          write(lungfo,*)' '
          write(lungfo,*)'      SOUTINALL:'
          write(lungfo,*)' '
          write(lungfo,*)'      Horizontal lineare transfer-matrix:'
          write(lungfo,*)' '
          write(lungfo,*)'      ',tfmh(1,1),tfmh(1,2)
          write(lungfo,*)'      ',tfmh(2,1),tfmh(2,2)
          write(lungfo,*)' '
          write(lungfo,*)'      Vertical lineare transfer-matrix:'
          write(lungfo,*)' '
          write(lungfo,*)'      ',tfmv(1,1),tfmv(1,2)
          write(lungfo,*)'      ',tfmv(2,1),tfmv(2,2)
          write(lungfo,*)' '

        else if(ilintra.lt.0) then

          open(unit=99,file='wave_lintra.dat',status='unknown',
     &      recl=256)
          write(99,*)'* ',icode,code(1:len_trim(code))
          write(99,*)tfmh(1,1),tfmh(1,2)
          write(99,*)tfmh(2,1),tfmh(2,2)
          write(99,*)tfmv(1,1),tfmv(1,2)
          write(99,*)tfmv(2,1),tfmv(2,2)
          write(99,*)wbetasub(1:16,1)
          write(99,*)wbetasub(1:16,2)
          write(99,*)wbetasub(1:16,3)
          close(99)

        endif !(ilintra.eq.0) then

        if (ilintra.ne.0) then !(ilintra.ne.0) then

          w22=tfmh

          dum22(1,1)=1.0d0
          dum22(1,2)=0.0d0
          dum22(2,1)=0.0d0
          dum22(2,2)=1.0d0

          call deqinv(2,w22,2,iw2,ifail,2,dum22)

          if (ifail.ne.0) then
            write(6,*)'*** Error in SOUINTALL_OMP: Matrix invertation failed'
            write(6,*)'Please check horizontal beta functions.'
            write(lungfo,*)'*** Error in SOUINTALL_OMP: Matrix invertation failed'
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
            write(6,*)'*** Error in SOUINTALL_OMP: Matrix invertation failed'
            write(6,*)'Please check vertical beta functions.'
            write(lungfo,*)'*** Error in SOUINTALL_OMP: Matrix invertation failed'
            write(lungfo,*)'Please check vertical beta functions.'
            stop '*** Program WAVE aborted ***'
          endif

          tfmvi=w22

          tfmdehi=tfmhi
          tfmdevi=tfmvi

        endif !(ilintra.eq.0) then

      endif !ibunch.ne.0

      if (mpinr.ne.0) then
        allocate(phaserphi(nobsvrphi));
        allocate(expom1rphi(nobsvrphi));
        allocate(afferphi(6,nobsvrphi*nfreq))
        allocate(unphrphi(6,nobsvrphi*nfreq))
        phaserphi=0.0d0
        expom1rphi=(0.0d0,0.0d0)
        afferphi=(0.0d0,0.0d0)
        unphrphi=(0.0d0,0.0d0)
      else
        allocate(affe(6,nobsv*nfreq))
        affe=(0.0d0,0.0d0)
      endif

      allocate(ampz(nfreq))
      ampz=0.0d0
      allocate(azcos(nfreq))
      allocate(azsin(nfreq))
      allocate(phrnrn(nfreq))
      allocate(phexp(nfreq))
      azcos=1.0d0
      azsin=0.0d0
      phrnrn=0.0d0
      phexp=0.0d0

      JPIN=IPIN
      if (ipin.eq.3) ipin=0
      JFOLD=IFOLD
      JEFOLD=IEFOLD

      IF (IAMPLI.LT.0) THEN
        IF (ISOUR.EQ.1) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      repetition of amplitude activated:'
          WRITE(LUNGFO,*)
     &      '      A -> A * (1 + exp(i*phi) + exp(i*2*phi)...'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      IAMPLI:',IAMPLI
          WRITE(LUNGFO,*)
        ENDIF !ISOUR
        IF (AMPRAN.NE.0.D0) THEN
          IF (ISOUR.EQ.1) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'      phase errors for repetition activated:'
            WRITE(LUNGFO,*)
     &        '      A -> A * (1 + exp(i*phi*xran1) + exp(i*2*phi*xran2)...'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'      IAMPSEED: ',IAMPSEED
            WRITE(LUNGFO,*)
          ENDIF !ISOUR
c          ALLOCATE(XRANA(-IAMPLI))
          CALL RNORML(XRANA,-IAMPLI)
          RMS=0.D0
          DO IREP=1,-IAMPLI
            XRANA(IREP)=AMPRAN*XRANA(IREP)
            RMS=RMS+XRANA(IREP)**2
          ENDDO
          RMS=SQRT(RMS/(-IAMPLI))
          WRITE(LUNGFO,*)
     &      '      rel. rms phase AMPRAN (input): ',SNGL(AMPRAN)
          WRITE(LUNGFO,*)
     &      '      rel. rms phase error 1. source (from generated errors): '
     &      ,SNGL(RMS)
        ENDIF   !(AMPRAN.NE.0.D0)

        IF (ABS(IMAMPLI).EQ.3) THEN
C--- WRITE FILES FOR PROGRAM PHASE OF JOHANNES BAHRDT
          CALL PHASE_BAHRDT
        ENDIF

      ENDIF   !IAMPLI

      JDX10=NOBSV/10
      JX10=JDX10

      IF (JDX10.LT.1) JDX10=1

      IF (ISOUR.EQ.1.AND.NOBSV.GT.1) THEN
        WRITE(6,*)' '
        WRITE(6,*)
     &    ' counting from 1 to 10 for first source to show progress:'
        WRITE(6,*)' '
        WRITE(6,*)' '
      ENDIF

      IF (ISPECMODE.EQ.1) CALL TRASOU(ISOUR)

      IF (JPIN.NE.2) THEN

C CALCULATE DISTRIBUTION IN PINHOLE {

        IX10=1
        KINSIDE=0

        ielec=0
        phrnrn=0.0
        phexp=cdexp(dcmplx(0.0d0,dble(phrnrn*twopi1)))

        egamma=sourceg(1,1,isour)

        do ibun=1,nbunch

          do isub=1,neinbunch

            ielec=ielec+1
            ampz=0.0d0

            if (ibunch.ne.0) then

c5.2.2013 bug?              if (iubunch.eq.0.or.iubunch.eq.1.or.ibunch.eq.1) then
              if (iubunch.eq.0.or.iubunch.eq.1) then

                if (bunchlen.eq.0.0d0) then
                  CALL GRNDMm(phrnrn,1)
                  phrnrn=phrnrn(1)
c                  CALL GRNDMm(phrnrn,nfreq)
c                  phexp=cdexp(dcmplx(0.0d0,dble(phrnrn*twopi1)))
                  phexp=cdexp(dcmplx(0.0d0,dble(phrnrn(1)*twopi1)))
                else if (bunchlen.gt.0.0d0)  then !(bunchlen.eq.0.0d0) then
c                  CALL RNORML(phrnrn,nfreq)
                  CALL RNORML(phrnrn,1)
                  phrnrn=phrnrn(1)
                  do kfreq=1,nfreq
                    schwingungen = bunchlen*1.0d9 /
c     &                (WTOE1/freq(kfreq)) * phrnrn(kfreq)
     &                (WTOE1/freq(kfreq)) * phrnrn(1)
                    schwingungen=mod(schwingungen,1.0d0)
                    phexp(kfreq)=cdexp(dcmplx(0.0d0,schwingungen*twopi1))
                  enddo
                else ! bunchlen
                  phrnrn=0.0d0
                  phexp=(1.0d0,0.0d0)
                endif !(bunchlen.eq.0.0d0) then

                bunchx=phrnrn(1)*bunchlen !long. position in bunch

              else if (iubunch.eq.2.or.iubunch.eq.3.or.iubunch.eq.4) then

c here for all frequencies have the same phase
                if (nsource.gt.1) then
                  write(6,*)'*** WARNING IN SOUINTALL_OMP: USE OF IBUNCH IS RATHER'
                  write(6,*)'*** TRICKY SINCE MULTIPLE SOURCES ARE NOT YET FULLY'
                  write(6,*)'*** IMPLEMENTED! ELECTRONS ON FILE wave_phasespace.dat'
                  write(6,*)'*** ARE READ FOR EACH SOURCE AND ALL VARIABLES MUST'
                  write(6,*)'*** MUST MATCH THE SOURCE, GOOD LUCK...'
                  write(6,*)'*** TYPICAL PROBLEMS: NUMBER OF STEPS IN SOURCES SEEMS TO BE EXCEEDED'
                  write(6,*)'*** OR END OF FILE FOR wave_phasespace.dat IS REACHED AND SO ON!'
                endif
                call bunch(tphase)
                bunchx=tphase*clight1 !long. position in bunch

                do kfreq=1,nfreq
                  phexp(kfreq)=cdexp(dcmplx(0.0d0,tphase*freq(kfreq)/hbarev1))
                enddo !kfreq

              else if (iubunch.eq.-2) then
                write(6,*)
     &            '*** Error in SOUINTALL_OMP: IUBUNCH=-2 is not supported anymore ***'
                write(6,*)'*** Program WAVE aborted ***'
                write(lungfo,*)
     &            '*** Error in SOUINTALL_OMP: IUBUNCH=-2 is not supported anymore ***'
                write(lungfo,*)'*** Program WAVE aborted ***'
                stop

              else if (iubunch.eq.-1) then

                call ubunch(xub,yub,zub,ypub,zpub,gammaub,tphase)
                bunchx=tphase*clight1 !long. position in bunch
                egamma=gammaub
                xelec=xub
                yelec=yub
                zelec=zub
                ypelec=ypub
                zpelec=zpub

                do kfreq=1,nfreq
                  schwingungen = bunchlen*1.0d9*clight1 /
     &              (freq(kfreq)/hbarev1) * tphase
                  schwingungen=mod(schwingungen,1.0d0)
                  phexp(kfreq)=cdexp(dcmplx(0.0d0,schwingungen*twopi1))
                enddo !kfreq

              else !if (iubunch.eq.0) then

                write(LUNGFO,*)
     &            '*** Error in SOUINTALL_OMP: Bad value of IUBUNCH (not -1,0,1,2,3,4)!  ***'
                write(6,*)
     &            '*** Error in SOUINTALL_OMP: Bad value of IUBUNCH (not -1,0,1,2,3,4)!  ***'
                stop '*** Program WAVE aborted ***'

              endif !iubunch

              if (iubunch.ne.3.and.iubunch.ne.-1) then

                xelec=SOURCEAO(1,1,ISOUR)
                yelec=SOURCEAO(2,1,ISOUR)
                zelec=SOURCEAO(3,1,ISOUR)

                vx1=SOURCEAO(1,2,ISOUR)
                vy1=SOURCEAO(2,2,ISOUR)
                vz1=SOURCEAO(3,2,ISOUR)

                ypelec=vy1/VX1
                zpelec=vz1/VX1
                egamma=sourceg(1,1,isour)

                xco=xelec
                yco=yelec
                zco=zelec
                ypco=ypelec
                zpco=zpelec
                egammaco=egamma

              endif !(iubunch.ne.3.and.iubunch.ne.1) then

              if (iubunch.ne.3.and.iubunch.ne.4.and.iubunch.ne.-1) then

                CALL RNORML(XRAN,5)

                if (espread.gt.0.0d0.and.iefold.eq.0) then
                  dpp=espread*xran(1)
                  egamma=egamma*(1.0d0+dpp)
                endif

                if (iubunch.ne.1) then

                  if (ifold.eq.0) then
                    if (bsigy(isour).gt.0.0d0) then
                      yy=bsigy(isour)*xran(2)
                    else
                      yy=0.0d0
                    endif
                    if (bsigyp(isour).gt.0.0d0) then
                      yyp=bsigyp(isour)*xran(3)
                    else
                      yyp=0.0d0
                    endif
                    if (bsigz(isour).gt.0.0d0) then
                      zz=bsigz(isour)*xran(4)
                    else
                      zz=0.0d0
                    endif
                    if (bsigzp(isour).gt.0.0d0) then
                      zzp=bsigzp(isour)*xran(5)
                    else
                      zzp=0.0d0
                    endif

                    zelec= tfmhi(1,1)*zz+tfmhi(1,2)*zzp
                    zpelec=tfmhi(2,1)*zz+tfmhi(2,2)*zzp
                    yelec= tfmvi(1,1)*yy+tfmvi(1,2)*yyp
                    ypelec=tfmvi(2,1)*yy+tfmvi(2,2)*yyp

                  endif !ifold

c simple treatment of closed orbit, assume small angles
                  zelec=zz+zco
                  zpelec=zzp+zpco
                  yelec=yy+yco
                  ypelec=yyp+ypco

                else !iubunch.ne.1

c phasespace ellipse: gammah*z**2+2*alphah*z*zp+betah*zp**2=eps0h

                  CALL RNORML(XRAN,5)

                  if (espread.gt.0.0d0.and.iefold.eq.0) then
                    egamma=egamma*(1.0d0+espread*xran(1))
                  endif

                  if (ifold.eq.0) then

                    ! assume beta(s)=beta0(s)+s**2/beta(0) and alpha0=-s/beta(0)
                    ! and a drift transfer-matrix ((1,s),(1,0))

                    alphah=-betaph/2.0d0
                    gammah=(1.0d0+alphah**2)/betah
                    beta0h=1.0d0/gammah
                    s0h=alphah/gammah

                    zz=sqrt(eps0h*beta0h)*xran(2)
                    zzp=-sqrt(eps0h/beta0h)*xran(3)
                    zelec=zz-s0h*zzp !inverse transformation
                    zpelec=zzp

                    alphav=-betapv/2.0d0
                    gammav=(1.0d0+alphav**2)/betav
                    beta0v=1.0d0/gammav
                    s0v=alphav/gammah

                    yy=sqrt(eps0v*beta0v)*xran(4)
                    yyp=-sqrt(eps0v/beta0v)*xran(5)
                    yelec=yy-s0v*yyp
                    ypelec=yyp

c simple treatment of closed orbit, assume small angles

                    xelec=sourceao(1,1,isour)
                    yelec=yelec+sourceao(2,1,isour)
                    zelec=zelec+sourceao(3,1,isour)
                    ypelec=ypelec+sourceao(2,2,isour)/sourceao(1,2,isour)
                    zpelec=zpelec+sourceao(3,2,isour)/sourceao(1,2,isour)

                  else !if (ifold.eq.0) then
                    write(6,*)
     &                '*** Error in SOUINTALL_OMP: IFOLD.ne.0.and.IUBUNCH.eq.1 are incompatible'
                    write(6,*)'*** Program WAVE aborted ***'
                    write(lungfo,*)
     &                '*** Error in SOUINTALL_OMP: IFOLD.ne.0.and.IUBUNCH.eq.1 are incompatible'
                    write(lungfo,*)'*** Program WAVE aborted ***'
                    stop
                  endif !(ifold.eq.0) then

                endif !(iubunch.ne.1) then

              endif !iubunch.ne.3.and.iubunch.ne.4.and.iubunch.ne.-1

              if (ielec.eq.1.and.ibunch.eq.-1) then
                egamma=egammafirst
                xelec=xfirst
                yelec=yfirst
                zelec=zfirst
                ypelec=ypfirst
                zpelec=zpfirst
              endif

              if (iwbunch.gt.0) then
                write(22,'(7e20.10)')egamma,bunchx,
     &            xelec,yelec,zelec,ypelec,zpelec
              endif

            endif !ibunch

            zelec=zelec+dpp*disp0
            zpelec=zpelec+dpp*ddisp0

            vn=clight1*dsqrt((1.0d0-1.0d0/egamma)*(1.0d0+1.0d0/egamma))

            vxelec=vn/sqrt(1.0d0+ypelec**2+zpelec**2)
            vyelec=vxelec*ypelec
            vzelec=vxelec*zpelec

            kobsv=nobsv

            jobunch=icbrill
            if (iobunch.ne.-9999) then
              jobunch=iobunch
            endif

            if (nobsv.eq.1.and.ielec.eq.1) then
              write(6,*)' '
              write(6,*)
     &          ' counting from 1 to 10 for first source to show progress:'
              write(6,*)' '
              write(6,*)' '
            endif

            if (ielec.eq.1) then
              iizaehl=0 !total number of steps in souintana_omp
            endif

            call tracks_omp(isour)

            if (mpinr.eq.0) then

              inside=1
              if (jpin.eq.3) inside=-3
              call souintana_ini_omp(1,jobunch,inside)
              call souintana_omp(isour,jobunch,inside)
              inside=1
              if (ielec.eq.1) then
c                NZAEHL10=KZAEHL*nelec*nobsv/10
                NZAEHL10=nelec*nobsv/10
                MZAEHL=NZAEHL10
                IX10=1
                lzaehl=0
              endif

              if (jpin.gt.0) call souintana_omp(isour,1,inside)

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(inside)
!$OMP& SHARED(lzaehl,mthreads,isour,ielec,jpin,jobunch,kobsv,ampzmax,kobs,mpinr,mzaehl,ix10,nzaehl10,iizaehl,kzaehl)

!$OMP DO
              do iobsv=2,kobsv-1

                lzaehl=lzaehl+1

                if (iobsv.eq.jobunch) cycle

                call souintana_omp(isour,iobsv,inside)

                if (inside.eq.0) then
                  kinside=1
                  write(lungfo,*)
                  write(lungfo,*)
     &              '      found observation point outside radiation cone;'
                  write(lungfo,*)
     &              '      number of source, number and (X,Y,Z) of observation point:'
                  write(lungfo,*)
     &              '      ',isour,iobsv,sngl(obsv(1,iobsv)),
     &              sngl(obsv(2,iobsv)),sngl(obsv(3,iobsv))
                  write(lungfo,*)
                endif   !inside


                IF (ISOUR.eq.1.and.lZAEHL.GE.MZAEHL.and.ix10.lt.10) THEN
                  CALL date_and_time(dtday,dttime,dtzone,idatetime)
                  WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
                  IX10=IX10+1
                  MZAEHL=MZAEHL+NZAEHL10
                ENDIF

              enddo !nobsv

!$OMP END DO
!$OMP END PARALLEL

              if (kobsv.gt.1) then

                iobsv=kobsv

                call souintana_omp(isour,iobsv,inside)

                if (inside.eq.0) then

                  kinside=1

                  write(lungfo,*)
                  write(lungfo,*)
     &              '      found observation point outside radiation cone;'
                  write(lungfo,*)
     &              '      number of source, number and (X,Y,Z) of observation point:'
                  write(lungfo,*)
     &              '      ',isour,iobsv,sngl(obsv(1,iobsv)),
     &              sngl(obsv(2,iobsv)),sngl(obsv(3,iobsv))
                  write(lungfo,*)

                endif   !inside


              endif !(kobsv.gt.1) then

            else

              call souintrphi_omp(isour,inside)

            endif !(mpinr.eq.0) then

          enddo !isub=1,neinbunch
        enddo !ibun=1,nbunch

C CALCULATE DISTRIBUTION IN PINHOLE }

      ELSE !JPIN.NE.2

        if (mpinr.gt.0)
     &    stop '*** Error in SOUINTALL_OMP: MPINR NOT ALLOWED FOR IPIN.EQ.2'

      ENDIF !JPIN.NE.2

      IF (ISPECMODE.EQ.1.AND.ISOUR.EQ.NSOURCE) THEN
        DEALLOCATE(DWT)
        DEALLOCATE(DWX)
        DEALLOCATE(DWX2P)
        DEALLOCATE(DWB)
        DEALLOCATE(DWB2P)
        DEALLOCATE(DWY)
        DEALLOCATE(DWY2P)
        DEALLOCATE(DWZ)
        DEALLOCATE(DWZ2P)
      ENDIF

      IF (AMPRAN.NE.0.D0.AND.IAMPLI.LT.0) THEN
        DEALLOCATE(XRANA)
      ENDIF   !(AMPRAN.NE.0.D0)

      IF (ISOUR.EQ.NSOURCE.AND.KINSIDE.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN SOUINTALL_OMP:'
        WRITE(6,*)'there are observation points outside radiation cone'
        WRITE(6,*)'maybe WGWINFC, WBL0CUT... or collimators not suitable'
        WRITE(6,*)'see output file WAVE.OUT for details'
        WRITE(6,*)
      ENDIF !KINSIDE

      IPIN=JPIN
      IFOLD=JFOLD
      IEFOLD=JEFOLD

      deallocate(phrnrn)
      deallocate(phexp)
      if (mpinr.ne.0) then
        deallocate(phaserphi);
        deallocate(expom1rphi);
        deallocate(afferphi)
        deallocate(unphrphi)
      else
        deallocate(affe)
      endif
      deallocate(ampz)
      deallocate(azcos)
      deallocate(azsin)

      if (ibunch.ne.0.and.isour.eq.nsource.and.ilintra.lt.0) then
        deallocate(wbeta)
        deallocate(wbetak)
      endif !ibunch.ne.0

      deallocate(ampzmax,kobs)

      RETURN
      END
