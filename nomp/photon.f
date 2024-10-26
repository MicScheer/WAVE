*CMZ :  4.00/15 27/04/2022  08.13.49  by  Michael Scheer
*CMZ :  4.00/13 29/11/2021  14.07.42  by  Michael Scheer
*CMZ :  4.00/07 27/04/2020  16.44.15  by  Michael Scheer
*CMZ :  4.00/04 17/05/2019  14.17.20  by  Michael Scheer
*CMZ :  3.06/00 19/02/2019  11.42.27  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.13.08  by  Michael Scheer
*CMZ :  3.05/04 28/06/2018  09.36.14  by  Michael Scheer
*CMZ :  3.03/04 27/10/2017  18.09.30  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.13.38  by  Michael Scheer
*CMZ :  3.02/00 15/10/2014  09.27.04  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  16.46.42  by  Michael Scheer
*CMZ :  2.70/07 14/01/2013  13.39.56  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.47.54  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  12.38.24  by  Michael Scheer
*CMZ :  2.68/02 08/06/2012  09.54.11  by  Michael Scheer
*CMZ :  2.67/06 24/05/2012  14.15.39  by  Michael Scheer
*CMZ :  2.67/05 16/05/2012  14.19.25  by  Michael Scheer
*CMZ :  2.67/04 15/05/2012  11.27.15  by  Michael Scheer
*CMZ :  2.67/03 09/05/2012  16.24.45  by  Michael Scheer
*-- Author :    Michael Scheer   09/05/2012
      subroutine photon(x,y,z,veln,gamma,bx,by,bz,dgamma,dtim,mode)
*KEEP,gplhint.
*KEND.

      use wave_g1

      implicit none

      integer mode,ical,i,ncy,ieof,iwarn

      double precision x,y,z,veln(3),bmag(3),bx,by,bz,ebeam,elmom,gamma,
     &  bparn,bper(3),bpern,epho,eec,bpervn(3),
     &  dgamma,b2per,
     &  pdum,dtim,ec,photons,de,deecg1,eecg1,g1,yrnint10,sigv,dum

      real rnrn(2),hrndm1m,xran(1),rr(2)
      real*8 fill(100)
      double precision, dimension (:), allocatable ::
     &  xrn,yrn,yrnint,coef,work1,work2,work3,work4,
     &  coefpsi,eecpsilog,sigpsi,cylog

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,photon.
      include 'photon.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      data ical/0/,iwarn/0/,ncy/0/

      save

      sigv=0.0d0
      nbing1=max(nbing1,10)

      if (mode.eq.-1) then

        if (ical.eq.0) then

          if (iphmode.lt.0) then
            ncy=11 !27.4.2020
!            ncy=0
!            open(unit=99,file='wave_cy.dat',status='old')
!1           call util_skip_comment_end(99,ieof)
!            if (ieof.ne.0) goto 9
!            read(99,*,end=9) dum
!            ncy=ncy+1
!            goto 1
!9           rewind(99)
            allocate(coefpsi(max(nbing1,ncy)))
            allocate(eecpsilog(max(nbing1,ncy)))
            allocate(cylog(max(nbing1,ncy)))
            allocate(sigpsi(max(nbing1,ncy)))
            do i=1,ncy
!              call util_skip_comment_end(99,ieof)
!              read(99,*)
              eecpsilog(i)=eecpsilog1(i)
              sigpsi(i)=sigpsi1(i)
              cylog(i)=cylog1(i)
            enddo
            eecpsilog(1:ncy)=log(eecpsilog(1:ncy))
            cylog(1:ncy)=log(cylog(1:ncy)/1000.0d0) !mrad -> rad
!            close(99)
          endif !iphmode.lt.0

          qfrms=0.0d0
          qfmean=0.0d0
          nqfphotons=0

          if (ihisini_c.eq.0.and.ihphotons.ne.0) then
            ihisini_c=-2
            call hisini
          endif !ihisini_c

          write(lungfo,*)' '
          write(lungfo,*)'     Subroutine photon (quantum fluctuations):'
          write(lungfo,*)' '

          if (iphmode.eq.1) then

            write(lungfo,*)
     &        '     *** Error in photon:'
            write(lungfo,*)
     &        '     *** IPHMODE=1 requires HBOOK package of CERN, which is'
            write(lungfo,*)
     &        '     *** not available in this version of WAVE'
            write(lungfo,*)
     &        '     --- Program WAVE aborted ---'
            write(6,*)
     &        '     *** Error in photon:'
            write(6,*)
     &        '     *** IPHMODE=1 requires HBOOK package of CERN, which is'
            write(6,*)
     &        '     *** not available in this version of WAVE'
            write(6,*)
     &        '     --- Program WAVE aborted ---'

            stop


          else !iphmode

            write(lungfo,*)'      number of sampling bins (NBING1)',nbing1
            write(lungfo,*)
     &        '      max. photon energy normalized to characteristic energy (EECMAXG1):'
            write(lungfo,*)
     &        '      ',eecmaxg1
            write(lungfo,*)' '

            if (iphmode.lt.0) then
              write(lungfo,*)'     transversal momentum is taken into account'
              write(lungfo,*)' '
            endif

            if (nbing1.le.2) then
              write(6,*)'*** Error in photon: NBING1 must be greater then 1'
              write(lungfo,*)'*** Error in photon: NBING1 must be greater then 1'
              stop
            endif

            allocate(xrn(max(nbing1,ncy)))
            allocate(yrn(max(nbing1,ncy)))
            allocate(yrnint(max(nbing1,ncy)))
            allocate(coef(max(nbing1,ncy)))
            allocate(work1(max(nbing1,ncy)))
            allocate(work2(max(nbing1,ncy)))
            allocate(work3(max(nbing1,ncy)))
            allocate(work4(max(nbing1,ncy)))

            deecg1=eecmaxg1/(nbing1-1)

            eecg1=0.0d0
            do i=1,10
              eecg1=eecg1+deecg1/10.0d0
              call util_g1_static(eecg1,g1)
              xrn(i)=eecg1
              yrn(i)=g1/eecg1
            enddo

            yrnint(1)=0.0d0
            do i=2,10
              yrnint(i)=yrnint(i-1)+(yrn(i)+yrn(i-1))/2.0d0*(xrn(i)-xrn(i-1))
            enddo
            yrnint10=yrnint(10)

            eecg1=0.0d0
            do i=10,nbing1
              eecg1=eecg1+deecg1
              call util_g1_static(eecg1,g1)
              xrn(i)=eecg1
              yrn(i)=g1/eecg1
            enddo

            call util_spline_running_integral(
     &        xrn(10:nbing1),yrn(10:nbing1),nbing1-10+1,yrnint(10:nbing1),
     &        coef,work1,work2,work3,work4)

            yrnint(10)=yrnint10
            yrnint(11:nbing1)=yrnint(11:nbing1)+yrnint10
            yrnint=yrnint/yrnint(nbing1)

            do i=2,nbing1
              if (yrnint(i).le.yrnint(i-1)) then
                stop '*** Error in photon: Bad integration of G1 ***'
              endif
            enddo

            call util_spline_coef(
     &        yrnint,xrn,nbing1,0.0d0,0.0d0,
     &        coef,work1,work2,work3,work4)

            if (iphmode.lt.0) then
              call util_spline_coef(
     &          eecpsilog,cylog,ncy,0.0d0,0.0d0,
     &          coefpsi,work1,work2,work3,work4)
            endif

          endif !iphmode

          !dgamma=-pdum*gamma**2*b**2*dt
          pdum=cgam1/2.0d0/pi1*clight1*(clight1/1.0d9)**2*emassg1

          ical=1

        endif !ical

        bmag(1)=bx
        bmag(2)=by
        bmag(3)=bz

        elmom=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
        ebeam=emassg1*gamma !GeV

        bparn=(bmag(1)*veln(1)+bmag(2)*veln(2)+bmag(3)*veln(3))
        bper=bmag-bparn*veln
        b2per=bper(1)**2+bper(2)**2+bper(3)**2
        bpern=sqrt(b2per)
        bpervn=bper/bpern

        ec=ecdipkev1*bpern*ebeam**2*1.0d-6 !GeV

        !dgamma = pdum * gamma**2 * b2per * dtim
        !dN = 15*sqrt(3)/8 * dE/Ec = 3.2476 * de/ec

        de=pdum*b2per*gamma*ebeam*dtim !GeV

        if (ec.ne.0.0d0) then
          photons=3.2476d0*de/ec !number of photons
        else
          photons=0.0d0
        endif

        call util_random(2,rnrn)  !S. 39

        if (photons.ge.1.0d0.and.iwarn.eq.0) then
          write(6,*)'*** Warning in PHOTON: Step size to large, ***'
          write(6,*)'*** i.e. probabilty to generate photon is greater than one! ***'
          write(6,*)'*** Check NLPOI ***'
          write(lungfo,*)'*** Warning in PHOTON: Step size to large, ***'
          write(lungfo,*)'*** i.e. probabilty to generate photon is greater than one! ***'
          write(lungfo,*)'*** Check NLPOI ***'
          iwarn=1
        endif

        if(rnrn(1).le.photons) then

          if (iphmode.eq.1) then
          else !iphmode
            call util_spline_inter(yrnint,xrn,coef,nbing1,
     &        dble(rnrn(2)),eec,-1)
            if (eec.lt.0.0d0) then
              print*,'*** Warning in PHOTON: Negative photon energy occured ***'
              print*,'rnrn:',rnrn
              print*,'setting Epho/Ec = 1.e-6'
              eec=1.0d-6
            endif
            if (iphmode.lt.0) then
              call util_spline_inter(eecpsilog,cylog,coefpsi,ncy,log(eec),
     &          sigv,-1)
              sigv=0.408*exp(sigv)/ebeam !rad
              call util_random_gauss(1,xran,rr)
              sigv=xran(1)*sigv
            else
              sigv=0.0d0
            endif !iphmode
          endif !iphmode

          epho=eec*ec
          dgamma=-epho/ebeam*gamma
          dpphoton=-(veln+sigv*bpervn)*epho
c          dpphoton(2)=dpphoton(2)*1.0d10

          nqfphotons=nqfphotons+1
          qfmean=qfmean+epho
          qfrms=qfrms+epho**2

          if (ihphotons.ne.0) then
            fill(1)=nutrack
            fill(2)=nustep
            fill(3)=x
            fill(4)=y
            fill(5)=z
            fill(6)=elmom*veln(1)*1.0d9
            fill(7)=elmom*veln(2)*1.0d9
            fill(8)=elmom*veln(3)*1.0d9
            fill(9)=ebeam*1.0d9
            fill(10)=-dpphoton(1)*1.0d9
            fill(11)=-dpphoton(2)*1.0d9
            fill(12)=-dpphoton(3)*1.0d9
            fill(13)=epho*1.0d9 !eV
            fill(14)=ec*1.0d9 !eV
            call hfm(nidphotons,fill)
          endif
        else
          dpphoton=0.0d0
          dgamma=0.0d0
        endif !(rnrn.le.wrad)

      else !mode
        write(6,*)'*** Error in photon: Invalid mode:', mode,'  ***'
        stop
      endif !mode

      return
      end
