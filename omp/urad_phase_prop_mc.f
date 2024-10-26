*CMZ :          06/08/2024  14.45.08  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  21.58.40  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_mc(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double complex, dimension(:), allocatable :: expom,dexpom,phshift
      double complex :: apolh,apolr,apoll,apol45

      double precision, dimension(:,:), allocatable :: obsv
      double precision dx,dx2,dy,dyph,dzph,dz,y,z,omc,domc,phlowz,phlowy,dzy2,eps(6),
     &  dr,drred,da,x,xobs,yobs,zobs,rlambda1,ans,stok1,stok2,stok3,stok4,zmin,ymin,fsum

      integer :: mthreads,iy,iz,n,jy,jz,iobs,ieps,ifrq,iobfr,jobs,jobfr,nobsv,iepho

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      nobsv=(npinzo_u-1)*(npinyo_u-1)
      nobsvprop_u=npinyprop_u*npinzprop_u

      allocate(expom(nobsv*nepho_u),dexpom(nobsv*nepho_u),phshift(nobsv),obsv(3,nobsv))

      ymin=pincen_u(2)-pinh_u/2.0d0
      zmin=pincen_u(3)-pinw_u/2.0d0

      dy=pinh_u/max(1,npinyo_u-1)
      dz=pinw_u/max(1,npinzo_u-1)

      iobs=0
      do iy=1,npinyo_u-1
        do iz=1,npinzo_u-1
          iobs=iobs+1
          obsv(1,iobs)=pincen_u(1)
          obsv(2,iobs)=ymin+(dble(iy)-0.5)*dy
          obsv(3,iobs)=zmin+(dble(iz)-0.5)*dz
        enddo
      enddo

      aradprop_u=(0.0d0,0.0d0)

      if (npinyprop_u.gt.1) then
        dyph=pinhprop_u/1000.0d0/dble(npinyprop_u-1)
        phlowy=-pinhprop_u/2000.0d0
      else
        dyph=pinhprop_u/1000.0d0
        phlowy=-pinhprop_u/1000.0d0
      endif

      if (npinzprop_u.gt.1) then
        dzph=pinwprop_u/1000.0d0/dble(npinzprop_u-1)
        phlowz=-pinwprop_u/2000.0d0
      else
        dzph=pinwprop_u/1000.0d0
        phlowz=-pinwprop_u/1000.0d0
      endif

      da=dz*dy

      n=0

      x=pinxprop_u/1000.0d0
      y=phlowy-dyph

      do iy=1,npinyprop_u
        y=y+dyph
        if (abs(y).lt.1.0d-12) y=0.0d0
        z=phlowz-dzph
        do iz=1,npinzprop_u
          n=n+1
          z=z+dzph
          if (abs(z).lt.1.0d-12) z=0.0d0
          obsvprop_u(1:3,n)=[x,y,z]
        enddo
      enddo

      omc=epho_u(1)/(hbarev1*clight1)
      if(nepho_u.gt.1) then
        domc=(epho_u(2)-epho_u(1))/(hbarev1*clight1)
      endif

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(domc,omc,da,obsvprop_u,obsv,npinzo_u,npinyo_u,nobsv,nobsvprop_u,epho_u,nepho_u,aradprop_u,fieldbunch)

!$OMP DO

      do jobs=1,nobsvprop_u
c        ith=OMP_GET_THREAD_NUM()+1

        x=obsvprop_u(1,jobs)
        y=obsvprop_u(2,jobs)
        z=obsvprop_u(3,jobs)

        DO IOBS=1,nobsv

          XOBS=obsv(1,IOBS)
          YOBS=obsv(2,IOBS)
          ZOBS=obsv(3,IOBS)

          dx=xobs-x
          dx2=dx*dx
          DY=YOBS-y
          DZ=ZOBS-z
          DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

          IF (DZY2.GT.0.01D0*dx2) THEN
            WRITE(6,*)'*** ERROR IN URAD_phase_prop_mc ***'
            WRITE(6,*)'CHECK INPUT FILE AND INCREASE PinX'
            WRITE(6,*)'*** PROGRAM ABORTED ***'
            STOP
          ENDIF

          EPS(1)=DZY2/dx2
          DO IEPS=2,6
            EPS(IEPS)=EPS(IEPS-1)*EPS(1)
          ENDDO !IEPS

c      TAYLOR-EXPANSION DONE WITH REDUCE
c     IN "WTAY1.RED";
c     on rounded;
c     on numval;
c     precision 13;
c     F:=SQRT(1+EPS);
c     DR:=TAY1(F,EPS,6);
c     ON FORT;
c     OUT "RED.FOR";
c     DR;
c     SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

          ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &      -0.0390625D0*eps(4)+
     &      0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

          DR=DABS(dx*(ANS+1.0D0))
          DRRED=-DABS(dx*ANS)

          IF (DR.NE.0.0d0) THEN
            EXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*OMC))/DR
          ELSE
            EXPOM(IOBS)=1.0D0
          ENDIF

          if (nepho_u.gt.1) then
            DEXPOM(IOBS)=CDEXP(DCMPLX(0.0d0,DRRED*DOMC))
          endif
c            print*,ith,iobs,expom(iobs)
c+seq,dum2.
        ENDDO   !NOBS

        DO ifrq=1,nepho_u

          jOBFR=jOBS+NOBSVprop_u*(ifrq-1)

          RLAMBDA1=epho_u(ifrq)/WTOE1*1.0D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)

          iobs=0
          do iy=1,npinyo_u-1
            do iz=1,npinzo_u-1
              iobs=iobs+1

              IOBFR=IOBS+nobsv*(ifrq-1)

              IF (ifrq.EQ.1) THEN
                PHSHIFT(IOBS)=EXPOM(IOBFR)
              ELSE
                PHSHIFT(IOBS)=PHSHIFT(IOBS)*DEXPOM(IOBS)
              ENDIF   !(ifrq.EQ.1)

              if (dx.gt.0.0d0) then
                aradprop_u(1:6,jobfr)=aradprop_u(1:6,jobfr)+
     &            fieldbunch(1:6,iz,iy,ifrq)*PHSHIFT(IOBS)*da*rlambda1
              else
                aradprop_u(1:6,jobfr)=aradprop_u(1:6,jobfr)+
     &            dconjg(fieldbunch(1:6,iz,iy,ifrq))*PHSHIFT(IOBS)*da*rlambda1
              endif
            ENDDO   !NFREQ

          ENDDO  !NOBSVz
        ENDDO  !NOBSVy

      ENDDO !nobsvprop_u

!$OMP END DO

!$OMP END PARALLEL

      if (ifieldprop_u.ne.2) then

        do ifrq=1,nepho_u
          do jobs=1,nobsvprop_u

            jobfr=jobs+nobsvprop_u*(ifrq-1)

            apolh=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(1,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(1,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(1,3))

            apolr=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(2,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(2,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(2,3))

            apoll=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(3,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(3,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(3,3))

            apol45=
     &        aradprop_u(1,jobfr)*dconjg(vstokes(4,1))+
     &        aradprop_u(2,jobfr)*dconjg(vstokes(4,2))+
     &        aradprop_u(3,jobfr)*dconjg(vstokes(4,3))

            stok1=dreal(
     &        apolr*conjg(apolr)+
     &        apoll*conjg(apoll))

            stok2=-stok1+
     &        dreal(2.*apolh*conjg(apolh))

            stok3=
     &        dreal(2.*apol45*conjg(apol45))-
     &        stok1

            stok4=dreal(
     &        apolr*conjg(apolr)-
     &        apoll*conjg(apoll))

            stokesprop_u(1,jobfr)=stok1
            stokesprop_u(2,jobfr)=stok2
            stokesprop_u(3,jobfr)=stok3
            stokesprop_u(4,jobfr)=stok4

          enddo
        enddo

      endif !(ifieldprop_u.ne.2)

      obsvprop_u=obsvprop_u*1000.0d0

      return
      end
