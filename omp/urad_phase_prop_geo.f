*CMZ :          05/09/2024  16.27.42  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  09.37.27  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_geo

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0),a3(3)

      double complex :: expom,apolh,apolr,apoll,apol45

      double precision dx,dx2,dy,dyph,dzph,dz,y,z,omc,domc,phlowz,phlowy,dzy2,eps(6),
     &  dr,drred,da,x,xobs,yobs,zobs,rlambda1,ans,stok1,stok2,stok3,stok4,rn(3)

      integer :: ifrq,iobph,iobs,ieps,iobfr

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
c+seq,uservar.

      nobsvprop_u=npinyprop_u*npinzprop_u

      x=pinxprop_u/1000.0d0

      do ifrq=1,nepho_u

        omc=epho_u(ifrq)/(hbarev1*clight1)

        do iobs=1,nobsv_u

          iobph=iobs+nobsv_u*(ifrq-1)

          xobs=obsv_u(1,iobs)/1000.0d0
          yobs=obsv_u(2,iobs)/1000.0d0
          zobs=obsv_u(3,iobs)/1000.0d0

          dx=xobs-x
          dx2=dx*dx

          rn(1)=real(arad_u(2,iobph)*conjg(arad_u(6,iobph))-arad_u(3,iobph)*conjg(arad_u(5,iobph)))
          rn(2)=real(arad_u(3,iobph)*conjg(arad_u(4,iobph))-arad_u(1,iobph)*conjg(arad_u(6,iobph)))
          rn(3)=real(arad_u(1,iobph)*conjg(arad_u(5,iobph))-arad_u(2,iobph)*conjg(arad_u(4,iobph)))
          rn=rn/norm2(rn)

          y=yobs-rn(2)/rn(1)*dx
          z=zobs-rn(3)/rn(1)*dx

          obsvprop_u(1:3,iobph)=[x,y,z]
          obsvprop_u(4:6,iobph)=rn

          dy=yobs-y
          dz=zobs-z
          dzy2=dz*dz+dy*dy
          da=dx2+dzy2 !convert to solid angle

          if (dzy2.le.0.01d0*dx2) THEN

            eps(1)=dzy2/dx2
            do ieps=2,6
              eps(ieps)=eps(ieps-1)*eps(1)
            enddo !ieps

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
     &        -0.0390625D0*eps(4)+
     &        0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

            dr=dabs(dx*(ans+1.0d0))
            drred=-dabs(dx*ans)

          else
            dr=abs(dx)*sqrt(1.0d0+eps(1))
            drred=dx-dr
          endif

          rlambda1=epho_u(ifrq)/wtoe1*1.0d9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)
          expom=cdexp(dcmplx(0.0d0,drred*omc))/dr

          if (dx.gt.0.0d0) then
            aradprop_u(1:6,iobph)=arad_u(1:6,iobph)*da*expom*rlambda1
          else
            aradprop_u(1:6,iobph)=dconjg(arad_u(1:6,iobph))*da*expom*rlambda1
          endif

        ENDDO   !NFREQ

      ENDDO   !NOBS

      do ifrq=1,nepho_u
        do iobs=1,nobsvprop_u

          iobfr=iobs+nobsvprop_u*(ifrq-1)

          apolh=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(1,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(1,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(1,3))

          apolr=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(2,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(2,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(2,3))

          apoll=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(3,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(3,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(3,3))

          apol45=
     &      aradprop_u(1,iobfr)*dconjg(vstokes(4,1))+
     &      aradprop_u(2,iobfr)*dconjg(vstokes(4,2))+
     &      aradprop_u(3,iobfr)*dconjg(vstokes(4,3))

          stok1=dreal(
     &      apolr*conjg(apolr)+
     &      apoll*conjg(apoll))

          stok2=-stok1+
     &      dreal(2.*apolh*conjg(apolh))

          stok3=
     &      dreal(2.*apol45*conjg(apol45))-
     &      stok1

          stok4=dreal(
     &      apolr*conjg(apolr)-
     &      apoll*conjg(apoll))

          stokesprop_u(1,iobfr)=stok1
          stokesprop_u(2,iobfr)=stok2
          stokesprop_u(3,iobfr)=stok3
          stokesprop_u(4,iobfr)=stok4

        enddo
      enddo

      obsvprop_u=obsvprop_u*1000.0d0

      return
      end
