*CMZ :          01/08/2024  10.22.17  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_wigner

      use omp_lib
      use uradphasemod
      use wignermod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      complex*16 :: wig2d(nzwig,nywig,nzthewig,nythewig)
      complex*16 :: wkern(npinz_u,npiny_u,npinz_u,npiny_u),
     &  fkt(npinz_u,npiny_u),fkr(npinz_u,npiny_u)

      real*8 ::
     &  er2d(npinzprop_u,npinyprop_u),
     &  ei2d(npinzprop_u,npinyprop_u)

      real*8 ::
     &  ekr2d(npinz_u,npiny_u),
     &  eki2d(npinz_u,npiny_u),
     &  thez2(npinz_u),they2(npiny_u)

      real*8 ::
     &  z(max(npinz_u,npinzprop_u)),y(max(npinz_u,npinyprop_u)),
     &  t(max(nythewig,nzthewig)),
     &  aradr(2,npinzprop_u,npinyprop_u,nepho_u),
     &  aradi(2,npinzprop_u,npinyprop_u,nepho_u)

      real*8 om,dthe

      integer :: ixy,kx,ky,ktime=1,iepho,istat,nx,nt,iz,iwy,iwz,i1,i2,iwt,iobsv,
     &  i,iy,iyw,iytw,izw,iztw,iobs,iobph
     &  ,izleft,izright,iylow,iyhigh,ifail,itz,ity,lz,ly

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase_wigner',1)

        allocate(
     &    wigkr(npinz_u,npiny_u,npinz_u,npiny_u,nepho_u,4),
     &    wigki(npinz_u,npiny_u,npinz_u,npiny_u,nepho_u,4)
     &    )

      wigkr=0.0d0
      wigki=0.0d0

      allocate(
     &  wigr(nzwig,nywig,nzthewig,nythewig,nepho_u,4),
     &  wigi(nzwig,nywig,nzthewig,nythewig,nepho_u,4))

      wigr=0.0d0
      wigi=0.0d0

      if (ifieldprop_u.eq.0) then

        dthe=pinwwig/1000.0d0/(max(1,nzwig-1))
        z(1)=-pinwwig/2.0d0
        do iz=2,nzwig
          z(iz)=z(iz-1)+dthe
          if (abs(z(iz)).lt.1.0d-12) z(iz)=0.0d0
        enddo

        dthe=pinhwig/(max(1,nywig-1))
        y(1)=-pinhwig/2.0d0
        do iy=2,nywig
          y(iy)=y(iy-1)+dthe
          if (abs(y(iy)).lt.1.0d-12) y(iy)=0.0d0
        enddo

        dthe=thezwig/(max(1,nzthewig-1))
        thez2(1)=-thezwig/2.0d0
        do iz=2,nzthewig
          thez2(iz)=thez2(iz-1)+dthe
          if (abs(thez2(iz)).lt.1.0d-12) thez2(iz)=0.0d0
        enddo

        dthe=theywig/(max(1,nythewig-1))
        they2(1)=-theywig/2.0d0
        do iy=2,nythewig
          they2(iy)=they2(iy-1)+dthe
          if (abs(they2(iy)).lt.1.0d-12) they2(iy)=0.0d0
        enddo

        do iepho=1,nepho_u

          iobs=0
          do iz=1,npinz_u
            do iy=1,npiny_u

              iobs=iobs+1
              iobph=iobs+nobsv_u*(iepho-1)

              ekr2d(iz,iy)=dreal(arad_u(3,iobph))
              eki2d(iz,iy)=dimag(arad_u(3,iobph))

            enddo
          enddo

          call util_wigner_2d_kernel(npinz_u,npiny_u,ekr2d,eki2d,wkern,istat)

          wigkr(1:npinz_u,1:npiny_u,1:npinz_u,1:npiny_u,iepho,1)=
     &      dreal(wkern(1:npinz_u,1:npiny_u,1:npinz_u,1:npiny_u))
          wigki(1:npinz_u,1:npiny_u,1:npinz_u,1:npiny_u,iepho,1)=
     &      dimag(wkern(1:npinz_u,1:npiny_u,1:npinz_u,1:npiny_u))

          do itz=1,nzthewig
            do ity=1,nythewig

              fkt(1:npinz_u,1:npiny_u)=wkern(itz,ity,1:npinz_u,1:npiny_u)

c              if (itz.eq.nzthewig/2+1.and.ity.eq.nythewig/2+1) then
c                print*,'hallo'
c              endif

              call util_fourier_linear_complex_2d(npinz_u,npiny_u,thez2,they2,
     &          fkt,npinz_u,npiny_u,z,y,fkr,ifail)

c              if (itz.eq.nzthewig/2+1.and.ity.eq.nythewig/2+1) then
c                ly=npiny_u/2+1
c                do lz=1,npinz_u
cc                  do ly=1,npiny_u
cc                  write(77,*)lz,thez2(lz),z(lz),dreal(fkt(lz,ly)),
c                  write(77,*)z(lz),dreal(fkt(lz,ly)),
c     &              dimag(fkt(lz,ly)),dreal(fkr(lz,ly)),dimag(fkr(lz,ly))
cc                  enddo
c                enddo
c                stop
c              endif

              wigr(1:npinz_u,1:npiny_u,itz,ity,iepho,1)=dreal(fkr(1:npinz_u,1:npiny_u))*4.0d0
              wigi(1:npinz_u,1:npiny_u,itz,ity,iepho,1)=dimag(fkr(1:npinz_u,1:npiny_u))*4.0d0

            enddo
          enddo

        enddo

        goto 9999
      endif

      iobsv=0
      do iepho=1,nepho_u
        iy=1
        iz=0
        do i=1,nobsvprop_u
          iz=iz+1
          if (iz.gt.npinzprop_u) then
            iz=1
            iy=iy+1
          endif
          iobsv=iobsv+1
          aradr(1:2,iz,iy,iepho)=dreal(aradprop_u(2:3,iobsv))
          aradi(1:2,iz,iy,iepho)=dimag(aradprop_u(2:3,iobsv))
        enddo
      enddo

      z(1:npinzprop_u)=obsvzprop_u(1:npinzprop_u)
      y(1:npinyprop_u)=obsvyprop_u(1:npinyprop_u)

      do iepho=1,nepho_u

        om=epho_u(iepho)/hbarev1

        er2d(:,:)=aradr(2,:,:,iepho)
        ei2d(:,:)=aradi(2,:,:,iepho)

        call util_wigner_2d(
     &    npinzprop_u,nzwig,nzthewig,z,wigthez,
     &    npinyprop_u,nywig,nythewig,y,wigthey,
     &    om,er2d,ei2d,wig2d,mthreads_u,istat)
        wigr(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,1)=
     &    dreal(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))
        wigi(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,1)=
     &    dimag(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))

        er2d(:,:)=aradr(2,:,:,iepho)
        ei2d(:,:)=aradi(1,:,:,iepho)

        call util_wigner_2d(
     &    npinzprop_u,nzwig,nzthewig,z,wigthez,
     &    npinyprop_u,nywig,nythewig,y,wigthey,
     &    om,er2d,ei2d,wig2d,mthreads_u,istat)
        wigr(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,2)=
     &    dreal(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))
        wigi(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,2)=
     &    dimag(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))

        er2d(:,:)=aradr(1,:,:,iepho)
        ei2d(:,:)=aradi(2,:,:,iepho)

        call util_wigner_2d(
     &    npinzprop_u,nzwig,nzthewig,z,wigthez,
     &    npinyprop_u,nywig,nythewig,y,wigthey,
     &    om,er2d,ei2d,wig2d,mthreads_u,istat)
        wigr(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,3)=
     &    dreal(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))
        wigi(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,3)=
     &    dimag(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))

        er2d(:,:)=aradr(1,:,:,iepho)
        ei2d(:,:)=aradi(1,:,:,iepho)

        call util_wigner_2d(
     &    npinzprop_u,nzwig,nzthewig,z,wigthez,
     &    npinyprop_u,nywig,nythewig,y,wigthey,
     &    om,er2d,ei2d,wig2d,mthreads_u,istat)
        wigr(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,4)=
     &    dreal(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))
        wigi(1:nzwig,1:nywig,1:nzthewig,1:nythewig,iepho,4)=
     &    dimag(wig2d(1:nzwig,1:nywig,1:nzthewig,1:nythewig))

      enddo !nepho

9999  if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase_wigner',0)

      return
      end
