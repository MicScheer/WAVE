*CMZ :          23/08/2026  09.58.03  by  Michael Scheer
*-- Author :    Michael Scheer   05/01/2026
        subroutine urad_phase_amp_genpho(zi,yi,ny,nz,obsv,
     &    moderan,nelec,noranone,npho,nepho,
     &    ephomin,ephomax,
     &    sigz,sigzp,sigy,sigyp,
     &    ndimpho,photons,ndimele,electrons,
     &    arad,specnor_si)

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      double complex arad(6,nz*ny*nepho),ara(6),e(3),b(3),crn(3)

      double complex :: apol,amp0(6),damp(6),amp(6),ampn(6),zexp,
     &  apolh,apolr,apoll,apol45,cero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)
     &  ,cjvsto(4,3)

      real*8 :: zi,yi,obsv(3,nz*ny),rn(3),
     &  z(nz),y(ny),
     &  exr(nz,ny,nepho),
     &  exi(nz,ny,nepho),
     &  eyr(nz,ny,nepho),
     &  eyi(nz,ny,nepho),
     &  ezr(nz,ny,nepho),
     &  ezi(nz,ny,nepho),
     &  bxr(nz,ny,nepho),
     &  bxi(nz,ny,nepho),
     &  byr(nz,ny,nepho),
     &  byi(nz,ny,nepho),
     &  bzr(nz,ny,nepho),
     &  bzi(nz,ny,nepho),
     &  erx,eix,ery,eiy,erz,eiz,
     &  brx,bix,bry,biy,brz,biz,
     &  zel,yel,zpel,ypel,
     &  zph,yph,specnor_si,
     &  c=3.0d8,stok1,stok2,stok3,stok4

      integer :: ny,nz,iy,iz,ly,lz,iepho,moderan,nelec,noranone,npho,nepho,nobsv,ngam,nele,
     &  n2z,n2y,n1z,n1y,kel,iel,iph,l,kpho,ndimpho,ndimele,ndimapho=9

      real eran(nelec*4),rr(2),phran(npho*2),depho,epho(nepho),
     &  sigz,sigzp,sigy,sigyp,
     &  ephomin,ephomax,zl,yl,tzl,tyl,ephol,
     &  photons(ndimpho),
     &  electrons(ndimele),
     &  aradr(6,nz*ny*nepho),
     &  a(16),p(9),
     &  z2(2),y2(2),ty2(2),tz2(2)

      callutil_break

      cjvsto=dconjg(vstokes)

      l=0
      do iepho=1,nepho
        do iy=1,ny
          do iz=1,nz
            l=l+1
            exr(iz,iy,iepho)=dreal(arad(1,l))*c
            exi(iz,iy,iepho)=dimag(arad(1,l))*c
            eyr(iz,iy,iepho)=dreal(arad(2,l))*c
            eyi(iz,iy,iepho)=dimag(arad(2,l))*c
            ezr(iz,iy,iepho)=dreal(arad(3,l))*c
            ezi(iz,iy,iepho)=dimag(arad(3,l))*c
            bxr(iz,iy,iepho)=dreal(arad(4,l))*c
            bxi(iz,iy,iepho)=dimag(arad(4,l))*c
            byr(iz,iy,iepho)=dreal(arad(5,l))*c
            byi(iz,iy,iepho)=dimag(arad(5,l))*c
            bzr(iz,iy,iepho)=dreal(arad(6,l))*c
            bzi(iz,iy,iepho)=dimag(arad(6,l))*c
          enddo
        enddo
      enddo

      nobsv=nz*ny

      z(1:nz)=obsv(3,1:nz)

      ly=0
      do iy=1,ny*nz,nz
        ly=ly+1
        y(ly)=obsv(2,iy)
      enddo

      depho=0.0
      if (nepho.gt.1) then
        ephol=ephomax-ephomin
        depho=ephol/(nepho-1)
        epho(1)=ephomin
        do iepho=2,nepho
          epho(iepho)=epho(iepho-1)+depho
        enddo
      else
        epho(1)=(ephomax+ephomin)/2.0
      endif

      !allutil_break

      call util_random_gauss(nelec*4,eran,rr)

      nele=0
      ngam=0
      kel=0

      do iel=1,nelec

        if (noranone.ne.0.and.iel.gt.1) then
          kel=kel+1
          zpel=sigzp*eran(kel)
          kel=kel+1
c         zel=zi+sigz*eran(kel)+zpel*obsv(1,nobsv/2+1)
          zel=sigz*eran(kel)+zpel*obsv(1,nobsv/2+1)
          kel=kel+1
          ypel=sigyp*eran(kel)
          kel=kel+1
c          yel=yi+sigy*eran(kel)+ypel*obsv(1,nobsv/2+1)
          yel=sigy*eran(kel)+ypel*obsv(1,nobsv/2+1)
        else
c          zel=zi
          zel=0.0
          zpel=0.0
c          yel=yi
          yel=0.0
          ypel=0.0
        endif

        nele=nele+1
        electrons(nele)=zel
        nele=nele+1
        electrons(nele)=yel
        nele=nele+1
        electrons(nele)=zpel
        nele=nele+1
        electrons(nele)=ypel

      enddo

      if (moderan.eq.0) then

        nele=0

        do iel=1,nelec

          nele=nele+1
          zel=electrons(nele)
          nele=nele+1
          yel=electrons(nele)
          nele=nele+1
          zpel=electrons(nele)
          nele=nele+1
          ypel=electrons(nele)

          do iy=1,ny

            yph=y(iy)

            do iz=1,nz

              zph=z(iz)

              do iepho=1,nepho

                erx=exr(iz,iy,iepho)
                eix=exi(iz,iy,iepho)
                ery=eyr(iz,iy,iepho)
                eiy=eyi(iz,iy,iepho)
                erz=ezr(iz,iy,iepho)
                eiz=ezi(iz,iy,iepho)

                brx=bxr(iz,iy,iepho)
                bix=bxi(iz,iy,iepho)
                bry=byr(iz,iy,iepho)
                biy=byi(iz,iy,iepho)
                brz=bzr(iz,iy,iepho)
                biz=bzi(iz,iy,iepho)

                ara(1)=dcmplx(erx,eix)
                ara(2)=dcmplx(ery,eiy)
                ara(3)=dcmplx(erz,eiz)
                ara(4)=dcmplx(brx,bix)
                ara(5)=dcmplx(bry,biy)
                ara(6)=dcmplx(brz,biz)

                rn(1)=real(ara(2)*conjg(ara(6))-ara(3)*conjg(ara(5)))
                rn(2)=real(ara(3)*conjg(ara(4))-ara(1)*conjg(ara(6)))
                rn(3)=real(ara(1)*conjg(ara(5))-ara(2)*conjg(ara(4)))

                if (rn(1).eq.0.0d0) then
                  ngam=ngam+ndimapho
                  cycle
                endif


c                rn=rn/norm2(rn)
c                rn(2)=rn(2)/rn(1)
c                rn(3)=rn(3)/rn(1)
c                rn(1)=sqrt(1.0d0-(rn(2)**2+rn(3)**2))

                ara=ara/c

                apolh=
     &            ara(1)*cjvsto(1,1)
     &            +ara(2)*cjvsto(1,2)
     &            +ara(3)*cjvsto(1,3)

                apolr=
     &            ara(1)*cjvsto(2,1)
     &            +ara(2)*cjvsto(2,2)
     &            +ara(3)*cjvsto(2,3)

                apoll=
     &            ara(1)*cjvsto(3,1)
     &            +ara(2)*cjvsto(3,2)
     &            +ara(3)*cjvsto(3,3)

                apol45=
     &            ara(1)*cjvsto(4,1)
     &            +ara(2)*cjvsto(4,2)
     &            +ara(3)*cjvsto(4,3)

                stok1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
                stok2=dreal(-stok1+2.0d0*apolh*conjg(apolh))
                stok3=dreal(2.0d0*apol45*conjg(apol45)-stok1)
                stok4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

                ngam=ngam+1
                photons(ngam)=epho(iepho)
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                ngam=ngam+1
                photons(ngam)=zph+zel
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                print*,'z',iz,iy,ngam,zph
                ngam=ngam+1
                photons(ngam)=yph+yel
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                print*,'y',iz,iy,ngam,yph
                ngam=ngam+1
                photons(ngam)=rn(3)/rn(1)+zpel
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                ngam=ngam+1
                photons(ngam)=rn(2)/rn(1)+ypel
                ngam=ngam+1
                photons(ngam)=stok1*specnor_si
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                ngam=ngam+1
                photons(ngam)=stok2*specnor_si
                ngam=ngam+1
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                photons(ngam)=stok3*specnor_si
                if (photons(ngam).ne.photons(ngam)) print*,ngam
                ngam=ngam+1
                photons(ngam)=stok4*specnor_si
                if (photons(ngam).ne.photons(ngam)) print*,ngam

              enddo !nepho

            enddo !nz
          enddo !ny

        enddo !nelec

      else

        nele=0

        do iel=1,nelec

          nele=nele+1
          zel=electrons(nele)
          nele=nele+1
          yel=electrons(nele)
          nele=nele+1
          zpel=electrons(nele)
          nele=nele+1
          ypel=electrons(nele)

          call util_random(npho*2,phran)
c          phran(1:2)=0.0

          kpho=0

          do iph=1,npho

            kpho=kpho+1
            iz=int(phran(kpho)*nz)+1
            kpho=kpho+1
            iy=int(phran(kpho)*ny)+1

c            if (iph.eq.1) then
c              iz=nz/2+1
c              iy=ny/2+1
c            endif

            yph=y(iy)
            zph=z(iz)

            do iepho=1,nepho

              erx=exr(iz,iy,iepho)
              eix=exi(iz,iy,iepho)
              ery=eyr(iz,iy,iepho)
              eiy=eyi(iz,iy,iepho)
              erz=ezr(iz,iy,iepho)
              eiz=ezi(iz,iy,iepho)

              brx=bxr(iz,iy,iepho)
              bix=bxi(iz,iy,iepho)
              bry=byr(iz,iy,iepho)
              biy=byi(iz,iy,iepho)
              brz=bzr(iz,iy,iepho)
              biz=bzi(iz,iy,iepho)

              ara(1)=dcmplx(erx,eix)
              ara(2)=dcmplx(ery,eiy)
              ara(3)=dcmplx(erz,eiz)
              ara(4)=dcmplx(brx,bix)
              ara(5)=dcmplx(bry,biy)
              ara(6)=dcmplx(brz,biz)

              rn(1)=real(ara(2)*conjg(ara(6))-ara(3)*conjg(ara(5)))
              rn(2)=real(ara(3)*conjg(ara(4))-ara(1)*conjg(ara(6)))
              rn(3)=real(ara(1)*conjg(ara(5))-ara(2)*conjg(ara(4)))


              ara=ara/c

c                rn=rn/norm2(rn)
c                rn(2)=rn(2)/rn(1)
c                rn(3)=rn(3)/rn(1)
c                rn(1)=sqrt(1.0d0-(rn(2)**2+rn(3)**2))

              apolh=
     &          ara(1)*cjvsto(1,1)
     &          +ara(2)*cjvsto(1,2)
     &          +ara(3)*cjvsto(1,3)

              apolr=
     &          ara(1)*cjvsto(2,1)
     &          +ara(2)*cjvsto(2,2)
     &          +ara(3)*cjvsto(2,3)

              apoll=
     &          ara(1)*cjvsto(3,1)
     &          +ara(2)*cjvsto(3,2)
     &          +ara(3)*cjvsto(3,3)

              apol45=
     &          ara(1)*cjvsto(4,1)
     &          +ara(2)*cjvsto(4,2)
     &          +ara(3)*cjvsto(4,3)

              stok1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
              stok2=dreal(-stok1+2.0d0*apolh*conjg(apolh))
              stok3=dreal(2.0d0*apol45*conjg(apol45)-stok1)
              stok4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

              ngam=ngam+1
              photons(ngam)=epho(iepho)
              ngam=ngam+1
              photons(ngam)=zph+zel
              ngam=ngam+1
              photons(ngam)=yph+yel
              ngam=ngam+1
              photons(ngam)=rn(3)/rn(1)+zpel
              ngam=ngam+1
              photons(ngam)=rn(2)/rn(1)+ypel
c                ngam=ngam+1
c                photons(ngam)=(erx**2+eix**2+ery**2+eiy**2+erz**2+eiz**2)*specnor_si
              ngam=ngam+1
              photons(ngam)=stok1*specnor_si
              ngam=ngam+1
              photons(ngam)=stok2*specnor_si
              ngam=ngam+1
              photons(ngam)=stok3*specnor_si
              ngam=ngam+1
              photons(ngam)=stok4*specnor_si

            enddo !nepho

          enddo !npho

        enddo !nelec

      endif !moderan

      callutil_break
      end
