*CMZ :  4.02/01 15/07/2026  11.29.35  by  Michael Scheer
*-- Author :    Michael Scheer   05/01/2026
      subroutine urad_phase_genpho(nthreads,ith,modegrid,nelec,noranone,
     &  npho,npola,mz,my,mtz,mty,nepho,
     &  zmin,zmax,ymin,ymax,
     &  tzmin,tzmax,tymin,tymax,
     &  ephomin,ephomax,
     &  sigz,sigzp,sigy,sigyp,
     &  photons,electrons,wigner)

      implicit none

      integer :: ith,nelec,noranone,npho,npola(5),mz,my,mtz,mty,nz,ny,ntz,nty,kel,
     &  n1,n2,n1z,n2z,n1y,n2y,n1tz,n2tz,n1ty,n2ty,nele,
     &  iel,iph,iepho,lz,ltz,ly,lty,lpola,ngam,modegrid,nepho,nthreads,na(9),istat

      real eran(nelec*4),rr(2),phran(npho*5),depho,epho(nepho),
     &  z(abs(mz)),y(abs(my)),tz(abs(mtz)),ty(abs(mty)),
     &  dz,dy,dtz,dty,zel,yel,zpel,ypel,zcen,ycen,tzcen,tycen,
     &  zmin,zmax,ymin,ymax,zph,yph,tzph,typh,
     &  tzmin,tzmax,tymin,tymax,
     &  sigz,sigzp,sigy,sigyp,
     &  ephomin,ephomax,zl,yl,tzl,tyl,ephol,
     &  wigner(npola(5),abs(mz),abs(my),abs(mtz),abs(mty),nepho),
     &  photons(7*nepho*npho*nelec*nthreads),
     &  electrons(4*nepho*nelec*nthreads),
     &  a(16),p(9),wint(2,2,2,2),
     &  z2(2),y2(2),ty2(2),tz2(2),w

      nz=iabs(mz)
      ny=iabs(my)
      ntz=iabs(mtz)
      nty=iabs(mty)

      dz=0.0
      zl=zmax-zmin
      if (zl.gt.0.0) then
        dz=zl/(nz-1)
        z(1)=zmin
        do lz=2,nz
          z(lz)=z(lz-1)+dz
        enddo
      else
        z(1)=(zmax+zmin)/2.0
      endif

      zcen=(zmax+zmin)/2.0

      dtz=0.0
      tzl=tzmax-tzmin
      if (tzl.gt.0) then
        dtz=tzl/(ntz-1)
        tz(1)=tzmin
        do ltz=2,ntz
          tz(ltz)=tz(ltz-1)+dtz
        enddo
      else
        tz(1)=(tzmax+tzmin)/2.0
      endif

      tzcen=(tzmax+tzmin)/2.0

      dy=0.0
      yl=ymax-ymin
      if (yl.gt.0.0) then
        dy=yl/(ny-1)
        y(1)=ymin
        do ly=2,ny
          y(ly)=y(ly-1)+dy
        enddo
      else
        y(1)=(ymax+ymin)/2.0
      endif

      ycen=(ymax+ymin)/2.0

      dty=0.0
      tyl=tymax-tymin
      if (tyl.gt.0) then
        dty=tyl/(nty-1)
        ty(1)=tymin
        do lty=2,nty
          ty(lty)=ty(lty-1)+dty
        enddo
      else
        ty(1)=(tymax+tymin)/2.0
      endif

      tycen=(tymax+tymin)/2.0

      depho=0.0
      ephol=ephomax-ephomin
      if (ephol.gt.0) then
        depho=ephol/(nepho-1)
        epho(1)=ephomin
        do iepho=2,nepho
          epho(iepho)=epho(iepho-1)+depho
        enddo
      else
        epho(1)=(ephomax+ephomin)/2.0
      endif

      call util_random_gauss(nelec*4,eran,rr)

      if (modegrid.eq.1) then

        ngam=(ith-1)*7*npho*nelec*nepho
        nele=(ith-1)*4*nelec
        kel=0
        do iel=1,nelec
          zel=sigz*  eran(1+kel)
          zpel=sigzp*eran(2+kel)
          yel=sigy*  eran(3+kel)
          ypel=sigyp*eran(4+kel)
          kel=kel+4
          nele=nele+1
          electrons(nele)=zel
          nele=nele+1
          electrons(nele)=yel
          nele=nele+1
          electrons(nele)=zpel
          nele=nele+1
          electrons(nele)=ypel
          call util_random(npho*5,phran)
          do iph=1,npho
            lz=int(phran(1+(iph-1)*5)*nz)+1
            ly=int(phran(2+(iph-1)*5)*ny)+1
            ltz=int(phran(3+(iph-1)*5)*ntz)+1
            lty=int(phran(4+(iph-1)*5)*nty)+1
            lpola=int(phran(5+(iph-1)*5)*npola(5))+1
            do iepho=1,nepho
              ngam=ngam+1
              photons(ngam)=lpola
              ngam=ngam+1
              photons(ngam)=epho(iepho)
              ngam=ngam+1
              photons(ngam)=z(lz)+zel
              ngam=ngam+1
              photons(ngam)=y(ly)+yel
              ngam=ngam+1
              photons(ngam)=tz(ltz)+zpel
              ngam=ngam+1
              photons(ngam)=ty(lty)+ypel
              ngam=ngam+1
              photons(ngam)=wigner(lpola,lz,ly,ltz,lty,iepho)
            enddo !nepho
          enddo !npho
        enddo !nselec
      else
        ngam=(ith-1)*7*nelec*npho*nepho
        nele=(ith-1)*4*nelec
        kel=0
        na(1:4)=2
        do iel=1,nelec

          if (noranone.ne.0.and.iel.gt.1) then
            zel=sigz*  eran(1+kel)
            zpel=sigzp*eran(2+kel)
            yel=sigy*  eran(3+kel)
            ypel=sigyp*eran(4+kel)
          else
            zel=0.0
            zpel=0.0
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

          kel=kel+4

          call util_random(npho*5,phran)
          !all util_break

          do iph=1,npho

            if (mz.gt.0) then
              zph=phran(1+(iph-1)*5)*zl
            else
              zph=zcen+zl/2.0
            endif

            if (nz.gt.1) then
              n1z=int(zph/dz)+1
              if (n1z.lt.nz) then
                n2z=n1z+1
              else
                n1z=nz-1
                n2z=n1z
              endif
              z2=[z(n1z),z(n2z)]
            else
              n1z=1
              n2z=n1z
              z2=[z(n1z),z(n2z)+1.0]
            endif

            zph=zmin+zph

            if (my.gt.0) then
              yph=phran(2+(iph-1)*5)*yl
            else
              yph=ycen+yl/2.0
            endif

            if (ny.gt.1) then
              n1y=int(yph/dy)+1
              if (n1y.lt.ny) then
                n2y=n1y+1
              else
                n1y=ny-1
                n2y=n1y
              endif
              y2=[y(n1y),y(n2y)]
            else
              n1y=1
              n2y=n1y
              y2=[y(n1y),y(n2y)+1.0]
            endif

            yph=ymin+yph

            if (mtz.gt.0) then
              tzph=phran(3+(iph-1)*5)*tzl
            else
              tzph=tzcen+tzl/2.0
            endif

            if (ntz.gt.1) then
              n1tz=int(tzph/dtz)+1
              if (n1tz.lt.ntz) then
                n2tz=n1tz+1
              else
                n1tz=ntz-1
                n2tz=n1tz
              endif
              tz2=[tz(n1tz),tz(n2tz)]
            else
              n1tz=1
              n2tz=n1tz
              tz2=[tz(n1tz),tz(n2tz)+1.0]
            endif

            tzph=tzmin+tzph

            if (mty.gt.0) then
              typh=phran(4+(iph-1)*5)*tyl
            else
              typh=tycen+tyl/2.0
            endif

            if (nty.gt.1) then
              n1ty=int(typh/dty)+1
              if (n1ty.lt.nty) then
                n2ty=n1ty+1
              else
                n1ty=nty-1
                n2ty=n1ty
              endif
              ty2=[ty(n1ty),ty(n2ty)]
            else
              n1ty=1
              n2ty=n1ty
              ty2=[ty(n1ty),ty(n2ty)+1.0]
            endif

            typh=tymin+typh

            lpola=int(phran(5+(iph-1)*5)*npola(5))+1

            do iepho=1,nepho

              p(1:4)=[zph,yph,tzph,typh]

              wint(1,1,1,1)=wigner(lpola,n1z,n1y,n1tz,n1ty,iepho)
              wint(2,1,1,1)=wigner(lpola,n2z,n1y,n1tz,n1ty,iepho)
              wint(1,2,1,1)=wigner(lpola,n1z,n2y,n1tz,n1ty,iepho)
              wint(2,2,1,1)=wigner(lpola,n2z,n2y,n1tz,n1ty,iepho)

              wint(1,1,2,1)=wigner(lpola,n1z,n1y,n2tz,n1ty,iepho)
              wint(2,1,2,1)=wigner(lpola,n2z,n1y,n2tz,n1ty,iepho)
              wint(1,2,2,1)=wigner(lpola,n1z,n2y,n2tz,n1ty,iepho)
              wint(2,2,2,1)=wigner(lpola,n2z,n2y,n2tz,n1ty,iepho)

              wint(1,1,1,2)=wigner(lpola,n1z,n1y,n1tz,n2ty,iepho)
              wint(2,1,1,2)=wigner(lpola,n2z,n1y,n1tz,n2ty,iepho)
              wint(1,2,1,2)=wigner(lpola,n1z,n2y,n1tz,n2ty,iepho)
              wint(2,2,1,2)=wigner(lpola,n2z,n2y,n1tz,n2ty,iepho)

              wint(1,1,2,2)=wigner(lpola,n1z,n1y,n2tz,n2ty,iepho)
              wint(2,1,2,2)=wigner(lpola,n2z,n1y,n2tz,n2ty,iepho)
              wint(1,2,2,2)=wigner(lpola,n1z,n2y,n2tz,n2ty,iepho)
              wint(2,2,2,2)=wigner(lpola,n2z,n2y,n2tz,n2ty,iepho)

              call util_linear_inter_4d_real(2,2,2,2,z2,y2,tz2,ty2,wint,zph,yph,tzph,typh,
     &          w,istat)

              ngam=ngam+1
              photons(ngam)=lpola
              ngam=ngam+1
              photons(ngam)=epho(iepho)
              ngam=ngam+1
              photons(ngam)=zph+zel
              ngam=ngam+1
              photons(ngam)=yph+yel
              ngam=ngam+1
              photons(ngam)=tzph+zpel
              ngam=ngam+1
              photons(ngam)=typh+ypel
              ngam=ngam+1
              photons(ngam)=w

            enddo !nepho
          enddo !npho
        enddo !nelec
      endif

      end
