*CMZ :  2.57/04 01/02/2006  10.29.42  by  Michael Scheer
*CMZ :  2.57/00 03/11/2005  16.29.31  by  Michael Scheer
*CMZ :  2.56/00 11/10/2005  15.29.20  by  Michael Scheer
*-- Author :    Michael Scheer   10/10/2005
      module wlsf90m

      double precision
     &  px1,px2,b0,xl,zl,freqlow,freqhig,ebeam,brho,
     &  x01,x02,xp01,xp02,pz,ec,b0x,epspx,aperl,aperr,aperphi0

      integer nfreq,nx0,nxp0,icirc

      namelist /wlsn/
     &  nfreq,nx0,nxp0,icirc,
     &  px1,px2,pz,b0,xl,zl,freqlow,freqhig,ebeam,
     &  x01,x02,xp01,xp02,epspx,aperl,aperr,aperphi0

      end module

      program wlscosb

      use wlsf90m

      implicit none

      include '/home/scheer/util/util_phycon_incl.f'

      double precision
     &  dx0,dxp0,
     &  x,xp,z,b,fluxden,xp1,xp2,
     &  xk,zk,x0,xp0,x1,x2,
     &  z1,z2,f1,f2,fdeni,freq,dfreq,fluxdeni1,fluxdeni2,b1,b2,
     &  px1true,px2true,pxtrue,pxeps,zeps,pxfinal,fluxdeni0,fluxdeni00,
     &  corraper,corraper1,corraper2,
     &  efilt(10000),filt(10000),facfilt,
     &  filt2p(10000),
     &  ws2(10000),
     &  ws3(10000),
     &  ws4(10000),
     &  ws1(10000)

      integer ix0,ixp0,ifreq,iz,icount,luno,ifilt,nfilt,mode

      character(256) comfilt

      data luno/20/

      open(unit=99,file='wlscosb.filter',status='old')
      read(99,'(A)')comfilt
      read(99,*)nfilt
      do ifilt=1,nfilt
        read(99,*)efilt(ifilt),filt(ifilt)
      enddo
      close(99)

      call UTIL_SPLINE_COEF(efilt,filt,nfilt,0.0d0,0.0d0,filt2p,ws1,ws2,ws3,ws4)

      open(unit=luno,file='wlscosb.out',status='unknown')
      open(unit=99,file='wlscosb.nam',status='old',readonly)
      read(99,wlsn)
      close(99)

      if (xl.ne.0.0d0) then
        xk=2.0d0*pi/xl
      else
        xk=0.0d0
      endif

      if (zl.ne.0.0d0) then
        zk=2.0d0*pi/zl
      else
        stop '*** Error in WLSCOSB: ZL must not be zero!'
      endif

      nx0=abs(nx0)
      nxp0=abs(nxp0)
      nfreq=abs(nfreq)
      freqlow=abs(freqlow)
      freqhig=abs(freqhig)

      if (nfreq.le.1) then
        stop '**** Error in WLSCOSB: nfreq.le.1'
      endif

      if (freqhig.le.freqlow) then
        stop '**** Error in WLSCOSB: freqhig.le.freqlow'
      endif

      if (nfreq.gt.1) then
        dfreq=(freqhig-freqlow)/(nfreq-1)
      else
        dfreq=0.0d0
      endif

      if (nx0.gt.1) then
        dx0=(x02-x01)/(nx0-1)
      else
        dx0=0.0d0
      endif

      if (nxp0.gt.1) then
        dxp0=(xp02-xp01)/(nxp0-1)
      else
        dxp0=0.0d0
      endif

      brho=-ebeam*1.0d9/clight
      ec=665.*ebeam**2*b0

      print*
      print*
      print*,'Program WLSCOSB'
      print*
      print*,'Comment on filter file:'
      print*,comfilt
      print*
      print*,'z longitudinal, x horizontal coords'
      print*,'Approximation: B0(x) = B(x0)'
      print*
      print*,'EBEAM, B0: ',ebeam,b0
      print*,'Lx, Lz:', xl, zl
      print*,'kx, kz:', xk, zk
      print*,''
      print*,'Pz:',pz
      print*,'Px1, Px2:',px1,px2
      print*,''
      print*,'icirc:',icirc
      if (icirc.eq.0) then
        print*,'AperL,AperHW:',aperl,aperr
      else
        print*,'AperL,AperR:',aperl,aperr
      endif
      print*,'AperPhi0:',aperphi0
      print*,''
      print*,'Ec: ',ec
      print*,'Eg1, Eg2: ',freqlow/ec,freqhig/ec
      print*,''

      write(luno,*)
      write(luno,*)
      write(luno,*)'Program WLSCOSB'
      write(luno,*)
      write(luno,*)'Approximation: B0(x) = B(x0)'
      write(luno,*)
      write(luno,*)'EBEAM, B0: ',ebeam,b0
      write(luno,*)'Lx, Lz:', xl, zl
      write(luno,*)'kx, kz:', xk, zk
      write(luno,*)''
      write(luno,*)'Pz:',pz
      write(luno,*)'Px1, Px2:',px1,px2
      write(luno,*)''
      write(luno,*)'AperL,AperR:',aperl,aperr
      write(luno,*)'AperPhi0:',aperphi0
      write(luno,*)''
      write(luno,*)'Ec: ',ec
      write(luno,*)'Eg1, Eg2: ',freqlow/ec,freqhig/ec
      write(luno,*)''

      do iz=1,2

        do ix0=1,nx0

          x0=x01+(ix0-1)*dx0

          b0x=b0*cos(xk*x0)
          ec=665.*ebeam**2*b0x

          mode=-1
          fdeni=0.0d0
          freq=freqlow-dfreq
          do ifreq=1,nfreq

            freq=freq+dfreq

            call photons(b0,freq,fluxden)
            call UTIL_SPLINE_INTER(efilt,filt,filt2p,nfilt,freq,facfilt,mode)

            if (ifreq.eq.1) then
              f2=fluxden*facfilt
            else
              f1=f2
              f2=fluxden*facfilt
              fdeni=fdeni+(f2+f1)/2.0d0*dfreq
            endif

          enddo !ifreq

          fluxdeni00=fdeni

          mode=-1
          fdeni=0.0d0
          freq=freqlow-dfreq
          do ifreq=1,nfreq

            freq=freq+dfreq

            call photons(b0x,freq,fluxden)
            call UTIL_SPLINE_INTER(efilt,filt,filt2p,nfilt,freq,facfilt,mode)

            if (ifreq.eq.1) then
              f2=fluxden*facfilt
            else
              f1=f2
              f2=fluxden*facfilt
              fdeni=fdeni+(f2+f1)/2.0d0*dfreq
            endif

          enddo !ifreq

          fluxdeni0=fdeni

          do ixp0=1,nxp0

            xp0=xp01+(ixp0-1)*dxp0

            !Näherungen: L=L+pz, sin(zkz)=kz*z, x(z)=x0

            z1=((px1-x0)/pz-xp0)*brho/b0x
            z2=((px2-x0)/pz-xp0)*brho/b0x

            if (iz.eq.1) then
              z=z1
              pxtrue=px1
              pxfinal=px1
            else
              z=z2
              pxtrue=px2
              pxfinal=px2
            endif

            icount=0
            zeps=z+epspx
            call traj(zk,x0,xp0,pxeps,zeps,x,xp,b)

1           call traj(zk,x0,xp0,pxtrue,z,x,xp,b)

            if (abs(pxfinal-pxtrue).gt.epspx) then
              icount=icount+1
              call util_adjust(pxfinal,zeps,pxeps,z,pxtrue)
              if (icount.le.1000) then
                goto 1
              else
                print*,'*** Warning: icout.gt.1000'
                write(luno,*),'*** Warning: icout.gt.1000'
              endif
            endif

            call caper(icirc,aperl,aperr,aperphi0,pz,pxtrue,corraper)

            if (iz.eq.1) then
              xp1=xp
              x1=x
              z1=z
              b1=b
              px1true=pxtrue
              corraper1=corraper
            else
              xp2=xp
              x2=x
              z2=z
              b2=b
              px2true=pxtrue
              corraper2=corraper
            endif

            fdeni=0.0d0

            mode=-1
            freq=freqlow-dfreq
            do ifreq=1,nfreq

              freq=freq+dfreq

              call photons(b,freq,fluxden)
              call UTIL_SPLINE_INTER(efilt,filt,filt2p,nfilt,freq,facfilt,mode)

              if (ifreq.eq.1) then
                f2=fluxden*facfilt
              else
                f1=f2
                f2=fluxden*facfilt
                fdeni=fdeni+(f2+f1)/2.0d0*dfreq
              endif

            enddo !ifreq

            if (iz.eq.1) then
              fluxdeni1=fdeni*corraper1
            else
              fluxdeni2=fdeni*corraper2
              print*
              print*,'x0,xp0:',x0,xp0
              print*
              print*,'z1, z2:'
              print*,z1, z2
              print*
              print*,'px1true,px2true:'
              print*,px1true,px2true
              print*
              print*,'corraper1,corraper2:'
              print*,corraper1,corraper2
              print*
              print*,'(b0x-b0)/b0:',
     &          (b0x-b0)/b0
              print*
              print*,'db1/b0:',
     &          (b1-b0)/b0
              print*,'db2/b0:',
     &          (b2-b0)/b0
              print*,'db/b0:',
     &          (b2-b1)/b0
              print*
              print*,'db1/b0x:',
     &          (b1-b0x)/b0x
              print*,'db2/b0x:',
     &          (b2-b0x)/b0x
              print*,'db/b0x:',
     &          (b2-b1)/b0x
              print*
              print*,'(f0-f00)/f00:',
     &          (fluxdeni0-fluxdeni00)/fluxdeni00
              print*
              print*,'df1/f00:',
     &          (fluxdeni1-fluxdeni00)/fluxdeni00
              print*,'df2/f00:',
     &          (fluxdeni2-fluxdeni00)/fluxdeni00
              print*
              print*,'df1/f0:',
     &          (fluxdeni1-fluxdeni0)/fluxdeni0
              print*,'df2/f0:',
     &          (fluxdeni2-fluxdeni0)/fluxdeni0
              write(90,'(7e15.5)')px2true,z2,x2,b2,xp2,corraper2,
     &          (fluxdeni2-fluxdeni0)/fluxdeni0
              print*,'df/f:',
     &          (fluxdeni2-fluxdeni1)/(fluxdeni2+fluxdeni1)*2.0d0
              print*
              if (px2true.ne.px1true) then
                print*,'slope/mm:',
     &            (fluxdeni2-fluxdeni1)/(fluxdeni2+fluxdeni1)*2.0d0/
     &            (px2true-px1true)/1000.
              else
                print*,'slope/mm: 0.0'
              endif
              print*
              write(luno,*)
              write(luno,*)'x0,xp0:',x0,xp0
              write(luno,*)
              write(luno,*),'z1, z2:'
              write(luno,*),z1, z2
              write(luno,*)
              write(luno,*)'px1true,px2true:'
              write(luno,*)px1true,px2true
              write(luno,*)
              write(luno,*)'(b0x-b0)/b0:',
     &          (b0x-b0)/b0
              write(luno,*)
              write(luno,*)'db1/b0:',
     &          (b1-b0)/b0
              write(luno,*)'db2/b0:',
     &          (b2-b0)/b0
              write(luno,*)'db/b0:',
     &          (b2-b1)/b0
              write(luno,*)
              write(luno,*)'db1/b0x:',
     &          (b1-b0x)/b0x
              write(luno,*)'db2/b0x:',
     &          (b2-b0x)/b0x
              write(luno,*)'db/b0x:',
     &          (b2-b1)/b0x
              write(luno,*)
              write(luno,*)'(f0-f00)/f00:',
     &          (fluxdeni0-fluxdeni00)/fluxdeni00
              write(luno,*)
              write(luno,*)'df1/f00:',
     &          (fluxdeni1-fluxdeni00)/fluxdeni00
              write(luno,*)'df2/f0:',
     &          (fluxdeni2-fluxdeni00)/fluxdeni00
              write(luno,*)
              write(luno,*)'df1/f0:',
     &          (fluxdeni1-fluxdeni0)/fluxdeni0
              write(luno,*)'df2/f0:',
     &          (fluxdeni2-fluxdeni0)/fluxdeni0
              write(luno,*)'df/f:',
     &          (fluxdeni2-fluxdeni1)/(fluxdeni2+fluxdeni1)*2.0d0
              write(luno,*)
              if (px2true.ne.px1true) then
                write(luno,*)'slope/mm:',
     &            (fluxdeni2-fluxdeni1)/(fluxdeni2+fluxdeni1)*2.0d0/
     &            (px2true-px1true)/1000.
              else
                write(luno,*)'slope/mm: 0.0'
              endif
            endif

          enddo !ixp0
        enddo !ix0

      enddo !z1,z2

      close(luno)

      stop
      end

      subroutine traj(zk,x0,xp0,px,z,x,xp,b)

C Berechnet ausgehend von z und px, die tatsächlichen x,xp,b unter der
c Näherung, dass B transversal unabhängig von z x also
c      B(x)=B(x0)=B(0,0)cos(kx*x) ist.

      use wlsf90m

      implicit none

      include '/home/scheer/util/util_phycon_incl.f'

      double precision zk,x0,xp0,px,z,x,xp,b
      double precision zkz,szkz,czkz,br

      zkz=zk*z

      szkz=sin(zkz)
      czkz=cos(zkz)

      br=b0x/brho

      b=b0x*czkz

      xp=br/zk*szkz+xp0
      x=br/zk/zk*(1.0d0-czkz)+xp0*z+x0

      px=x+(pz-z)*xp

      return
      end

      subroutine photons(b,freq,f)

      use wlsf90m

      implicit none

      double precision b,freq,f,y

      if (b.ne.0.0d0) then
        y=freq/(665.*ebeam**2*abs(b))
        f=y*exp(-y)
      endif
calt      f=1.0d0+freq/ec*(b-b0x)/b0x

      return
      end

      subroutine caper(icirc,al,ar,alpha0,pz,px,corr)

      implicit none

      double precision pi
      parameter (PI=3.14159265359d0)

      double precision al,ar,pz,px,corr,
     &  alpha,ys,xs,f,f0,phi,alpha0

      integer icirc

      alpha=abs(atan(px/pz)+alpha0)

      if (icirc.ne.0) then

C Berechnet effektive Öffnung durch einen Kollimator
C Ansatz: Schnittfläche zweier verschobener Kreise

        xs=alpha*al/2.0d0
        ys=sqrt(ar**2-xs**2)
        phi=2.0d0*acos(xs/ar)

        f0=pi*ar**2
        f=2.0d0*(f0*phi/(2.0d0*pi)-ys*xs)

        corr=f/f0

      else  !(icirc.ne.0) then

c Nährung: Der cos wird vernachlässigt, d.h. z.B. bei al=0 hängt die Apertur nicht
c         vom Winkel ab, was falsch ist, da sie mit cos(alpha) kleiner wird.

        corr=max(1.0d0-tan(alpha)*al/(2.0d0*ar),0.0d0)

      endif  !(icirc.ne.0) then

      return
      end

      include '/home/scheer/util/util_adjust.f'
