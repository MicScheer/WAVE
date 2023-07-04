*CMZ :  3.06/00 15/02/2019  14.44.39  by  Michael Scheer
*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  3.02/03 10/11/2014  10.47.09  by  Michael Scheer
*CMZ :  3.02/00 09/10/2014  14.54.58  by  Michael Scheer
*CMZ :  3.01/05 12/06/2014  09.31.15  by  Michael Scheer
*CMZ :  3.01/02 28/01/2014  17.00.03  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  08.31.37  by  Michael Scheer
*CMZ :  2.63/03 18/02/2009  13.15.35  by  Michael Scheer
*CMZ :  2.54/01 01/03/2005  09.47.58  by  Michael Scheer
*CMZ :  2.53/05 25/02/2005  11.55.05  by  Michael Scheer
*CMZ :  2.52/05 17/08/2004  08.54.30  by  Michael Scheer
*CMZ :  1.01/01 11/08/2004  13.30.53  by  Michael Scheer
*CMZ :  1.01/00 02/03/2004  17.00.13  by  Michael Scheer
*CMZ :  1.00/01 27/02/2004  14.29.35  by  Michael Scheer
*CMZ :  1.00/00 26/02/2004  17.21.29  by  Michael Scheer
*CMZ :  0.99/13 26/02/2004  16.14.57  by  Michael Scheer
*CMZ :  0.99/12 26/02/2004  12.02.34  by  Michael Scheer
*CMZ :  0.99/11 25/02/2004  15.21.06  by  Michael Scheer
*CMZ :  0.99/10 25/02/2004  13.42.35  by  Michael Scheer
*CMZ :  0.99/09 20/02/2004  17.26.48  by  Michael Scheer
*CMZ :  0.99/08 20/02/2004  16.32.55  by  Michael Scheer
*CMZ :  0.99/07 16/02/2004  15.21.29  by  Michael Scheer
*CMZ :  0.99/03 12/02/2004  13.55.05  by  Michael Scheer
*CMZ :  0.99/00 26/01/2004  17.03.49  by  Michael Scheer
*CMZ :  0.00/08 23/01/2004  12.52.23  by  Michael Scheer
*CMZ :  0.00/07 16/01/2004  11.05.44  by  Michael Scheer
*CMZ :  0.00/06 09/01/2004  15.55.17  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  14.52.54  by  Michael Scheer
*CMZ :  0.00/04 23/12/2003  10.15.07  by  Michael Scheer
*CMZ :  0.00/02 15/12/2003  12.43.34  by  Michael Scheer
*CMZ :  0.00/01 10/12/2003  17.56.52  by  Michael Scheer
*-- Author :    Michael Scheer   02/12/2003
      subroutine bpolyplot(iplot,xmin,xmax,ymin,ymax,zmin,zmax,theta,phi,
     &  usercom)
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      implicit none

*KEEP,mshplt.
      include 'mshplt.cmn'
*KEND.

      real, dimension (:), allocatable :: xpl,ypl,zpl,zplm,xmpl,ympl,zmpl

      real
     &  xplb(2),yplb(2),zplb(2),
     &  xplbo(2),yplbo(2),zplbo(2),eps

      real xmin,xmax,ymin,ymax,zmin,zmax,theta,phi,
     &  x,y,z,bx,by,bz,dx,dy,dz,bxo,byo,bzo,bo,
     &  xmn,xmx,ymn,ymx,zmn,zmx,
     &  xmmn,xmmx,ymmn,ymmx,zmmn,zmmx,
     &  xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,
     &  rmag,rcol,rcolo,rplan,rcorn,
     &  xc,yc,zc,xmc(1),ymc(1),zmc(1),dot0,circ0,pscal,vn,vnx,vny,vnz,
     &  rmtyp20,rmtyp24,rmtyp31

      integer i,iplot,iplot1,iplot10,iplot100,ibatch,idev,
     &  imag,icol,iplan,icorn,
     &  iplano,ncorno,iline,nline,
     &  ncorn,nplanmax,ncornmax,idx,igird,imago,impl,izero
     &  ,lunbase

      character(20) cdx
      character(23) cdxmm
      character(256) usercom

*KEEP,pawcmn.
*KEND.

      data dot0/25./
      data circ0/5./
      data rmtyp20/20./
      data rmtyp24/-9999./
      data rmtyp31/31./

      data eps/0.01/

      nplanmax=0
      ncornmax=0

      iplano=0
      nline=0

      xmn=1.e10
      xmx=-1.e10
      ymn=1.e10
      ymx=-1.e10
      zmn=1.e10
      zmx=-1.e10

      open(unit=99,file='polymag.mag',status='old')

1     read(99,*,end=9) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
      if (bx**2+by**2+bz**2.eq.0.0d0) goto 1

      nline=nline+1

      imag=rmag
      icol=rcol
      iplan=rplan
      icorn=rcorn

      if (iplan.gt.nplanmax) nplanmax=iplan
      if (abs(icorn).gt.ncornmax) ncornmax=abs(icorn)

      if (x.lt.xmn) xmn=x
      if (x.gt.xmx) xmx=x
      if (y.lt.ymn) ymn=y
      if (y.gt.ymx) ymx=y
      if (z.lt.zmn) zmn=z
      if (z.gt.zmx) zmx=z

      goto 1

9     close (99)

      allocate(xpl(ncornmax))
      allocate(ypl(ncornmax))
      allocate(zpl(ncornmax))
      allocate(zplm(ncornmax))

      allocate(xmpl(ncornmax*nplanmax))
      allocate(ympl(ncornmax*nplanmax))
      allocate(zmpl(ncornmax*nplanmax))

c initialize plotting

      call util_test_batch(ibatch)

      if (ibatch.eq.0) then
        idev=1   !x11
      else
        idev=0
      endif

c      call mplint(idev)
      lunbase=300000
      call util_get_free_lun(lunbase)
      print*,' '
      print*,
     &  '--- BPOLYPLOT: File polymag*.eps will be generated according plotting option in polymag.in'
      print*,' '
      call mshplt_init(lunbase,-20.,-20.,0,0,800,800,
     &  'polymag.eps',
     &  '',
     &  '',
c     &  'polymag_viewer.sh',
c     &  'polymag_kill_viewer.sh',
     &  0.0)
      call mplopt('DATE',1)
      call mplset('YGTI',1.0)
      call mplset('GSIZ',0.5)

      iplot100=abs(iplot)/100
      iplot10=(abs(iplot)-iplot100*100)/10
      iplot1=abs(iplot)-iplot100*100-iplot10*10

c--- Open plotfiles {

      open(unit=99,file='polymag.mag',status='old')

c--- Open plotfiles }

      iplano=1

      if (xmin.eq.9999.) then
        xplmin=xmn
      else
        xplmin=xmin
      endif

      if (xmax.eq.9999.) then
        xplmax=xmx
      else
        xplmax=xmax
      endif

      if (ymin.eq.9999.) then
        yplmin=ymn
      else
        yplmin=ymin
      endif

      if (ymax.eq.9999.) then
        yplmax=ymx
      else
        yplmax=ymax
      endif

      if (zmin.eq.9999.) then
        zplmin=zmn
      else
        zplmin=zmin
      endif

      if (zmax.eq.9999.) then
        zplmax=zmx
      else
        zplmax=zmax
      endif

      dx=(xplmax-xplmin)*0.05
      if (xmin.eq.9999.) then
        xplmin=xplmin-dx
      endif

      if (xmax.eq.9999.) then
        xplmax=xplmax+dx
      endif

      dy=(yplmax-yplmin)*0.05
      if (ymin.eq.9999.) then
        yplmin=yplmin-dy
      endif

      if (ymax.eq.9999.) then
        yplmax=yplmax+dy
      endif

      dz=(zplmax-zplmin)*0.05
      if (zmin.eq.9999.) then
        zplmin=zplmin-dz
      endif

      if (zmax.eq.9999.) then
        zplmax=zplmax+dz
      endif

      if (xplmax.le.xplmin.or.zplmax.le.zplmin.or.zplmax.le.zplmin) then
        print *,'*** Warning in BPOLYPLOT: Bad coordinate system for plotting '
        goto 9999
      endif

      pscal=sqrt(10000./((xplmax-xplmin)*(zplmax-zplmin)))

c--- Only 3D plot {

      if (abs(iplot100).eq.1) then

c y is vertical (WAVE-system)

        call mplfr3(xplmin,xplmax,-zplmax,-zplmin,yplmin,yplmax,theta,phi,'W')

        do iline=1,nline

11        read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 11

          imag=rmag
          icol=rcol
          iplan=rplan
          icorn=rcorn

          ncorn=abs(icorn)

          if (iline.eq.nline) then

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol

          endif !(iline.eq.nline

          if(iplan.ne.iplano.or.iline.eq.nline) then

            do i=1,ncorno
              if (
     &          xpl(i).lt.xplmin.or.
     &          xpl(i).gt.xplmax.or.
     &          ypl(i).lt.yplmin.or.
     &          ypl(i).gt.yplmax.or.
     &          zpl(i).lt.zplmin.or.
     &          zpl(i).gt.zplmax
     &          ) goto 7
            enddo

            do i=1,ncorno
              zplm(i)=-zpl(i)
            enddo

            call mgset('PLCI',1.)
            call mpl3(ncorno,xpl,zplm,ypl)

            xc=0.
            yc=0.
            zc=0.

            do i=1,ncorno-1
              xc=xc+xpl(i)
              yc=yc+ypl(i)
              zc=zc+zpl(i)
            enddo

            xc=xc/(ncorno-1)
            yc=yc/(ncorno-1)
            zc=zc/(ncorno-1)

            do i=1,ncorno
              xpl(i)=xpl(i)+(xc-xpl(i))*0.015*rcolo
              ypl(i)=ypl(i)+(yc-ypl(i))*0.015*rcolo
              zpl(i)=zpl(i)+(zc-zpl(i))*0.015*rcolo
              zplm(i)=-zpl(i)
            enddo

            call mgset('PLCI',rcolo)
            call mpl3(ncorno,xpl,zplm,ypl)
 7        continue

        endif !iplano

        xpl(ncorn)=x
        ypl(ncorn)=y
        zpl(ncorn)=z

        iplano=iplan
        ncorno=ncorn
        rcolo=rcol

      enddo !nline

      dx=xplmax-xplmin
      idx=nint(log10(dx)-1)
      dx=10.**idx

      xpl(1)=xplmin
      xpl(2)=xpl(1)+dx
      ypl(1)=yplmin-(yplmax-yplmin)*0.03
      ypl(2)=ypl(1)
      zpl(1)=-zplmax-(-zplmin+zplmax)*0.03
      zpl(2)=zpl(1)

      call mgset('PLCI',1.)
      call mpl3(2,xpl,zpl,ypl)

      xpl(2)=xpl(1)
      zpl(2)=zpl(1)+(-zplmin+zplmax)*0.03

      call mgset('PLCI',1.)
      call mpl3(2,xpl,zpl,ypl)

      xpl(1)=xpl(1)+dx
      xpl(2)=xpl(1)
      zpl(2)=zpl(1)+(-zplmin+zplmax)*0.03

      call mgset('PLCI',1.)
      call mpl3(2,xpl,zpl,ypl)
      call miztoc(nint(dx),cdx)

      dx=(xplmax-xplmin)/10.
      dy=(yplmax-yplmin)/10.
      dz=(zplmin-zplmax)/10.

      xpl(1)=xplmin-dx*7.
      xpl(2)=xpl(1)+dx*3.
      ypl(1)=yplmin+dy*3.
      ypl(2)=ypl(1)
      zpl(1)=-zplmax-dz*0.03
      zpl(2)=zpl(1)

      call mpl3(2,xpl,zpl,ypl)

      xpl(1)=xplmin-dx*7.
      xpl(2)=xpl(1)
      ypl(1)=yplmin+dy*3.
      ypl(2)=ypl(1)+dy*3.
      zpl(1)=-zplmax-dz*0.03
      zpl(2)=zpl(1)

      call mpl3(2,xpl,zpl,ypl)

      xpl(1)=xplmin-dx*7.
      xpl(2)=xpl(1)
      ypl(1)=yplmin+dy*3.
      ypl(2)=ypl(1)
      zpl(1)=-zplmax+dz*0.
      zpl(2)=zpl(1)+dz*3.

      call mpl3(2,xpl,zpl,ypl)

      call mplfra(0.,10.,0.,10.,'ABS')
      do i=1,idx+1
        cdxmm(i:i)=cdx(i:i)
      enddo
      cdxmm(idx+2:idx+4)=' mm'
      do i=idx+5,23
        cdxmm(i:i)=' '
      enddo

      call mtx(2.,9.5,usercom)
      call mtx(2.0,1.05,'x')
      call mtx(0.05,2.5,'y')
      call mtx(0.95,0.05,'z')
      call mtx(4.5,0.05,'scale ')
      call mtx(6.,0.05,cdxmm)

      call muwk(0,0)

      if (idev.ne.0) then
      endif

      else  if (abs(iplot100).eq.2) then

c z is vertical

        call mplfr3(xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,theta,phi,'W')

        do iline=1,nline

21        read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 21

          imag=rmag
          icol=rcol
          iplan=rplan
          icorn=rcorn

          ncorn=abs(icorn)

          if (iline.eq.nline) then

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol

          endif !(iline.eq.nline

          if(iplan.ne.iplano.or.iline.eq.nline) then

            do i=1,ncorno
              if (
     &            xpl(i).lt.xplmin.or.
     &            xpl(i).gt.xplmax.or.
     &            ypl(i).lt.yplmin.or.
     &            ypl(i).gt.yplmax.or.
     &            zpl(i).lt.zplmin.or.
     &            zpl(i).gt.zplmax
     &          ) goto 17
            enddo

            call mgset('PLCI',1.)
            call mpl3(ncorno,xpl,ypl,zpl)

            xc=0.
            yc=0.
            zc=0.

            do i=1,ncorno-1
              xc=xc+xpl(i)
              yc=yc+ypl(i)
              zc=zc+zpl(i)
            enddo

            xc=xc/(ncorno-1)
            yc=yc/(ncorno-1)
            zc=zc/(ncorno-1)

            do i=1,ncorno
              xpl(i)=xpl(i)+(xc-xpl(i))*0.015*rcolo
              ypl(i)=ypl(i)+(yc-ypl(i))*0.015*rcolo
              zpl(i)=zpl(i)+(zc-zpl(i))*0.015*rcolo
              zplm(i)=-zpl(i)
            enddo

            call mgset('PLCI',rcolo)
            call mpl3(ncorno,xpl,ypl,zpl)
 17         continue

          endif !iplano

          xpl(ncorn)=x
          ypl(ncorn)=y
          zpl(ncorn)=z

          iplano=iplan
          ncorno=ncorn
          rcolo=rcol

        enddo !nline

        dx=xplmax-xplmin
        idx=nint(log10(dx)-1)
        dx=10.**idx

        xpl(1)=xplmin
        xpl(2)=xpl(1)+dx
        ypl(1)=yplmin-(yplmax-yplmin)*0.03
        ypl(2)=ypl(1)
        zpl(1)=zplmin-(zplmax-zplmin)*0.03
        zpl(2)=zpl(1)
        zplm(2)=-zpl(1)
        zplm(1)=-zpl(2)

        call mgset('PLCI',1.)
        call mpl3(2,xpl,ypl,zpl)

        xpl(2)=xpl(1)
        ypl(2)=ypl(1)+(yplmax-yplmin)*0.03

        call mgset('PLCI',1.)
        call mpl3(2,xpl,ypl,zpl)

        xpl(1)=xpl(1)+dx
        xpl(2)=xpl(1)
        ypl(2)=ypl(1)+(yplmax-yplmin)*0.03

        call mgset('PLCI',1.)
        call mpl3(2,xpl,ypl,zpl)

        dx=(xplmax-xplmin)/10.
        dy=(yplmax-yplmin)/10.
        dz=(zplmax-zplmin)/10.

        xpl(1)=xplmin-dx*5.5
        xpl(2)=xpl(1)+dx*3.
        ypl(1)=yplmin+dy*1.
        ypl(2)=ypl(1)
        zpl(1)=zplmin+dz*0.03
        zpl(2)=zpl(1)

        call mpl3(2,xpl,zpl,ypl)

        xpl(1)=xplmin-dx*5.5
        xpl(2)=xpl(1)
        ypl(1)=yplmin+dy*1.
        ypl(2)=ypl(1)+dy*3.
        zpl(1)=zplmin+dz*0.03
        zpl(2)=zpl(1)

        call mpl3(2,xpl,zpl,ypl)

        xpl(1)=xplmin-dx*5.5
        xpl(2)=xpl(1)
        ypl(1)=yplmin+dy*1.
        ypl(2)=ypl(1)
        zpl(1)=zplmin+dz*0.
        zpl(2)=zpl(1)+dz*3.

        call mpl3(2,xpl,zpl,ypl)

        dx=xplmax-xplmin
        idx=nint(log10(dx)-1)
        dx=10.**idx

        call miztoc(nint(dx),cdx)

        call mplfra(0.,10.,0.,10.,'ABS')
        do i=1,idx+1
          cdxmm(i:i)=cdx(i:i)
        enddo
        cdxmm(idx+2:idx+4)=' mm'
        do i=idx+5,23
          cdxmm(i:i)=' '
        enddo

        call mtx(2.3,0.2,'x')
        call mtx(0.8,1.5,'z')
        call mtx(0.0,0.6,'y')
        call mtx(2.,9.5,usercom)
        call mtx(4.5,0.05,'scale ')
        call mtx(6.,0.05,cdxmm)

        call muwk(0,0)

        if (idev.ne.0) then
        endif

      endif !iplot100

c--- Only 3D plot }

c--- 3D, top and side views {

      if (iplot10.ne.0) then

        call mtitle(usercom)

        call mplset('YMGL',0.5)
        call mplzon(1,1,1,' ')
        call mplfra(0.,10.,0.,10.,'AB')
        call mgset('CHHE',0.4)
        call mtx(3.5,5.3,'upper magnets')
        call mtx(3.5,0.15,'lower magnets')
        call mgset('CHHE',0.25)
        call mplset('YMGL',2.)

        call mplzon(2,2,1,'S')

        rewind(99)
        iplano=1

c--- 3D {

        if (abs(iplot10).eq.1) then

c y is vertical (WAVE-system)

          call mplfr3(xplmin,xplmax,-zplmax,-zplmin,yplmin,yplmax,theta,phi,'W')

          do iline=1,nline

31          read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
            if (bx**2+by**2+bz**2.eq.0.0d0) goto 31

            imag=rmag
            icol=rcol
            iplan=rplan
            icorn=rcorn

            ncorn=abs(icorn)

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if(iplan.ne.iplano.or.iline.eq.nline) then

              do i=1,ncorno
                if (
     &            xpl(i).lt.xplmin.or.
     &            xpl(i).gt.xplmax.or.
     &            ypl(i).lt.yplmin.or.
     &            ypl(i).gt.yplmax.or.
     &            zpl(i).lt.zplmin.or.
     &            zpl(i).gt.zplmax
     &            ) goto 8
              enddo

              do i=1,ncorno
                zplm(i)=-zpl(i)
              enddo

              call mgset('PLCI',1.)
              call mpl3(ncorno,xpl,zplm,ypl)

              xc=0.
              yc=0.
              zc=0.

              do i=1,ncorno-1
                xc=xc+xpl(i)
                yc=yc+ypl(i)
                zc=zc+zpl(i)
              enddo

              xc=xc/(ncorno-1)
              yc=yc/(ncorno-1)
              zc=zc/(ncorno-1)

              do i=1,ncorno
                xpl(i)=xpl(i)+(xc-xpl(i))*0.015*rcolo
                ypl(i)=ypl(i)+(yc-ypl(i))*0.015*rcolo
                zpl(i)=zpl(i)+(zc-zpl(i))*0.015*rcolo
                zplm(i)=-zpl(i)
              enddo

              call mgset('PLCI',rcolo)
              call mpl3(ncorno,xpl,zplm,ypl)
 8            continue

            endif !iplano

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol

          enddo !nline

          dx=xplmax-xplmin
          idx=nint(log10(dx)-1)
          dx=10.**idx

          xpl(1)=xplmin
          xpl(2)=xpl(1)+dx
          ypl(1)=yplmin-(yplmax-yplmin)*0.03
          ypl(2)=ypl(1)
          zpl(1)=-zplmax-(-zplmin+zplmax)*0.03
          zpl(2)=zpl(1)

          call mgset('PLCI',1.)
          call mpl3(2,xpl,zpl,ypl)

          xpl(2)=xpl(1)
          zpl(2)=zpl(1)+(-zplmin+zplmax)*0.03

          call mgset('PLCI',1.)
          call mpl3(2,xpl,zpl,ypl)

          xpl(1)=xpl(1)+dx
          xpl(2)=xpl(1)
          zpl(2)=zpl(1)+(-zplmin+zplmax)*0.03

          call mgset('PLCI',1.)
          call mpl3(2,xpl,zpl,ypl)
          call miztoc(nint(dx),cdx)

          dx=(xplmax-xplmin)/10.
          dy=(yplmax-yplmin)/10.
          dz=(zplmin-zplmax)/10.

          xpl(1)=xplmin-dx*7.
          xpl(2)=xpl(1)+dx*3.
          ypl(1)=yplmin+dy*3.
          ypl(2)=ypl(1)
          zpl(1)=-zplmax-dz*0.03
          zpl(2)=zpl(1)

          call mpl3(2,xpl,zpl,ypl)

          xpl(1)=xplmin-dx*7.
          xpl(2)=xpl(1)
          ypl(1)=yplmin+dy*3.
          ypl(2)=ypl(1)+dy*3.
          zpl(1)=-zplmax-dz*0.03
          zpl(2)=zpl(1)

          call mpl3(2,xpl,zpl,ypl)

          xpl(1)=xplmin-dx*7.
          xpl(2)=xpl(1)
          ypl(1)=yplmin+dy*3.
          ypl(2)=ypl(1)
          zpl(1)=-zplmax+dz*0.
          zpl(2)=zpl(1)+dz*3.

          call mpl3(2,xpl,zpl,ypl)

          call mplfra(0.,10.,0.,10.,'ABS')
          do i=1,idx+1
            cdxmm(i:i)=cdx(i:i)
          enddo
          cdxmm(idx+2:idx+4)=' mm'
          do i=idx+5,23
            cdxmm(i:i)=' '
          enddo

c          call mtx(2.,9.5,usercom)
          call mtx(2.0,1.05,'x')
          call mtx(-0.25,2.5,'y')
          call mtx(0.95,0.05,'z')
          call mtx(4.5,0.05,'scale ')
          call mtx(6.,0.05,cdxmm)

        else  if (abs(iplot10).eq.2) then

c z is vertical

          call mplfr3(xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,theta,phi,'W')

          do iline=1,nline

41          read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
            if (bx**2+by**2+bz**2.eq.0.0d0) goto 41

            imag=rmag
            icol=rcol
            iplan=rplan
            icorn=rcorn

            ncorn=abs(icorn)

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if(iplan.ne.iplano.or.iline.eq.nline) then

              do i=1,ncorno
                if (
     &            xpl(i).lt.xplmin.or.
     &            xpl(i).gt.xplmax.or.
     &            ypl(i).lt.yplmin.or.
     &            ypl(i).gt.yplmax.or.
     &            zpl(i).lt.zplmin.or.
     &            zpl(i).gt.zplmax
     &            ) goto 18
              enddo

              call mgset('PLCI',1.)
              call mpl3(ncorno,xpl,ypl,zpl)

              xc=0.
              yc=0.
              zc=0.

              do i=1,ncorno-1
                xc=xc+xpl(i)
                yc=yc+ypl(i)
                zc=zc+zpl(i)
              enddo

              xc=xc/(ncorno-1)
              yc=yc/(ncorno-1)
              zc=zc/(ncorno-1)

              do i=1,ncorno
                xpl(i)=xpl(i)+(xc-xpl(i))*0.015*rcolo
                ypl(i)=ypl(i)+(yc-ypl(i))*0.015*rcolo
                zpl(i)=zpl(i)+(zc-zpl(i))*0.015*rcolo
                zplm(i)=-zpl(i)
              enddo

              call mgset('PLCI',rcolo)
              call mpl3(ncorno,xpl,ypl,zpl)
 18           continue

            endif !iplano

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol

          enddo !nline

          dx=xplmax-xplmin
          idx=nint(log10(dx)-1)
          dx=10.**idx

          xpl(1)=xplmin
          xpl(2)=xpl(1)+dx
          ypl(1)=yplmin-(yplmax-yplmin)*0.03
          ypl(2)=ypl(1)
          zpl(1)=zplmin-(zplmax-zplmin)*0.03
          zpl(2)=zpl(1)
          zplm(2)=-zpl(1)
          zplm(1)=-zpl(2)

          call mgset('PLCI',1.)
          call mpl3(2,xpl,ypl,zpl)

          xpl(2)=xpl(1)
          ypl(2)=ypl(1)+(yplmax-yplmin)*0.03

          call mgset('PLCI',1.)
          call mpl3(2,xpl,ypl,zpl)

          xpl(1)=xpl(1)+dx
          xpl(2)=xpl(1)
          ypl(2)=ypl(1)+(yplmax-yplmin)*0.03

          call mgset('PLCI',1.)
          call mpl3(2,xpl,ypl,zpl)

          dx=(xplmax-xplmin)/10.
          dy=(yplmax-yplmin)/10.
          dz=(zplmax-zplmin)/10.

          xpl(1)=xplmin-dx*5.5
          xpl(2)=xpl(1)+dx*3.
          ypl(1)=yplmin+dy*1.
          ypl(2)=ypl(1)
          zpl(1)=zplmin+dz*0.03
          zpl(2)=zpl(1)

          call mpl3(2,xpl,zpl,ypl)

          xpl(1)=xplmin-dx*5.5
          xpl(2)=xpl(1)
          ypl(1)=yplmin+dy*1.
          ypl(2)=ypl(1)+dy*3.
          zpl(1)=zplmin+dz*0.03
          zpl(2)=zpl(1)

          call mpl3(2,xpl,zpl,ypl)

          xpl(1)=xplmin-dx*5.5
          xpl(2)=xpl(1)
          ypl(1)=yplmin+dy*1.
          ypl(2)=ypl(1)
          zpl(1)=zplmin+dz*0.
          zpl(2)=zpl(1)+dz*3.

          call mpl3(2,xpl,zpl,ypl)

          dx=xplmax-xplmin
          idx=nint(log10(dx)-1)
          dx=10.**idx

          call miztoc(nint(dx),cdx)

          call mplfra(0.,10.,0.,10.,'ABS')
          do i=1,idx+1
            cdxmm(i:i)=cdx(i:i)
          enddo
          cdxmm(idx+2:idx+4)=' mm'
          do i=idx+5,23
            cdxmm(i:i)=' '
          enddo

          call mtx(2.3,0.2,'x')
          call mtx(0.8,1.5,'z')
          call mtx(0.0,0.6,'y')
c          call mtx(2.,9.5,usercom)
          call mtx(4.5,0.05,'scale ')
          call mtx(6.,0.05,cdxmm)

        endif !iplot10

c--- 3D }

c--- y vs z or z vs y {

        if (iplot10.eq.1) then
          call mplfra(zplmin,zplmax,yplmin,yplmax,' ')
          call mplax('z (mm)', 'y (mm)')
          call bpolypl2(forzpl,forypl,forcol,23)
        else if (iplot10.eq.2) then
          call mplfra(yplmin,yplmax,zplmin,zplmax,' ')
          call mplax('y (mm)', 'z (mm)')
          call bpolypl2(forypl,forzpl,forcol,23)
        endif

        rewind(99)
        iplano=1

        do iline=1,nline

51        read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 51

          imag=rmag
          icol=rcol
          iplan=rplan
          icorn=rcorn

          ncorn=abs(icorn)

          if (iline.eq.nline) then

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol

          endif !(iline.eq.nline

          if(iplan.ne.iplano.or.iline.eq.nline) then

            do i=1,ncorno
              zplm(i)=-zpl(i)
            enddo

            call mgset('PLCI',1.)

            if (iplot10.eq.1) then
              call mpl(ncorno,zpl,ypl)
            else if (iplot10.eq.2) then
              call mpl(ncorno,ypl,zpl)
            endif

            xc=0.
            yc=0.
            zc=0.

            do i=1,ncorno-1
              xc=xc+xpl(i)
              yc=yc+ypl(i)
              zc=zc+zpl(i)
            enddo

            xc=xc/(ncorno-1)
            yc=yc/(ncorno-1)
            zc=zc/(ncorno-1)

            izero=0
            do i=1,ncorno
              if (
     &          abs(yc-ypl(i)).gt.1.0e-6 .and. abs(zc-zpl(i)).gt.1.0e-6
     &          ) izero=1
              xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
              ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
              zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
              zplm(i)=-zpl(i)
            enddo

            call mgset('PLCI',rcolo)

            if (iplot10.eq.1) then
              if (izero.ne.0) call mpl(ncorno,zpl,ypl)
            else if (iplot10.eq.2) then
              if (izero.ne.0) call mpl(ncorno,ypl,zpl)
            endif

          endif !iplano

          xpl(ncorn)=x
          ypl(ncorn)=y
          zpl(ncorn)=z

          iplano=iplan
          ncorno=ncorn
          rcolo=rcol

        enddo !nline

        if (iplot10.eq.1) then
          call bpolypl2(forzpl,forypl,forcol,23)
        else if (iplot10.eq.2) then
          call bpolypl2(forypl,forzpl,forcol,23)
        endif

c--- y vs z or z vs y }

c--- top views of girder {

        call mplzon(1,4,3,'S')

        do igird=1,2

c--- z vs x, y is vertical coordinate {

          call mplfra(xplmin,xplmax,zplmin,zplmax,' ')
          call mplax('x (mm)', 'z (mm)')
          call bpolypl2(forxpl,forzpl,forcol,13)

          rewind(99)
61        read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 61
          backspace(99)
          imago=rmag
          iplano=1
          impl=0

          do iline=1,nline

71          read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
            if (bx**2+by**2+bz**2.eq.0.0d0) goto 71

            imag=rmag
            icol=rcol
            iplan=rplan
            icorn=rcorn

            ncorn=abs(icorn)

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if (imag.ne.imago.or.iline.eq.nline) then

              xmc(1)=0.
              ymc(1)=0.
              zmc(1)=0.

              xmmx=-1.0e30
              xmmn= 1.0e30
              ymmx=-1.0e30
              ymmn= 1.0e30
              zmmx=-1.0e30
              zmmn= 1.0e30

              do i=1,impl-1
                xmc(1)=xmc(1)+xmpl(i)
                ymc(1)=ymc(1)+ympl(i)
                zmc(1)=zmc(1)+zmpl(i)
                if (xmpl(i).gt.xmmx) xmmx=xmpl(i)
                if (xmpl(i).lt.xmmn) xmmn=xmpl(i)
                if (ympl(i).gt.ymmx) ymmx=ympl(i)
                if (ympl(i).lt.ymmn) ymmn=ympl(i)
                if (zmpl(i).gt.zmmx) zmmx=zmpl(i)
                if (zmpl(i).lt.zmmn) zmmn=zmpl(i)
              enddo

              xmc(1)=xmc(1)/(impl-1)
              ymc(1)=ymc(1)/(impl-1)
              zmc(1)=zmc(1)/(impl-1)

              dx=xmmx-xmmn
              dy=ymmx-ymmn
              dz=zmmx-zmmn

              impl=0

            endif !imag.ne.imago

            impl=impl+1

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              xmpl(impl)=x
              ympl(impl)=y
              zmpl(impl)=z

              bxo=bx
              byo=by
              bzo=bz

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if(iplan.ne.iplano.or.iline.eq.nline) then

              xc=0.
              yc=0.
              zc=0.

              do i=1,ncorno-1
                xc=xc+xpl(i)
                yc=yc+ypl(i)
                zc=zc+zpl(i)
              enddo

              xc=xc/(ncorno-1)
              yc=yc/(ncorno-1)
              zc=zc/(ncorno-1)

              if (igird.eq.1.and.yc.ge.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,zpl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(xc-xpl(i)).gt.1.0e-6 .and. abs(zc-zpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,zpl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,zplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(zplbo(2)-zplbo(1))**2)

                  if (vn.ne.0.0d0) then

                    vnx=(xplbo(2)-xplbo(1))/vn
                    vnz=(zplbo(2)-zplbo(1))/vn

                    xplb(1)=xplbo(2)+vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)-vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)

                    xplb(1)=xplbo(2)-vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)+vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)

                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  endif

                endif !imago

              else if (igird.eq.2.and.yc.le.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,zpl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(xc-xpl(i)).gt.1.0e-6 .and. abs(zc-zpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,zpl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,zplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(zplbo(2)-zplbo(1))**2)

                  if (vn.ne.0.0d0) then

                    vnx=(xplbo(2)-xplbo(1))/vn
                    vnz=(zplbo(2)-zplbo(1))/vn

                    xplb(1)=xplbo(2)+vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)-vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)

                    xplb(1)=xplbo(2)-vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)+vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)

                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  endif

                endif !imago

              endif !yc

            endif !iplano

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            xmpl(impl)=x
            ympl(impl)=y
            zmpl(impl)=z

            bxo=bx
            byo=by
            bzo=bz

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol
            imago=imag

          enddo !nline

          call bpolypl2(forxpl,forzpl,forcol,13)

        enddo !igird

c--- y vs x }

        call muwk(0,0)

        if (idev.ne.0) then
        endif

      else if (iplot10.eq.2) then

        call mplset('YMGL',0.5)
        call mplzon(1,1,1,' ')
        call mplfra(0.,10.,0.,10.,'AB')
        call mgset('CHHE',0.4)
        call mtx(3.5,5.3,'upper magnets')
        call mtx(3.5,0.15,'lower magnets')
        call mgset('CHHE',0.25)
        call mplset('YMGL',2.)

        call mplzon(1,4,3,'S')

        do igird=1,2

c--- z vs x, z is vertical coordinate {

          call mplfra(xplmin,xplmax,yplmin,yplmax,' ')
          call mplax('x (mm)', 'y (mm)')
          call bpolypl2(forxpl,forypl,forcol,12)

          rewind(99)
81        read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 81
          backspace(99)
          imago=rmag
          iplano=1
          impl=0

          do iline=1,nline

91          read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
            if (bx**2+by**2+bz**2.eq.0.0d0) goto 91

            imag=rmag
            icol=rcol
            iplan=rplan
            icorn=rcorn

            ncorn=abs(icorn)

            if (imag.ne.imago.or.iline.eq.nline) then

              xmc(1)=0.
              ymc(1)=0.
              zmc(1)=0.

              xmmx=-1.0e30
              xmmn= 1.0e30
              ymmx=-1.0e30
              ymmn= 1.0e30
              zmmx=-1.0e30
              zmmn= 1.0e30

              do i=1,impl-1
                xmc(1)=xmc(1)+xmpl(i)
                ymc(1)=ymc(1)+ympl(i)
                zmc(1)=zmc(1)+zmpl(i)
                if (xmpl(i).gt.xmmx) xmmx=xmpl(i)
                if (xmpl(i).lt.xmmn) xmmn=xmpl(i)
                if (ympl(i).gt.ymmx) ymmx=ympl(i)
                if (ympl(i).lt.ymmn) ymmn=ympl(i)
                if (zmpl(i).gt.zmmx) zmmx=zmpl(i)
                if (zmpl(i).lt.zmmn) zmmn=zmpl(i)
              enddo

              xmc(1)=xmc(1)/(impl-1)
              ymc(1)=ymc(1)/(impl-1)
              zmc(1)=zmc(1)/(impl-1)

              dx=xmmx-xmmn
              dy=ymmx-ymmn
              dz=zmmx-zmmn

              impl=0

            endif !imag.ne.imago

            impl=impl+1

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              xmpl(impl)=x
              ympl(impl)=y
              zmpl(impl)=z

              bxo=bx
              byo=by
              bzo=bz

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if(iplan.ne.iplano.or.iline.eq.nline) then

              xc=0.
              yc=0.
              zc=0.

              do i=1,ncorno-1
                xc=xc+xpl(i)
                yc=yc+ypl(i)
                zc=zc+zpl(i)
              enddo

              xc=xc/(ncorno-1)
              yc=yc/(ncorno-1)
              zc=zc/(ncorno-1)

              if (igird.eq.1.and.zc.ge.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,ypl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(yc-ypl(i)).gt.1.0e-6 .and. abs(xc-xpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,ypl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,yplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(yplbo(2)-yplbo(1))**2)

                  if (vn.ne.0.0d0) then
                    vnx=(xplbo(2)-xplbo(1))/vn
                    vny=(yplbo(2)-yplbo(1))/vn

                    xplb(1)=xplbo(2)+vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)-vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)

                    xplb(1)=xplbo(2)-vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)+vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)
                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  endif

                endif !imago

              else if (igird.eq.2.and.zc.le.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,ypl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(yc-ypl(i)).gt.1.0e-6 .and. abs(xc-xpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,ypl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,yplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(yplbo(2)-yplbo(1))**2)

                  if (vn.ne.0.0d0) then
                    vnx=(xplbo(2)-xplbo(1))/vn
                    vny=(yplbo(2)-yplbo(1))/vn

                    xplb(1)=xplbo(2)+vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)-vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)

                    xplb(1)=xplbo(2)-vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)+vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)
                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  endif

                endif !imago

              endif !yc

            endif !iplano

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            xmpl(impl)=x
            ympl(impl)=y
            zmpl(impl)=z

            bxo=bx
            byo=by
            bzo=bz

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol
            imago=imag

          enddo !nline

          call bpolypl2(forxpl,forypl,forcol,12)

        enddo !igird

c--- y vs x }

        call muwk(0,0)

        if (idev.ne.0) then
        endif

        call muwk(0,0)

        if (idev.ne.0) then
        endif

      endif !iplot10

c--- top views of girder}

c--- 3D, top and side views }

c--- top views of girder {

      if (iplot1.eq.1) then

        call mtitle(usercom)

        call mplset('YMGL',0.5)
        call mplzon(1,1,1,' ')
        call mplfra(0.,10.,0.,10.,'AB')
        call mgset('CHHE',0.4)
        call mtx(3.5,5.3,'upper magnets')
        call mtx(3.5,0.15,'lower magnets')
        call mgset('CHHE',0.25)
        call mplset('YMGL',2.)

        call mplzon(1,2,1,'S')

        do igird=1,2

c--- z vs x, y is vertical coordinate {

          call mplfra(xplmin,xplmax,zplmin,zplmax,' ')
          call mplax('x (mm)', 'z (mm)')
          call bpolypl2(forxpl,forzpl,forcol,13)

          rewind(99)
101       read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 101
          backspace(99)
          imago=rmag
          iplano=1
          impl=0

          do iline=1,nline

102         read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
            if (bx**2+by**2+bz**2.eq.0.0d0) goto 102

            imag=rmag
            icol=rcol
            iplan=rplan
            icorn=rcorn

            ncorn=abs(icorn)

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if (imag.ne.imago.or.iline.eq.nline) then

              xmc(1)=0.
              ymc(1)=0.
              zmc(1)=0.

              xmmx=-1.0e30
              xmmn= 1.0e30
              ymmx=-1.0e30
              ymmn= 1.0e30
              zmmx=-1.0e30
              zmmn= 1.0e30

              do i=1,impl-1
                xmc(1)=xmc(1)+xmpl(i)
                ymc(1)=ymc(1)+ympl(i)
                zmc(1)=zmc(1)+zmpl(i)
                if (xmpl(i).gt.xmmx) xmmx=xmpl(i)
                if (xmpl(i).lt.xmmn) xmmn=xmpl(i)
                if (ympl(i).gt.ymmx) ymmx=ympl(i)
                if (ympl(i).lt.ymmn) ymmn=ympl(i)
                if (zmpl(i).gt.zmmx) zmmx=zmpl(i)
                if (zmpl(i).lt.zmmn) zmmn=zmpl(i)
              enddo

              xmc(1)=xmc(1)/(impl-1)
              ymc(1)=ymc(1)/(impl-1)
              zmc(1)=zmc(1)/(impl-1)

              dx=xmmx-xmmn
              dy=ymmx-ymmn
              dz=zmmx-zmmn

              impl=0

            endif !imag.ne.imago

            impl=impl+1

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              xmpl(impl)=x
              ympl(impl)=y
              zmpl(impl)=z

              bxo=bx
              byo=by
              bzo=bz

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if(iplan.ne.iplano.or.iline.eq.nline) then

              xc=0.
              yc=0.
              zc=0.

              do i=1,ncorno-1
                xc=xc+xpl(i)
                yc=yc+ypl(i)
                zc=zc+zpl(i)
              enddo

              xc=xc/(ncorno-1)
              yc=yc/(ncorno-1)
              zc=zc/(ncorno-1)

              if (igird.eq.1.and.yc.ge.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,zpl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(xc-xpl(i)).gt.1.0e-6 .and. abs(zc-zpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,zpl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,zplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(zplbo(2)-zplbo(1))**2)

                  if (vn.ne.0.0d0) then
                    vnx=(xplbo(2)-xplbo(1))/vn
                    vnz=(zplbo(2)-zplbo(1))/vn

                    xplb(1)=xplbo(2)+vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)-vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)

                    xplb(1)=xplbo(2)-vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)+vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)
                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  endif

                endif !imago

              else if (igird.eq.2.and.yc.le.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,zpl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(xc-xpl(i)).gt.1.0e-6 .and. abs(zc-zpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,zpl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,zplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(zplbo(2)-zplbo(1))**2)
                  if (vn.ne.0.0d0) then
                    vnx=(xplbo(2)-xplbo(1))/vn
                    vnz=(zplbo(2)-zplbo(1))/vn

                    xplb(1)=xplbo(2)+vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)-vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)

                    xplb(1)=xplbo(2)-vnz*dx/10.0-vnx*dx/10.0
                    zplb(1)=zplbo(2)+vnx*dz/10.0-vnz*dz/10.0

                    call mpl(2,xplb,zplb)
                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),zmc(1))
                  endif

                endif !imago

              endif !yc

            endif !iplano

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            xmpl(impl)=x
            ympl(impl)=y
            zmpl(impl)=z

            bxo=bx
            byo=by
            bzo=bz

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol
            imago=imag

          enddo !nline

          call bpolypl2(forxpl,forzpl,forcol,13)

        enddo !igird

c--- y vs x }

        call muwk(0,0)

        if (idev.ne.0) then
        endif

      else if (iplot1.eq.2) then

        call mtitle(usercom)

        call mplset('YMGL',0.5)
        call mplzon(1,1,1,' ')
        call mplfra(0.,10.,0.,10.,'AB')
        call mgset('CHHE',0.4)
        call mtx(3.5,5.3,'upper magnets')
        call mtx(3.5,0.15,'lower magnets')
        call mgset('CHHE',0.25)
        call mplset('YMGL',2.)

        call mplzon(1,2,1,'S')

        do igird=1,2

c--- z vs x, z is vertical coordinate {

          call mplfra(xplmin,xplmax,yplmin,yplmax,' ')
          call mplax('x (mm)', 'y (mm)')
          call bpolypl2(forxpl,forypl,forcol,12)

          rewind(99)
103       read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
          if (bx**2+by**2+bz**2.eq.0.0d0) goto 103
          backspace(99)
          imago=rmag
          iplano=1
          impl=0

          do iline=1,nline

104         read(99,*) rmag,rcol,rplan,rcorn,x,y,z,bx,by,bz
            if (bx**2+by**2+bz**2.eq.0.0d0) goto 104

            imag=rmag
            icol=rcol
            iplan=rplan
            icorn=rcorn

            ncorn=abs(icorn)

            if (imag.ne.imago.or.iline.eq.nline) then

              xmc(1)=0.
              ymc(1)=0.
              zmc(1)=0.

              xmmx=-1.0e30
              xmmn= 1.0e30
              ymmx=-1.0e30
              ymmn= 1.0e30
              zmmx=-1.0e30
              zmmn= 1.0e30

              do i=1,impl-1
                xmc(1)=xmc(1)+xmpl(i)
                ymc(1)=ymc(1)+ympl(i)
                zmc(1)=zmc(1)+zmpl(i)
                if (xmpl(i).gt.xmmx) xmmx=xmpl(i)
                if (xmpl(i).lt.xmmn) xmmn=xmpl(i)
                if (ympl(i).gt.ymmx) ymmx=ympl(i)
                if (ympl(i).lt.ymmn) ymmn=ympl(i)
                if (zmpl(i).gt.zmmx) zmmx=zmpl(i)
                if (zmpl(i).lt.zmmn) zmmn=zmpl(i)
              enddo

              xmc(1)=xmc(1)/(impl-1)
              ymc(1)=ymc(1)/(impl-1)
              zmc(1)=zmc(1)/(impl-1)

              dx=xmmx-xmmn
              dy=ymmx-ymmn
              dz=zmmx-zmmn

              impl=0

            endif !imag.ne.imago

            impl=impl+1

            if (iline.eq.nline) then

              xpl(ncorn)=x
              ypl(ncorn)=y
              zpl(ncorn)=z

              xmpl(impl)=x
              ympl(impl)=y
              zmpl(impl)=z

              bxo=bx
              byo=by
              bzo=bz

              iplano=iplan
              ncorno=ncorn
              rcolo=rcol

            endif !(iline.eq.nline

            if(iplan.ne.iplano.or.iline.eq.nline) then

              xc=0.
              yc=0.
              zc=0.

              do i=1,ncorno-1
                xc=xc+xpl(i)
                yc=yc+ypl(i)
                zc=zc+zpl(i)
              enddo

              xc=xc/(ncorno-1)
              yc=yc/(ncorno-1)
              zc=zc/(ncorno-1)

              if (igird.eq.1.and.zc.ge.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,ypl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(yc-ypl(i)).gt.1.0e-6 .and. abs(xc-xpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,ypl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,yplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(yplbo(2)-yplbo(1))**2)
                  if (vn.ne.0.0d0) then
                    vnx=(xplbo(2)-xplbo(1))/vn
                    vny=(yplbo(2)-yplbo(1))/vn

                    xplb(1)=xplbo(2)+vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)-vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)

                    xplb(1)=xplbo(2)-vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)+vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)
                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  endif

                endif !imago

              else if (igird.eq.2.and.zc.le.0.0) then

                do i=1,ncorno
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',1.)
                call mpl(ncorno,xpl,ypl)

                izero=0
                do i=1,ncorno
                  if (
     &              abs(yc-ypl(i)).gt.1.0e-6 .and. abs(xc-xpl(i)).gt.1.0e-6
     &              ) izero=1
                  xpl(i)=xpl(i)+(xc-xpl(i))*0.03*rcolo
                  ypl(i)=ypl(i)+(yc-ypl(i))*0.03*rcolo
                  zpl(i)=zpl(i)+(zc-zpl(i))*0.03*rcolo
                  zplm(i)=-zpl(i)
                enddo

                call mgset('PLCI',rcolo)
                if (izero.ne.0) call mpl(ncorno,xpl,ypl)

                if (imag.ne.imago.or.iline.eq.nline) then

                  bo=sqrt(bxo*bxo+byo*byo+bzo*bzo)

                  if (abs(bxo).lt.bo*eps) bxo=0.0
                  if (abs(byo).lt.bo*eps) byo=0.0
                  if (abs(bzo).lt.bo*eps) bzo=0.0

                  xplb(1)=xmc(1)-2.*bxo/bo*dx/8.
                  xplb(2)=xmc(1)+2.*bxo/bo*dx/8.
                  yplb(1)=ymc(1)-2.*byo/bo*dy/8.
                  yplb(2)=ymc(1)+2.*byo/bo*dy/8.
                  zplb(1)=zmc(1)-2.*bzo/bo*dz/8.
                  zplb(2)=zmc(1)+2.*bzo/bo*dz/8.

                  xplbo(1)=xplb(1)
                  xplbo(2)=xplb(2)
                  yplbo(1)=yplb(1)
                  yplbo(2)=yplb(2)
                  zplbo(1)=zplb(1)
                  zplbo(2)=zplb(2)

                  call mpl(2,xplbo,yplbo)

                  vn=sqrt((xplbo(2)-xplbo(1))**2+(yplbo(2)-yplbo(1))**2)
                  if (vn.ne.0.0d0) then
                    vnx=(xplbo(2)-xplbo(1))/vn
                    vny=(yplbo(2)-yplbo(1))/vn

                    xplb(1)=xplbo(2)+vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)-vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)

                    xplb(1)=xplbo(2)-vny*dx/10.0-vnx*dx/10.0
                    yplb(1)=yplbo(2)+vnx*dy/10.0-vny*dy/10.0

                    call mpl(2,xplb,yplb)
                  endif !vn

                  if (byo.gt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp20)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  else if (byo.lt.0.0) then
                    call mgset('MTYP',rmtyp24)
                    call mgset('MSCF',circ0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                    call mgset('MTYP',rmtyp31)
                    call mgset('MSCF',dot0*pscal/5.)
                    call mpm(1,xmc(1),ymc(1))
                  endif

                endif !imago

              endif !yc

            endif !iplano

            xpl(ncorn)=x
            ypl(ncorn)=y
            zpl(ncorn)=z

            xmpl(impl)=x
            ympl(impl)=y
            zmpl(impl)=z

            bxo=bx
            byo=by
            bzo=bz

            iplano=iplan
            ncorno=ncorn
            rcolo=rcol
            imago=imag

          enddo !nline

          call bpolypl2(forxpl,forypl,forcol,12)

        enddo !igird

c--- y vs x }

        call muwk(0,0)

        if (idev.ne.0) then
        endif

      endif !iplot1

c--- top views of girder}

9999  close (99)

      deallocate(xpl)
      deallocate(ypl)
      deallocate(zpl)
      deallocate(zplm)
      deallocate(xmpl)
      deallocate(ympl)
      deallocate(zmpl)

      call mplend

      return
      end
