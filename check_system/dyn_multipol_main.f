*CMZ :  2.62/01 17/04/2007  13.58.17  by  Michael Scheer
*CMZ :  2.61/04 29/03/2007  15.30.35  by  Michael Scheer
*CMZ :  2.61/02 20/03/2007  16.55.13  by  Michael Scheer
*CMZ :  2.57/05 08/01/2007  18.02.51  by  Michael Scheer
*-- Author :    Michael Scheer   13/12/2006
      program dyn_multipol_main

c---  Calcutes dynamic multipoles and optionally applies shiming results.
c     Effects of L-shims are substracted from dynamic multipoles.
c
c usage:

c        ue112_dyn_multipole_s3_sh_0.dat and
c        ue112_dyn_multipole_s3_sh_56.dat contain distributions without L-shims
c        for shift=0 and shift=56


c        uname=0: just calculation of dynamic multipoles
c        uname=1: apply fit results
c        fit_dyn_multipol.out contains therot-values to be taken
c        fit_dyn_multipol.in must contain correct number of slices

      implicit none

*KEEP,dyn_multipol.
      integer nzp,nzf0,nzfpi,nzfpih,nz,iplot,mode,itherot,nslice,
     &  iuheli,irecu,imodutype,ITHEMSYM,ITHESYML,ITHESYMD,
     &  iustep,moduname,modustep,nperi

      parameter (nzp=1000)

      integer nparp
      parameter (nparp=100)
      integer nvarp
      parameter (nvarp=100)

      double precision z0(nzp),
     &  f0(nzp),
     &  zpi(nzp),
     &  fpi(nzp),
     &  zpih(nzp),
     &  fpih(nzp),
     &  y2f0(nzp),
     &  y2fpi(nzp),
     &  y2fpih(nzp),
     &  ws1(nzp),
     &  ws2(nzp),
     &  ws3(nzp),
     &  ws4(nzp),
     &  phi,z,f,shift,radin,bcr,
     &  paropt(nparp),var(3,nvarp)
     &  ,zmin,zmax,dz
     &  ,xstart,xstop,xianf,xiend,urshift,period,xiustep,xfustep
     &  ,thegrotu,thegrotl,scalmod

      real xplmin,xplmax,yplmin,yplmax

      character(50) comment

      common/dynmulc/
     &  moduname,modustep,
     &  z0,f0,zpi,fpi,zpih,fpih,phi,z,f,shift,radin,bcr,
     &  y2f0,
     &  y2fpi,
     &  y2fpih,ws1,ws2,ws3,ws4,paropt,var,
     &  xplmin,xplmax,yplmin,yplmax,zmin,zmax,dz,
     &  nzf0,nzfpi,nzfpih,nz,iplot,mode,itherot,nslice,
     &  iuheli,irecu,imodutype,ITHEMSYM,ITHESYML,ITHESYMD,
     &  iustep,nperi
     &  ,xstart,xstop,xianf,xiend,urshift,period,xiustep,xfustep,
     &  thegrotu,thegrotl,scalmod,
     &  comment

      namelist/dynmuln/
     &  moduname,modustep,xstart,xianf,xiend,xstop,mode,
     &  iplot,xplmin,xplmax,dz,yplmin,yplmax,
     &  irecu,iuheli,urshift,period,imodutype,ITHEMSYM,ITHESYML,ITHESYMD,
     &  itherot,nperi,
     &  nslice,radin,bcr,xiustep,xfustep,
     &  thegrotu,thegrotl,scalmod,
     &  comment
*KEND.

      double precision thestart(nparp),dthe,chi2min,
     &  zi,zf0,zpf0,yf0,ypf0,chi2,
     &  dyn,twopi,ff0,ffpi,cs22,sn22,sn2,brho,
     &  vdyn(4,1000),xf0,
     &  zianf,yianf,zpianf,ypianf,
     &  ziend,yiend,zpiend,ypiend,zpout,ypout,
     &  xx(3),yy(3),yp(3),xopt,aopt,a(3)

      real xpl(1000),ypl(1000),dx,dy

      integer ivar,nvar,ipar,istat1,istat2,modunameo,nmagmod

      REAL*4 RPAW(50000)
      COMMON/PAWC/RPAW

      character c1

      data chi2min/1.0d30/
      data twopi/6.283185307180d0/

      open(unit=99,file='fit_dyn_multipol.in',status='old')

      read(99,dynmuln)

      if (xianf.eq.9999.)   xianf=xstart
      if (xiend.eq.9999.)   xiend=xstop
      if (xiustep.eq.9999.) xiustep=-(nperi+2)*period/2.0d0/1000.0d0
      if (xfustep.eq.9999.) xfustep=+(nperi+2)*period/2.0d0/1000.0d0

      if (thegrotu.eq.9999.)
     &  thegrotu=acos(sin(urshift/period*twopi)**2)/twopi*360.
      if (thegrotl.eq.9999.)
     &  thegrotl=acos(sin(urshift/period*twopi)**2)/twopi*360.

      modunameo=moduname

      if (moduname.eq.1) then

        iustep=0

      else if (moduname.eq.2) then

        iustep=1
        moduname=1

      else if (moduname.eq.9999) then
c noch nicht aktiviert
        call util_skip_comment(99)
        read(99,*)dthe

        call util_skip_comment(99)
        read(99,*)thestart(1)

c L-Shims initialization{

        nzf0=0
        open(unit=99,file='ue112_dyn_multipole_s3_sh_0.dat',status='old')
1       continue

        read(99,*,end=91)z,f

        nzf0=nzf0+1
        z0(nzf0)=z
        f0(nzf0)=f

        goto 1

91      close(99)

        nzfpi=0
        open(unit=99,file='ue112_dyn_multipole_s3_sh_56.dat',status='old')
2       continue

        read(99,*,end=92)z,f

        nzfpi=nzfpi+1
        zpi(nzfpi)=z
        fpi(nzfpi)=f

        goto 2
92      close(99)

        call util_spline_coef(z0,f0,nzf0,0.0d0,0.0d0,y2f0,ws1,ws2,ws3,ws4)
        call util_spline_coef(zpi,fpi,nzfpi,0.0d0,0.0d0,y2fpi,ws1,ws2,ws3,ws4)

c L-Shims init}

      else !moduname

        stop '*** Error in dyn_multipol_main: Bad mode for UNAME/UOUT'

      endif !moduname

      close(99)

      zmin=xplmin
      zmax=xplmax
      xplmin=1.e30
      xplmax=-1.e30
      yplmin=1.e30
      yplmax=-1.e30

c -------------------------------------------------

c actual loop{

      open(unit=99,file='dyn_multipol_pawtit.txt')
      write(99,*)comment
      close(99)

      if (modunameo.eq.2) then
        call system
     &    ("mv ue112_dyn_multipole_lshimmed.zibint ue112_dyn_multipole_lshimmed.zibint.bck")
        call system
     &    ("mv ue112_dyn_multipole_lshimmed.zfbint ue112_dyn_multipole_lshimmed.zfbint.bck")
        open(unit=71,file='ue112_dyn_multipole_lshimmed.zibint',status='new')
        open(unit=72,file='ue112_dyn_multipole_lshimmed.zfbint',status='new')
      endif

      call system("mv dyn_multipole.loop dyn_multipole.loop.bck")
      call system("mv dyn_multipole_vert.loop dyn_multipole_vert.loop.bck")
      open(unit=88,file='dyn_multipole.loop',status='new',recl=256)
      open(unit=89,file='dyn_multipole_vert.loop',status='new',recl=256)

      write(88,*)xianf,xiend
      write(89,*)xianf,xiend

      z=zmin
      nvar=0

      do while(z.le.zmax)

        open(unit=80,file='uname.fit',status='unknown')

        if (moduname.eq.1) then

          write(80,'(a)')comment

          write(80,*)moduname,' !moduname'
          write(80,*)modustep,' !modustep'

          write(80,*)z/1000.0d0,' !z'

          write(80,*)iustep,' !iustep'
          write(80,*)xiustep,xfustep,' !xiustep,xfustep'

          write(80,*)xstart,' !xstart'
          write(80,*)xianf,' !xianf'
          write(80,*)xiend,' !xiend'
          write(80,*)xstop,' !xstop'

          write(80,*)mode,' !mode'

          write(80,*)irecu,' !irecu'
          write(80,*)iuheli,' !iuheli'
          write(80,*)urshift,' !urshift, bzw. shellana*L_z'
          write(80,*)period,' !period-length'
          write(80,*)nperi,' !krecper bzw. nperella'

          write(80,*)imodutype,' !BELLANA oder REC-Undulator'
          write(80,*)itherot,' !Korrekturmagnete'
          write(80,*)scalmod,' !SCALMOD'
          write(80,*)thegrotu,' !THEGROTU'
          write(80,*)thegrotl,' !THEGROTL'
          write(80,*)ithemsym,' !ithemsym'
          write(80,*)ithesyml,' !ithesyml'
          write(80,*)ithesymd,' !ithesymd'
          write(80,*)nslice,' !nslice'
          write(80,*)radin,' !radin'
          write(80,*)bcr,' !bcr'

        else if (moduname.eq.9999) then

          open(unit=99,file='fit_dyn_multipol.out',status='old')
          read(99,*)chi2min
          do ipar=1,nmagmod
            read(99,*)thestart(ipar)
          enddo
          close(99)

          nvar=nmagmod

          write(80,*)mode
          write(80,*)nmagmod
          write(80,*)nslice
          write(80,*)shift
          write(80,*)radin
          write(80,*)bcr
          do ipar=1,nmagmod
            write(80,*)thestart(ipar)
          enddo

        endif !moduname

        close(80) !uname.fit

        call system('cp wave.in.ue112.9999 wave.in')
        call system('/home/scheer/wav/bin/wave.exe')

        if (moduname.eq.1) then

          open(unit=80,file='uout.fit',status='old',err=9123)

          read(80,*)brho

          read(80,*)xianf
          read(80,*)xiend
          read(80,*)zianf,yianf
          read(80,*)zpianf,ypianf
          read(80,*)ziend,yiend
          read(80,*)zpiend,ypiend

          read(80,*)xf0
          read(80,*)zf0,yf0
          read(80,*)zpf0,ypf0

          close(80)

          zianf=zianf*1000.0d0
          yianf=yianf*1000.0d0
          zpianf=zpianf*1000.0d0
          ypianf=ypianf*1000.0d0

          ziend=ziend*1000.0d0
          yiend=yiend*1000.0d0
          zpiend=zpiend*1000.0d0
          ypiend=ypiend*1000.0d0

          zf0=zf0*1000.0d0
          yf0=yf0*1000.0d0
          zpf0=zpf0*1000.0d0
          ypf0=ypf0*1000.0d0

          zpout=zpiend-zpianf
          ypout=ypiend-ypianf

          write(88,*)z,zpout,zpf0
          write(89,*)z,ypout,ypf0

          if (modunameo.eq.2) then
            write(71,*)  z,zpf0*brho,ypf0*brho
            write(72,*)zf0,zpf0*brho,ypf0*brho
          endif

          if (iplot.ne.0) then
            nvar=nvar+1
            xpl(nvar)=z
            ypl(nvar)=zpout
            if (xpl(nvar).gt.xplmax) xplmax=xpl(nvar)
            if (xpl(nvar).lt.xplmin) xplmin=xpl(nvar)
            if (ypl(nvar).gt.yplmax) yplmax=ypl(nvar)
            if (ypl(nvar).lt.yplmin) yplmin=ypl(nvar)
          endif

        endif !moduname

9123    z=z+dz

      enddo !z

      close(88)
      close(89)

      if (modunameo.eq.2) then
        close(71)
        close(72)
      endif

      if (moduname.eq.1) then

        if (iplot.ne.0) then
          dx=(xplmax-xplmin)*0.1
          xplmax=xplmax+dx
          xplmin=xplmin-dx
          dy=(yplmax-yplmin)*0.1
          yplmax=yplmax+dy
          yplmin=yplmin-dy
          call hlimit(50000)
          CALL HPLINT(1)
          CALL system('rm -f dyn_multipol_fit.eps')
          CALL IGMETA(40,-113)
          CALL HTITLE(comment)
          CALL HPLOPT('DATE',1)
          CALL HPLOPT('GRID',1)
          call igset('MTYP',24.)
          call hplfra(xplmin,xplmax,yplmin,yplmax,' ')
          call hplax(
     &      'horizontal diplacement "M#mm"N#','horizontal kick "M#mrad"N#')
          call ipl(nvar,xpl,ypl)
          call ipm(nvar,xpl,ypl)
          call iuwk(0,0)
          call igmeta(0,0)
          call system("mv fort.40 dyn_multipole_fit.eps")
          print*,'--- Enter RETURN key to continue!'
          read(5,'(a)')c1
          call hplend
          close(40)
        endif

        stop

      endif !moduname

c actual loop}


      chi2=0.0
      phi=shift/period*twopi
      sn2=sin(phi)**2
      cs22=cos(phi/2.0d0)**2
      sn22=1.0d0-cs22

      open(unit=80,file='ue112_dyn_multipole.loop',status='old')

      nvar=0
3     continue

      read(80,*,end=93)zi,zf0,zpf0,yf0,ypf0
c      ,ziia,zia0,zpia0,yia0,ypia0

      if (zi.eq.-9999.) goto 3

      if (zi.lt.z0(1).and.zi.ge.z0(1)-(z0(2)-z0(1))) then

        xx(1:3)=z0(1:3)
        yy(1:3)=f0(1:3)
        call util_parabel(xx,yy,a,yp,xopt,aopt,istat1)
        if (istat1.ne.0) then
          print*,'*** Error in UTIL_PARABEL ***'
          goto 3
        endif
        ff0=a(1)+a(2)*zi+a(3)*zi*zi

      else if (zi.gt.z0(nzf0).and.zi.le.z0(nzf0)+(z0(nzf0)-z0(nzf0-1))) then

        xx(1:3)=z0(nzf0-2:nzf0)
        yy(1:3)=f0(nzf0-2:nzf0)
        call util_parabel(xx,yy,a,yp,xopt,aopt,istat1)
        if (istat1.ne.0) then
          print*,'*** Error in UTIL_PARABEL ***'
          goto 3
        endif
        ff0=a(1)+a(2)*zi+a(3)*zi*zi

      else
        call util_spline_inter_status(z0,f0,y2f0,nzf0,zi,ff0,-1,istat1)
        if (istat1.ne.0) goto 3
      endif

      if (zi.lt.zpi(1).and.zi.ge.zpi(1)-(zpi(2)-zpi(1))) then

        xx(1:3)=zpi(1:3)
        yy(1:3)=fpi(1:3)
        call util_parabel(xx,yy,a,yp,xopt,aopt,istat1)
        if (istat1.ne.0) then
          print*,'*** Error in UTIL_PARABEL ***'
          goto 3
        endif
        ffpi=a(1)+a(2)*zi+a(3)*zi*zi

      else if (zi.gt.zpi(nzf0).and.zi.le.zpi(nzf0)+(zpi(nzf0)-zpi(nzf0-1))) then

        xx(1:3)=zpi(nzfpi-2:nzfpi)
        yy(1:3)=fpi(nzfpi-2:nzfpi)
        call util_parabel(xx,yy,a,yp,xopt,aopt,istat1)
        if (istat1.ne.0) then
          print*,'*** Error in UTIL_PARABEL ***'
          goto 3
        endif
        ffpi=a(1)+a(2)*zi+a(3)*zi*zi

      else
        call util_spline_inter_status(zpi,fpi,y2fpi,nzfpi,zi,ffpi,0,istat2)
        if (istat2.ne.0) goto 3
      endif

      nvar=nvar+1

      dyn=ff0*cs22+ffpi*sn22

      vdyn(1,nvar)=zi
      vdyn(2,nvar)=zpf0
      if (sn2.ne.0.0d0) then
        vdyn(3,nvar)=dyn-zpf0/sn2
      else
        vdyn(3,nvar)=dyn-zpf0
      endif
      vdyn(4,nvar)=zf0

      if (iplot.ne.0) then
        xpl(nvar)=vdyn(1,nvar)
        ypl(nvar)=vdyn(3,nvar)*brho
        if (xpl(nvar).gt.xplmax) xplmax=xpl(nvar)
        if (xpl(nvar).lt.xplmin) xplmin=xpl(nvar)
        if (ypl(nvar).gt.yplmax) yplmax=ypl(nvar)
        if (ypl(nvar).lt.yplmin) yplmin=ypl(nvar)
      endif

      chi2=chi2+vdyn(3,nvar)**2

      goto 3

93    close(80)

      call system("mv dyn_multipol.out dyn_multipol.out.bck")
      call system
     &("mv ue112_dyn_multipole_lshimmed.zibint ue112_dyn_multipole_lshimmed.zibint.bck")
      call system
     &("mv ue112_dyn_multipole_lshimmed.zfbint ue112_dyn_multipole_lshimmed.zfbint.bck")

      open(unit=70,file='dyn_multipol.out',status='new')
      open(unit=71,file='ue112_dyn_multipole_lshimmed.zibint',status='new')
      open(unit=72,file='ue112_dyn_multipole_lshimmed.zfbint',status='new')

      write(70,*) chi2

      do ipar=1,nmagmod
        write(70,*)thestart(ipar)
      enddo

      do ivar=1,nvar
        write(70,*)vdyn(1,ivar),vdyn(2,ivar),vdyn(3,ivar)
        write(71,*)vdyn(1,ivar),vdyn(3,ivar)*brho,'-9999.'
        write(72,*)vdyn(4,ivar),vdyn(3,ivar)*brho,'-9999.'
      enddo

      close(70)
      close(71)
      close(72)

      if (iplot.ne.0) then
        dx=(xplmax-xplmin)*0.1
        xplmax=xplmax+dx
        xplmin=xplmin-dx
        dy=(yplmax-yplmin)*0.1
        yplmax=yplmax+dy
        yplmin=yplmin-dy
        call hlimit(50000)
        CALL HPLINT(1)
        CALL system('rm -f dyn_multipol_fit.eps')
        CALL IGMETA(40,-113)
        CALL HTITLE(comment)
        CALL HPLOPT('DATE',1)
        call igset('MTYP',24.)
        call hplfra(xplmin,xplmax,yplmin,yplmax,' ')
        call ipl(nvar,xpl,ypl)
        call ipm(nvar,xpl,ypl)
        call iuwk(0,0)
        call igmeta(0,0)
        call system("mv fort.40 dyn_multipole_fit.eps")
        print*,'--- Enter RETURN key to continue!'
        read(5,'(a)')c1
        call hplend
        close(40)
      endif

      stop
      end
