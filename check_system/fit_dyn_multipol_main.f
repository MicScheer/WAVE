*CMZ :  2.61/02 07/03/2007  13.08.18  by  Michael Scheer
*CMZ :  2.57/05 13/12/2006  16.03.04  by  Michael Scheer
*CMZ :  2.33/07 04/05/2001  15.03.35  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.34.14  by  Michael Scheer
*CMZ :  1.00/00 04/08/97  12.38.30  by  Michael Scheer
*-- Author :  Michael Scheer
      program fit_dyn_multipol_main

c---  main program to use minuit in FORTRAN-driven mode

      implicit none

*KEEP,DYN_MULTIPOL.
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

      external fcn
      external futil
        double precision futil

        double precision thestart(nparp),dthe,chi2min,
     &    theamx,theamn,dum

        integer lread,lwrite,lsave,
     &    ivar,nvar,npar,icond,ipar

        REAL*4 RPAW(50000)
        COMMON/PAWC/RPAW

        character(50) c50
        character(4) c4
        character(10) c10

c     unit 80 für inputfile reserviert

      data lread/5/
      data lwrite/6/
        data lsave/7/

      call mninit(lread,lwrite,lsave)

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

        open(unit=99,file='fit_dyn_multipol.in',status='old')

        read(99,'(a)')c50
        read(99,*)mode
        read(99,*)nslice
        read(99,*)iplot,xplmin,xplmax,dz,yplmin,yplmax
        read(99,*)shift
        read(99,*)npar
        read(99,*)radin
        read(99,*)bcr
        read(99,*)dthe
        read(99,*)thestart(1)

        zmin=xplmin
        zmax=xplmax

        if (iplot.ne.0) then
          call hlimit(50000)
          CALL HPLINT(1)
          CALL system('rm -f fit_dyn_multipol_fit.eps')
          CALL IGMETA(40,-113)
          CALL HTITLE(c50)
          CALL HPLOPT('DATE',1)
          call hplfra(xplmin,xplmax,yplmin,yplmax,' ')
        endif

      if (thestart(1).ne.9999.0d0) then
          do ipar=2,npar
            read(99,*)thestart(ipar)
          enddo
        endif

        close(99)

        if (thestart(1).eq.9999.0d0) then
          open(unit=99,file='fit_dyn_multipol.out',status='old')
          read(99,*)chi2min
          do ipar=1,npar
            read(99,*)thestart(ipar)
          enddo
          close(99)
        endif

        if (dthe.eq.9999.0d0) then
          theamx=-1.0d30
          theamn=1.0d30
          do ipar=1,npar
            if (thestart(ipar).gt.theamx) theamx=abs(thestart(ipar))
            if (thestart(ipar).lt.theamn) theamn=abs(thestart(ipar))
          enddo
          if (theamn.eq.0.0d0.and.theamx.eq.0.0d0) then
            print*,'*** Error in fit_dyn_multipol_main:'
            print*,'    No default value for step size defined, since'
            print*,'    all parameters are zero.'
            print*,'    Set value in file fit_dyn_multipol.in'
          else if (theamn.ne.0.0d0) then
              dthe=theamn/100.0d0
          else
              dthe=theamx/100.0d0
          endif
        endif

        call system("rm fit_dyn_multipol_fcn.term")
        call system("touch fit_dyn_multipol_fcn.term")

        call system("mv fit_dyn_multipol_fcn.fit fit_dyn_multipol_fcn.fit.bck")
      open(unit=70,file='fit_dyn_multipol_fcn.fit',status='new',recl=256)
      close(70)

        call mnseti(c50)

        nvar=npar

        do ivar=1,nvar
          write(c4,'(i4)')ivar
          c10='TheRot'//c4
          call mnparm(ivar,c10,thestart(ivar),dthe,0.0d0,0.0d0,icond)
          if (icond.ne.0) then
            print*, '*** Error for Parameter ',ivar
            stop
          endif
        enddo !nvar

        call mncomd(fcn,'MINIMIZE',icond,futil)

        write(6,*)
        write(6,*)'Minuit status after MINIMIZE'
        write(6,*)icond
        write(6,*)

      close(lsave)

        call fcn(npar,dum,dum,dum,3,futil)

        stop
      end

      include '/home/scheer/wav/cmz/fit_dyn_multipol_fcn.f'
      include '/home/scheer/wav/cmz/fcn_futil.f'
