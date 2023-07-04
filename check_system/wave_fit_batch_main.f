*CMZ :  4.00/11 30/05/2021  14.27.36  by  Michael Scheer
*CMZ :  2.66/13 02/06/2010  14.32.47  by  Michael Scheer
*CMZ :  2.66/07 10/12/2009  16.14.26  by  Michael Scheer
*-- Author :    Michael Scheer   10/12/2009
      program wave_fit_batch_main

c---  main program to use minuit in FORTRAN-driven mode

      implicit none

*KEEP,WAVE_FIT_BATCH.

      integer nvarp,nvar,maxloop
      parameter (nvarp=100)

      double precision
     &  chi2,chi2min,
     &  varfin(nvarp),var(nvarp),paropt(nvarp),weight(nvarp)

      common/fncc/
     &  chi2,chi2min,varfin,var,paropt,weight,nvar,maxloop
*KEND.

      integer nparp
      parameter (nparp=100)

      external fcn
      external function futil
      double precision futil

      double precision parstart(nparp),parstep(nparp),dum(nparp)

      integer lread,lwrite,lsave,
     &  ivar,npar,icond,ipar

      character(50) c50
      character(4) c4
      character(10) c10

c     unit 80 fuer inputfile reserviert

      data lread/5/
      data lwrite/6/
      data lsave/7/

      call mninit(lread,lwrite,lsave)

      call system("touch wave_fit_batch.term")

      open(unit=99,file='wave_fit_batch.in',status='old')

      call util_skip_comment(99)
      read(99,'(a)')c50

      call util_skip_comment(99)
      read(99,*)npar,maxloop
      do ipar=1,npar
        call util_skip_comment(99)
        read(99,*)parstart(ipar),parstep(ipar)
      enddo

      call util_skip_comment(99)
      read(99,*)nvar
      do ivar=1,nvar
        call util_skip_comment(99)
        read(99,*)varfin(ivar),weight(ivar)
      enddo

      close(99)

      open(unit=70,file='wave_fit_batch.fit',status='new',recl=256)
      close(70)

      call mnseti(c50)

      do ipar=1,npar
        write(c4,'(i4)')ipar
        c10='param'//c4
        call mnparm(ipar,c10,parstart(ipar),parstep(ipar),0.0d0,0.0d0,icond)
        if (icond.ne.0) then
          print*, '*** Error for Parameter ',ipar
          stop
        endif
      enddo

      call mncomd(fcn,'MINIMIZE',icond,futil)

      write(6,*)
      write(6,*)'Minuit status after MINIMIZE'
      write(6,*)icond
      write(6,*)

      close(lsave)

      call fcn(npar,dum,dum(1),dum,3,futil)

      stop
      end

      include 'fcn_batch.f'
      include 'futil.f'
