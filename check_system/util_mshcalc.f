*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ :  1.10/00 09/11/2016  13.57.10  by  Michael Scheer
*-- Author :    Michael Scheer   13/10/2016
      subroutine  util_mshcalc(clinein,res)

      implicit none

      integer ndimp,nfunp
      parameter (ndimp=1000,nfunp=18)

      double precision res,xopy(ndimp)
      double precision :: pi=3.141592653589793d0

      integer ndim,istat
      integer lenc,jc,ic1,ic,kbo,kbc,iw,nfirst,nlast,ianf,iend,kfun
      integer :: nres=0,ifun=0

      character(*) clinein
      character(2048) cline
      character(4096) cbrack,cwork
      character(32) c32,cfun(nfunp)
      character c1

      equivalence(c1,ic1)
      ic1=0

      cfun(1)='sind'
      cfun(2)='cosd'
      cfun(3)='tand'
      cfun(4)='cotd'
      cfun(5)='sin'
      cfun(6)='cos'
      cfun(7)='tan'
      cfun(8)='cot'
      cfun(9)='asind'
      cfun(10)='acosd'
      cfun(11)='atand'
      cfun(12)='acotd'
      cfun(13)='asin'
      cfun(14)='acos'
      cfun(15)='atan'
      cfun(16)='acot'
      cfun(17)='sqrt'
      cfun(18)='int'

      ndim=ndimp

      call util_string_trim(clinein,nfirst,nlast)
      cline=''
      cline=clinein(nfirst:nlast)

1     call util_string_trim(cline,nfirst,nlast)
      lenc=nlast-nfirst+1
      if (lenc.eq.0) goto 9999

      cwork=''
      lenc=nlast-nfirst+1
      cwork=cline(nfirst:nlast)
      cline=''
      iw=0
      ic=0
      do jc=1,lenc
        ic=ic+1
        c1=cwork(ic:ic)
        if (c1.ne.' ') then
          iw=iw+1
          cline(iw:iw)=c1
        endif
      enddo
      call util_string_trim(cline,nfirst,nlast)
      lenc=nlast-nfirst+1

      kfun=0

      do ifun=1,nfunp
        ianf=1
        iend=lenc
        call util_string_substring(cline,trim(cfun(ifun)),ianf,iend,istat)
        if (istat.eq.0) then
          kfun=ifun
          exit
        endif
      enddo !all known functions

      kbo=0
      kbc=0

      call util_string_trim(cline,nfirst,nlast)
      lenc=nlast-nfirst+1
      do ic=1,lenc
        c1=cline(ic:ic)
        if (c1.eq.'[') then
          kbo=ic+1
        else if (c1.eq.']') then
          kbc=ic-1
          exit
        endif
      enddo !lenc

      if (kbo.gt.0.and.kfun.eq.0) then
        print*,"*** Error in util_mshcalc: Unknown function in line"
        call util_string_trim(cline,nfirst,nlast)
        print*,cline(nfirst:nlast)
        stop
      endif

      if (kbo.gt.0.and.kbc.le.kbo) then
        print*,"*** Error in util_mshcalc: Missing bracket in line"
        call util_string_trim(cline,nfirst,nlast)
        print*,cline(nfirst:nlast)
        stop
      endif

      if (kbo.gt.0) then
        cbrack=cline(kbo:kbc)
      else
        cbrack=cline
      endif

      call util_mshcalc_basics(cbrack,xopy,res)
      nres=nres+1
      xopy(nres)=res

      if (kbo.eq.0.or.kfun.eq.0) goto 9999

      write(c32,*)nres
      call  util_string_trim(c32,nfirst,nlast)

      if (kfun.eq.1) then
        xopy(nres)=sin(res*pi/180.0d0)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.2) then
        xopy(nres)=cos(res*pi/180.0d0)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.3) then
        xopy(nres)=tan(res*pi/180.0d0)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.4) then
        xopy(nres)=1.0d0/tan(res*pi/180.0d0)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.5) then
        xopy(nres)=sin(res)
        cline(kbo-4:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.6) then
        xopy(nres)=cos(res)
        cline(kbo-4:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.7) then
        xopy(nres)=tan(res)
        cline(kbo-4:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.8) then
        xopy(nres)=1.0d0/tan(res)
        cline(kbo-4:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.9) then
        xopy(nres)=asin(res)/pi*180.0d0
        cline(kbo-6:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.10) then
        xopy(nres)=acos(res)/pi*180.0d0
        cline(kbo-6:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.11) then
        xopy(nres)=atan(res)/pi*180.0d0
        cline(kbo-6:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.12) then
        xopy(nres)=atan(1.0d0/res)/pi*180.0d0
        cline(kbo-6:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.13) then
        xopy(nres)=asin(res)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.14) then
        xopy(nres)=acos(res)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.15) then
        xopy(nres)=atan(res)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.16) then
        xopy(nres)=atan(1.0d0/res)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.17) then
        xopy(nres)=sqrt(res)
        cline(kbo-5:kbc+1)='x'//c32(nfirst:nlast)
      else if (kfun.eq.18) then
        xopy(nres)=int(res)
        cline(kbo-4:kbc+1)='x'//c32(nfirst:nlast)
      endif

c      call util_string_trim(cline,nfirst,nlast)
c      print*,cline(nfirst:nlast)

      goto 1

9999  continue

      res=xopy(nres)

      return
      end
