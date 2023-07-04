*CMZ :  4.01/02 07/05/2023  11.33.00  by  Michael Scheer
*CMZ :  4.00/16 08/08/2022  14.33.41  by  Michael Scheer
*CMZ :  4.00/15 05/07/2022  14.01.32  by  Michael Scheer
*CMZ :  4.00/09 14/08/2020  07.07.05  by  Michael Scheer
*CMZ :          10/08/2020  10.48.44  by  Michael Scheer
      subroutine undumagwav(xin,y,z,bx,by,bz,ifail)

      implicit none
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,undumagc.
      include 'undumagc.cmn'
*KEND.

      double precision xin,x,y,z,bx,by,bz,xr,yr,zr
      double precision, save :: xper=0.0d0,xpermin=0.0d0,xpermax=0.0d0
      real g(3)
      integer ifail

      integer, save :: iwarnu=0,ical=0

      if (ical.eq.0) then
        if (kbundumag.eq.1.or.kbundumag.eq.2) then
          if (nunduperw.gt.nunduclc) then
            xper=perlenclc*(nunduperw-nunduclc)
          endif
        else if (kbundumag.eq.3.or.kbundumag.eq.4) then
          if (nperiodw_h.gt.nunduclc) then
            xper=perlenclc*(nperiodw_h-nunduclc)
          endif
        endif
        xpermin=-xper/2.0d0
        xpermax=-xpermin
        ical=1
      endif

      if (xper.gt.0.0d0) then
        if (xin.gt.xpermin.and.xin.lt.xpermax) then
          x=mod(xin-xpermin,perlenclc)
        else if (xin.le.xpermax) then
          x=xin-xpermin
        else
          x=xin-xpermax
        endif
      else
        x=xin
      endif

      if (kbundumap_c.eq.0) then
        call util_random(3,g)
        g=g-0.5
        if (kbundumag.eq.1.or.kbundumag.eq.2) then
          xr=x+urandox*g(1)
          yr=y+urandoy*g(2)
          zr=z+urandoz*g(3)
        else if (kbundumag.eq.3.or.kbundumag.eq.4) then
          xr=x+urandox_h*g(1)
          yr=y+urandoy_h*g(2)
          zr=z+urandoz_h*g(3)
        else
          xr=x+(urandox+urandox_h)/2.*g(1)
          yr=y+(urandoy+urandoy_h)/2.*g(2)
          zr=z+(urandoz+urandoz_h)/2.*g(3)
        endif
      else
        xr=x
        yr=y
        zr=z
      endif

      call undumag_field(xr,yr,zr,bx,by,bz,ifail)

      if (ifail.ne.0.and.ifail.ne.-1) then
        if (iwarnu.eq.0) then
          print*," "
          print*,"WARNING IN UNDUMAGWAV: BAD RETURN FROM UNDUMAG_FIELD"
          print*,sngl(x),sngl(y),sngl(z),sngl(bx),sngl(by),sngl(bz)
          print*,"*** Further warnings suppressed ***"
          write(lungfo,*)" "
          write(lungfo,*)"*** WARNING IN UNDUMAGWAV: BAD RETURNS FROM UNDUMAG_FIELD ***"
          iwarnu=1
        endif
      endif

      return
      end
