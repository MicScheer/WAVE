*CMZ :          15/11/2018  15.14.34  by  Michael Scheer
*CMZ : 00.00/15 12/10/2013  12.22.24  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  17.37.24  by  Michael Scheer
*-- Author :    Michael Scheer   06/07/2007
      subroutine util_spline_adjust(ifail)

c Tool to find best x-value meeting given y-value of a function

      implicit none

      integer ifail,npoi,i

      double precision xx,yy,yf,xf,yold
      double precision, dimension (:), allocatable :: x,y,x2p,w1,w2,w3,w4

      ifail=0

      open(unit=99,file='util_spline_adjust.dat',status='old')

      call util_skip_comment_end(99,ifail)
      if (ifail.ne.0) then
        print*,
     &    '*** Error 1 in UTIL_SPLINE_ADJUST: No data on file util_spline_adjust.dat'
        ifail=1
        return
      endif
      read(99,*)yf

      do i=1,3
        call util_skip_comment_end(99,ifail)
        if (ifail.ne.0) then
          print*,
     &      '*** Error 2 in UTIL_SPLINE_ADJUST: No data on file util_spline_adjust.dat'
          ifail=2
          return
        endif
        read(99,*)xx,yy
      enddo

      npoi=3

1     call util_skip_comment_end(99,ifail)
      if (ifail.eq.0) then
        read(99,*)xx,yy
        npoi=npoi+1
        goto 1
      endif

      allocate(x(npoi))
      allocate(y(npoi))
      allocate(x2p(npoi))
      allocate(w1(npoi))
      allocate(w2(npoi))
      allocate(w3(npoi))
      allocate(w4(npoi))

      rewind(99)
      call util_skip_comment_end(99,ifail)
      if (ifail.ne.0)
     &  stop '*** Error in UTIL_SPLINE_ADJUST: No data on file util_spline_adjust.dat'
      read(99,*)yf

      do i=1,npoi
        call util_skip_comment_end(99,ifail)
        read(99,*)x(i),y(i)
      enddo

      close(99)

      yold=y(1)
      do i=2,npoi
        if (y(i).eq.yold) then
          print*,
     &      '*** Error 3 in UTIL_SPLINE_ADJUST: Data not monoton on file util_spline_adjust.dat'
          ifail=3
          return
        endif
        yold=y(i)
      enddo

      call util_spline_coef(y,x,npoi,9999.0d0,9999.0d0,x2p,w1,w2,w3,w4)
      call util_spline_inter(y,x,x2p,npoi,yf,xf,-1)

      print*,'result of util_spline_adjust:',xf

      return
      end
