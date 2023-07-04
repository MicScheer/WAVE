*CMZ : 00.00/15 12/04/2012  15.10.13  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/06 09/07/2007  08.27.20  by  Michael Scheer
*-- Author :    Michael Scheer   06/07/2007
      subroutine util_spline_adjust_minmax(ifail)

c Tool to find best x-value meeting given y-value of a function

      implicit none

      integer ifail,npoi,i,jfail

      double precision xx,yy,xpf,fopt
      double precision, dimension (:), allocatable :: x,y,yp,y2p,w1,w2,w3,w4

      ifail=0

      open(unit=99,file='util_spline_adjust_minmax.dat',status='old')

      do i=1,3
        call util_skip_comment_end(99,ifail)
        if (ifail.ne.0) then
          print*,
     &'*** Error 2 in UTIL_SPLINE_adjust_minmax: Not enough data on file util_spline_adjust_minmax.dat'
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
      allocate(yp(npoi))
      allocate(y2p(npoi))
      allocate(w1(npoi))
      allocate(w2(npoi))
      allocate(w3(npoi))
      allocate(w4(npoi))

      rewind(99)

      do i=1,npoi
        call util_skip_comment_end(99,ifail)
        read(99,*)x(i),y(i)
      enddo

      close(99)

      call util_max_parabel(npoi,x,y,xpf,fopt,w1,w2,jfail)

      print*,'Result of util_max_parabel:'
      print*,xpf,fopt

      if (ifail.ne.0.and.jfail.eq.0) ifail=2

      call util_spline_coef_deriv(
     &  x,y,npoi,-9999.0d0,-9999.0d0,yp,y2p,w1,w2,w3,w4)
      call util_sort_func(npoi,yp,x)
      call util_spline_interpolation(npoi,yp,x,0.0d0,xpf,y2p,w1,w2,w3,w4,-1)
      call util_sort_func(npoi,x,y)
      call util_spline_interpolation_stat(npoi,x,y,xpf,fopt,y2p,w1,w2,w3,w4,-1
     &  ,ifail)

      if (ifail.ne.0) then
        print*,'*** ERROR: UTIL_SPLINE_ADJUST_MINMAX FAILED ***'
      else
        print*,'Result of util_spline_adjust_minmax:'
        print*,xpf,fopt
      endif

      return
      end
