*CMZ :          05/09/2020  08.51.27  by  Michael Scheer
*CMZ : 00.00/15 12/10/2013  12.17.44  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  17.37.24  by  Michael Scheer
*-- Author :    Michael Scheer   06/07/2007
      subroutine util_spline_adjust_batch(ifail)

c Tool to find best x-value meeting given y-value of a function

      implicit none

      integer ifail,npoi,i,maxiter,maxspline,ibad,iter,nstart

      double precision xfirst,xstep,xf,yf,yerror,errmax,ydum
      double precision, dimension (:), allocatable :: x,y,x2p,w1,w2,w3,w4

      ifail=0

      open(unit=99,file='util_spline_adjust_batch.in',status='old')

      call util_skip_comment_end(99,ifail)
      if (ifail.ne.0) then
        print*,
     &    '*** Error 1 in UTIL_spline_adjust_batch: No data on file util_spline_adjust_batch.dat'
        ifail=1
        return
      endif

      read(99,*)xfirst,xstep
      read(99,*)yf
      read(99,*)yerror,maxiter,maxspline

      allocate(x(maxspline))
      allocate(y(maxspline))
      allocate(x2p(maxspline))
      allocate(w1(maxspline))
      allocate(w2(maxspline))
      allocate(w3(maxspline))
      allocate(w4(maxspline))

      nstart=0
11    read(99,*,end=99)xf,ydum
      nstart=nstart+1
      x(nstart)=xf
      y(nstart)=ydum
      goto 11

99    close(99)

      if (nstart.gt.0) call util_sort_func(nstart,x,y)

      npoi=3

      if (nstart.eq.0) then
        nstart=1
        x(1)=xfirst
        x(2)=x(1)+xstep
        x(3)=x(2)+xstep
      else if (nstart.eq.1) then
        nstart=2
        x(2)=x(1)+xstep
        x(3)=x(2)+xstep
      else if (nstart.eq.2) then
        nstart=3
        x(3)=x(2)+xstep
      else
        npoi=nstart
        nstart=npoi+1
      endif

      call util_sort_func(npoi,x,y)

      errmax=-1.0d30
      print*,' '
      do i=nstart,npoi

        open(unit=99,file='util_spline_adjust_batch-prog-in.dat',status='unknown')
        write(99,*)x(i)
        close(99)

        call system("util_spline_adjust_batch.sh")

        open(unit=99,file='util_spline_adjust_batch-prog-out.dat',status='old')
        read(99,*)y(i)
        close(99)

        print*,'x:',x(i),y(i)
        write(20,*)x(i),y(i)

      enddo

      open(unit=20,file='util_spline_adjust_batch.log',status='unknown')
      call util_sort_func(npoi,x,y)
      do i=1,npoi
        write(20,*)x(i),y(i)
      enddo

      do iter=1,maxiter

        call util_sort_func(npoi,x,y)

        print*,' '
        print*,'Iteration:',iter
        print*,' '

        do i=1,npoi
          print*, x(i),y(i)
        enddo
        print*,' '

        if (y(2).ge.y(1)) then
          do i=2,npoi
            if (y(i).le.y(i-1)) stop '*** y-values not strictly monotonic ***'
          enddo
        else
          do i=2,npoi
            if (y(i).ge.y(i-1)) stop '*** y-values not strictly monotonic ***'
          enddo
        endif

        call util_sort_func(npoi,y,x)
        call util_spline_coef(y,x,npoi,9999.0d0,9999.0d0,x2p,w1,w2,w3,w4)
        call util_spline_inter(y,x,x2p,npoi,yf,xf,-1)

        open(unit=99,file='util_spline_adjust_batch-prog-in.dat',status='unknown')
        write(99,*)xf
        close(99)

        print*,' '
        print*,'x:',xf
        print*,' '
        call system("util_spline_adjust_batch.sh")

        open(unit=99,file='util_spline_adjust_batch-prog-out.dat',status='old')
        read(99,*)ydum
        close(99)

        if (npoi.lt.maxspline) then
          npoi=npoi+1
          x(npoi)=xf
          y(npoi)=ydum
        else
          errmax=-1.0d30
          do i=1,npoi
            if (abs(y(i)-yf).gt.errmax) then
              ibad=i
              errmax=abs(y(i)-yf)
            endif
          enddo
          x(ibad)=xf
          y(ibad)=ydum
        endif

        write(20,*)xf,ydum

        if (abs(ydum-yf).le.yerror) then
          goto 9
        endif

      enddo !iter

9     continue

      call util_sort_func(npoi,x,y)

      print*,' '
      print*,'Iteration:',iter
      print*,' '

      do i=1,npoi
        print*, x(i),y(i)
      enddo
      print*,' '

      print*,' '
      if (abs(ydum-yf).le.yerror) then
        print*,'util_spline_adjust_batch converged'
      else
        print*,'util_spline_adjust_batch did NOT converge'
      endif
      print*,' '
      do i=1,npoi
        if(x(i).eq.xf) print*,'xf, yf:',x(i),y(i)
      enddo

      open(unit=99,file='util_spline_adjust_batch.out',status='unknown')
      do i=1,npoi
        write(99,*)x(i),y(i)
      enddo
      do i=1,npoi
        if(x(i).eq.xf) write(99,*)'* xf, yf:',x(i),y(i)
      enddo
      close(99)

      return
      end
