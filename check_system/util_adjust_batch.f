*CMZ :          05/09/2020  08.57.27  by  Michael Scheer
*CMZ : 00.00/15 12/04/2012  15.10.13  by  Michael Scheer
*CMZ : 00.00/07 04/06/2010  14.31.08  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  17.37.24  by  Michael Scheer
*-- Author :    Michael Scheer   06/07/2007
      subroutine util_adjust_batch(ifail)

c Tool to find best x-value meeting given y-value of a function

      implicit none

      integer ifail,i,maxiter,ibad,iter

      double precision xfirst,xstep,xf,yf,yerror,ydum,errmax,xbad
      double precision x(2),y(2),a,b

      character(1024)c1024
      character(2048)c2048

      ifail=0

      call system("pwd > util_adjust_batch.pwd")
      open(unit=99,file="util_adjust_batch.pwd")
      read(99,'(a)')c1024
      close(99)

      open(unit=99,file='util_adjust_batch.in',status='old')

      call util_skip_comment_end(99,ifail)
      if (ifail.ne.0) then
        print*,
     &    '*** Error 1 in util_adjust_batch: No data on file util_adjust_batch.in'
        ifail=1
        return
      endif

      read(99,*)xfirst,xstep
      read(99,*)yf
      read(99,*)yerror,maxiter

      x(1)=xfirst
      x(2)=xfirst+xstep

      open(unit=20,file='util_adjust_batch.log',status='unknown')

      errmax=-1.0d30
      do i=1,2

        open(unit=99,file='util_adjust_batch-prog-in.dat',status='unknown')
        write(99,*)x(i)
        close(99)

        print*,' '
        print*,'x:',x(i)
        print*,' '
        c2048=". "//trim(c1024)//"/util_adjust_batch.sh"
        call system(trim(c2048))

        open(unit=99,file='util_adjust_batch-prog-out.dat',status='old')
        read(99,*)y(i)
        close(99)

        write(20,*)x(i),y(i)

      enddo

      if (maxiter.lt.1) then

        do i=1,2
          write(6,*)x(i),y(i)
        enddo

      endif

      do iter=1,maxiter

        print*,' '
        print*,'Iteration:',iter
        print*,' '

        do i=1,2
          print*, x(i),y(i)
        enddo
        print*,' '

        a=(y(2)-y(1))/(x(2)-x(1))

        if (a.eq.0.0d0) then
          xf=(xbad+x(1))/2.0d0
        else
          b=y(2)-a*x(2)
          xf=(yf-b)/a
        endif

        open(unit=99,file='util_adjust_batch-prog-in.dat',status='unknown')
        write(99,*)xf
        close(99)

        print*,' '
        print*,'x:',xf
        print*,' '
        call system(trim(c2048))

        open(unit=99,file='util_adjust_batch-prog-out.dat',status='old')
        read(99,*)ydum
        close(99)

        errmax=-1.0d30
        do i=1,2
          if (abs(y(i)-yf).gt.errmax) then
            ibad=i
            errmax=abs(y(i)-yf)
          endif
        enddo
        xbad=x(ibad)
        x(ibad)=xf
        y(ibad)=ydum

        write(20,*)xf,ydum

        if (abs(ydum-yf).le.yerror) then
          goto 9
        endif

      enddo !iter

9     continue

      print*,' '
      print*,'Iteration:',iter
      print*,' '

      do i=1,2
        print*, x(i),y(i)
      enddo
      print*,' '

      print*,' '
      if (abs(ydum-yf).le.yerror) then
        print*,'util_adjust_batch converged'
      else
        print*,'util_adjust_batch did NOT converge'
      endif

      print*,' '
      do i=1,2
        if(x(i).eq.xf) print*,'xf, yf:',x(i),y(i)
      enddo

      open(unit=99,file='util_adjust_batch.out',status='unknown')
      do i=1,2
        write(99,*)x(i),y(i)
      enddo
      do i=1,2
        if(x(i).eq.xf) write(99,*)'* xf, yf:',x(i),y(i)
      enddo
      close(99)

      return
      end
