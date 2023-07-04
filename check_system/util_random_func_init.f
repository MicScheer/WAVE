*CMZ :          23/11/2021  15.41.27  by  Michael Scheer
*CMZ : 00.00/16 03/11/2014  14.34.42  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2014
      subroutine util_random_func_init(func,n,x,fint,coef)

c Initializes the arrays fint and coef for the external
c function util_random_func

      implicit none

      real*8 fint(*),x(*),coef(*)
      real*8, dimension (:), allocatable :: w1,w2,w3,w4,fwork

      external func
      real func,ff

      integer i,n

      allocate(fwork(n))
      allocate(w1(n))
      allocate(w2(n))
      allocate(w3(n))
      allocate(w4(n))

      do i=1,n
        ff=func(sngl(x(i)))
        fwork(i)=dble(ff)
      enddo

      call util_spline_running_integral(x,fwork,n,fint,coef,w1,w2,w3,w4)
      fint(1:n)=fint(1:n)/fint(n)

      do i=2,n
        if (fint(i).le.fint(i-1)) then
          print*,"*** Error in util_random_init: Integral of function not monoton, can not be reversed ***"
          stop
        endif
      enddo

      call util_spline_coef(fint,x,n,0.0d0,0.0d0,coef,w1,w2,w3,w4)

      deallocate(fwork)
      deallocate(w1)
      deallocate(w2)
      deallocate(w3)
      deallocate(w4)

      return
      end
