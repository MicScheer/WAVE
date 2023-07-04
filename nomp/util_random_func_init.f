*CMZ :  3.02/03 03/11/2014  14.19.56  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2014
      subroutine util_random_func_init(func,n,x,fint,coef)

c Initializes the array fwork for the external function func to be used
c by subroutine util_random_func

      implicit none

      real*8 fint(*),x(*),coef(*)
      real*8, dimension (:), allocatable :: w1,w2,w3,w4,fwork

      external function func
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
      call util_spline_coef(fint,x,n,0.0d0,0.0d0,coef,w1,w2,w3,w4)

      deallocate(fwork)
      deallocate(w1)
      deallocate(w2)
      deallocate(w3)
      deallocate(w4)

      return
      end
