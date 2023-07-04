*CMZ :  4.00/15 27/04/2022  08.09.29  by  Michael Scheer
*CMZ :  3.06/00 14/01/2019  17.27.04  by  Michael Scheer
*CMZ :  3.03/02 18/01/2016  13.02.27  by  Michael Scheer
*CMZ :  3.02/03 30/10/2014  17.17.28  by  Michael Scheer
*-- Author :    Michael Scheer   05/09/2014
      subroutine util_random_gauss(n,g,rr)

      implicit none

c Based on textbook "Numerical Recipies"
c The subroutine util_random is called, which can be initialized
c and controlled by util_random_init, util_random_set_seed, and
c util_random_get_seed.

      real g(n),r,v1,v2,fac,rr(2)
      integer n,i

      if (n.eq.1) then
1       call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 1
        fac=sqrt(-2.*log(r)/r)
        g(1)=v1*fac
      else
        do i=1,(n/2*2),2
11      call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 11
        fac=sqrt(-2.*log(r)/r)
        g(i)=v1*fac
        g(i+1)=v2*fac
        enddo
      endif

      if (mod(n,2).ne.0) then
12      call util_random(2,rr)
        v1=2.*rr(1)-1.
        v2=2.*rr(2)-1.
        r=v1*v1+v2*v2
        if (r.ge.1..or.r.eq.0.0) goto 12
        fac=sqrt(-2.*log(r)/r)
        g(n)=v1*fac
      endif

      return
      end
