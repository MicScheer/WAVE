*CMZ :  3.02/03 27/10/2014  17.36.49  by  Michael Scheer
*-- Author :    Michael Scheer   27/10/2014
      subroutine util_random_func(n,x,fint,coef,nran,ran)

      implicit none

      real*8 x(*),fint(*),coef(*),r8
      real ran(*),rn(1000)

      integer i,n,loop,l,nran,l1

      loop=nran/1000

      if (loop.gt.0) then
        do l=1,loop
          call util_random(1000,rn)
          l1=(l-1)*1000
          do i=1,1000
            call util_spline_inter(fint,x,
     &        coef,n,dble(rn(i)),r8,-1)
            ran(l1+i)=sngl(r8)
          enddo
        enddo
      endif

      call util_random(nran-loop*1000,rn)
      l1=loop*1000
      do i=1,nran-loop*1000
        call util_spline_inter(fint,x,coef,n,dble(rn(i)),r8,-1)
        ran(l1+i)=sngl(r8)
      enddo

      return
      end
