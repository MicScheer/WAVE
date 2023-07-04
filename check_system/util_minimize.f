*CMZ :          12/05/2017  16.59.03  by  Michael Scheer
*-- Author :    Michael Scheer   12/05/2017
      subroutine util_minimize(xmin,fmin,step,xtol,ftol,maxiter,fcn,ifail)

      implicit none

      external fcn

      double precision xmin,fmin,step,xtol,ftol,
     &  x(3),y(3),yp(3),a(3),dfmax

      integer ifail,ix,maxiter,i,iter

      x(1)=xmin
      x(2)=x(1)+step
      x(3)=x(2)+step

      call fcn(x(1),y(1),ifail)
      call fcn(x(2),y(2),ifail)
      call fcn(x(3),y(3),ifail)

      do iter=1,maxiter

        call util_parabel(x,y,a,yp,xmin,fmin,ifail)

        if (ifail.ne.0) return

        if (Abs(fmin).le.ftol.or.abs(xmin).le.xtol) then
          ifail=0
          maxiter=iter
          return
        endif

        dfmax=-1.0d30
        do i=1,3
          if (abs(y(i)-fmin).gt.dfmax) then
            dfmax=abs(y(i)-fmin)
            ix=i
          endif
        enddo

        x(ix)=xmin
        call fcn(x(ix),y(ix),ifail)

      enddo

      maxiter=iter

      return
      end
