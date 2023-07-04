*CMZ :          29/04/2017  12.15.05  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_newton(n,x,f,g,step,xtol,ftol,maxiter,xold,
     &  fcn,ifail)

! This routine tries to find x(1:n) such, that f=0

! Input:
!        n: Number of x-values
!    x(n): Start value
! step(n):  Initial relative step size
!     ftol: Accuracy of f
!     xtol: Accuracy of x
!  maxiter: Max. number of iterations
!  xold(n): Workingspace

! external function f=fcn(...,x,...)
! call fcn(n,x,f,g,ifail)


! Output:
!       x: Solution
!       g: gradiant df/dx(i)
!       f: Residual
! maxiter: Last iteration
!   ifail: Status

      implicit none

      double precision fg,ftol,xtol,x(n),f,fold,g(n),xold(n),step(n)
      integer i,iter,n,maxiter,ifail,iconv,jconv

      external fcn

      ifail=0

      call fcn(n,x,f,g,ifail)
      if (ifail.ne.0) then
        ifail=-1
        goto 9999
      endif

      do iter=1,maxiter

        iconv=1
        jconv=1

        xold=x
        fold=f

        do i=1,n
          if (abs(g(i)).gt.1.0d-15) x(i)=(g(i)*xold(i)-fold)/g(i)
          if (abs(x(i)-xold(i)).gt.xtol) jconv=0
        enddo !n

        call fcn(n,x,f,g,ifail)
        print*,iter,x,f

        if (ifail.ne.0) then
          ifail=-2
          goto 9999
        endif

        if (abs(f).gt.ftol) iconv=0

        if (iconv.ne.0.or.jconv.ne.0) then
          maxiter=iter
          exit
        endif

      enddo !maxiter

      f=f-fold
      step=x-xold

9999  continue

      return
      end
