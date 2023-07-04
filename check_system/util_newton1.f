*CMZ :          05/09/2020  08.29.39  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_newton1(x,f,g,step,xtol,ftol,maxiter,xold,
     &  fcn,ifail)

! This routine tries to find x such, that f=0

! Input:
!    x: Start value
! step:  Initial relative step size
!     ftol: Accuracy of f
!     xtol: Accuracy of x
!  maxiter: Max. number of iterations
!  xold: Workingspace

! external function f=fcn(...,x,...)
! call fcn(x,f,g,ifail)


! Output:
!       x: Solution
!       g: gradiant df/dx
!       f: Residual
! maxiter: Last iteration
!   ifail: Status

      implicit none

      double precision ftol,xtol,x,f,fold,g,xold,step
      integer iter,n,maxiter,ifail,iconv,jconv

      external fcn

      ifail=0

      call fcn(x,f,g,ifail)
      if (ifail.ne.0) then
        ifail=-1
        goto 9999
      endif

      do iter=1,maxiter

        iconv=1
        jconv=1

        xold=x
        fold=f

        if (abs(g).gt.1.0d-15) x=(g*xold-fold)/g
        if (abs(x-xold).gt.xtol) jconv=0

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
