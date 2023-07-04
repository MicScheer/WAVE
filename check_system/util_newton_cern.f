*CMZ :          28/04/2017  11.03.47  by  Michael Scheer
*-- Author :
      subroutine util_newton_cern(n,x,f,ftol,xtol,maxf,iprt,info,sub,w)

      implicit none

! Input:
!      n: Number of equations and variables
!      x(1:n): Start values
!      ftol: Accuracy parameter for test 1
!      xtol: Accuracy parameter for test 2
!      maxf: Max. number of iterations, something like 50 * (n+3)
!      iprt: Print flag; 0: no printing, i: print x(i) at each iteration
!      w: workingspace of size n*(n+3)

! Output:
!      x(1:n): Solution
!      f(1:n): Residuals
!      info: 0: Bad input
!            1: Test 1 ok, i.e. max|Fi|<= FTOL
!            2: Test 2 ok, i.e. max|Xi-Xi'|<= XTOL * max|Xi|
!            3: Tests 1 and 2 ok
!            4: max. number of iterations reached or exceeded
!            5: Jacobinan matrix is singular
!            6: Not making good progress
!            7: Iterations are diverging
!            8: Convergence, but problems with too small a xtol, or Jacobian
!               nearly singular or variable badly scaled

      double precision ftol,xtol,x(*),f(*),w(*)
      integer n,maxf,iprt,info

      external sub

      call dsnleq(n,x,f,ftol,xtol,maxf,iprt,info,sub,w) !CERN C201

      return
      end
