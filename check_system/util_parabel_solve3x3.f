*CMZ :          31/10/2022  17.16.03  by  Michael Scheer
*CMZ : 00.00/16 19/03/2014  12.30.26  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.27.13  by  Michael Scheer
*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_PARABEL_solve3x3(Xin,Yin,A,YP,XOPT,yopt,IFAIL)

C--- CALCULATES A(1),A(2),A(3), THE DERIVATIVES YP(X(1)),YP(X(2)),YP(X(3)),
C    AND THE EXTREMUM (XOPT,A(XOPT)) OF PARABOLA Y=A1+A2*X+A3*X**2
C    FROM COORDINATES OF THE THREE POINTS (X(1),Y(1)),(X(2),Y(2)),(X(3),Y(3))
C

      IMPLICIT NONE

      INTEGER IFAIL,i

      REAL*8 a(3),a33(3,3),xin(3),yin(3),yp(3),xopt,yopt,x(3),y(3)

      IFAIL=0

      x=xin
      y=yin
      call util_sort_func(3,x,y)

      do i=1,3
        a33(i,1)=1.0d0
        a33(i,2)=x(i)
        a33(i,3)=x(i)**2
        a(i)=y(i)
        !print*,a33(i,:),a(i)
        !print*
      enddo

      call util_solve_3x3(a33,a,ifail)

      if (ifail.ne.0) return

      yp(:)=a(2)+2.0d0*a(3)*x(:)

      if (a(3).ne.0.0d0) then
        xopt=-a(2)/2.0d0/a(3)
        yopt=a(1)+a(2)*xopt+a(3)*xopt**2
      endif

      if (xopt.lt.x(1).or.xopt.gt.x(3)) ifail=2

      RETURN
      END
