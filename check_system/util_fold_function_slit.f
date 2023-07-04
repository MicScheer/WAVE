*CMZ :          23/04/2023  16.09.55  by  Michael Scheer
*CMZ : 00.00/15 12/10/2013  12.22.24  by  Michael Scheer
*CMZ : 00.00/07 02/05/2008  13.10.35  by  Michael Scheer
*CMZ : 00.00/02 12/07/2004  16.15.58  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.24.55  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE util_fold_function_slit(n,x,F,S,FG)

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG(X)=INT{F(x)*C(x-X),Dx}

C--   INPUT:

C-       n:   NUMBER OF x,F-VALUES
C-       x:   ARRAY OF X-VALUES (MUST BE IN ASCENDING ORDER)
C-       F: ARRAY OF FUNCTION-VALUES
c-       S: Width of slit

C--   OUTPUT:

C-       FG:   FG(x(i)) IS CALCULATED points kl.le.i.le.kh

      IMPLICIT NONE

      INTEGER n,I,n2,ifail,modeini

      double precision, dimension(:), allocatable :: coef,w1,w2,w3,w4,ff,xg
      REAL*8 x(n),F(n),S,FG(n),s2,x3(3),y3(3),a(3),yp(3),xopt,yopt,dx

C- CHECK ASCENDING ORDER

      fg=f*0.0d0
      if (n.lt.3) return

      if (s.eq.0.0d0) then
        fg=f
        return
      endif

      DO I=2,n
        IF (x(I).LE.x(I-1))
     &    STOP '*** ERROR SR util_fold_function_slit:
     &    ARRAY x NOT IN ASCENDING ORDER ***'
      ENDDO

      n2=n+2
      allocate(ff(n2),coef(n2),xg(n2),w1(n2),w2(n2),w3(n2),w4(n2))

      s2=s/2.0d0

      x3=x(1:3)
      y3=f(1:3)
      dx=x(2)-x(1)
      xg(1)=x(1)-s2
      call util_parabel(x3,y3,a,yp,xopt,yopt,IFAIL)
      ff(1)=a(1)+a(2)*xg(1)+a(3)*xg(1)**2

      x3=x(n-2:n)
      y3=f(n-2:n)
      dx=x(n)-x(n-1)
      xg(n2)=x(n)+s2
      call util_parabel(x3,y3,a,yp,xopt,yopt,IFAIL)
      ff(n2)=a(1)+a(2)*xg(n2)+a(3)*xg(n2)**2

C- SPLINES OF FUNCTION F

      xg(2:n)=x(1:n)
      ff(2:n)=f(1:n)

      do i=1,n
        if (i.eq.1) then
          modeini=-1
        else
          modeini=0
        endif
        call util_spline_integral_window(xg,ff,n2,x(i)-s2,x(i)+s2,fg(i)
     &    ,coef,w1,w2,w3,w4,modeini,ifail)
        fg(i)=fg(i)/s
      enddo

      deallocate(ff,coef,xg,w1,w2,w3,w4)

      return
      end
