*CMZ :  3.07/00 08/03/2019  19.50.00  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.56.50  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.67/05 16/05/2012  12.45.37  by  Michael Scheer
*CMZ : 00.00/07 12/10/2009  12.17.45  by  Michael Scheer
*CMZ : 00.00/02 14/04/2003  12.46.09  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.48  by  Michael Scheer
*-- Author : Michael Scheer
c      SUBROUTINE UTIL_SPLINE_COEF_omp(X,Y,N,YP1,YPN,Y2,AA,BB,CC,C)
      SUBROUTINE UTIL_SPLINE_COEF_omp(N,YP1,YPN)
*KEEP,gplhint.
*KEND.

C--- CALCULATES SPLINE COEFFICIENTS

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       X: ARRAY OF X-VALUES
C-       Y: ARRAY OF Y-VALUES
C-       YP1:  SECOND DERIVATIVE AT FIRST X-VALUE
C-       YPN:  SECOND DERIVATIVE AT LAST X-VALUE

C--   OUPUT:

C-       Y2:   SPLINE-COEFFICIENTS

C--   WORKINGSPACE: AA(N),BB(N),CC(N),C(N)

      use wobsvmod

      IMPLICIT NONE

      INTEGER N,J
c      REAL*8  X(N),Y(N),y2(N),AA(N),BB(N),CC(N),C(N)

      REAL*8 YP1,YPN

      double precision xx(3),yy(3),a(3),yp(3),xopt,yopt
      INTEGER ifail

      IF (N.LT.3) then
        if (abs(yp1).eq.9999.0d0) then
          wobsv3_th(1)=0.0d0
        else
          wobsv3_th(1)=yp1
        endif
        if (abs(ypn).eq.9999.0d0) then
          wobsv3_th(n)=0.0d0
        else
          wobsv3_th(n)=ypn
        endif
        RETURN
      endif

c      print*,"N:",n

      if (abs(yp1).eq.9999.0d0) then
        xx=x_th(1:3)
        yy=wobsv1_th(1:3)
        call UTIL_parabel_omp(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          wobsv3_th(1)=2.0d0*a(3)
        else
          wobsv3_th(1)=0.0d0
        endif
      else
        wobsv3_th(1)=YP1
      endif

      if (abs(ypn).eq.9999.0d0) then
        xx=x_th(n-2:n)
        yy=wobsv1_th(n-2:n)
        call UTIL_parabel_omp(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          wobsv3_th(n)=2.0d0*a(3)
        else
          wobsv3_th(N)=0.0d0
        endif
      else
        wobsv3_th(N)=YPN
      endif

      wobsv7_th(1)=wobsv3_th(1)
      wobsv7_th(N)=wobsv3_th(n)

      wobsv5_th(1)=1.D0
      wobsv6_th(1)=0.D0
      wobsv6_th(N)=1.D0

      DO J=2,N-1
        if(x_th(j+1).eq.x_th(j)) then
          write(6,*)
          write(6,*)
     &      '*** Error in util_spline_coef_omp: Intervall of zero length'
          write(6,*)'j, x_th(j), x_th(j+1):',j,x_th(j),x_th(j+1)
          write(6,*)
          stop
        endif
        wobsv4_th(J)=(x_th(J  )-x_th(J-1))/6.D0
        wobsv5_th(J)=(x_th(J+1)-x_th(J-1))/3.D0
        wobsv6_th(J)=(x_th(J+1)-x_th(J  ))/6.D0
        wobsv7_th(J)=(wobsv1_th(J+1)-wobsv1_th(J  ))/(x_th(J+1)-x_th(J  ))
     &    -(wobsv1_th(J  )-wobsv1_th(J-1))/(x_th(J  )-x_th(J-1))
      ENDDO !J

      DO J=2,N-1

        wobsv5_th(J)=wobsv5_th(J)-wobsv4_th(J)*wobsv6_th(J-1)
        wobsv7_th(J)= wobsv7_th(J)-wobsv4_th(J)* wobsv7_th(J-1)
C          wobsv4_th(J)=wobsv4_th(J)-wobsv4_th(J)*wobsv5_th(J-1)

        wobsv6_th(J)=wobsv6_th(J)/wobsv5_th(J)
        wobsv7_th(J)= wobsv7_th(J)/wobsv5_th(J)
        wobsv5_th(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
         wobsv3_th(J)=wobsv7_th(J)-wobsv6_th(J)*wobsv3_th(J+1)
      ENDDO

      RETURN
      END
