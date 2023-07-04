*CMZ :  3.05/13 30/08/2018  16.19.30  by  Michael Scheer
*CMZ : 00.00/02 27/06/2005  19.03.55  by  Michael Scheer
*CMZ : 00.00/01 05/06/96  16.08.05  by  Michael Scheer
*-- Author :    Michael Scheer   03/06/96
      SUBROUTINE UTIL_DETERMINANTE_4(A,DET,ifail)

C CALCULATES DETERMINANT OF MATRIX A(3,3)

      IMPLICIT NONE

      DOUBLE PRECISION A(4,4),DET
      integer ifail

      ifail=0

      det=
     &   a(1,1)*a(2,2)*a(3,3)*a(4,4)
     &  -a(1,1)*a(2,2)*a(3,4)*a(4,3)
     &  -a(1,1)*a(2,3)*a(3,2)*a(4,4)
     &  +a(1,1)*a(2,3)*a(3,4)*a(4,2)
     &  +a(1,1)*a(2,4)*a(3,2)*a(4,3)
     &  -a(1,1)*a(2,4)*a(3,3)*a(4,2)
     &  -a(1,2)*a(2,1)*a(3,3)*a(4,4)
     &  +a(1,2)*a(2,1)*a(3,4)*a(4,3)
     &  +a(1,2)*a(2,3)*a(3,1)*a(4,4)
     &  -a(1,2)*a(2,3)*a(3,4)*a(4,1)
     &  -a(1,2)*a(2,4)*a(3,1)*a(4,3)
     &  +a(1,2)*a(2,4)*a(3,3)*a(4,1)
     &  +a(1,3)*a(2,1)*a(3,2)*a(4,4)
     &  -a(1,3)*a(2,1)*a(3,4)*a(4,2)
     &  -a(1,3)*a(2,2)*a(3,1)*a(4,4)
     &  +a(1,3)*a(2,2)*a(3,4)*a(4,1)
     &  +a(1,3)*a(2,4)*a(3,1)*a(4,2)
     &  -a(1,3)*a(2,4)*a(3,2)*a(4,1)
     &  -a(1,4)*a(2,1)*a(3,2)*a(4,3)
     &  +a(1,4)*a(2,1)*a(3,3)*a(4,2)
     &  +a(1,4)*a(2,2)*a(3,1)*a(4,3)
     &  -a(1,4)*a(2,2)*a(3,3)*a(4,1)
     &  -a(1,4)*a(2,3)*a(3,1)*a(4,2)
     &  +a(1,4)*a(2,3)*a(3,2)*a(4,1)

      if (det.ne.det) then
        ifail=1
      else if (1.0d0/det.eq.0.0d0) then
        ifail=2
      endif

      RETURN
      END
