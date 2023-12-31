*CMZ :          25/08/2014  16.15.57  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  14.20.46  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
*
* $Id: rfstft.F,v 1.1 1996/04/17 12:32:04 mclareni Exp $
*
* $Log: rfstft.F,v $
* Revision 1.1  1996/04/17 12:32:04  mclareni
* Add d/rfstft.F (D705) and to Imakefile. cfstft.F becomes D706.
* In tests, add d705m.F for rfstft and d706m.F for cfstft and the corresponding
* additions to main.F and Imakefile.
*
*
      SUBROUTINE RFSTFT(MS,A)

      COMPLEX A(0:*),T,T1,T2,U,W

cmsh      PARAMETER (PI = 3.14159 26535 89793D0)
      PARAMETER (PI = 3.141592653589793)

      IF(MS .EQ. 0) THEN
       A(0)=REAL(A(0))
       RETURN
      ENDIF
      M=ABS(MS)-1
      N=2**M
      U=(0.,1.)
      IF(MS .LT. 0) THEN
       CALL CFSTFT(-M,A)
       F=0.25/N
       DO 1 I = 0,N-1
    1  A(I)=F*A(I)
       A(N)=A(0)
       U=CONJG(U)
      ENDIF
cmsh    2 PHI=PI/SIGN(N,MS)
      PHI=PI/SIGN(N,MS)
      W=CMPLX(COS(PHI),SIN(PHI))
      DO 3 J = 0,N/2
      T=CONJG(A(N-J))
      T1=A(J)+T
      T2=(A(J)-T)*U
      A(J)=T1+T2
      A(N-J)=CONJG(T1-T2)
    3 U=U*W
      IF(MS .GT. 0) CALL CFSTFT(M,A)
      RETURN
      END
