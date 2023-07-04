*CMZ :  4.00/07 10/07/2020  08.22.44  by  Michael Scheer
*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/07 07/05/2008  13.09.13  by  Michael Scheer
*CMZ : 00.00/02 27/06/2003  09.36.13  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_INTER_HUNT(XA,YA,Y2A,N,X,Y)

C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       YA2:  ARRAY SPLINE COEFFICIENTS
C-       X: Y(X) IS CALCULATED

C--   OUTPUT:

C-       Y: Y(X) IS CALCULATED

      IMPLICIT NONE

      INTEGER N,KLO,KHI,KLOLD,K,KD

      REAL*8 Y,X,H,A,B

      REAL*8 XA(N),YA(N),Y2A(N)

      DATA KLOLD/1/

      save

      IF(     XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &    .OR.
     &    XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1))) THEN
        WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
        WRITE(6,*)'X:',X
        WRITE(6,*)'X(1),X(N):',XA(1),XA(N)
        STOP '*** ERROR IN UTIL_SPLINE_INTER_HUNT: X OUT OF RANGE ***'
      ENDIF

      IF (KLOLD.GT.N) THEN
        KLO=1
      ELSE
        KLO=KLOLD
      ENDIF

      IF (X.GE.XA(KLO)) THEN

C HUNT UP

        KD=1
1       KHI=MIN(KLO+KD,N)
        IF (X.GT.XA(KHI)) THEN
          KD=2*KD
          KLO=KHI
          GOTO 1
        ENDIF


      ELSE  !(X.GE.XA(KLO))

C HUNT DOWN
        KD=1
        KHI=KLO
2       KLO=MAX(KHI-KD,1)
        IF (X.LT.XA(KLO)) THEN
          KD=2*KD
          KHI=KLO
          GOTO 2
        ENDIF

      ENDIF

3     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 3
      ENDIF

      H=XA(KHI)-XA(KLO)

      IF (H.le.0.D0) STOP '*** ERROR SR UTIL_SPLINE_INTER: BAD INPUT ***'

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H

      Y=A*YA(KLO)+B*YA(KHI)+
     &  (A*(A*A-1.0d0)*Y2A(KLO)+B*(B*B-1.0d0)*Y2A(KHI))*(H**2)/6.D0

      KLOLD=KLO

      RETURN
      END
