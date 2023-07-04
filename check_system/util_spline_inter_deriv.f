*CMZ : 00.00/07 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/02 27/06/2003  15.56.18  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_INTER_DERIV(XA,YA,Y2A,N,X,Y,YP,MODE)

C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       YA2:  ARRAY SPLINE COEFFICIENTS
C-       X: Y(X) IS CALCULATED
C-       MODE: CONTROL FLAG:
C-             MODE.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             MODE.LT.0: NEW INITIALIZATION

C--   OUTPUT:

C-       Y:  Y(X) IS CALCULATED
C-      YP: DY/DX IS CALCULATED

      IMPLICIT NONE

      INTEGER NOLD,N,KLO,KHI,KLOLD,K,MODE

      REAL*8 Y,X,XA1OLD,XANOLD,H,A,B,YP,H26
      REAL*8 A2,A3,B2,B3,XL,YL,XH,YH,Y2L,Y2H

      REAL*8 XA(N),YA(N),Y2A(N)

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      IF(     XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &    .OR.
     &    XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1))) THEN
        WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
        WRITE(6,*)'X:',X
        STOP '***SR UTIL_SPLINE_INTER_DERIV: X OUT OF RANGE ***'
      ENDIF

      IF (MODE.LT.0.OR.KLOLD.GE.N) THEN
        KLO=1
      ELSE IF(NOLD.EQ.N
     &    .AND. XA(1).EQ.XA1OLD
     &    .AND. XA(N).EQ.XANOLD
     &    .AND. X.GT.XA(KLOLD)
     &    ) THEN
        KLO=KLOLD
      ELSE
        KLO=1
      ENDIF


      IF (X.LT.XA(KLO+1)) THEN
      KHI=KLO+1
      GOTO 2
      ENDIF

      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

2     XL=XA(KLO)
      XH=XA(KHI)
      YL=YA(KLO)
      YH=YA(KHI)

      H=XH-XL

      IF (H.EQ.0.D0) STOP '*** ERROR SR UTIL_SPLINE_INTER_DERIV: BAD INPUT ***'

      A=(XH-X)/H
      B=(X-XL)/H
      A2=A*A
      A3=A2*A
      B2=B*B
      B3=B2*B
      Y2L=Y2A(KLO)
      Y2H=Y2A(KHI)
      H26=H*H/6.D0

      Y=A*YL+B*YH + (A*(A-1.D0)*(A+1.D0)*Y2L + B*(B-1.D0)*(B+1.D0)*Y2H)*H26

      YP=(YH-YL)/H
     & +((3.D0*B2-1.D0)*Y2H-(3.D0*A2-1.D0)*Y2L)/6.D0*H

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)

      RETURN
      END
