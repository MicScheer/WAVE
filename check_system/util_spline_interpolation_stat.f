*CMZ : 00.00/15 12/10/2013  12.17.18  by  Michael Scheer
*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  13.46.51  by  Michael Scheer
*CMZ : 00.00/06 06/07/2007  16.55.15  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.16.22  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_INTERPOLATION_STAT(N,XA,YA,X,Y,Y2A,W1,W2,W3,W4,MODE,ISTAT)

C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       X: Y(X) IS CALCULATED
C-       MODE: CONTROL FLAG:
C-             MODE.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             MODE.LT.0: NEW INITIALIZATION

C--   OUTPUT:

C-       Y: Y(X) IS CALCULATED

      IMPLICIT NONE

      double precision Y,X,XA1OLD,XANOLD,H,A,B
      double precision XA(N),YA(N),Y2A(N),
     &  W1(N),W2(N),W3(N),W4(N)

      INTEGER NOLD,N,KLO,KHI,KLOLD,K,MODE,ISTAT

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      istat=0

      IF(N.lt.3) THEN
        print*, '*** ERROR IN UTIL_SPLINE_INTERPOLATION_STAT: LESS THAN 3 POINTS GIVEN***'
        istat=-1
        return
      ENDIF

      IF(XA(1).GE.XA(N)) then
        print*, '*** ERROR IN UTIL_SPLINE_INTERPOLATION_STAT: XA(1) .GE. XA(N) ***'
        istat=-1
        return
      ENDIF

      IF(X.LT.XA(1)) then
        WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
        WRITE(6,*)'X:',X
        print*, '*** ERROR IN UTIL_SPLINE_INTER: X OUT OF RANGE ***'
        istat=-1
        return
      ENDIF

      IF(X.GT.XA(N)) then
        WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
        WRITE(6,*)'X:',X
        print*, '*** ERROR IN UTIL_SPLINE_INTER: X OUT OF RANGE ***'
        istat=-1
        return
      ENDIF

      IF (MODE.LT.0.OR.KLOLD.GE.N) THEN
        CALL UTIL_SPLINE_COEF(XA,YA,N,9999.0d0,9999.0d0,Y2A,W1,W2,W3,W4)
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

2     H=XA(KHI)-XA(KLO)

      IF (H.le.0.D0) THEN
        PRINT*,'*** ERROR IN UTIL_SPLINE_INTER: BAD INPUT ***'
        istat=-1
        return
      ENDIF

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     &  (A*(A+1.D0)*(A-1.D0)*Y2A(KLO)+B*(B+1.D0)*
     &  (B-1.D0)*Y2A(KHI))*(H**2)/6.D0

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)

      RETURN
      END
