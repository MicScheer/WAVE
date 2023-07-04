*CMZ : 00.00/15 12/10/2013  12.22.24  by  Michael Scheer
*CMZ : 00.00/07 02/05/2008  13.10.35  by  Michael Scheer
*CMZ : 00.00/02 21/06/2004  11.46.09  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_FOLD_FOURIER
     &  (X,Y,NX,NFOLD,AK,NKX,YFOLD,COEF,WORK1,WORK2,WORK3,WORK4,IFAIL)

C---  CALCULATES FOLDING SUM(INTEGRAL(F(x)*a(k)*COS((k-1)*kx0*(x-x0)),dx0),k)
C---  ATTENTION: AK(1) = a(k=0)
C---  X IS AN EQUALLY SPACED GRID
C---  2*NFOLD DEFINES WIDTH OF FOLDING FUNCTION

      IMPLICIT NONE

      INTEGER NX,NX1,NX2,IFAIL,NKX,NKXO,NFOLDO,IAKO

      INTEGER K,I,NFOLD,IX0,J,NLOCALP
      PARAMETER (NLOCALP=1024)

      INTEGER IH,IL,KL,KH,KLA,KHA,K1L,K1H,K1LA,K1HA

      REAL*8 X(NX),Y(NX),YFOLD(NX),ADDCOS,AK(NKX),AK1H24
      REAL*8 COEF(NX),WORK1(NX),WORK2(NX),WORK3(NX),WORK4(NX)

      REAL*8 XL,XH,YL,YH,YPL,YPH,HXK0,CS(0:NLOCALP),SN(0:NLOCALP)
      REAL*8 H,XL0,XK0,PI,X0,XKK,XKK2,XKK4,AKK(NLOCALP),HO,XK0O,AKO(NLOCALP)
      REAL*8 CSKIXL, CSKIXH,CS1X0L,CS1X0H
      REAL*8 SN1X0L,SN1X0H

      REAL*8 HK3XK03(NLOCALP),H2K2XK02(NLOCALP),HKXK0(NLOCALP),H2
     &        ,XK02K2(NLOCALP),YPL3,YPL6,YPH3,YPH6

      DATA PI/3.141592653589793D0/
      DATA NKXO/0/,HO/0.D0/,XK0O/0.D0/,NFOLDO/0/,AKK/NLOCALP*0.D0/

      IFAIL=0

        IF (NFOLD*NKX.GE.NLOCALP) THEN
          WRITE(6,*)'ERROR IN UTIL_FOLD_FOURIER: DIMENSION NLOCALP EXCEEDED'
          IFAIL=1
          RETURN
        ENDIF

        IF (NFOLD.LT.1) THEN
          WRITE(6,*)
     &    'ERROR IN UTIL_FOLD_FOURIER: NUMBER OF MASHES FOR FOLDING FUNCTION'
          WRITE(6,*)
     &    'LOWER THAN ONE'
          IFAIL=1
          RETURN
        ENDIF

        IF (NFOLD.GT.NX/2) THEN
          WRITE(6,*)
     &    'ERROR IN UTIL_FOLD_FOURIER: NUMBER OF MASHES FOR FOLDING FUNCTION'
          WRITE(6,*)
     &    'TOO LARGE'
          IFAIL=1
          RETURN
        ENDIF

        CALL UTIL_SPLINE_COEF(X,Y,NX,9999.0d0,9999.0d0,COEF,WORK1,WORK2,WORK3,WORK4)

C--- ADJUST TO MASH

        H=x(2)-x(1)
        H2=H*H
        XK0=PI/(NFOLD*ABS(H))
        XL0=2.D0*PI/XK0
        HXK0=H*XK0

      IAKO=0
      DO K=1,NKX
          IF (AK(K).NE.AKO(K)) THEN
              IAKO=1
              GOTO 10
          ENDIF
      ENDDO

10    IF (IAKO.NE.0.OR.NFOLDO.NE.NFOLD.OR.NKXO.NE.NKX.OR.
     &      HO.NE.H.OR.XK0.NE.XK0O) THEN

              XKK=0.D0
         DO K=1,NKX-1
                  XKK=XKK+XK0
                  XKK2=XKK*XKK
                  XKK4=XKK2*XKK2
                  AKK(K+1)=AK(K+1)/(h*XKK4)
              ENDDO
              AK1H24=AK(1)*H/24.D0

              XK02K2(1)=XK0
              HKXK0(1)=HXK0

              DO K=2,NFOLD*NKX
                  XK02K2(K)=XK02K2(K-1)+XK0
                  HKXK0(K)=HKXK0(K-1)+HXK0
         ENDDO

              DO K=1,NFOLD*NKX
                  XK02K2(K)=XK02K2(K)*XK02K2(K)
                  H2K2XK02(K)=HKXK0(K)*HKXK0(K)
         ENDDO

              DO K=1,NFOLD*NKX
                  HK3XK03(K)=H2K2XK02(K)*HKXK0(K)/H2
         ENDDO

              CS(0)=1.D0
              SN(0)=0.D0
              CS(1)=COS(HKXK0(1))
              SN(1)=SIN(HKXK0(1))
              DO K=2,NFOLD*NKX
                  CS(K)=CS(K-1)*CS(1)-SN(K-1)*SN(1)
                  SN(K)=CS(K-1)*SN(1)+SN(K-1)*CS(1)
         ENDDO

      ENDIF   !OLD-VALUES

        DO I=1,NX
        YFOLD(I)=0.0D0
      ENDDO

      NX1=1+NFOLD
      NX2=NX-NFOLD

        DO IX0=NX1,NX2 !LOOP OVER FOLDING REGION

          DO J=1,2*NFOLD   !INTEGRATE FULL INTERVALS

            I=IX0-NFOLD+J-1

          IL=I
            IH=I+1

          X0=X(IX0)
          XH=X(IH)
          XL=X(IL)
          YH=Y(IH)
          YL=Y(IL)
          YPL=COEF(IL)
          YPH=COEF(IH)
          YPL3=YPL/3.D0
          YPL6=YPL/6.D0
          YPH3=YPH/3.D0
          YPH6=YPH/6.D0

              YFOLD(IX0)=YFOLD(IX0)
     &       +aK1H24*(-H2*(yph+ypl)+12.0*(yh+yl))

              KL=0
              KH=0
         DO K=1,NKX-1

              KL=KL+IX0-IL
              KH=KH+IX0-IH
              K1L=K-KL
              K1H=K+KH
              KLA=ABS(KL)
              KHA=ABS(KH)
              K1HA=ABS(K1H)
              K1LA=ABS(K1L)

              CSKIXL=CS(KLA)
              CSKIXH=CS(KHA)
              CS1X0L=CS(K1LA)
              CS1X0H=CS(K1HA)

              IF (K1L.LT.0) THEN
                 SN1X0L=SN(K1LA)
              ELSE
                 SN1X0L=-SN(K1LA)
              ENDIF

              IF (K1H.LT.0) THEN
                 SN1X0H=SN(K1HA)
              ELSE
                 SN1X0H=-SN(K1HA)
              ENDIF

      addcos=akk(K+1)*(
     & CS1X0L*(H2K2XK02(K)*yph3+XK02K2(K)*yh-yph)
     &+CS1X0H*(H2K2XK02(K)*ypl3+XK02K2(K)*yl-ypl)
     &+CSKIXH*(H2K2XK02(K)*ypl6-XK02K2(K)*yl+ypl)
     &+CSKIXL*(H2K2XK02(K)*yph6-XK02K2(K)*yh+yph)
     &+SN1X0L*(HK3XK03(K)*yh-HKXK0(K)*yph)+SN1X0H*(HK3XK03(K)*yl-HKXK0(K)*ypl))

                  YFOLD(IX0)=YFOLD(IX0)+ADDCOS

              ENDDO   !K

          ENDDO !J
          ENDDO !IX0

      NKXO=NKX
      NFOLDO=NFOLD
      HO=H
      XK0O=XK0

      DO K=1,NKX
          AKO(K)=AK(K)
      ENDDO

      RETURN
      END
