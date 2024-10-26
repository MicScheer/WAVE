*CMZ :  3.07/00 09/03/2019  12.05.11  by  Michael Scheer
*CMZ :  3.06/00 18/02/2019  18.55.48  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  12.31.20  by  Michael Scheer
*CMZ :  2.63/03 02/05/2008  14.41.00  by  Michael Scheer
*CMZ :  2.52/04 12/07/2004  16.15.58  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.24.55  by  Michael Scheer
*-- Author : Michael Scheer
c      SUBROUTINE UTIL_FOLD_FUNCTION_GAUSS_omp(NF,XF,F,SIGMA,RNSIGMA,FG,
c     &  COEF,WS1,WS2,WS3,WS4)
      SUBROUTINE UTIL_FOLD_FUNCTION_GAUSS_omp(NF,SIGMA,RNSIGMA)
*KEEP,gplhint.
*KEND.

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG(X)=INT{F(XF)*G(XF-X),DXF}

C--   INPUT:

C-       NF:   NUMBER OF XF,F-VALUES
C-       XF:   ARRAY OF X-VALUES (MUST BE IN ASCENDING ORDER)
C-       F: ARRAY OF FUNCTION-VALUES
C-       SIGMA:  SIGMA OF GAUSSIAN
C-       RNSIGMA: NUMBER OF SIGMAS TO BE CONSIDERED

C--   OUTPUT:

C-       FG:   FG(X0) IS CALCULATED
      use wobsvmod

      IMPLICIT NONE

      EXTERNAL FUNCTION DERF

      INTEGER NF,IH,IL,I

      REAL*8 SIGMA,X0,RNSIGMA
c      REAL*8 XF(NF),F(NF),SIGMA,X0,FG(NF),RNSIGMA
c      REAL*8 COEF(NF)
c      REAL*8 WS1(NF),WS2(NF),WS3(NF),WS4(NF)

      REAL*8 CH,CL,CH2,CL2,CHCL,CH2CL,CHCL2,XL,XH,YL,YH,H,H61,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,EPS

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

C- CHECK ASCENDING ORDER

      DO I=2,NF
        IF (x_th(I).LE.x_th(I-1)) then
          print*,""
          print*,'*** ERROR SR UTIL_FOLD_FUNCTION_GAUSS_omp:'
          print*,'ARRAY XF NOT IN ASCENDING ORDER ***'
          print*,""
          stop
        endif
      ENDDO

      EPS=(x_th(NF)-x_th(1))*1.0D-10

C- SPLINES OF FUNCTION F

      CALL UTIL_SPLINE_COEF_omp(NF,-9999.0d0,-9999.0d0)

      DO IL=1,NF-1

        IH=IL+1

        XL=x_th(IL)
        XH=x_th(IH)
        YL=wobsv1_th(IL)
        YH=wobsv1_th(IH)
        CL=wobsv3_th(IL)
        CH=wobsv3_th(IH)

        CL2=2.0D0*CL
        CH2=2.0D0*CH
        CHCL=CH-CL
        CHCL2=CH+CL2
        CH2CL=CH2+CL

        XHXL=XH*XL
        XH2=XH*XH
        XL2=XL*XL

        H=XH-XL

        IF (H.LE.0.0D0) THEN
          PRINT*,
     &      '*** ERROR SR UTIL_FOLD_GAUSS_omp:'
          PRINT*,
     &      '*** ARRAY XF NOT IN ASCENDING ORDER'
          STOP
        ENDIF

        H61=1.0D0/(6.0D0*H)

        wobsv4_th(IL)=((CHCL2*XH-CH2CL*XL)*XHXL+6.0D0*(XH*YL-XL*YH))*H61
        wobsv5_th(IL)=((CH2-CL2)*XHXL-CHCL2*XH2+CH2CL*XL2+6.0D0*(YH-YL))*H61
        wobsv6_th(IL)=(-CH*XL+CL*XH)/(2.0D0*H)
        wobsv7_th(IL)=CHCL*H61

        wobsv2_th(IL)=0.0D0

      ENDDO !NF-1

      wobsv2_th(NF)=0.0D0

      SN=RNSIGMA*SIGMA
      S2=SIGMA*SIGMA
      S22=2.0D0*S2
      SNR21=1.0D0/(ROOT2*SIGMA)
      SR2PI1=R2PI1/SIGMA

      DO I=1,NF

        X0=x_th(I)

        X02=X0*X0
        X02S2=S2+X02
        X023S2=S22+X02S2

        IF (X0-SN.GE.x_th(1)-EPS.AND.X0+SN.LE.x_th(NF)+EPS) THEN

C UPPER BRANCH

          DO IL=I,NF-1

            IH=IL+1

            XL=x_th(IL)
            XH=x_th(IH)

            IF (XL-X0.LE.SN) THEN

              IF (XH-X0.GT.SN) XH=X0+SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          wobsv5_th(IL)+wobsv6_th(IL)*(XH+X0)+wobsv7_th(IL)*(S22+XH**2+XH*X0+X02))
     &          +SQPI2*SIGMA*(
     &          wobsv4_th(IL)+wobsv6_th(IL)*X02S2
     &          +X0*(wobsv5_th(IL)+wobsv7_th(IL)*X023S2))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          wobsv5_th(IL)+wobsv6_th(IL)*(XL+X0)+wobsv7_th(IL)*(S22+XL**2+XL*X0+X02))
     &          +SQPI2*SIGMA*(
     &          wobsv4_th(IL)+wobsv6_th(IL)*X02S2
     &          +X0*(wobsv5_th(IL)+wobsv7_th(IL)*X023S2))*
     &          DERF(DX*SNR21))

              wobsv2_th(I)=wobsv2_th(I)+FGH-FGL

            ELSE
              GOTO 81
            ENDIF ! (X-SN.GE.x_th(1).AND.X+SN.LE.x_th(NF))

          ENDDO !IL

 81       CONTINUE

C LOWER BRANCH

          DO IH=I,2,-1

            IL=IH-1

            XL=x_th(IL)
            XH=x_th(IH)

            IF (X0-XH.LE.SN) THEN

              IF (X0-XL.GT.SN) XL=X0-SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          wobsv5_th(IL)+wobsv6_th(IL)*(XH+X0)+wobsv7_th(IL)*(S22+XH**2+XH*X0+X02))
     &          +SQPI2*SIGMA*(
     &          wobsv4_th(IL)+wobsv6_th(IL)*X02S2
     &          +X0*(wobsv5_th(IL)+wobsv7_th(IL)*X023S2))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          wobsv5_th(IL)+wobsv6_th(IL)*(XL+X0)+wobsv7_th(IL)*(S22+XL**2+XL*X0+X02))
     &          +SQPI2*SIGMA*(
     &          wobsv4_th(IL)+wobsv6_th(IL)*X02S2
     &          +X0*(wobsv5_th(IL)+wobsv7_th(IL)*X023S2))*
     &          DERF(DX*SNR21))

              wobsv2_th(I)=wobsv2_th(I)+FGH-FGL

            ELSE
              GOTO 82
            ENDIF ! (X-SN.GE.x_th(1).AND.X+SN.LE.x_th(NF))

          ENDDO !IH

 82       CONTINUE

        ELSE IF (X0+SN.GT.x_th(NF)) THEN

          GOTO 88

        ENDIF ! (X-SN.GE.x_th(1).AND.X+SN.LE.x_th(NF))

      ENDDO !NF

 88   CONTINUE

      RETURN
      END
