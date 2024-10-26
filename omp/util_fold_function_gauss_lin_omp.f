*CMZ :  4.01/03 09/06/2023  14.38.30  by  Michael Scheer
*CMZ :  4.00/07 06/04/2020  08.51.14  by  Michael Scheer
*CMZ :  3.06/00 18/02/2019  19.28.56  by  Michael Scheer
*CMZ :  2.70/02 14/12/2012  10.35.01  by  Michael Scheer
*-- Author : Michael Scheer
      !SUBROUTINE util_fold_function_gauss_lin(NF,XF,F,SIGMA,DNSIGMA,
      !&  FG,WS1,WS2)
      SUBROUTINE util_fold_function_gauss_lin_omp(NF,SIGMA,DNSIGMA)
*KEEP,gplhint.
*KEND.

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG(X)=INT{F(XF)*G(XF-X),DXF}

      use wobsvmod

      IMPLICIT NONE

      EXTERNAL FUNCTION DERF

      INTEGER NF,IH,IL,I

      REAL*8 SIGMA,X0,DNSIGMA
      !REAL*8 XF(NF),F(NF),SIGMA,X0,FG(NF),DNSIGMA
      !REAL*8 WS1(NF),WS2(NF)

      REAL*8 XL,XH,YL,YH,H,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,EPS

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

C- CHECK ASCENDING ORDER

      DO I=2,NF
        IF (x_th(I).LE.x_th(I-1))
     &    STOP '*** ERROR SR UTIL_FOLD_FUNCTION_GAUSS_LIN:
     &    ARRAY XF NOT IN ASCENDING ORDER ***'
      ENDDO

      EPS=(x_th(NF)-x_th(1))*1.0D-10

      DO IL=1,NF-1

        IH=IL+1

        XL=x_th(IL)
        XH=x_th(IH)
        YL=wobsv1_th(IL)
        YH=wobsv1_th(IH)

        XHXL=XH*XL
        XH2=XH*XH
        XL2=XL*XL

        H=XH-XL

        IF (H.LE.0.0D0) THEN
          PRINT*,
     &      '*** ERROR SR UTIL_FOLD_GAUSS_LIN:'
          PRINT*,
     &      '*** ARRAY XF NOT IN ASCENDING ORDER'
          STOP
        ENDIF

        wobsv3_th(IL)=(XH*YL-XL*YH)/h
        wobsv4_th(IL)=(YH-YL)/h

        wobsv2_th(IL)=0.0D0

      ENDDO !NF-1

      wobsv2_th(NF)=0.0D0

      SN=DNSIGMA*SIGMA
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
     &          SR2PI1*(-EXP(-DX**2/S22)*S2*wobsv4_th(IL)
     &          +SQPI2*SIGMA*(wobsv3_th(IL)+X0*wobsv4_th(IL))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(-EXP(-DX**2/S22)*S2*wobsv4_th(IL)
     &          +SQPI2*SIGMA*(wobsv3_th(IL)+X0*wobsv4_th(IL))*
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
     &          wobsv4_th(IL))
     &          +SQPI2*SIGMA*(
     &          wobsv3_th(IL)
     &          +X0*(wobsv4_th(IL)))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          wobsv4_th(IL))
     &          +SQPI2*SIGMA*(
     &          wobsv3_th(IL)
     &          +X0*(wobsv4_th(IL)))*
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
