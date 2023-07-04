*CMZ :  3.06/00 18/02/2019  18.55.48  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  12.31.20  by  Michael Scheer
*CMZ :  2.63/03 02/05/2008  14.41.00  by  Michael Scheer
*CMZ :  2.52/04 12/07/2004  16.15.58  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.24.55  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_FOLD_FUNCTION_GAUSS(NF,XF,F,SIGMA,RNSIGMA,FG,
     &  COEF,WS1,WS2,WS3,WS4)
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
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

      IMPLICIT NONE

      EXTERNAL FUNCTION DERF

      INTEGER NF,IH,IL,I

      REAL*8 XF(NF),F(NF),SIGMA,X0,FG(NF),RNSIGMA
      REAL*8 COEF(NF)
      REAL*8 WS1(NF),WS2(NF),WS3(NF),WS4(NF)

      REAL*8 CH,CL,CH2,CL2,CHCL,CH2CL,CHCL2,XL,XH,YL,YH,H,H61,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,EPS

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

C- CHECK ASCENDING ORDER

      DO I=2,NF
        IF (XF(I).LE.XF(I-1))
     &    STOP '*** ERROR SR UTIL_FOLD_FUNCTION_GAUSS:
     &    ARRAY XF NOT IN ASCENDING ORDER ***'
      ENDDO

      EPS=(XF(NF)-XF(1))*1.0D-10

C- SPLINES OF FUNCTION F

      CALL UTIL_SPLINE_COEF(XF,F,NF,-9999.0d0,-9999.0d0,COEF,WS1,WS2,WS3,WS4)

      DO IL=1,NF-1

        IH=IL+1

        XL=XF(IL)
        XH=XF(IH)
        YL=F(IL)
        YH=F(IH)
        CL=COEF(IL)
        CH=COEF(IH)

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
     &      '*** ERROR SR UTIL_FOLD_GAUSS:'
          PRINT*,
     &      '*** ARRAY XF NOT IN ASCENDING ORDER'
          STOP
        ENDIF

        H61=1.0D0/(6.0D0*H)

        WS1(IL)=((CHCL2*XH-CH2CL*XL)*XHXL+6.0D0*(XH*YL-XL*YH))*H61
        WS2(IL)=((CH2-CL2)*XHXL-CHCL2*XH2+CH2CL*XL2+6.0D0*(YH-YL))*H61
        WS3(IL)=(-CH*XL+CL*XH)/(2.0D0*H)
        WS4(IL)=CHCL*H61

        FG(IL)=0.0D0

      ENDDO !NF-1

      FG(NF)=0.0D0

      SN=RNSIGMA*SIGMA
      S2=SIGMA*SIGMA
      S22=2.0D0*S2
      SNR21=1.0D0/(ROOT2*SIGMA)
      SR2PI1=R2PI1/SIGMA

      DO I=1,NF

        X0=XF(I)

        X02=X0*X0
        X02S2=S2+X02
        X023S2=S22+X02S2

        IF (X0-SN.GE.XF(1)-EPS.AND.X0+SN.LE.XF(NF)+EPS) THEN

C UPPER BRANCH

          DO IL=I,NF-1

            IH=IL+1

            XL=XF(IL)
            XH=XF(IH)

            IF (XL-X0.LE.SN) THEN

              IF (XH-X0.GT.SN) XH=X0+SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XH+X0)+WS4(IL)*(S22+XH**2+XH*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XL+X0)+WS4(IL)*(S22+XL**2+XL*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              FG(I)=FG(I)+FGH-FGL

            ELSE
              GOTO 81
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IL

 81       CONTINUE

C LOWER BRANCH

          DO IH=I,2,-1

            IL=IH-1

            XL=XF(IL)
            XH=XF(IH)

            IF (X0-XH.LE.SN) THEN

              IF (X0-XL.GT.SN) XL=X0-SN

              DX=XH-X0

              FGH=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XH+X0)+WS4(IL)*(S22+XH**2+XH*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL)+WS3(IL)*(XL+X0)+WS4(IL)*(S22+XL**2+XL*X0+X02))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)+WS3(IL)*X02S2
     &          +X0*(WS2(IL)+WS4(IL)*X023S2))*
     &          DERF(DX*SNR21))

              FG(I)=FG(I)+FGH-FGL

            ELSE
              GOTO 82
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IH

 82       CONTINUE

        ELSE IF (X0+SN.GT.XF(NF)) THEN

          GOTO 88

        ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

      ENDDO !NF

 88   CONTINUE

      RETURN
      END
