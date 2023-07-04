*CMZ :  4.00/07 06/04/2020  14.51.16  by  Michael Scheer
*CMZ :  3.06/00 18/02/2019  19.28.56  by  Michael Scheer
*CMZ :  2.70/02 14/12/2012  10.35.01  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE util_fold_gauss_lin(NF,XF,F,SIGMA,DNSIGMA,IX,FG,WS1,WS2)
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

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG=INT{F(XF)*G(XF-X),DXF}

C--   INPUT:

C-       NF:   NUMBER OF XF,F-VALUES
C-       XF:   ARRAY OF X-VALUES (MUST BE IN ASCENDING ORDER)
C-       F: ARRAY OF FUNCTION-VALUES
C-       SIGMA:  SIGMA OF GAUSSIAN
C-       DNSIGMA: NUMBER OF SIGMAS TO BE CONSIDERED
C-       RESULT FG REFERES TO XF(IX)

C--   OUTPUT:

C-       FG:   FG(X0) IS CALCULATED

      IMPLICIT NONE

      EXTERNAL FUNCTION DERF

      INTEGER NF,IH,IL,I,ix

      REAL*8 XF(NF),F(NF),SIGMA,X0,FG,DNSIGMA
      REAL*8 WS1(NF),WS2(NF)

      REAL*8 XL,XH,YL,YH,H,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,EPS

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

C- CHECK ASCENDING ORDER

      DO I=2,NF
        IF (XF(I).LE.XF(I-1))
     &    STOP '*** ERROR IN UTIL_FOLD_GAUSS_LIN:
     &    ARRAY XF NOT IN ASCENDING ORDER ***'
      ENDDO

      EPS=(XF(NF)-XF(1))*1.0D-10

      DO IL=1,NF-1

        IH=IL+1

        XL=XF(IL)
        XH=XF(IH)
        YL=F(IL)
        YH=F(IH)

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

        WS1(IL)=(XH*YL-XL*YH)/h
        WS2(IL)=(YH-YL)/h

      ENDDO !NF-1

      FG=0.0D0

      SN=DNSIGMA*SIGMA
      S2=SIGMA*SIGMA
      S22=2.0D0*S2
      SNR21=1.0D0/(ROOT2*SIGMA)
      SR2PI1=R2PI1/SIGMA

c      DO I=1,NF

        I=IX

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
     &          SR2PI1*(-EXP(-DX**2/S22)*S2*WS2(IL)
     &          +SQPI2*SIGMA*(WS1(IL)+X0*WS2(IL))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(-EXP(-DX**2/S22)*S2*WS2(IL)
     &          +SQPI2*SIGMA*(WS1(IL)+X0*WS2(IL))*
     &          DERF(DX*SNR21))

              FG=FG+FGH-FGL

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
     &          WS2(IL))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)
     &          +X0*(WS2(IL)))*
     &          DERF(DX*SNR21))

              DX=XL-X0

              FGL=
     &          SR2PI1*(
     &          -EXP(-DX**2/S22)*S2*(
     &          WS2(IL))
     &          +SQPI2*SIGMA*(
     &          WS1(IL)
     &          +X0*(WS2(IL)))*
     &          DERF(DX*SNR21))

              FG=FG+FGH-FGL

            ELSE
              GOTO 82
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IH

 82       CONTINUE

        ELSE IF (X0+SN.GT.XF(NF)) THEN

          GOTO 88

        ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

c      ENDDO !NF

 88   CONTINUE

      RETURN
      END
