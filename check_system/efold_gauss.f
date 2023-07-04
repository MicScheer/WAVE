*CMZ :  4.00/07 05/04/2020  19.23.35  by  Michael Scheer
*CMZ :  4.00/04 14/05/2019  10.42.12  by  Michael Scheer
*CMZ :  3.03/02 03/12/2015  13.58.36  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.63/03 02/05/2008  14.41.00  by  Michael Scheer
*CMZ :  2.52/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.51/03 23/06/2004  12.05.57  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.45  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.09.58  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.57.23  by  Michael Scheer
*CMZ :  1.00/00 19/08/97  16.17.52  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  10.12.13  by  Michael Scheer
*CMZ : 00.01/06 14/02/95  10.17.13  by  Michael Scheer
*CMZU: 00.01/04 18/01/95  18.16.43  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.21.01  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.43  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.46  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE EFOLD_GAUSS(NF,XF,F,SIGMA,NSIGMA,X0,kmode,FG)
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
C-       NSIGMA: NUMBER OF SIGMAS TO BE CONSIDERED
C-       X0:   FG(X0) IS CALCULATED
C-               X0-NSIGM*SIGMA MUST NOT BE LOWER THAN XF(1)
C-               X0+NSIGM*SIGMA MUST NOT EXCEED XF(NF)
C-       kmode: CONTROL FLAG:
C-             kmode.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             kmode.LT.0: NEW INITIALIZATION

C--   OUTPUT:

C-       FG:   FG(X0) IS CALCULATED


      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.

      INTEGER NF,kmode,MSTEPP,NSTEP,ISTEP,NSIGMA,NLOW,NHIGH
      INTEGER ICOUNT,IL,IH,I,kl,kh

      DOUBLE PRECISION XF(NF),F(NF),SIGMA,X0,FG,XLOW,XHIGH,DXW,G,FWW
      DOUBLE PRECISION XW(NDFREQP),F2(NDFREQP),COEF(NDFREQP),FW(NDFREQP),DR2P1
      DOUBLE PRECISION VS1(NDFREQP),VS2(NDFREQP),VS3(NDFREQP),VS4(NDFREQP)
      DOUBLE PRECISION WS1(NDFREQP),WS2(NDFREQP),WS3(NDFREQP),WS4(NDFREQP)

      EXTERNAL FUNCTION DERF

       DOUBLE PRECISION  CH,CL,CH2,CL2,CHCL,CH2CL,CHCL2,XL,XH,YL,YH,H,H61,
     &  XHXL,XH2,XL2,SN,S2,ROOT2,SNR21,DERF,FGH,FGL,R2PI1,DX,X02,S22,
     &  SQPI2,X02S2,X023S2,SR2PI1,RNSIGMA

      DATA ROOT2/1.4142135623731D0/
      DATA R2PI1/0.398942280401433D0/
      DATA DR2P1/0.398942280401433D0/
      DATA SQPI2/1.2533141373155D0/

      save

      IF (NF.GT.NDFREQP) THEN
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED  ***'
      ENDIF

C- SPLINES OF FUNCTION F

      IF (kmode.LT.0) THEN
        CALL util_spline_coef(XF,F,NF,-9999.0d0,-9999.0d0,F2,WS1,WS2,WS3,WS4)
      ENDIF !kmode


      IF (IEFOLD.GT.0) THEN

        IF (kmode.LT.0) THEN
          DO IL=1,NF-1

            IH=IL+1

            XL=XF(IL)
            XH=XF(IH)
            YL=F(IL)
            YH=F(IH)
            CL=F2(IL)
            CH=F2(IH)

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
     &          '*** ERROR SR EFOLD_GAUSS:'
              PRINT*,
     &          '*** E-ARRAY NOT IN ASCENDING ORDER'
              STOP
            ENDIF

            H61=1.0D0/(6.0D0*H)

            WS1(IL)=((CHCL2*XH-CH2CL*XL)*XHXL+6.0D0*(XH*YL-XL*YH))*H61
            WS2(IL)=((CH2-CL2)*XHXL-CHCL2*XH2+CH2CL*XL2+6.0D0*(YH-YL))*H61
            WS3(IL)=(-CH*XL+CL*XH)/(2.0D0*H)
            WS4(IL)=CHCL*H61

          ENDDO !NF-1
        ENDIF !kmode

        if (kmode.lt.0) then
          nlow=1
        endif

        if (x0.le.xf(nlow)) then
          nlow=1
        endif

        I=-1
        DO IL=nlow,nf
          IF (X0.EQ.XF(IL)) THEN
            I=IL
            nlow=il
            GOTO 1
          ENDIF
        ENDDO

1       IF (I.EQ.-1) THEN
          PRINT*,
     &      '*** ERROR SR EFOLD_GAUSS:'
          PRINT*,
     &      '*** BAD ENERGY FOR FOLDING'
          PRINT*,X0
          STOP
        ENDIF

        FG=0.0D0

        RNSIGMA=NSIGMA
        SN=RNSIGMA*SIGMA
        S2=SIGMA*SIGMA
        S22=2.0D0*S2
        SNR21=1.0D0/(ROOT2*SIGMA)
        SR2PI1=R2PI1/SIGMA

        X02=X0*X0
        X02S2=S2+X02
        X023S2=S22+X02S2

        IF (X0-SN.GE.XF(1).AND.X0+SN.LE.XF(NF)) THEN

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

              FG=FG+FGH-FGL

            ELSE
              GOTO 82
            ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))

          ENDDO !IH
 82       CONTINUE

        ELSE IF (X0+SN.GT.XF(NF)) THEN
          GOTO 88
        ENDIF ! (X-SN.GE.XF(1).AND.X+SN.LE.XF(NF))
 88     CONTINUE

      ELSE !IEFOLD

C OLD VERSION

C-- SET UP FOLDING BUFFER

      XLOW =X0-NSIGMA*SIGMA
      XHIGH=X0+NSIGMA*SIGMA

      IF (DABS(XLOW-XHIGH)/XHIGH.LT.1.D-6) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(LUNGFO,*)'DABS(XLOW-XHIGH)/XHIGH.LT.1.D-LUNGFO'
        WRITE(LUNGFO,*)'ESPREAD TOO SMALL, CHECK INPUT FILE'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(6,*)'DABS(XLOW-XHIGH)/XHIGH.LT.1.D-6'
        WRITE(6,*)'ESPREAD TOO SMALL, CHECK INPUT FILE'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED  ***'
      ENDIF

      IF (XLOW.LT.XF(1)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(LUNGFO,*)'X0 IS LOWER THEN XF(1)'
        WRITE(LUNGFO,*)'CHECK INPUT TO ROUTINE'
        WRITE(LUNGFO,*)'X0,XF(1):',X0,XF(1)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(6,*)'X0 IS LOWER THEN XF(1)'
        WRITE(6,*)'CHECK INPUT TO ROUTINE'
        WRITE(6,*)'X0,XF(1):',X0,XF(1)
        STOP '*** PROGRAM WAVE ABORTED  ***'
      ENDIF

      IF (XHIGH.GT.XF(NF)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(LUNGFO,*)'X0 IS GREATER THEN XF(NF)'
        WRITE(LUNGFO,*)'CHECK INPUT TO ROUTINE'
        WRITE(LUNGFO,*)'X0,XF(NF):',X0,XF(NF)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN EFOLD_GAUSS ***'
        WRITE(6,*)'X0 IS GREATER THEN XF(NF)'
        WRITE(6,*)'CHECK INPUT TO ROUTINE'
        WRITE(6,*)'X0,XF(NF):',X0,XF(NF)
        STOP '*** PROGRAM WAVE ABORTED  ***'
      ENDIF


C- NUMBER OF POINTS WITHIN FOLDING INTERVALL

      if (kmode.lt.0) then
        nlow=1
        nhigh=nf
      endif

      if (x0.le.xf(nlow)) then
        nlow=1
      endif

      kl=nlow

      DO ISTEP=kl,NF
        IF (XF(ISTEP).LT.XLOW) THEN
          NLOW=ISTEP
        ELSE
          GOTO 100
        ENDIF
      ENDDO
100   CONTINUE

      nhigh=nhigh+1
      if (nhigh.gt.nf) then
        nhigh=nf
      endif

      kh=nhigh
      DO ISTEP=kh,1,-1
        IF (XF(ISTEP).GT.XHIGH) THEN
          NHIGH=ISTEP
        ELSE
          GOTO 101
        ENDIF
      ENDDO
101   CONTINUE

      MSTEPP=(NSIGMA*5)/2*2+1
      NSTEP=MAX(NHIGH-NLOW,MSTEPP)

      XW(1)=XLOW
      DXW=(XHIGH-XLOW)/(NSTEP-1)
      DO ISTEP=2,NSTEP
        XW(ISTEP)=XW(ISTEP-1)+DXW
      ENDDO

CERROR 19.8.97      DO ISTEP=1,NSTEP-1
      DO ISTEP=1,NSTEP
        IF (XW(ISTEP).GE.XF(1).AND.XW(ISTEP).LE.XF(NF)) THEN
          CALL WAVE_SPLINE_INTER(XF,F,F2,NF,XW(ISTEP),FWW,-1,ICOUNT)
          G=DEXP(-0.5D0*((XW(ISTEP)-X0)/SIGMA)**2)*DR2P1/SIGMA
          FW(ISTEP)=FWW*G
C         WRITE(8,*)XW(ISTEP),FWW,G,FW(ISTEP)
        ELSE
          FW(ISTEP)=0.0d0
        ENDIF
      ENDDO

C- ACTUAL INTEGRATION

      CALL WAVE_SPLINE_INTEGRAL(XW,FW,NSTEP,FG,COEF,VS1,VS2,VS3,VS4)

C     CALL WAVE_SIMPSON(NSTEP,XW,FW,WS1(1))
C     IF ( (WS1(1)+FG).NE.0.D0
C     & .AND.DABS((WS1(1)-FG)/(WS1(1)+FG) ).GT.0.02
C     &) THEN
C        WRITE(6,*)
C        WRITE(6,*)'SPLINE RESULTAT:',FG
C        WRITE(6,*)'SIMPSON RESULTAT:',WS1(1)
C        WRITE(6,*)
C     ENDIF

      ENDIF !IEFOLD

      RETURN
      END
