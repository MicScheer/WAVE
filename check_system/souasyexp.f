*CMZ :  2.50/00 28/09/2009  13.02.09  by  Michael Scheer
*CMZ :  2.31/01 24/04/2001  17.56.57  by  Michael Scheer
*CMZ :  2.20/12 11/04/2001  16.40.07  by  Michael Scheer
*CMZ :  2.20/11 11/04/2001  15.48.38  by  Michael Scheer
*CMZ :  2.20/10 10/04/2001  11.24.36  by  Michael Scheer
*CMZ :  2.20/09 03/04/2001  10.29.02  by  Michael Scheer
*CMZ :  2.15/01 30/03/2001  19.35.06  by  Michael Scheer
*CMZ :  2.20/07 18/03/2001  16.50.02  by  Michael Scheer
*CMZ :  2.20/05 15/03/2001  14.46.11  by  Michael Scheer
*CMZ :  2.20/04 09/03/2001  16.47.40  by  Michael Scheer
*CMZ :  2.20/03 23/02/2001  15.04.13  by  Michael Scheer
*CMZ :  2.20/02 21/02/2001  11.30.46  by  Michael Scheer
*CMZ :  2.20/01 20/02/2001  14.18.37  by  Michael Scheer
*-- Author : Michael Scheer

      SUBROUTINE SOUASYEXP(IVELOFIELD,IROIASY
     &                     ,X2,Y2,Z2,VX2,VY2,VZ2,OBSVX,OBSVY,OBSVZ
     &                     ,BET1N,EXPOMR,EXPOMI,OM,BETP,BETPP,C,DMYGAMMA
     &                     ,A1R,A1I,A2R,A2I)
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

C--- EVALUATE ASYMPTOTIC EXPANSION UP TO SECOND ORDER

         IMPLICIT NONE

      DOUBLE PRECISION AI1R,AI1I,AI2R,AI2I,A1R(3),A1I(3),A2R(3),A2I(3)
        DOUBLE PRECISION RX,RY,RZ,RNX,RNY,RNZ,EXPOMR,EXPOMI,OM,C,DMYGAMMA
      DOUBLE PRECISION X2,Y2,Z2,VX2,VY2,VZ2,OBSVX,OBSVY,OBSVZ

      DOUBLE PRECISION RN(3),BET(3),BETP(3),BETPP(3),RP,RNBET(3),RNP(3),V(3)
     &                  ,VDUM1(3),VDUM2(3),VDUM3(3)
     &                  ,F(3),F1(3),F2,F3,F1P(3),F2P,F3P
     &                  ,FP(3),F23,F23P,RNPBETP(3),BET1N,BET1NP,FPPHIP(3)
     &                  ,R,R1,FV(3),FVP(3),FV1(3),FV1P(3),FVPPHIP(3)

      INTEGER ICOMP,IVELOFIELD,IROIASY



C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

          RX=OBSVX-X2
          RY=OBSVY-Y2
          RZ=OBSVZ-Z2

          R=SQRT(RX*RX+RY*RY+RZ*RZ)
          R1=1.D0/R

          RNX=RX*R1
          RNY=RY*R1
          RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

          RN(1)=RNX
          RN(2)=RNY
          RN(3)=RNZ

          V(1)=VX2
          V(2)=VY2
          V(3)=VZ2

          BET(1)=V(1)/C
          BET(2)=V(2)/C
          BET(3)=V(3)/C

          RP=-RN(1)*V(1)-RN(2)*V(2)-RN(3)*V(3)

          RNP(1)=(-V(1)-RN(1)*RP)/R;
          RNP(2)=(-V(2)-RN(2)*RP)/R;
          RNP(3)=(-V(3)-RN(3)*RP)/R;

          RNBET(1)=RN(1)-BET(1)
          RNBET(2)=RN(2)-BET(2)
          RNBET(3)=RN(3)-BET(3)

          RNPBETP(1)=RNP(1)-BETP(1)
          RNPBETP(2)=RNP(2)-BETP(2)
          RNPBETP(3)=RNP(3)-BETP(3)

          F2=R
          F3=BET1N*BET1N
          F23=F2*F3

          CALL UTIL_VCROSS_VCROSS(RN,RNBET,BETP,F1)

          F(1)=F1(1)/F23
          F(2)=F1(2)/F23
          F(3)=F1(3)/F23

          CALL UTIL_VCROSS_VCROSS(RNP,RNBET,BETP,VDUM1)
          CALL UTIL_VCROSS_VCROSS(RN,RNPBETP,BETP,VDUM2)
          CALL UTIL_VCROSS_VCROSS(RN,RNBET,BETPP,VDUM3)

          F1P(1)=VDUM1(1)+VDUM2(1)+VDUM3(1)
          F1P(2)=VDUM1(2)+VDUM2(2)+VDUM3(2)
          F1P(3)=VDUM1(3)+VDUM2(3)+VDUM3(3)

          F2P=RP
          BET1NP=
     &                 (-BETP(1)*RN(1)-BETP(2)*RN(2)-BETP(3)*RN(3)
     &                   -BET(1)*RNP(1)-BET(2)*RNP(2)-BET(3)*RNP(3))
          F3P=2.D0*BET1N*BET1NP
          F23P=F2P*F3+F2*F3P

          FP(1)=(F1P(1)*F23-F1(1)*F23P)/(F23*F23)
          FP(2)=(F1P(2)*F23-F1(2)*F23P)/(F23*F23)
          FP(3)=(F1P(3)*F23-F1(3)*F23P)/(F23*F23)

          FPPHIP(1)=(FP(1)*BET1N-F(1)*BET1NP)/F3
          FPPHIP(2)=(FP(2)*BET1N-F(2)*BET1NP)/F3
          FPPHIP(3)=(FP(3)*BET1N-F(3)*BET1NP)/F3

          FV1(1)=RNBET(1)*C/DMYGAMMA/DMYGAMMA/F2
          FV1(2)=RNBET(2)*C/DMYGAMMA/DMYGAMMA/F2
          FV1(3)=RNBET(3)*C/DMYGAMMA/DMYGAMMA/F2

          FV(1)=FV1(1)/F23
          FV(2)=FV1(2)/F23
          FV(3)=FV1(3)/F23

          FV1P(1)=C/DMYGAMMA/DMYGAMMA*
     &                      (RNPBETP(1)*F2-RNBET(1)*F2P)/F2/F2
          FV1P(2)=C/DMYGAMMA/DMYGAMMA*
     &                      (RNPBETP(2)*F2-RNBET(2)*F2P)/F2/F2
          FV1P(3)=C/DMYGAMMA/DMYGAMMA*
     &                      (RNPBETP(3)*F2-RNBET(3)*F2P)/F2/F2

          FVP(1)=(FV1P(1)*F23-FV1(1)*F23P)/(F23*F23)
          FVP(2)=(FV1P(2)*F23-FV1(2)*F23P)/(F23*F23)
          FVP(3)=(FV1P(3)*F23-FV1(3)*F23P)/(F23*F23)

          FVPPHIP(1)=(FVP(1)*BET1N-FV(1)*BET1NP)/F3
          FVPPHIP(2)=(FVP(2)*BET1N-FV(2)*BET1NP)/F3
          FVPPHIP(3)=(FVP(3)*BET1N-FV(3)*BET1NP)/F3

                 DO ICOMP=1,3

             IF (IVELOFIELD.EQ.0) THEN
              AI1R=(F(ICOMP)+FV(ICOMP))/BET1N*EXPOMI/OM
              AI1I=-(F(ICOMP)+FV(ICOMP))/BET1N*EXPOMR/OM
              AI2R=(FPPHIP(ICOMP)+FVPPHIP(ICOMP))/BET1N/OM*EXPOMR/OM
              AI2I=(FPPHIP(ICOMP)+FVPPHIP(ICOMP))/BET1N/OM*EXPOMI/OM
                   ELSE IF (IVELOFIELD.EQ.1) THEN
              AI1R=F(ICOMP)/BET1N*EXPOMI/OM
              AI1I=-F(ICOMP)/BET1N*EXPOMR/OM
              AI2R=FPPHIP(ICOMP)/BET1N/OM*EXPOMR/OM
              AI2I=FPPHIP(ICOMP)/BET1N/OM*EXPOMI/OM
                   ELSE IF (IVELOFIELD.LT.0) THEN
              AI1R=FV(ICOMP)/BET1N*EXPOMI/OM
              AI1I=-FV(ICOMP)/BET1N*EXPOMR/OM
              AI2R=FVPPHIP(ICOMP)/BET1N/OM*EXPOMR/OM
              AI2I=FVPPHIP(ICOMP)/BET1N/OM*EXPOMI/OM
             ELSE    !IVELOFIELD
                     WRITE(6,*)
     &                 '*** ERROR IN SOUASYEXP: BAD VALUE OF IVELOFIELD  ***'
                     WRITE(6,*) '*** PROGRAM WAVE ABORTED  ***'
                     STOP
                   ENDIF   !IVELOFIELD

                A1R(ICOMP)=AI1R
                A1I(ICOMP)=AI1I

             IF (ABS(IROIASY).GT.1) THEN
                  A2R(ICOMP)=AI2R
                  A2I(ICOMP)=AI2I
             ELSE
                  A2R(ICOMP)=0.D0
                  A2I(ICOMP)=0.D0
             ENDIF

           ENDDO  !ICOMP

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

         RETURN
      END
