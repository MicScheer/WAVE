*CMZ :  2.63/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.52/00 24/06/2004  15.55.45  by  Michael Scheer
*CMZ :  2.41/01 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.41/00 20/03/2002  19.16.55  by  Michael Scheer
*CMZ :  2.30/01 12/04/2001  15.28.59  by  Michael Scheer
*CMZ :  2.20/11 11/04/2001  15.24.21  by  Michael Scheer
*CMZ :  2.20/10 10/04/2001  11.16.57  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.14.32  by  Michael Scheer
*CMZ :  1.00/00 06/08/97  17.48.45  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  10.39.53  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  10.39.41  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.50  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WAVE_TRACK_INTER(
     &           T,X,Y,Z,VX,VY,VZ,VXP,VYP,VZP,BS,ICOUNT,GAMMA)
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,track.
      include 'track.cmn'
*KEND.

      INTEGER KLO,KHI,K,ICOUNT,MODE,I

      DOUBLE PRECISION T,X,Y,Z,VX,VY,VZ,VXP,VYP,VZP,BS
     &                ,H,H6,H26,A,A2,A21H6,A3AH26,B,B2,B21H6,B3BH26
     &                ,DT,DT10,GAMMA

      IF (ICOUNT.EQ.0) THEN
          MODE=0
          DT=(DWT(2)-DWT(1))
          DT10=DT*1.D-10
          DO I=2,MCO
         IF (ABS(DWT(I)-DWT(I-1)-DT).GT.DT10) THEN
           MODE=1
           GOTO 19
                ENDIF
          ENDDO
19        KLO=1
          KHI=MCO
          ICOUNT=1
        ENDIF

      IF (MODE.EQ.1) THEN

        IF (KLO.GE.MCO.OR.KLO.LT.1.OR.KHI.GT.MCO.OR.KHI.LT.2) THEN
          KLO=1
          KHI=MCO
        ENDIF

        IF (T.GE.DWT(KLO).AND.T.LT.DWT(KLO+1)) THEN
          KHI=KLO+1
          GOTO 2
        ELSE IF (T.LT.DWT(KLO).OR.T.GE.DWT(KHI)) THEN
        KLO=1
        KHI=MCO
        ENDIF

      K=1
11    K=K*2
      KHI=KLO+K
      IF (KHI.GE.MCO) GOTO 12
      IF (T.GT.DWT(KHI)) THEN
          KLO=KHI
          GOTO 11
      ELSE
          GOTO 1
      ENDIF

12    KHI=MCO

1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(DWT(K).GT.T)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF

      ELSE !MODE

      IF (T.GE.DWT(1).AND.T.LT.DWT(MCO)) THEN
         KLO=T/DT+1
         KHI=KLO+1
         IF (KHI.GT.MCO) THEN
             KHI=MCO
             KLO=KHI-1
         ENDIF
      ELSE IF (T.LT.DWT(1)) THEN
          KLO=1
          KHI=2
      ELSE IF (T.GE.DWT(MCO)) THEN
          KLO=MCO-1
          KHI=MCO
      ENDIF

      ENDIF !MODE

2     H=DWT(KHI)-DWT(KLO)

      IF (H.EQ.0.) THEN
        WRITE(6,*) '*** ERROR IN WAVE_TRACK_INTER: BAD INPUT ***'
        STOP
      ENDIF

      H6=H/6.D0
      H26=H6*H
      A=(DWT(KHI)-T)/H
      A2=A*A
      A3AH26=(A2-1.D0)*A*H26
      A21H6=(-3.D0*A2+1.D0)*H6
      B=(T-DWT(KLO))/H
      B2=B*B
      B21H6=(3.D0*B2-1.D0)*H6
      B3BH26=(B2-1.D0)*B*H26

      X=A*DWX(KLO)+B*DWX(KHI)+A3AH26*DWX2P(KLO)+B3BH26*DWX2P(KHI)
      VX=(-DWX(KLO)+DWX(KHI))/H+A21H6*DWX2P(KLO)+B21H6*DWX2P(KHI)
      VXP=A*DWX2P(KLO)+B*DWX2P(KHI)

      Y=A*DWY(KLO)+B*DWY(KHI)+A3AH26*DWY2P(KLO)+B3BH26*DWY2P(KHI)
      VY=(-DWY(KLO)+DWY(KHI))/H+A21H6*DWY2P(KLO)+B21H6*DWY2P(KHI)
      VYP=A*DWY2P(KLO)+B*DWY2P(KHI)

      Z=A*DWZ(KLO)+B*DWZ(KHI)+A3AH26*DWZ2P(KLO)+B3BH26*DWZ2P(KHI)
      VZ=(-DWZ(KLO)+DWZ(KHI))/H+A21H6*DWZ2P(KLO)+B21H6*DWZ2P(KHI)
      VZP=A*DWZ2P(KLO)+B*DWZ2P(KHI)

      BS=A*DWB(KLO)+B*DWB(KHI)+A3AH26*DWB2P(KLO)+B3BH26*DWB2P(KHI)

      GAMMA=(TRAGAM(KLO)+TRAGAM(KHI))/2.0D0

      RETURN
      END
