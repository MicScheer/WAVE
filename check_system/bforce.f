*CMZ :  3.02/04 25/11/2014  10.04.56  by  Michael Scheer
*CMZ :  3.02/00 10/10/2014  11.20.47  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.61/02 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.52/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.01/00 04/12/97  14.32.11  by  Michael Scheer
*-- Author :    Michael Scheer   26/11/97

      SUBROUTINE BFORCE

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

        implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,bforce.
      include 'bforce.cmn'
*KEND.

      INTEGER ICAL
      DATA ICAL/0/

      INTEGER NPOIX,NPOIY,NPOIZ,IY,IZ,IX,I,NBFDIMP

      DOUBLE PRECISION RMU0,BFINT(3),XX,YY,ZZ,XBASE,YBASE,ZBASE
     &  ,DX,DY,DZ,BX,BY,BZ,AX,AY,AZ,BNX,BNY,BNZ
     &  ,F(3),R(3),T(3),RESULT,
     &  RX,RY,RZ,TORQTOTX,TORQTOTY,TORQTOTZ,TORQROTX,TORQROTY,TORQROTZ

      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE ::
     &  XB,YB,WS1,WS2,WS3,WS4,COEF,ZB

      DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE :: BB,BBI,BT,BTI

      NBFDIMP=MAX(NBFORCX,NBFORCY,NBFORCZ)+1

      ALLOCATE(XB(NBFDIMP))
      ALLOCATE(YB(NBFDIMP))
      ALLOCATE(ZB(NBFDIMP))
      ALLOCATE(WS1(NBFDIMP))
      ALLOCATE(WS2(NBFDIMP))
      ALLOCATE(WS3(NBFDIMP))
      ALLOCATE(WS4(NBFDIMP))
      ALLOCATE(COEF(NBFDIMP))

      ALLOCATE(BB(NBFDIMP,3))
      ALLOCATE(BBI(NBFDIMP,3))
      ALLOCATE(BT(NBFDIMP,3))
      ALLOCATE(BTI(NBFDIMP,3))

      IF (ICAL.EQ.0) THEN

        RMU0=4.D0*PI1/1.D7

        NPOIX=NBFORCX+1
        NPOIY=NBFORCY+1
        NPOIZ=NBFORCZ+1

        IF (
     &      NPOIX.GT.NBFDIMP
     &      .OR.
     &      NPOIY.GT.NBFDIMP
     &      .OR.
     &      NPOIZ.GT.NBFDIMP
     &      ) THEN

          WRITE(LUNGFO,*)'*** ERROR IN BFORCE:'
          WRITE(LUNGFO,*)
     &      'DIMENSION NBFDIMP IN BFORCE.CMN EXCEEDED'
          WRITE(LUNGFO,*)
     &      'CURRENT VALUE IS:',NBFDIMP
          WRITE(LUNGFO,*)
     &      '(NUMBER OF POINTS IS NUMBER OF MASHES+1 IN EACH DIMENSION)'
          WRITE(6,*)'*** ERROR IN BFORCE:'
          WRITE(6,*)
     &      'DIMENSION NBFDIMP IN BFORCE.CMN EXCEEDED'
          WRITE(6,*)
     &      'CURRENT VALUE IS:',NBFDIMP
          WRITE(6,*)
     &      '(NUMBER OF POINTS IS NUMBER OF MASHES+1 IN EACH DIMENSION)'
          STOP

        ENDIF   !NBFDIMP

        IF (NBFORCX.NE.0) THEN
          DX=BFLENX/NBFORCX
        ELSE
          DX=0.0
        ENDIF

        IF (NBFORCY.NE.0) THEN
          DY=BFLENY/NBFORCY
        ELSE
          DY=0.0
        ENDIF

        IF (NBFORCZ.NE.0) THEN
          DZ=BFLENZ/NBFORCZ
        ELSE
          DZ=0.0
        ENDIF

        XBASE=BFCENX-BFLENX/2.D0
        YBASE=BFCENY-BFLENY/2.D0
        ZBASE=BFCENZ-BFLENZ/2.D0

        ICAL=1

      ENDIF !ICAL

C FIRST YZ-PLANE, NORMAL VECTOR IS (-1,0,0){

      BNX=-1.D0
      BNY=0.D0
      BNZ=0.D0

      XX=XBASE
      YY=YBASE-DY

      DO IY=1,NPOIY

        YY=YY+DY
        YB(IY)=YY
        ZZ=ZBASE-DZ

        DO IZ=1,NPOIZ

          ZZ=ZZ+DZ
          CALL MYBFELD(XX,YY,ZZ,BX,BY,BZ,AX,AY,AZ)

          ZB(IZ)=ZZ
          BB(IZ,1)=BX*(BX*BNX+BY*BNY+BZ*BNZ)
     &      -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNX)
          BB(IZ,2)=BY*(BX*BNX+BY*BNY+BZ*BNZ)
     &      -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNY)
          BB(IZ,3)=BZ*(BX*BNX+BY*BNY+BZ*BNZ)
     &      -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNZ)

          F(1)=BB(IZ,1)
          F(2)=BB(IZ,2)
          F(3)=BB(IZ,3)

          R(1)=XX-BFCENX
          R(2)=YY-BFCENY
          R(3)=ZZ-BFCENZ

          CALL UTIL_VCROSS(R,F,T)

          BT(IZ,1)=T(1)
          BT(IZ,2)=T(2)
          BT(IZ,3)=T(3)

        ENDDO  !IZ

        CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,1),NPOIZ,RESULT
     &    ,COEF,WS1,WS2,WS3,WS4)
        BBI(IY,1)=RESULT
        CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,2),NPOIZ,RESULT
     &    ,COEF,WS1,WS2,WS3,WS4)
        BBI(IY,2)=RESULT
        CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,3),NPOIZ,RESULT
     &    ,COEF,WS1,WS2,WS3,WS4)
        BBI(IY,3)=RESULT

        CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,1),NPOIZ,RESULT
     &    ,COEF,WS1,WS2,WS3,WS4)
        BTI(IY,1)=RESULT
        CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,2),NPOIZ,RESULT
     &    ,COEF,WS1,WS2,WS3,WS4)
        BTI(IY,2)=RESULT
        CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,3),NPOIZ,RESULT
     &    ,COEF,WS1,WS2,WS3,WS4)
        BTI(IY,3)=RESULT

      ENDDO !IY

      CALL WAVE_SPLINE_INTEGRAL(YB,BBI(1,1),NPOIY,RESULT
     &  ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BBI(1,2),NPOIY,RESULT
     &  ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BBI(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      BFORCX(1)=BFINT(1)/RMU0
      BFORCY(1)=BFINT(2)/RMU0
      BFORCZ(1)=BFINT(3)/RMU0

      CALL WAVE_SPLINE_INTEGRAL(YB,BTI(1,1),NPOIY,RESULT
     &  ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BTI(1,2),NPOIY,RESULT
     &  ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BTI(1,3),NPOIY,RESULT
     &  ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      TORQX(1)=BFINT(1)/RMU0
      TORQY(1)=BFINT(2)/RMU0
      TORQZ(1)=BFINT(3)/RMU0

C }FIRST YZ-PLANE, NORMAL VECTOR IS (-1,0,0)

C SECOND YZ-PLANE, NORMAL VECTOR IS (1,0,0){

      BNX=1.D0
      BNY=0.D0
      BNZ=0.D0

      XX=XBASE+BFLENX
      YY=YBASE-DY

      DO IY=1,NPOIY

        YY=YY+DY
        YB(IY)=YY
        ZZ=ZBASE-DZ

        DO IZ=1,NPOIZ

          ZZ=ZZ+DZ
          CALL MYBFELD(XX,YY,ZZ,BX,BY,BZ,AX,AY,AZ)

          ZB(IZ)=ZZ
          BB(IZ,1)=BX*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNX)
          BB(IZ,2)=BY*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNY)
          BB(IZ,3)=BZ*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNZ)

          F(1)=BB(IZ,1)
          F(2)=BB(IZ,2)
          F(3)=BB(IZ,3)

          R(1)=XX-BFCENX
          R(2)=YY-BFCENY
          R(3)=ZZ-BFCENZ

          CALL UTIL_VCROSS(R,F,T)

          BT(IZ,1)=T(1)
          BT(IZ,2)=T(2)
          BT(IZ,3)=T(3)

        ENDDO  !IZ

      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,1),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IY,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,2),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IY,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,3),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IY,3)=RESULT

      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,1),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IY,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,2),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IY,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,3),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IY,3)=RESULT

      ENDDO !IY

      CALL WAVE_SPLINE_INTEGRAL(YB,BBI(1,1),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BBI(1,2),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BBI(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      BFORCX(3)=BFINT(1)/RMU0
      BFORCY(3)=BFINT(2)/RMU0
      BFORCZ(3)=BFINT(3)/RMU0

      CALL WAVE_SPLINE_INTEGRAL(YB,BTI(1,1),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BTI(1,2),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BTI(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      TORQX(3)=BFINT(1)/RMU0
      TORQY(3)=BFINT(2)/RMU0
      TORQZ(3)=BFINT(3)/RMU0

C } SECOND YZ-PLANE, NORMAL VECTOR IS (1,0,0)

C FIRST XZ-PLANE, NORMAL VECTOR IS (0,-1,0){

      BNX=0.D0
      BNY=-1.D0
      BNZ=0.D0

      XX=XBASE-DX
      YY=YBASE

      DO IX=1,NPOIX

        XX=XX+DX
        XB(IX)=XX
        ZZ=ZBASE-DZ

        DO IZ=1,NPOIZ

          ZZ=ZZ+DZ
          CALL MYBFELD(XX,YY,ZZ,BX,BY,BZ,AX,AY,AZ)

          ZB(IZ)=ZZ
          BB(IZ,1)=BX*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNX)
          BB(IZ,2)=BY*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNY)
          BB(IZ,3)=BZ*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNZ)

          F(1)=BB(IZ,1)
          F(2)=BB(IZ,2)
          F(3)=BB(IZ,3)

          R(1)=XX-BFCENX
          R(2)=YY-BFCENY
          R(3)=ZZ-BFCENZ

          CALL UTIL_VCROSS(R,F,T)

          BT(IZ,1)=T(1)
          BT(IZ,2)=T(2)
          BT(IZ,3)=T(3)

        ENDDO  !IZ

      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,1),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,2),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,3),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,3)=RESULT

      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,1),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,2),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,3),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,3)=RESULT

      ENDDO !IX

      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      BFORCX(2)=BFINT(1)/RMU0
      BFORCY(2)=BFINT(2)/RMU0
      BFORCZ(2)=BFINT(3)/RMU0

      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      TORQX(2)=BFINT(1)/RMU0
      TORQY(2)=BFINT(2)/RMU0
      TORQZ(2)=BFINT(3)/RMU0

C }FIRST XZ-PLANE, NORMAL VECTOR IS (0,-1,0)

C SECOND XZ-PLANE, NORMAL VECTOR IS (0,+1,0){

      BNX=0.D0
      BNY=1.D0
      BNZ=0.D0

      XX=XBASE-DX
      YY=YBASE+BFLENY

      DO IX=1,NPOIX

        XX=XX+DX
        XB(IX)=XX
        ZZ=ZBASE-DZ

        DO IZ=1,NPOIZ

          ZZ=ZZ+DZ
          CALL MYBFELD(XX,YY,ZZ,BX,BY,BZ,AX,AY,AZ)

          ZB(IZ)=ZZ
          BB(IZ,1)=BX*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNX)
          BB(IZ,2)=BY*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNY)
          BB(IZ,3)=BZ*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNZ)

          F(1)=BB(IZ,1)
          F(2)=BB(IZ,2)
          F(3)=BB(IZ,3)

          R(1)=XX-BFCENX
          R(2)=YY-BFCENY
          R(3)=ZZ-BFCENZ

          CALL UTIL_VCROSS(R,F,T)

          BT(IZ,1)=T(1)
          BT(IZ,2)=T(2)
          BT(IZ,3)=T(3)

        ENDDO  !IZ

      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,1),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,2),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BB(1,3),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,3)=RESULT

      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,1),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,2),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(ZB,BT(1,3),NPOIZ,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,3)=RESULT

      ENDDO !IX

      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      BFORCX(4)=BFINT(1)/RMU0
      BFORCY(4)=BFINT(2)/RMU0
      BFORCZ(4)=BFINT(3)/RMU0

      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      TORQX(4)=BFINT(1)/RMU0
      TORQY(4)=BFINT(2)/RMU0
      TORQZ(4)=BFINT(3)/RMU0

C }SECOND XZ-PLANE, NORMAL VECTOR IS (0,+1,0)

C FIRST XY-PLANE, NORMAL VECTOR IS (0,0,-1){

      BNX=0.D0
      BNY=0.D0
      BNZ=-1.D0

      XX=XBASE-DX
      ZZ=ZBASE

      DO IX=1,NPOIX

        XX=XX+DX
        XB(IX)=XX
        YY=YBASE-DY

        DO IY=1,NPOIY

          YY=YY+DY
          CALL MYBFELD(XX,YY,ZZ,BX,BY,BZ,AX,AY,AZ)

          YB(IY)=YY
          BB(IY,1)=BX*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNX)
          BB(IY,2)=BY*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNY)
          BB(IY,3)=BZ*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNZ)

          F(1)=BB(IY,1)
          F(2)=BB(IY,2)
          F(3)=BB(IY,3)

          R(1)=XX-BFCENX
          R(2)=YY-BFCENY
          R(3)=ZZ-BFCENZ

          CALL UTIL_VCROSS(R,F,T)

          BT(IY,1)=T(1)
          BT(IY,2)=T(2)
          BT(IY,3)=T(3)

        ENDDO  !IY

      CALL WAVE_SPLINE_INTEGRAL(YB,BB(1,1),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BB(1,2),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BB(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,3)=RESULT

      CALL WAVE_SPLINE_INTEGRAL(YB,BT(1,1),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BT(1,2),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BT(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,3)=RESULT

      ENDDO !IX

      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      BFORCX(5)=BFINT(1)/RMU0
      BFORCY(5)=BFINT(2)/RMU0
      BFORCZ(5)=BFINT(3)/RMU0

      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      TORQX(5)=BFINT(1)/RMU0
      TORQY(5)=BFINT(2)/RMU0
      TORQZ(5)=BFINT(3)/RMU0

C }FIRST XY-PLANE, NORMAL VECTOR IS (0,0,-1)

C SECOND XY-PLANE, NORMAL VECTOR IS (0,0,+1){

      BNX=0.D0
      BNY=0.D0
      BNZ=1.D0

      XX=XBASE-DX
      ZZ=ZBASE+BFLENZ

      DO IX=1,NPOIX

        XX=XX+DX
        XB(IX)=XX
        YY=YBASE-DY

        DO IY=1,NPOIY

          YY=YY+DY
          CALL MYBFELD(XX,YY,ZZ,BX,BY,BZ,AX,AY,AZ)

          YB(IY)=YY
          BB(IY,1)=BX*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNX)
          BB(IY,2)=BY*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNY)
          BB(IY,3)=BZ*(BX*BNX+BY*BNY+BZ*BNZ)
     &              -0.5D0*((BX*BX+BY*BY+BZ*BZ)*BNZ)

          F(1)=BB(IY,1)
          F(2)=BB(IY,2)
          F(3)=BB(IY,3)

          R(1)=XX-BFCENX
          R(2)=YY-BFCENY
          R(3)=ZZ-BFCENZ

          CALL UTIL_VCROSS(R,F,T)

          BT(IY,1)=T(1)
          BT(IY,2)=T(2)
          BT(IY,3)=T(3)

        ENDDO  !IY

      CALL WAVE_SPLINE_INTEGRAL(YB,BB(1,1),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BB(1,2),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BB(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BBI(IX,3)=RESULT

      CALL WAVE_SPLINE_INTEGRAL(YB,BT(1,1),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BT(1,2),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(YB,BT(1,3),NPOIY,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BTI(IX,3)=RESULT

      ENDDO !IX

      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BBI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      BFORCX(6)=BFINT(1)/RMU0
      BFORCY(6)=BFINT(2)/RMU0
      BFORCZ(6)=BFINT(3)/RMU0

      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,1),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(1)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,2),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(2)=RESULT
      CALL WAVE_SPLINE_INTEGRAL(XB,BTI(1,3),NPOIX,RESULT
     & ,COEF,WS1,WS2,WS3,WS4)
      BFINT(3)=RESULT

      TORQX(6)=BFINT(1)/RMU0
      TORQY(6)=BFINT(2)/RMU0
      TORQZ(6)=BFINT(3)/RMU0

C }SECOND XY-PLANE, NORMAL VECTOR IS (0,0,+1)


      BFORCX(7)=0.D0
      BFORCY(7)=0.D0
      BFORCZ(7)=0.D0
      TORQX(7)=0.D0
      TORQY(7)=0.D0
      TORQZ(7)=0.D0
      DO I=1,6
          BFORCX(7)=BFORCX(7)+BFORCX(I)
          BFORCY(7)=BFORCY(7)+BFORCY(I)
          BFORCZ(7)=BFORCZ(7)+BFORCZ(I)
          TORQX(7)=TORQX(7)+TORQX(I)
          TORQY(7)=TORQY(7)+TORQY(I)
          TORQZ(7)=TORQZ(7)+TORQZ(I)
      ENDDO

        RX=BFCENX-TORQCENX
        RY=BFCENY-TORQCENY
        RZ=BFCENZ-TORQCENZ

        TORQROTX=RY*BFORCZ(7)-RZ*BFORCY(7)
        TORQROTY=RZ*BFORCX(7)-RX*BFORCZ(7)
        TORQROTZ=RX*BFORCY(7)-RY*BFORCX(7)

        TORQTOTX=TORQROTX+TORQX(7)
        TORQTOTY=TORQROTY+TORQY(7)
        TORQTOTZ=TORQROTZ+TORQZ(7)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR BFORCE:'
      WRITE(LUNGFO,*)

        IF (KBREC.NE.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'IF CONSIDERED  VOLUME IS INSIDE OR PARTIALLY INSIDE'
          WRITE(LUNGFO,*)'MAGNET, RESULTS ARE WRONG!!'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN BFORCE:'
          WRITE(6,*)
          WRITE(6,*)'IF CONSIDERED  VOLUME IS INSIDE OR PARTIALLY INSIDE'
          WRITE(6,*)'MAGNET, RESULTS ARE WRONG!!'
          WRITE(6,*)
        ENDIF

      WRITE(LUNGFO,*)'     Box parameters:'
      WRITE(LUNGFO,*)'     BFCENX,BFCENY,BFCENZ:'
     &                      ,SNGL(BFCENX),SNGL(BFCENY),SNGL(BFCENZ)
      WRITE(LUNGFO,*)'     BFLENX,BFLENY,BFLENZ:'
     &                      ,SNGL(BFLENX),SNGL(BFLENY),SNGL(BFLENZ)
      WRITE(LUNGFO,*)'     Reference point for total torque:'
      WRITE(LUNGFO,*)'     '
     &                      ,SNGL(TORQCENX),SNGL(TORQCENY),SNGL(TORQCENZ)
      WRITE(LUNGFO,*)'     NBFORCX,NBFORCY,NBFORCZ:'
     &                      ,NBFORCX,NBFORCY,NBFORCZ
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     Plane            Force [N]      Torque [Nm] (with resp. to box center):'

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     low X, X-comp. ',SNGL(BFORCX(1)),SNGL(TORQX(1))
      WRITE(LUNGFO,*)
     &'     low X, Y-comp. ',SNGL(BFORCY(1)),SNGL(TORQY(1))
      WRITE(LUNGFO,*)
     &'     low X, Z-comp. ',SNGL(BFORCZ(1)),SNGL(TORQZ(1))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     high X, X-comp.',SNGL(BFORCX(3)),SNGL(TORQX(3))
      WRITE(LUNGFO,*)
     &'     high X, Y-comp.',SNGL(BFORCY(3)),SNGL(TORQY(3))
      WRITE(LUNGFO,*)
     &'     high X, Z-comp.',SNGL(BFORCZ(3)),SNGL(TORQZ(3))

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     low Y, X-comp. ',SNGL(BFORCX(2)),SNGL(TORQX(2))
      WRITE(LUNGFO,*)
     &'     low Y, Y-comp. ',SNGL(BFORCY(2)),SNGL(TORQY(2))
      WRITE(LUNGFO,*)
     &'     low Y, Z-comp. ',SNGL(BFORCZ(2)),SNGL(TORQZ(2))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     high Y, X-comp.',SNGL(BFORCX(4)),SNGL(TORQX(4))
      WRITE(LUNGFO,*)
     &'     high Y, Y-comp.',SNGL(BFORCY(4)),SNGL(TORQY(4))
      WRITE(LUNGFO,*)
     &'     high Y, Z-comp.',SNGL(BFORCZ(4)),SNGL(TORQZ(4))

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     low Z, X-comp. ',SNGL(BFORCX(5)),SNGL(TORQX(5))
      WRITE(LUNGFO,*)
     &'     low Z, Y-comp. ',SNGL(BFORCY(5)),SNGL(TORQY(5))
      WRITE(LUNGFO,*)
     &'     low Z, Z-comp. ',SNGL(BFORCZ(5)),SNGL(TORQZ(5))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     high Z, X-comp.',SNGL(BFORCX(6)),SNGL(TORQX(6))
      WRITE(LUNGFO,*)
     &'     high Z, Y-comp.',SNGL(BFORCY(6)),SNGL(TORQY(6))
      WRITE(LUNGFO,*)
     &'     high Z, Z-comp.',SNGL(BFORCZ(6)),SNGL(TORQZ(6))

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     Sum X:         ',SNGL(BFORCX(7)),SNGL(TORQX(7))
      WRITE(LUNGFO,*)
     &'     Sum Y:         ',SNGL(BFORCY(7)),SNGL(TORQY(7))
      WRITE(LUNGFO,*)
     &'     Sum Z:         ',SNGL(BFORCZ(7)),SNGL(TORQZ(7))
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     Torque due to r x F:         '
      WRITE(LUNGFO,*)
     &'     ',SNGL(TORQROTX),SNGL(TORQROTY),SNGL(TORQROTZ)
      WRITE(LUNGFO,*)
     &'     Total torque:         '
      WRITE(LUNGFO,*)
     &'     ',SNGL(TORQTOTX),SNGL(TORQTOTY),SNGL(TORQTOTZ)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

        DEALLOCATE(XB)
        DEALLOCATE(YB)
        DEALLOCATE(ZB)
        DEALLOCATE(WS1)
        DEALLOCATE(WS2)
        DEALLOCATE(WS3)
        DEALLOCATE(WS4)
        DEALLOCATE(COEF)

        DEALLOCATE(BB)
        DEALLOCATE(BBI)
        DEALLOCATE(BT)
        DEALLOCATE(BTI)

      RETURN
      END
