*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.66/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.63/02 06/03/2008  15.55.08  by  Michael Scheer
*CMZ :  2.62/04 20/11/2007  09.47.54  by  Michael Scheer
*CMZ :  2.62/02 04/06/2007  13.30.31  by  Michael Scheer
*CMZ :  2.61/06 12/04/2007  09.54.44  by  Michael Scheer
*CMZ :  2.61/02 19/03/2007  15.22.59  by  Michael Scheer
*CMZ :  2.61/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.61/00 29/01/2007  14.25.16  by  Michael Scheer
*CMZ :  2.60/00 26/01/2007  10.42.26  by  Michael Scheer
*CMZ :  2.59/02 25/01/2007  14.58.31  by  Michael Scheer
*CMZ :  2.59/00 23/01/2007  15.38.35  by  Michael Scheer
*CMZ :  2.58/01 23/01/2007  13.16.29  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  11.26.00  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/12 04/10/96  17.37.39  by  Michael Scheer
*CMZ : 00.01/10 31/05/96  14.14.47  by  Michael Scheer
*CMZ : 00.01/09 20/05/96  14.42.08  by  Michael Scheer
*-- Author :    Michael Scheer   06/05/96

      SUBROUTINE BELLANA(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

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

c{ 7.4.2010:
c Adapted to version of POLYMAG
C NOT CLEAR, IF ALL COMMENTS ARE STILL CORRECT...
c} 7.4.2010:


C     V1=
C         Bi/2/4/KY * COS(KX*(X-X0)) * EXP(KY*Y) *  COS(KZ*Z-PHI/2)
C        +B0/1/4/KZ *                  EXP(KZ*Y) *  COS(KZ*Z-PHI/2)
C     V2=
C         Bi/2/4/KY * COS(KX*(X+X0)) * EXP(KY*Y) *  COS(KZ*Z+PHI/2)
C        +B0/1/4/KZ *                  EXP(KZ*Y) *  COS(KZ*Z+PHI/2)
C     V3=
C        -Bi/2/4/KY * COS(KX*(X+X0)) * EXP(-KY*Y) *  COS(KZ*Z-PHI/2)
C        -B0/1/4/KZ *                  EXP(-KZ*Y) *  COS(KZ*Z-PHI/2)
C     V4=
C        -Bi/2/4/KY * COS(KX*(X-X0)) * EXP(-KY*Y) *  COS(KZ*Z+PHI/2)
C        -B0/1/4/KZ *                  EXP(-KZ*Y) *  COS(KZ*Z+PHI/2)
C
C      X0 is distance of magnet center from device axis

C     V=V1+V2+V3+V4
C
C Noch berichtigen (Gode)
C     AX=-INT(DF(V,Y),Z)=INT(BY,Z)
C     AY=+INT(DF(V,X),Z)=-INT(BX,Z)
C     AZ=0
C
C
C     KY=SQRT(KX**2+KZ**2)
C
C     BX=-DF(V,X)
C     BY=-DF(V,Y)
C     BZ=-DF(V,Z)
C
C     BX1=
C        +KX/KY*Bi/2/4 * SIN(KX*(X-X0)) * EXP(KY*Y) *  COS(KZ*Z-PHI/2)
C     BY1=
C        -   Bi/2/4 * COS(KX*(X-X0)) * EXP(KY*Y) *  COS(KZ*Z-PHI/2)
C        -   B0/1/4 *                  EXP(KZ*Y) *  COS(KZ*Z-PHI/2)
C     BZ1=
C        +KZ/KY*Bi/2/4 * COS(KX*(X-X0)) * EXP(KY*Y) *  SIN(KZ*Z-PHI/2)
C        +      B0/1/4 *                  EXP(KZ*Y) *  SIN(KZ*Z-PHI/2)
C
C     BX2=
C        +KX/KY*Bi/2/4 * SIN(KX*(X+X0)) * EXP(KY*Y) *  COS(KZ*Z+PHI/2)
C     BY2=
C        -      Bi/2/4 * COS(KX*(X+X0)) * EXP(KY*Y) *  COS(KZ*Z+PHI/2)
C        -      B0/1/4 *                  EXP(KZ*Y) *  COS(KZ*Z+PHI/2)
C     BZ2=
C        +KZ/KY*Bi/2/4 * COS(KX*(X+X0)) * EXP(KY*Y) *  SIN(KZ*Z+PHI/2)
C        +      B0/1/4 *                  EXP(KZ*Y) *  SIN(KZ*Z+PHI/2)
C
C     BX3=
C        -KX/KY*Bi/2/4 * SIN(KX*(X+X0)) * EXP(-KY*Y) *  COS(KZ*Z-PHI/2)
C     BY3=
C        -      Bi/2/4 * COS(KX*(X+X0)) * EXP(-KY*Y) *  COS(KZ*Z-PHI/2)
C        -      B0/1/4 *                  EXP(-KZ*Y) *  COS(KZ*Z-PHI/2)
C     BZ3=
C        -KZ/KY*Bi/2/4 * COS(KX*(X+X0)) * EXP(-KY*Y) *  SIN(KZ*Z-PHI/2)
C        -      B0/1/4 *                  EXP(-KZ*Y) *  SIN(KZ*Z-PHI/2)
C
C     BX4=
C        -KX/KY*Bi/2/4 * SIN(KX*(X-X0)) * EXP(-KY*Y) *  COS(KZ*Z+PHI/2)
C     BY4=
C        -      Bi/2/4 * COS(KX*(X-X0)) * EXP(-KY*Y) *  COS(KZ*Z+PHI/2)
C        -      B0/1/4 *                  EXP(-KZ*Y) *  COS(KZ*Z+PHI/2)
C     BZ4=
C        -KZ/KY*Bi/2/4 * COS(KX*(X-X0)) * EXP(-KY*Y) *  SIN(KZ*Z+PHI/2)
C        -      B0/1/4 *                  EXP(-KZ*Y) *  SIN(KZ*Z+PHI/2)
C
C     AY1=
C        +KX/KY*Bi/2/4 * SIN(KX*(X-X0)) * EXP(KY*Y) * - SIN(KZ*Z-PHI/2)/KZ
C     AX1=
C        -   Bi/2/4 * COS(KX*(X-X0)) * EXP(KY*Y) *  SIN(KZ*Z-PHI/2)/KZ
C        -   B0/1/4 *                  EXP(KZ*Y) *  SIN(KZ*Z-PHI/2)/KZ
C
C     AY2=
C        +KX/KY*Bi/2/4 * SIN(KX*(X+X0)) * EXP(KY*Y) * - SIN(KZ*Z+PHI/2)/KZ
C     AX2=
C        -      Bi/2/4 * COS(KX*(X+X0)) * EXP(KY*Y) *  SIN(KZ*Z+PHI/2)/KZ
C        -      B0/1/4 *                  EXP(KZ*Y) *  SIN(KZ*Z+PHI/2)/KZ
C
C     AY3=
C        -KX/KY*Bi/2/4 * SIN(KX*(X+X0)) * EXP(-KY*Y) * - SIN(KZ*Z-PHI/2)/KZ
C
C     AX3=
C        -      Bi/2/4 * COS(KX*(X+X0)) * EXP(-KY*Y) *  SIN(KZ*Z-PHI/2)/KZ
C        -      B0/1/4 *                  EXP(-KZ*Y) *  SIN(KZ*Z-PHI/2)/KZ
C
C     AY4=
C        -KX/KY*Bi/2/4 * SIN(KX*(X-X0)) * EXP(-KY*Y) * - SIN(KZ*Z+PHI/2)/KZ
C     AX4=
C        -      Bi/2/4 * COS(KX*(X-X0)) * EXP(-KY*Y) *  SIN(KZ*Z+PHI/2)/KZ
C        -      B0/1/4 *                  EXP(-KZ*Y) *  SIN(KZ*Z+PHI/2)/KZ

      IMPLICIT NONE

      DOUBLE PRECISION BX,BY,BZ,AX,AY,AZ,X,Y,Z
     &  ,XIN,YIN,ZIN
     &  ,BXOUT,BYOUT,BZOUT
     &  ,AXOUT,AYOUT,AZOUT

      DOUBLE PRECISION XKXP,XKXM,YKY,PHI,PHI2,PHIROW,PHIROW2,BI2,ZKY
     &  ,BX1,BY1,BZ1
     &  ,BX2,BY2,BZ2
     &  ,BX3,BY3,BZ3
     &  ,BX4,BY4,BZ4
     &  ,AX1,AY1
     &  ,AX2,AY2
     &  ,AX3,AY3
     &  ,AX4,AY4
     &  ,COSXP,COSXM,EXPY,EXPZ
     &  ,SINXP,SINXM,EXPYM,EXPZM,ZKYK,XKYK,BRAD,PEL
     &  ,ZKZP1,ZKZP3,COSZP1,COSZP3,SINZP1,SINZP3
     &  ,ZKZP2,ZKZP4,COSZP2,COSZP4,SINZP2,SINZP4
     &  ,XKXI,YKI,HLEN100,HLEN75,BCONST
     &  ,GAPCOEF,YGAP,DEXPYGAP

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,ellana.
      include 'ellana.cmn'
*KEND.

      INTEGER ICAL,ICOEF,IDUM

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

        IF (XLELLANA.NE.0.D0) THEN
          XKELLANA=2.D0*PI1/XLELLANA
        ELSE
          XKELLANA=0.D0
        ENDIF

        IF (ZLELLANA.NE.0.D0) THEN
          ZKELLANA=2.D0*PI1/ZLELLANA
        ELSE
          ZKELLANA=0.D0
        ENDIF

        YKELLANA=DSQRT(XKELLANA**2+ZKELLANA**2)

        IF (YKELLANA.NE.0.D0) THEN
          YLELLANA=2.D0*PI1/YKELLANA
        ELSE
          YLELLANA=0.D0
        ENDIF

        PHI=SHELLANA*2.0D0*PI1
        PHI2=PHI/2.0D0
        PHIROW=ROWSHELLA*2.0D0*PI1
        PHIROW2=PHIROW/2.0D0
        PEL=EMASSG1*DSQRT( (DMYGAMMA+1.D0)*(DMYGAMMA-1.D0) )

        IF (B0ELLANA.NE.0.0D0) THEN
          BRAD=PEL/CLIGHT1*1.D9/B0ELLANA
        ELSE
          BRAD=0.0D0
        ENDIF

        IF (IELLCOEF.GT.0) THEN

          IF (IELLCOEF.GT.NELLCOEFP) THEN
            WRITE(6,*)'*** Error in BELLANA: Too many coefficients'
            WRITE(6,*)'*** i.e. more than',NELLCOEFP
            WRITE(6,*)'*** Program WAVE aborted ***'
            WRITE(LUNGFO,*)'*** Error in BELLANA: Too many coefficients'
            WRITE(LUNGFO,*)'*** i.e. more than',NELLCOEFP
            WRITE(LUNGFO,*)'*** Program WAVE aborted ***'
            STOP
          ENDIF

          OPEN(UNIT=99,FILE='bellana.coef',STATUS='OLD')

          read(99,*)gapcoef
          gapcoef=gapcoef/1000.0d0

          if (abs(gapcoef-refgapell).gt.1.0e-6) then
            print*,'*** Error in BELLANA: Bad reference gap on file bellana.coef'
            stop '*** Program WAVE aborted ***'
          endif

          DO ICOEF=1,IELLCOEF

            READ(99,*)IDUM,ELLCOEF(ICOEF)
            ELLCOEF(ICOEF)=ELLCOEF(ICOEF)*4.0D0
            IF (IDUM.NE.ICOEF-1) THEN
              WRITE(6,*)'*** Error in BELLANA: Bad numbering of coefficients'
              WRITE(6,*)'*** Program WAVE aborted ***'
              WRITE(LUNGFO,*)
     &          '*** Error in BELLANA: Bad numbering of coefficients'
              WRITE(LUNGFO,*)
     &          '*** Program WAVE aborted ***'
              STOP
            ENDIF
          ENDDO

          CLOSE(99)

        ELSE !IF (IELLCOEF.LE.0)

          ELLCOEF(1)=0.5D0
          ELLCOEF(2)=1.0D0

        ENDIF

        YGAP=(GAPELL-REFGAPELL)/2.0D0

        TLENELL=ZLELLANA*(NPERELLA+2)
        HLENELL=TLENELL/2.0D0
        HLEN75=HLENELL-ZLELLANA/2.0D0
        HLEN100=HLENELL-ZLELLANA

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     SR BELLANA, Parameter of elliptical undulator:'
        WRITE(LUNGFO,*)'     Field amplitude [T], bending radius[m], X0 [m]:  '
        WRITE(LUNGFO,*)'     ',SNGL(B0ELLANA),SNGL(BRAD),SNGL(X0ELLANA)
        WRITE(LUNGFO,*)'     Number of periods:  ',NPERELLA
        WRITE(LUNGFO,*)'     Half and total device length [m]:  '
        WRITE(LUNGFO,*)'     ',TLENELL,HLENELL
        WRITE(LUNGFO,*)'     Lx [m], Ly[m], Lz [m]:  '
        WRITE(LUNGFO,*)'     ',SNGL(XLELLANA),SNGL(YLELLANA),SNGL(ZLELLANA)
        WRITE(LUNGFO,*)'     Kx, Ky, Kz:  '
        WRITE(LUNGFO,*)'     ',SNGL(XKELLANA),SNGL(YKELLANA),SNGL(ZKELLANA)
        WRITE(LUNGFO,*)'     Shift and additional row shift in units of Lz:  '
        WRITE(LUNGFO,*)'     ',SNGL(SHELLANA),SNGL(ROWSHELLA)
        WRITE(LUNGFO,*)'     Gap and reference gap [m]:  '
        WRITE(LUNGFO,*)'     ',SNGL(GAPELL),SNGL(REFGAPELL)
        WRITE(LUNGFO,*)'     IELLS2S3:  ',IELLS2S3
        IF (IELLCOEF.LE.0) THEN
          WRITE(LUNGFO,*)'     Fourier coefficients:'
          WRITE(LUNGFO,*)'     0 0.125'
          WRITE(LUNGFO,*)'     1 0.25'
          IELLCOEF=2
        ELSE
          WRITE(LUNGFO,*)'     Fourier coefficients from file bellana.coef:'
          DO ICOEF=1,IELLCOEF
            WRITE(LUNGFO,*)'     ',icoef,ELLCOEF(ICOEF)
          ENDDO
        ENDIF
        WRITE(LUNGFO,*)

        WRITE(6,*)
        WRITE(6,*)'*** Warning in BELLANA: Vector potential not yet implemented'
        WRITE(6,*)

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** Warning in BELLANA: Vector potential not yet implemented'
        WRITE(LUNGFO,*)

        ICAL=1

      ENDIF !ICAL

      X=-ZIN
      Y=YIN
      Z=XIN

      IF (ABS(Z).GT.HLENELL) THEN
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        AXOUT=0.0D0
        AYOUT=0.0D0
        AZOUT=0.0D0
        RETURN
      ELSE IF (ABS(Z).GT.HLEN75) THEN
        BCONST=0.25D0
      ELSE IF (ABS(Z).GT.HLEN100) THEN
        BCONST=0.75D0
      ELSE
        BCONST=1.0D0
      ENDIF

      YKY=YKELLANA*Y
      ZKYK=ZKELLANA/YKELLANA
      ZKY=ZKELLANA*Y

      IF (IELLS2S3.GE.0) THEN
        ZKZP1=ZKELLANA*Z+PHI2
        ZKZP2=ZKELLANA*Z-PHI2
        ZKZP3=ZKELLANA*Z+PHI2+PHIROW
        ZKZP4=ZKELLANA*Z-PHI2+PHIROW
      ELSE
        ZKZP1=ZKELLANA*Z-PHI-PHIROW2
        ZKZP2=ZKELLANA*Z    -PHIROW2
        ZKZP3=ZKELLANA*Z+PHI+PHIROW2
        ZKZP4=ZKELLANA*Z    +PHIROW2
      ENDIF

      DEXPYGAP=DEXP(-ZKELLANA*YGAP)

      EXPZ=DEXP(ZKY)
      EXPZM=1.0D0/EXPZ
      EXPZ=EXPZ*DEXPYGAP
      EXPZM=EXPZM*DEXPYGAP

      COSZP1=DCOS(ZKZP1)
      SINZP1=DSIN(ZKZP1)
      COSZP3=DCOS(ZKZP3)
      SINZP3=DSIN(ZKZP3)

      COSZP2=DCOS(ZKZP2)
      SINZP2=DSIN(ZKZP2)
      COSZP4=DCOS(ZKZP4)
      SINZP4=DSIN(ZKZP4)

C first Fourier coefficient

      !X0ELLANA is distance of magnet center from device axis

      XKXP=XKELLANA*(X+X0ELLANA)
      XKXM=XKELLANA*(X-X0ELLANA)
      XKYK=XKELLANA/YKELLANA

      COSXP=DCOS(XKXP)
      COSXM=DCOS(XKXM)
      SINXP=DSIN(XKXP)
      SINXM=DSIN(XKXM)

      BI2=B0ELLANA/4.D0*ELLCOEF(1)

      BX1=0.0d0
      BY1=-     BI2 *         EXPZ * COSZP1
      BZ1=+     BI2 *         EXPZ * SINZP1

      BX2=0.0d0
      BY2=-     BI2 *         EXPZ * COSZP2
      BZ2=+     BI2 *         EXPZ * SINZP2

      BX3=0.0d0
      BY3=-     BI2 *         EXPZM * COSZP3
      BZ3=-     BI2 *         EXPZM * SINZP3

      BX4=0.0d0
      BY4=-     BI2 *         EXPZM * COSZP4
      BZ4=-     BI2 *         EXPZM * SINZP4

      AZ=0.0D0

      AX1=
     &  -   BI2 *                  EXPZ * SINZP1/ZKELLANA


      AY1=0.0D0

      AX2=
     &    -      BI2 *                  EXPZ * SINZP2/ZKELLANA
      AY2=0.0D0

      AX3=
     &  -      BI2 *                  EXPZM * SINZP3/ZKELLANA
      AY3=0.0D0

      AY4=0.0D0
      AX4=
     &    -      BI2 *                  EXPZM * SINZP4/ZKELLANA

      DO ICOEF=2,IELLCOEF

        XKXI=(ICOEF-1)*XKELLANA

        XKXP=XKXI*(X+X0ELLANA)
        XKXM=XKXI*(X-X0ELLANA)

        YKI=DSQRT(XKXI**2+ZKELLANA**2)

        YKY=YKI*Y
        DEXPYGAP=DEXP(-YKI*YGAP)
        EXPY=DEXP(YKY)
        EXPYM=1.0D0/EXPY
        EXPY=EXPY*DEXPYGAP
        EXPYM=EXPYM*DEXPYGAP

        ZKYK=ZKELLANA/YKI
        XKYK=XKXI/YKI

        COSXP=DCOS(XKXP)
        COSXM=DCOS(XKXM)
        SINXP=DSIN(XKXP)
        SINXM=DSIN(XKXM)

        BI2=B0ELLANA/8.D0*ELLCOEF(ICOEF)

        BX1=BX1+XKYK*BI2 * SINXM * EXPY * COSZP1
        BY1=BY1-     BI2 * COSXM * EXPY * COSZP1
        BZ1=BZ1+ZKYK*BI2 * COSXM * EXPY * SINZP1

        BX2=BX2+XKYK*BI2 * SINXP * EXPY * COSZP2
        BY2=BY2-     BI2 * COSXP * EXPY * COSZP2
        BZ2=BZ2+ZKYK*BI2 * COSXP * EXPY * SINZP2

        BX3=BX3-XKYK*BI2 * SINXP * EXPYM * COSZP3
        BY3=BY3-     BI2 * COSXP * EXPYM * COSZP3
        BZ3=BZ3-ZKYK*BI2 * COSXP * EXPYM * SINZP3

        BX4=BX4-XKYK*BI2 * SINXM * EXPYM * COSZP4
        BY4=BY4-     BI2 * COSXM * EXPYM * COSZP4
        BZ4=BZ4-ZKYK*BI2 * COSXM * EXPYM * SINZP4

        AY1=AY1
     &    +XKYK*BI2 * SINXM * EXPY * (-SINZP1)/ZKELLANA
        AX1=AX1
     &    -   BI2 * COSXM * EXPY * SINZP1/ZKELLANA

        AY2=AY2
     &    +XKYK*BI2 * SINXP * EXPY * (-SINZP2)/ZKELLANA
        AX2=AX2
     &    -      BI2 * COSXP * EXPY * SINZP2/ZKELLANA

        AY3=AY3
     &    -XKYK*BI2 * SINXP * EXPYM * (-SINZP3)/ZKELLANA

        AX3=AX3
     &    -      BI2 * COSXP * EXPYM * SINZP3/ZKELLANA

        AY4=AY4
     &    -XKYK*BI2 * SINXM * EXPYM * (-SINZP4)/ZKELLANA
        AX4=AX4
     &    -      BI2 * COSXM * EXPYM * SINZP4/ZKELLANA

      ENDDO !IELLCOEF-1

      BX=BX1+BX2+BX3+BX4
      BY=BY1+BY2+BY3+BY4
      BZ=BZ1+BZ2+BZ3+BZ4

      AX=AX+AX1+AX2+AX3+AX4
      AY=AY+AY1+AY2+AY3+AY4

      BXOUT=BZ*BCONST
      BYOUT=BY*BCONST
      BZOUT=-BX*BCONST

C Noch berichtigen (Gode)
      AXOUT=AZ*BCONST
      AYOUT=AY*BCONST
      AZOUT=-AX*BCONST

C     AX=-INT(DF(V,Y),Z)=INT(BY,Z)
C     AY=+INT(DF(V,X),Z)=-INT(BX,Z)
        AXOUT=0.0D0
        AYOUT=0.0D0
        AZOUT=0.0D0

      RETURN
      END
