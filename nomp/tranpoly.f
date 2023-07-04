*CMZ :  4.00/15 27/04/2022  08.13.16  by  Michael Scheer
*CMZ :  4.00/13 29/11/2021  14.07.42  by  Michael Scheer
*CMZ :  4.00/04 27/08/2019  11.49.27  by  Michael Scheer
*CMZ :  3.04/00 23/01/2018  13.05.43  by  Michael Scheer
*CMZ :  3.03/04 19/12/2017  13.15.36  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.12.35  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.43.01  by  Michael Scheer
*CMZ :  2.65/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.54/07 16/06/2005  12.04.22  by  Michael Scheer
*CMZ :  2.52/14 10/02/2005  10.17.13  by  Michael Scheer
*CMZ :  2.52/09 14/10/2004  16.20.51  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.42/03 13/09/2002  18.43.25  by  Michael Scheer
*CMZ :  2.41/13 22/08/2002  17.19.24  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.41/07 30/07/2002  11.46.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.34/08 12/09/2001  10.21.17  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  17.36.52  by  Michael Scheer
*CMZ :  2.20/01 04/09/2001  16.01.39  by  Michael Scheer
*CMZ :  2.16/08 19/10/2000  15.08.24  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 04/05/2000  15.22.24  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  2.12/00 03/06/99  15.29.39  by  Michael Scheer
*CMZ :  1.03/06 02/10/98  15.34.39  by  Michael Scheer
*CMZ : 00.02/00 19/11/96  14.56.00  by  Michael Scheer
*CMZ : 00.01/12 24/09/96  14.03.14  by  Michael Scheer
*CMZ : 00.01/10 02/09/96  16.27.29  by  Michael Scheer
*CMZ : 00.01/08 29/05/95  18.19.37  by  Michael Scheer
*CMZ : 00.01/07 10/03/95  11.25.51  by  Michael Scheer
*CMZ : 00.01/04 09/12/94  12.35.12  by  Michael Scheer
*CMZ : 00.01/03 28/11/94  09.53.54  by  Michael Scheer
*CMZ : 00.01/02 23/11/94  18.00.19  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.16.15  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.55.31  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.50  by  Michael Scheer
*-- Author : Michael Scheer
C************************************************************************
      SUBROUTINE TRANPOLY
C************************************************************************
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

C     Calculates Taylor expansion of the generating function by
c     fitting the coefficients from a set of particle trajectories

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,tranpo.
      include 'tranpo.cmn'
*KEEP,tralin.
      include 'tralin.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.

      CHARACTER(55) CODEERZ,CODEREF

      INTEGER NKOEF,IE,JE,KE,LE
     &         ,LUN1,LUNAKO,LUNERZ,LUNREF,LUNTRA,LUNSTR
     &         ,IREAD,MAXTRAP,NORDNGP,MKOEFV,NTOT
     &         ,I,J,K,L,II,JJ,KK,LL,IP,JP,KP,LP,IIP,JJP,KKP,LLP
     &         ,N,IERZ,IREF,ITRANC,IFAIL1,M,NP,I1,J1,L1,K1
     &         ,I1P,J1P,K1P,L1P,IWARNE,IWARNR,ICODEERZ,ICODEREF,
     &           ITRANCE,ITRANCR,kcode

      INTEGER IDETFAIL,IDUM,KTRAJ,LTRAJ

*KEEP,genfun.
      include 'genfun.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DIMENSION
     &  IE(MKOEF),JE(MKOEF),KE(MKOEF),LE(MKOEF)

      REAL*4 xran(1),rr(2)

      DOUBLE PRECISION
     &  BXI(MAXTRA),BYI(MAXTRA),BXF(MAXTRA),BYF(MAXTRA)
     & ,BZI(MAXTRA),BZF(MAXTRA),AZI(MAXTRA),AZF(MAXTRA)
     & ,AXI(MAXTRA),AYI(MAXTRA),AXF(MAXTRA),AYF(MAXTRA)
     & ,W2I(MAXTRA)

      DOUBLE PRECISION
     1     XI(4,MAXTRA),XF(4,MAXTRA),
     2     AKOEFF(NORDNG,NORDNG,NORDNG,NORDNG),
     2     AKOEFFC(NORDNG,NORDNG,NORDNG,NORDNG),
     *     AERZ(NORDNG,NORDNG,NORDNG,NORDNG),
     *     AREF(NORDNG,NORDNG,NORDNG,NORDNG),
     3     BVECT(MKOEF),BVECTC(MKOEF),BVECTCC(MKOEF),
     6     GMAT(MKOEF,MKOEF),GMATC(MKOEF,MKOEF),WORK(MKOEF),
     8     XE(MAXTRA,NORDNG+1),YE(MAXTRA,NORDNG+1),
     9     PE(MAXTRA,NORDNG+1),QE(MAXTRA,NORDNG+1),
     *     TRALIN(4,4),QUADR(4,4),DETCOP
     1    ,SI(MAXTRA),SF(MAXTRA)

      DOUBLE PRECISION
     &          GAMMA,BRHOABS,PEL,DUM,dgam
     &          ,X0,Y0,Z0,ZP0,YP0,BX0,BY0,BZ0,AX0,AY0,AZ0
     &          ,XFo,YFo,ZFo,ZPFo,YPFo,BXFo,BYFo,BZFo,AXFo,AYFo,AZFo
     &          ,XF0,YF0,ZF0,ZPF0,YPF0,BXF0,BYF0,BZF0,AXF0,AYF0,AZF0
     &         ,XICAVE,XFCAVE,XFXCAVE,XFYCAVE,XIXCAVE,XIYCAVE
     &         ,DIIP,DJJP,DKKP,DLLP,DI,DJ,DK,DL
     &          ,RESB2,RESBM,RESBAV,RESAM,RESA2,RESAAV,RESAMR,RESA2R
     &         ,RES,RESA,RESAR
     &          ,A1100,A2000,A0200,A0011,A0020,A0002,DETTRA1,DETTRA2
     &         ,ZLHAL,FY,FX,R2HAL,SINZ2,ZETAZ,ZLENGE,RHOHAL,B0HAL
     &         ,XLHAL,YLHAL,T11,T33,SIHPHIX,PHIX,QUADX,QUADFX,SINPHIX
     &         ,QUADY,QUADFY,SINPHIY,PHIY
     &         ,FOCEX,FOCEY,SIHPHIY
     &                  ,RESXMK, RESXAVK, RESXK,
     &                   RESYMK, RESYAVK, RESYK,
     &                   RESPXMK,RESPXAVK,RESPXK,
     &                   RESPYMK,RESPYAVK,RESPYK,RESAKO,
     &                   RESXME, RESXAVE, RESXE,
     &                   RESYME, RESYAVE, RESYE,
     &                   RESPXME,RESPXAVE,RESPXE,
     &                   RESPYME,RESPYAVE,RESPYE,RESERZ,
     &                   RESXMR, RESXAVR, RESXR,
     &                   RESYMR, RESYAVR, RESYR,
     &                   RESPXMR,RESPXAVR,RESPXR,
     &                   RESPYMR,RESPYAVR,RESPYR,RESREF,
     &                   RSXMTK, RSXAVTK, RSXTK,
     &                   RSYMTK, RSYAVTK, RSYTK,
     &                   RSPXMTK,RSPXAVTK,RSPXTK,
     &                   RSPYMTK,RSPYAVTK,RSPYTK,RESLIN,
     &                   RSXMT, RSXAVT, RSXT,
     &                   RSYMT, RSYAVT, RSYT,
     &                   RSPXMT,RSPXAVT,RSPXT,
     &                   RSPYMT,RSPYAVT,RSPYT,RESTRA,
     &                   RSXMTQ, RSXAVTQ, RSXTQ,
     &                   RSYMTQ, RSYAVTQ, RSYTQ,
     &                   RSPXMTQ,RSPXAVTQ,RSPXTQ,
     &                   RSPYMTQ,RSPYAVTQ,RSPYTQ,RESQUAD,
     &                   RESB,RESAAVR,RESR,AK,AE,AR,AKE,AKR

C     DOUBLE PRECISION XCLOSE,PCLOSE
C     DOUBLE PRECISION PXI0,PYI0,PXF0,PYF0

      DOUBLE PRECISION V0,V0X,V0Y,V0Z,VF0X,VF0Y,VF0Z,EWS(3),EWY(3),EWZ(3),EN
     &        ,ASR0,AYR0,AXR0,ASRF0,AYRF0,AXRF0
     &        ,VX1,VY1,VZ1,VX2,VY2,VZ2
     &        ,PXIR0,PYIR0,PXFR0,PYFR0,W0S,W0Y,W0X,WF0S,WF0Y,WF0X
     &        ,W0,WS1,WY1,WX1,WS2,WY2,WX2,ASR1,AYR1,AXR1
     &        ,ASR2,AYR2,AXR2,XPR1,YPR1,XPR2,YPR2
     &        ,X1,Y1,Z1,ZP1,YP1,X2,Y2,Z2,ZP2,YP2,AX1,AY1,AZ1,AX2,AY2,AZ2
     &        ,EWSF(3),EWYF(3),EWZF(3),SR1,YR1,XR1,SR2,YR2,XR2,sint

      DOUBLE PRECISION SINZ,COSZ

      DOUBLE PRECISION XKHAL,YKHAL,ZKHAL
     &  ,XKHALR,YKHALR,ZKHALR,XKHALI,YKHALI,ZKHALI

      double precision vx0,vxf0,vy0,vyf0,vz0,vzf0

*KEEP,ttracks.
      include 'ttracks.cmn'
*KEND.

*KEEP,pawcmn.
*KEND.

      DATA LUNAKO/15/
      DATA LUNERZ/13/
      DATA LUNREF/17/
      DATA LUNTRA/19/
      DATA LUNSTR/18/

      DATA IREAD/0/

      MAXTRAP=MAXTRA
      NORDNGP=NORDNG

C     PI=4.D0*DATAN(1.D0)

C--- READ AND WRITE CURRENT RUN NUMBER OF TRANPOLY

      OPEN(UNIT=LUNTRA,
     &      FILE='wave_trancode.dat',STATUS='OLD',
     &      FORM='FORMATTED',ERR=66)
          READ(LUNTRA,*,ERR=66)ITRANC
      CLOSE(LUNTRA)
66    ITRANC=ITRANC+1

      OPEN(UNIT=LUNTRA,FILE='wave_trancode.dat',
     &     STATUS='UNKNOWN',FORM='FORMATTED')
           WRITE(LUNTRA,*)ITRANC
C           IF (IGFLOAT.NE.0) ITRANC=-ITRANC
      CLOSE(LUNTRA)

      MKOEFV=MKOEF ! WEGEN UEBERGABE NACH SR DEQN

      LUN1=LUNGFO

      IF (LUN1.NE.LUNGFO)
     &OPEN(UNIT=LUN1,FILE='wave_tranpoly.out',STATUS='NEW',FORM='FORMATTED')

      WRITE (LUN1,*)
      WRITE (LUN1,*)
c      WRITE (LUN1,6789)
c6789  FORMAT(1H1)
      WRITE (LUN1,*) '     SR TRANPOLY'
      WRITE (LUN1,*) '     ==========='
      WRITE (LUN1,*)
      WRITE (LUN1,*)


      NKOEF=0
      DO I=1,NORDNG
      DO J=1,NORDNG
      DO K=1,NORDNG
      DO L=1,NORDNG

C---    THE SUBSEQUENT INDICES I,J,K,L RUN FROM 1...5,
C       THE CORRESPONDING MATH. INDICES i,j,k,l from 0...4

          IF (I-1 + J-1 + K-1 + L-1 .LT. NORDNG
     &      .AND.
     &      I-1 + J-1 + K-1 + L-1 .GT. 0
     &          .AND.

C         !MID-PLANE SYMMETRY
     &          ((K-1+L-1)/2*2.EQ.(K-1+L-1).OR.ISYM.EQ.0)
     &          .AND.

C         !HORI. FOCAL LENGTH = 0
     &          (IA11A20.EQ.0 .OR.
     &          (I-1.NE.2.OR.J-1.NE.0.OR.K-1.NE.0.OR.L-1.NE.0))  !A2000=0
     &          .AND.

C         !NOTE ALSO OTHER CHANGES DUE TO A1100=1
     &          (IA11A20.EQ.0 .OR.
     &          (I-1.NE.1.OR.J-1.NE.1.OR.K-1.NE.0.OR.L-1.NE.0))   !A1100=1
     &          .AND.

     &      (IA1000.EQ.0.OR.
     &          I-1.NE.1.OR.J-1.NE.0.OR.K-1.NE.0.OR.L-1.NE.0)  !A1000=0
     &          .AND.

     &      (IA0100.EQ.0.OR.
     &           I-1.NE.0.OR.J-1.NE.1.OR.K-1.NE.0.OR.L-1.NE.0) !A0100=0
     &          .AND.

     &      (IA0010.EQ.0.OR.
     &           I-1.NE.0.OR.J-1.NE.0.OR.K-1.NE.1.OR.L-1.NE.0) !A0010=0
     &          .AND.

     &      (IA0001.EQ.0.OR.
     &           I-1.NE.0.OR.J-1.NE.0.OR.K-1.NE.0.OR.L-1.NE.1) !A0001=0
     &          .AND.

C         !WLS IS HORZ. DRIFT
     &      (I-1+J-1.LT.3.OR.IWLSHOR.EQ.0)
     &          .AND.

     &      (K-1+L-1.EQ.0.OR.I2DIM.EQ.0)   !I2DIM CASE

     &        ) THEN

         NKOEF=NKOEF+1

       IF (NKOEF.GT.MKOEF) STOP
     &       '*** SR TRANPOLY: NKOEF.GT.MKOEF ***'

         IE(NKOEF)=I
         JE(NKOEF)=J
         KE(NKOEF)=K
         LE(NKOEF)=L

          ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO

C      OPEN(UNIT=LUNGFI,FILE=FILEI,STATUS='OLD',FORM='FORMATTED')

C--- LIES FILE MIT TRAJEKTORIEN

      OPEN(UNIT=LUNO,FILE=FILEO,STATUS='OLD') !,FORM='UNFORMATTED')

      read(luno,*)CODE
      read(luno,*)kCODE
      read(luno,*)GAMMA,dgam
      read(luno,*)NTOT

      V0=CLIGHT1*DSQRT((1.0D0-1.0D0/GAMMA)*(1.0D0+1.0D0/GAMMA))

C     WRITE(6,*)
C     WRITE(6,*) 'ANZAHL DER TRAJEKTORIEN AUF DATEN-FILE:',NTOT
C     WRITE(6,*)

      MTRAJ=NTOT

      IF (NTRAJ.EQ.-9999.OR.NTRAJ.EQ.0) THEN
        NTRAJ=NTOT
      ELSE IF (NTRAJ.EQ.-1) THEN
        NTRAJ=NKOEF
      ENDIF

      KTRAJ=NTOT/NTRAJ

      IF (MTRAJ.GT.NTOT)STOP '*** SR TRANPOLY: MTRAJ.GT.NTOT ***'

      IF (MTRAJ.GT.MAXTRA) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN TRANPOLY ***'
          WRITE(LUNGFO,*)'TOO MANY TRACKS ON FILE FILEO'
          WRITE(LUNGFO,*)'INCREASE PARAMETER MAXTRA IN FILE GENFUN.CMN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN TRANPOLY ***'
          WRITE(6,*)'TOO MANY TRACKS ON FILE FILEO'
          WRITE(6,*)'INCREASE PARAMETER MAXTRA IN FILE GENFUN.CMN'
          WRITE(6,*)
          STOP
      ENDIF

      read(luno,*)X0,Y0,Z0,ZP0,YP0,BX0,BY0,BZ0,AX0,AY0,AZ0
c18.1.2019      read(luno,*)XF0,YF0,ZF0,ZPF0,YPF0,BXF0,BYF0,BZF0,AXF0,AYF0,AZF0
c18.1.2018{
      read(luno,*)XFo,YFo,ZFo,ZPFo,YPFo,BXFo,BYFo,BZFo,AXFo,AYFo,AZFo
      read(luno,*)xf0,yf0,zf0,zpf0,ypf0,axf0,ayf0,azf0
      vxf0=v0/sqrt(1.0d0+zpf0*zpf0+ypf0*ypf0)
      vyf0=ypf0*vxf0
      vzf0=zpf0*vxf0
c18.1.2018}
      read(luno,*)OPNX,OPNY,OPNZ
      read(luno,*)OPNFX,OPNFY,OPNFZ
      read(luno,*)ZAPERT,YAPERT,DLAPER,sint

      DO J=1,NTOT
            read(luno,*) SI(J),XI(1,J),XI(2,J),XI(3,J),XI(4,J),
     &        BXI(J),BYI(J),BZI(J),AXI(J),AYI(J),AZI(J),
c            read(luno,*)
     &        SF(J),XF(1,J),XF(2,J),XF(3,J),XF(4,J),
     &        BXF(J),BYF(J),BZF(J),AXF(J),AYF(J),AZF(J),
     &        W2I(J)
      ENDDO

      CLOSE(LUNO)

C--- SUPERIMPOSE NOISE TO SIMULATE TRACKING UNCERTAINTIES

      IF (DRAUSCHX.NE.0..OR.DRAUSCHY.NE.0.) THEN

         DO I=1,NTOT
             call util_random_gauss(1,xran,rr)
             XF(1,I)=XF(1,I)+xran(1)*DRAUSCHX
             call util_random_gauss(1,xran,rr)
             XF(2,I)=XF(2,I)+xran(1)*DRAUSCHX
             call util_random_gauss(1,xran,rr)
             XF(3,I)=XF(3,I)+xran(1)*DRAUSCHY
             call util_random_gauss(1,xran,rr)
             XF(4,I)=XF(4,I)+xran(1)*DRAUSCHY
         ENDDO

      ENDIF

      IF (I2DIM.NE.0) THEN

         DO I=1,NTOT
             IF(    XI(3,I).NE.0.
     &            .OR. XI(4,I).NE.0.
     &            .OR. XF(3,I).NE.0.
     &            .OR. XF(4,I).NE.0. ) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &              '*** ERROR IN TRANPOLY ***'
            WRITE(LUNGFO,*)
     &              'FLAG I2DIM IS SET BUT TRACKS ARE NOT PLANAR'
            WRITE(LUNGFO,*)
     &              'SOMETHING IS WRONG WITH THE INPUT FILE FILEO'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)
     &              '*** ERROR IN TRANPOLY ***'
            WRITE(6,*)
     &              'FLAG I2DIM IS SET BUT TRACKS ARE NOT PLANAR'
            WRITE(6,*)
     &              'SOMETHING IS WRONG WITH THE INPUT FILE FILEO'
            WRITE(6,*)
                        STOP
                  ENDIF
         ENDDO

      ENDIF ! I2DIM

C--------------------------------------------------------------

      W0=V0
      PEL=EMASSE1*DSQRT( (GAMMA+1.D0)*(GAMMA-1.D0) )
      BRHOABS=PEL/CLIGHT1  !ABSOLUTE VALUE

C--- REFERENCE SYSTEM AT ENTRANCE

C     IF THE ENTRANCE PLANE OF THE DEVICE IS PERPENDICULAR TO THE REFERENCE
C       ORBIT, THEN THE REFERENCE ORBIT HAS ZERO SLOPE
C     AND DISPLACEMENT IN THIS SYSTEM AND ALL VARIABLES ARE RELATIVE TO THE
C     REFERENCE ORBIT. IF IN THE OTHER CASE THE ENTRANCE PLANE IS PERPEN-
C     DICULAR TO THE LONGITUDINAL AXIS OF THE DEVICE, THE REFERENCE SYSTEM
C     IS IDENTICAL TO THE LAB.-SYSTEM. BOTH CASES ARE DESCRIBED BY THE SAME
C     FORMULAS IF WE SUBTRACT THE REFERENCE ORBIT FROM ALL VARIABLES.
C     THE CLOSED ORBIT IS THEN TAKEN OUT AND ALL CLOSED ORBIT TERMES SHOULD
C     VANISH. THE SAME IS VALID FOR THE EXIT PLANE.


C     EWS HAS DIFFERENT MEANING THAN IN SR OPTIC

      EWS(1)=OPNX
      EWS(2)=OPNY
      EWS(3)=OPNZ

      EN=1.D0/DSQRT(EWS(3)*EWS(3)+EWS(1)*EWS(1))

C        UNIT-VECTOR EWZ~[EWS,(0,1,0)] (CROSS-PRODUCT)
      EWZ(1)=-EWS(3)*EN
      EWZ(2)= 0.
      EWZ(3)= EWS(1)*EN

C     UNIT-VECTOR EWY=[EWZ,EWS]

      EWY(1)= EWZ(2)*EWS(3) - EWZ(3)*EWS(2)
      EWY(2)= EWZ(3)*EWS(1) - EWZ(1)*EWS(3)
      EWY(3)= EWZ(1)*EWS(2) - EWZ(2)*EWS(1)

      EWSF(1)=OPNFX
      EWSF(2)=OPNFY
      EWSF(3)=OPNFZ

      EN=1.D0/DSQRT(EWSF(3)*EWSF(3)+EWSF(1)*EWSF(1))
      EWZF(1)=-EWSF(3)*EN
      EWZF(2)= 0.
      EWZF(3)= EWSF(1)*EN

      EWYF(1)= EWZF(2)*EWSF(3) - EWZF(3)*EWSF(2)
      EWYF(2)= EWZF(3)*EWSF(1) - EWZF(1)*EWSF(3)
      EWYF(3)= EWZF(1)*EWSF(2) - EWZF(2)*EWSF(1)

C--- NORMALIZED VELOCITIES OF THE REFERENCE ORBIT

      V0X=1.D0/DSQRT(1.D0+ZP0*ZP0+YP0*YP0)
      V0Y=YP0*V0X
      V0Z=ZP0*V0X

      VF0X=1.D0/DSQRT(1.D0+ZPF0*ZPF0+YPF0*YPF0)
      VF0Y=YPF0*VF0X
      VF0Z=ZPF0*VF0X

C--- TRANSFORM VELOCITIES OF REFERENCE ORBIT INTO REFERENCE SYSTEM

C     4.9.2001: (W0X,W0Y,W0S)=(V0X,V0Y,V0Z) OR
C     (W0X,W0Y,W0S)=(0,0,1)
C     CODE NOT CHANGED, ONLY INTERPRETION HAS CHANGED

C     FOR THE TIME BEING (7.1.1992) (W0S,W0Y,W0X)=(V0X,V0Y,V0Z) OR
C     (W0S,W0Y,W0X)=(1,0,0)

      W0S=V0X*EWS(1)+V0Y*EWS(2)+V0Z*EWS(3)
      W0Y=V0X*EWY(1)+V0Y*EWY(2)+V0Z*EWY(3)
      W0X=V0X*EWZ(1)+V0Y*EWZ(2)+V0Z*EWZ(3)

      WF0S=VF0X*EWSF(1)+VF0Y*EWSF(2)+VF0Z*EWSF(3)
      WF0Y=VF0X*EWYF(1)+VF0Y*EWYF(2)+VF0Z*EWYF(3)
      WF0X=VF0X*EWZF(1)+VF0Y*EWZF(2)+VF0Z*EWZF(3)

C--- TRANSFORM VECTOR POTENTIAL OF REFERENCE ORBIT INTO REFERENCE SYSTEM

      ASR0=AX0*EWS(1)+AY0*EWS(2)+AZ0*EWS(3)
      AYR0=AX0*EWY(1)+AY0*EWY(2)+AZ0*EWY(3)
      AXR0=AX0*EWZ(1)+AY0*EWZ(2)+AZ0*EWZ(3)

      ASRF0=AXF0*EWSF(1)+AYF0*EWSF(2)+AZF0*EWSF(3)
      AYRF0=AXF0*EWYF(1)+AYF0*EWYF(2)+AZF0*EWYF(3)
      AXRF0=AXF0*EWZF(1)+AYF0*EWZF(2)+AZF0*EWZF(3)

C--- CANONICAL VARIABLES OF REFERENCE ORBIT

C     PXI0=-AZ0/BRHOABS   !COORDINATE CHANGE
C     PYI0=-AY0/BRHOABS

C     PXF0=-AZF0/BRHOABS
C     PYF0=-AYF0/BRHOABS

      IF (ICHARGE.GT.0) THEN
            PXIR0=+AXR0/BRHOABS+W0X !  THIS IS A VECTOR COMPONENT !
            PYIR0=+AYR0/BRHOABS+W0Y
            PXFR0=+AXRF0/BRHOABS+WF0X
            PYFR0=+AYRF0/BRHOABS+WF0Y
      ELSE
            PXIR0=-AXR0/BRHOABS+W0X !  THIS IS A VECTOR COMPONENT !
            PYIR0=-AYR0/BRHOABS+W0Y
            PXFR0=-AXRF0/BRHOABS+WF0X
            PYFR0=-AYRF0/BRHOABS+WF0Y
      ENDIF

C--- LOOP OVER TRAJECTORIES -------------------------------------------

      DO I=1,MTRAJ

C--- EVERYTHING IN THE LAB.SYSTEM

         X1= SI(I)
         Z1= XI(1,I)
         ZP1=XI(2,I)
         Y1= XI(3,I)
         YP1=XI(4,I)

         VX1=1.D0/DSQRT(1.D0+ZP1**2+YP1**2)
         VY1=YP1*VX1
         VZ1=ZP1*VX1

         AX1=AXI(I)
         AY1=AYI(I)
         AZ1=AZI(I)

         X2= SF(I)
         Z2= XF(1,I)
         ZP2=XF(2,I)
         Y2= XF(3,I)
         YP2=XF(4,I)

         VX2=1.D0/DSQRT(1.D0+ZP2**2+YP2**2)
         VY2=YP2*VX2
         VZ2=ZP2*VX2

         AX2=AXF(I)
         AY2=AYF(I)
         AZ2=AZF(I)

C--- TRANSFORM VELOCITIES INTO REFERENCE SYSTEM

         WS1=VX1*EWS(1)+VY1*EWS(2)+VZ1*EWS(3)
         WY1=VX1*EWY(1)+VY1*EWY(2)+VZ1*EWY(3)
         WX1=VX1*EWZ(1)+VY1*EWZ(2)+VZ1*EWZ(3)

         WS2=VX2*EWSF(1)+VY2*EWSF(2)+VZ2*EWSF(3)
         WY2=VX2*EWYF(1)+VY2*EWYF(2)+VZ2*EWYF(3)
         WX2=VX2*EWZF(1)+VY2*EWZF(2)+VZ2*EWZF(3)

C--- TRANSFORM VECTOR POTENTIALS INTO REFERENCE SYSTEM

         ASR1=AX1*EWS(1)+AY1*EWS(2)+AZ1*EWS(3)
         AYR1=AX1*EWY(1)+AY1*EWY(2)+AZ1*EWY(3)
         AXR1=AX1*EWZ(1)+AY1*EWZ(2)+AZ1*EWZ(3)

         ASR2=AX2*EWSF(1)+AY2*EWSF(2)+AZ2*EWSF(3)
         AYR2=AX2*EWYF(1)+AY2*EWYF(2)+AZ2*EWYF(3)
         AXR2=AX2*EWZF(1)+AY2*EWZF(2)+AZ2*EWZF(3)


C--- SLOPES IN THE REFERENCE SYSTEM

         XPR1=WX1/WS1
         YPR1=WY1/WS1

         XPR2=WX2/WS2
         YPR2=WY2/WS2

C--- TRANSFORM COORDINATES INTO REFERENCE SYSTEM

         SR1=(X1-X0)*EWS(1)+(Y1-Y0)*EWS(2)+(Z1-Z0)*EWS(3)
         YR1=(X1-X0)*EWY(1)+(Y1-Y0)*EWY(2)+(Z1-Z0)*EWY(3)
         XR1=(X1-X0)*EWZ(1)+(Y1-Y0)*EWZ(2)+(Z1-Z0)*EWZ(3)

         SR2=(X2-XF0)*EWSF(1)+(Y2-YF0)*EWSF(2)+(Z2-ZF0)*EWSF(3)
         YR2=(X2-XF0)*EWYF(1)+(Y2-YF0)*EWYF(2)+(Z2-ZF0)*EWYF(3)
         XR2=(X2-XF0)*EWZF(1)+(Y2-YF0)*EWZF(2)+(Z2-ZF0)*EWZF(3)

C--- CANONICAL VARIABLES  (REFERENCE ORBIT IS SUBTRACTED)

         XIC(I)=XR1
         YIC(I)=YR1
         XFC(I)=XR2
         YFC(I)=YR2

      IF (ICHARGE.GT.0) THEN
         PXI(I)=+AXR1/BRHOABS+WX1-PXIR0
         PYI(I)=+AYR1/BRHOABS+WY1-PYIR0
         PXF(I)=+AXR2/BRHOABS+WX2-PXFR0
         PYF(I)=+AYR2/BRHOABS+WY2-PYFR0
      ELSE
         PXI(I)=-AXR1/BRHOABS+WX1-PXIR0
         PYI(I)=-AYR1/BRHOABS+WY1-PYIR0
         PXF(I)=-AXR2/BRHOABS+WX2-PXFR0
         PYF(I)=-AYR2/BRHOABS+WY2-PYFR0
      ENDIF

      ENDDO

C--------------------------------------------------------------
C240991 MEAN LENGTH OF VECTORS

      XICAVE=0.D0
      XFCAVE=0.D0
      XIXCAVE=0.D0
      XIYCAVE=0.D0
      XFXCAVE=0.D0
      XFYCAVE=0.D0
      DO I=1,MTRAJ
         XICAVE=XICAVE+DSQRT(XIC(I)**2+YIC(I)**2+PXF(I)**2+PYF(I)**2)
         XFCAVE=XFCAVE+DSQRT(XFC(I)**2+YFC(I)**2+PXI(I)**2+PYI(I)**2)
         XIXCAVE=XIXCAVE+DSQRT(XIC(I)**2+PXF(I)**2)
         XIYCAVE=XIYCAVE+DSQRT(YIC(I)**2+PYF(I)**2)
         XFXCAVE=XFXCAVE+DSQRT(XFC(I)**2+PXI(I)**2)
         XFYCAVE=XFYCAVE+DSQRT(YFC(I)**2+PYI(I)**2)
      ENDDO
      XICAVE=XICAVE/MTRAJ
      XFCAVE=XFCAVE/MTRAJ
      XIXCAVE=XIXCAVE/MTRAJ
      XIYCAVE=XIYCAVE/MTRAJ
      XFXCAVE=XFXCAVE/MTRAJ
      XFYCAVE=XFYCAVE/MTRAJ

C--------------------------------------------------------------
C230891
      IF(IA11A20.NE.0) THEN
         DO M=1,MTRAJ         ! i.e. A(1,1,0,0) = 1
             XFC(M)=XFC(M)-XIC(M)   !      A(2,0,0,0) = 0
             PXI(M)=PXI(M)-PXF(M)
         ENDDO
      ENDIF
C--------------------------------------------------------------
C     CHANGE SCALE

          X0=X0*DSCALE
          Y0=Y0*DSCALE
          Z0=Z0*DSCALE
          XF0=XF0*DSCALE
          YF0=YF0*DSCALE
          ZF0=ZF0*DSCALE

      DO M=1,MTRAJ
          XIC(M)=XIC(M)*DSCALE
          XFC(M)=XFC(M)*DSCALE
          PXI(M)=PXI(M)*DSCALE
          PXF(M)=PXF(M)*DSCALE
          YIC(M)=YIC(M)*DSCALE
          YFC(M)=YFC(M)*DSCALE
          PYI(M)=PYI(M)*DSCALE
          PYF(M)=PYF(M)*DSCALE
      ENDDO

C--------------------------------------------------------------
C     SET UP SYSTEM OF EQUATIONS

C--- INDEX POINTER IE(K1,K2,K3,K4), JE(K1,....

C---
      IF (MTRAJ.LT.NKOEF) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** ERROR IN TRANPOLY ***'
         WRITE(LUNGFO,*)'NOT ENOUGH TRACKS TO FIT COEFFICIENTS'
         WRITE(LUNGFO,*)'NUMBER OF TRACKS:      ',MTRAJ
         WRITE(LUNGFO,*)'NUMBER OF COEFFICIENTS:',NKOEF
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** ERROR IN TRANPOLY ***'
         WRITE(6,*)'NOT ENOUGH TRACKS TO FIT COEFFICIENTS'
         WRITE(6,*)'NUMBER OF TRACKS:      ',MTRAJ
         WRITE(6,*)'NUMBER OF COEFFICIENTS:',NKOEF
         WRITE(6,*)
         STOP
      ENDIF

C--- POWER TERMS FOR MATRIX

      LTRAJ=0
      DO M=1,NTOT,KTRAJ

        LTRAJ=LTRAJ+1

        IF (M.GT.NTRAJ) GOTO 1122

         XE(M,1)=0.D0         ! x**(i-1) in math. formulas
         YE(M,1)=0.D0
         PE(M,1)=0.D0
         QE(M,1)=0.D0

         XE(M,2)=1.D0         ! x**0 in math. formulas
         YE(M,2)=1.D0
         PE(M,2)=1.D0
         QE(M,2)=1.D0

         DO I=3,NORDNG

            XE(M,I)=XE(M,I-1)*XIC(M)  ! x**i in math. formulas
            YE(M,I)=YE(M,I-1)*YIC(M)
            PE(M,I)=PE(M,I-1)*PXF(M)
            QE(M,I)=QE(M,I-1)*PYF(M)

         ENDDO

         XE(M,NORDNG+1)=0.D0  ! just formal to make GMAT loop run
         YE(M,NORDNG+1)=0.D0
         PE(M,NORDNG+1)=0.D0
         QE(M,NORDNG+1)=0.D0

1122     CONTINUE

      ENDDO

C--- MATRIX GMAT

      DO NP=1,NKOEF
      DO N=1,NKOEF

        I=IE(N)
        J=JE(N)
        K=KE(N)
        L=LE(N)

        IP=IE(NP)
        JP=JE(NP)
        KP=KE(NP)
        LP=LE(NP)

        DIIP=DFLOAT((I-1)*(IP-1))   ! 0*0 <= DIIP <= 4*4
        DJJP=DFLOAT((J-1)*(JP-1))
        DKKP=DFLOAT((K-1)*(KP-1))
        DLLP=DFLOAT((L-1)*(LP-1))

        II=I+1
        JJ=J+1
        KK=K+1
        LL=L+1

        I1=II-1
        J1=JJ-1
        K1=KK-1
        L1=LL-1

        IIP=IP+1
        JJP=JP+1
        KKP=KP+1
        LLP=LP+1

        I1P=IIP-1
        J1P=JJP-1
        K1P=KKP-1
        L1P=LLP-1

        DO M=1,MTRAJ

          GMAT(NP,N)=GMAT(NP,N)

     &       + DIIP * XE(M,I1)  * PE(M,JJ)  * YE(M,KK)  * QE(M,LL)
     &      * XE(M,I1P) * PE(M,JJP) * YE(M,KKP) * QE(M,LLP)
     &       + DJJP * XE(M,II)  * PE(M,J1)  * YE(M,KK)  * QE(M,LL)
     &         * XE(M,IIP) * PE(M,J1P) * YE(M,KKP) * QE(M,LLP)
     &       + DKKP * XE(M,II)  * PE(M,JJ)  * YE(M,K1)  * QE(M,LL)
     &         * XE(M,IIP) * PE(M,JJP) * YE(M,K1P) * QE(M,LLP)
     &       + DLLP * XE(M,II)  * PE(M,JJ)  * YE(M,KK)  * QE(M,L1)
     &         * XE(M,IIP) * PE(M,JJP) * YE(M,KKP) * QE(M,L1P)

        ENDDO

      ENDDO
      ENDDO

C--- INHOMOGENEITY BVECT OF EQUATION SYSTEM

      DO N=1,NKOEF

        I=IE(N)
        J=JE(N)
        K=KE(N)
        L=LE(N)

        DI=DFLOAT(I-1)
        DJ=DFLOAT(J-1)
        DK=DFLOAT(K-1)
        DL=DFLOAT(L-1)

        II=I+1
        JJ=J+1
        KK=K+1
        LL=L+1

        I1=II-1
        J1=JJ-1
        K1=KK-1
        L1=LL-1

        DO M=1,MTRAJ

          BVECT(N)=BVECT(N)
     &   + DI * PXI(M)  * XE(M,I1) * PE(M,JJ) * YE(M,KK) * QE(M,LL)
     &   + DJ * XFC(M)  * XE(M,II) * PE(M,J1) * YE(M,KK) * QE(M,LL)
     &   + DK * PYI(M)  * XE(M,II) * PE(M,JJ) * YE(M,K1) * QE(M,LL)
     &   + DL * YFC(M)  * XE(M,II) * PE(M,JJ) * YE(M,KK) * QE(M,L1)

        ENDDO

      ENDDO

C--- SAVE SYSTEM OF EQUATIONS

      DO I=1,MKOEFV
      DO J=1,MKOEFV
         GMATC(I,J)=GMAT(I,J)
      ENDDO
         BVECTC(I)=BVECT(I)
      ENDDO

C---  SOLVE SYSTEM OF EQUATION USING CERN ROUTINE F010

      CALL DEQN(NKOEF,GMAT,MKOEFV,WORK,IFAIL1,1,BVECT)
      IF (IFAIL1.NE.0) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** ERROR IN TRANPOLY ***'
         WRITE(LUNGFO,*)
     & 'CERN ROUTINE F010 TO SOLVE SYSTEM OF EQUATIONS FAILED'
         WRITE(LUNGFO,*)
     & 'PERHAPS TOO MANY COEFFICIENTS TO FIT ??'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)'*** ERROR IN TRANPOLY ***'
         WRITE(6,*)
     & 'CERN ROUTINE F010 TO SOLVE SYSTEM OF EQUATIONS FAILED'
         WRITE(6,*)
     & 'PERHAPS TOO MANY COEFFICIENTS TO FIT ??'
         WRITE(6,*)
         STOP
      ENDIF

C--- GET RESULTS

      DO N=1,NKOEF

      I=IE(N)
      J=JE(N)
      K=KE(N)
      L=LE(N)

      AKOEFF(I,J,K,L)=BVECT(N)

      ENDDO

C--- CHECK SOLUTION OF THE EQUATION SYSTEM

      RESB2=0.D0
      RESBM=0.D0
      RESBAV=0.D0

      DO I=1,NKOEF
         BVECTCC(I)=0.D0
         DO J=1,NKOEF
         BVECTCC(I)=BVECTCC(I) +
     &            GMATC(I,J) * BVECT(J)
         ENDDO
         RES=BVECTCC(I)-BVECTC(I)
         IF(DABS(RES).GT.DABS(RESBM)) RESBM = RES
         RESBAV=RESBAV+RES
         RESB2 = RESB2 + RES*RES
      ENDDO
      RESB=DSQRT(RESB2/NKOEF)
      RESBAV=RESBAV/NKOEF

C--------------------------------------------------------------
C     RESCALING

      DO I=0,NORDNG-1
      DO J=0,NORDNG-1
      DO K=0,NORDNG-1
      DO L=0,NORDNG-1
      IF (I+J+K+L.LE.NORDNG-1)
     &   AKOEFF(I+1,J+1,K+1,L+1)=
     &   AKOEFF(I+1,J+1,K+1,L+1)*DSCALE**(I+K-1+J+L-1)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

          X0=X0/DSCALE
          Y0=Y0/DSCALE
          Z0=Z0/DSCALE
          XF0=XF0/DSCALE
          YF0=YF0/DSCALE
          ZF0=ZF0/DSCALE

      DO M=1,MTRAJ
          XIC(M)=XIC(M)/DSCALE
          XFC(M)=XFC(M)/DSCALE
          PXI(M)=PXI(M)/DSCALE
          PXF(M)=PXF(M)/DSCALE
          YIC(M)=YIC(M)/DSCALE
          YFC(M)=YFC(M)/DSCALE
          PYI(M)=PYI(M)/DSCALE
          PYF(M)=PYF(M)/DSCALE
      ENDDO

C230891
      IF(IA11A20.NE.0) THEN
        AKOEFF(2,2,1,1)=1.0D0 ! ADD A1100=1
        DO M=1,MTRAJ
          XFC(M)=XFC(M)+AKOEFF(2,2,1,1)*XIC(M)
          PXI(M)=PXI(M)+AKOEFF(2,2,1,1)*PXF(M)
        ENDDO
        NKOEF=NKOEF+1
        IE(NKOEF)=2
        JE(NKOEF)=2
        KE(NKOEF)=1
        LE(NKOEF)=1
      ENDIF

C--- COMPARE CALCULATED AND GIVEN COEFFICIENTS A(I,J,K,L)

      OPEN(UNIT=LUNERZ,FILE='wave_erzfun.in',STATUS='OLD',ERR=99)

C10.6.93 READ(LUNERZ,'(I10,1A60)')ICODEERZ,CODEERZ

      READ(LUNERZ,'(2I5,1A55)')ICODEERZ,ITRANCE,CODEERZ
          read(lunerz,*)idum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum
          read(lunerz,*)dum

      DO IREAD=1,1000000

         READ (LUNERZ,*, END=11) I,J,K,L,DUM

         IF(I+J+K+L.LE.NORDNG-1) THEN
            AERZ(I+1,J+1,K+1,L+1)=DUM
         ELSE
            IWARNE=1
         ENDIF

      END DO

11    IERZ=IREAD-1
99    CONTINUE
      IREAD=0
C     IF (IERZ.NE.0) IKOEFF=1
      CLOSE(LUNERZ)

      OPEN(UNIT=LUNREF,FILE='wave_erzfun.ref',STATUS='OLD',ERR=88)

C10.6.93 READ(LUNREF,'(I10,1A60)')ICODEREF,CODEREF

      READ(LUNREF,'(2I5,1A55)')ICODEREF,ITRANCR,CODEREF
          read(lunref,*)idum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum
          read(lunref,*)dum

C260891  READ(LUNREF,*)

C10.6.93 DO IREAD=1,NKOEF

      DO IREAD=1,1000000

        READ (LUNREF,*, END=22) I,J,K,L,DUM

        IF(I+J+K+L.LE.NORDNG-1) THEN
          AREF(I+1,J+1,K+1,L+1)=DUM
        ELSE
          IWARNR=1
        ENDIF
      END DO

22    IREF=IREAD-1
88    CONTINUE
C     IF (IREF.NE.0) IKOEFF=1
      CLOSE(LUNREF)

      RESAM=0.D0
      RESA2=0.D0
      RESAAV=0.D0
      RESAMR=0.D0
      RESA2R=0.D0
      RESAAVR=0.D0


      IF (IERZ.GT.0) THEN
        DO I=1,NORDNG
          DO J=1,NORDNG
            DO K=1,NORDNG
              DO L=1,NORDNG

                RES=AKOEFF(I,J,K,L)-AERZ(I,J,K,L)
                IF (DABS(RES).GT.DABS(RESAM)) RESAM=RES
                RESA2=RESA2+RES*RES
                RESAAV=RESAAV+RES
                RESR=AKOEFF(I,J,K,L)-AREF(I,J,K,L)
                IF (DABS(RESR).GT.DABS(RESAMR)) RESAMR=RESR
                RESA2R=RESA2R+RESR*RESR
                RESAAVR=RESAAVR+RESR

              ENDDO
            ENDDO
          ENDDO
        ENDDO

        RESA=DSQRT(RESA2/NKOEF)
        RESAAV=RESAAV/NKOEF
        RESAR=DSQRT(RESA2R/NKOEF)
        RESAAVR=RESAAVR/NKOEF
      ENDIF

C--- CALCULATE TRACKS FROM FITTED GENERATING FUNCTION (AKOEFF) AND COMPARE
C    WITH ORIGINAL TRACKING RESULTS

      CALL ARESI(      AKOEFF,MAXTRAP,NORDNGP,
     &                   RESXMK, RESXAVK, RESXK,
     &                   RESYMK, RESYAVK, RESYK,
     &                   RESPXMK,RESPXAVK,RESPXK,
     &                   RESPYMK,RESPYAVK,RESPYK,RESAKO)

C--- CALCULATE TRACKS FROM GIVEN GENERATING FUNCTION (AERZ) AND COMPARE
C    WITH ORIGINAL TRACKING RESULTS

      IF (IERZ.GT.0)
     &  CALL ARESI(  AERZ,MAXTRAP,NORDNGP,
     &  RESXME, RESXAVE, RESXE,
     &  RESYME, RESYAVE, RESYE,
     &  RESPXME,RESPXAVE,RESPXE,
     &  RESPYME,RESPYAVE,RESPYE,RESERZ)

C--- CALCULATE TRACKS FROM GIVEN GENERATING FUNCTION (AREF) AND COMPARE
C    WITH ORIGINAL TRACKING RESULTS

      IF (IREF.GT.0)
     &   CALL ARESI(     AREF,MAXTRAP,NORDNGP,
     &  RESXMR, RESXAVR, RESXR,
     &  RESYMR, RESYAVR, RESYR,
     &  RESPXMR,RESPXAVR,RESPXR,
     &  RESPYMR,RESPYAVR,RESPYR,RESREF)


C--- CALCULATE TRACKS FROM LINEAR TERMS OF GENERATING FUNCTION AND COMPARE
C    WITH ORIGINAL TRACKING RESULTS

      DO I=1,NORDNG
        DO J=1,NORDNG
          DO K=1,NORDNG
            DO L=1,NORDNG

              AKOEFFC(I,J,K,L)=0.D0
              IF(I-1+J-1+K-1+L-1.LE.2) AKOEFFC(I,J,K,L)=AKOEFF(I,J,K,L)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      CALL ARESI(AKOEFFC,MAXTRAP,NORDNGP,
     &  RSXMTK, RSXAVTK, RSXTK,
     &  RSYMTK, RSYAVTK, RSYTK,
     &  RSPXMTK,RSPXAVTK,RSPXTK,
     &  RSPYMTK,RSPYAVTK,RSPYTK,RESLIN)

C--- LINEAR TRANSFER MATRIX FOR CANONICAL VARIABLES

      TRACOP(1,1)=(-4.D0*AKOEFF(3,1,1,1)*AKOEFF(1,3,
     & 1,1)*AKOEFF(1,1,2,2)+2.D0*AKOEFF(3,1,1,1)*
     & AKOEFF(1,2,2,1)*AKOEFF(1,2,1,2)+AKOEFF(2,2
     & ,1,1)**2.D0*AKOEFF(1,1,2,2)-AKOEFF(2,2,1,1)*
     & AKOEFF(2,1,2,1)*AKOEFF(1,2,1,2)-AKOEFF(2,2
     & ,1,1)*AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1)+2.D0*
     & AKOEFF(2,1,2,1)*AKOEFF(2,1,1,2)*AKOEFF(1,3
     & ,1,1))/(AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)-
     & AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
      TRACOP(1,2)=(2.D0*AKOEFF(1,3,1,1)*AKOEFF(1,1,2
     & ,2)-AKOEFF(1,2,2,1)*AKOEFF(1,2,1,2))/(AKOEFF
     & (2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,
     & 2)*AKOEFF(1,2,2,1))
      TRACOP(1,3)=(AKOEFF(2,2,1,1)*AKOEFF(1,2,2,1
     & )*AKOEFF(1,1,2,2)-2.D0*AKOEFF(2,2,1,1)*AKOEFF(
     & 1,2,1,2)*AKOEFF(1,1,3,1)-2.D0*AKOEFF(2,1,
     & 2,1)*AKOEFF(1,3,1,1)*AKOEFF(1,1,2,2)+
     & AKOEFF(2,1,2,1)*AKOEFF(1,2,2,1)*AKOEFF(1,2
     & ,1,2)+4.D0*AKOEFF(2,1,1,2)*AKOEFF(1,3,1,1)*
     & AKOEFF(1,1,3,1)-AKOEFF(2,1,1,2)*AKOEFF(1,2
     & ,2,1)**2)/(AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     & -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
      TRACOP(1,4)=(AKOEFF(2,2,1,1)*AKOEFF(1,2,1,2
     & )-2.D0*AKOEFF(2,1,1,2)*AKOEFF(1,3,1,1))/(AKOEFF
     & (2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,
     & 2)*AKOEFF(1,2,2,1))
      TRACOP(2,1)=(-2.D0*AKOEFF(3,1,1,1)*AKOEFF(1,1,
     & 2,2)+AKOEFF(2,1,2,1)*AKOEFF(2,1,1,2))/(
     & AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1
     & ,1,2)*AKOEFF(1,2,2,1))
      TRACOP(2,2)=AKOEFF(1,1,2,2)/(AKOEFF(2,2,1,1
     & )*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)*AKOEFF(1,
     & 2,2,1))
      TRACOP(2,3)=(-AKOEFF(2,1,2,1)*AKOEFF(1,1,2,
     & 2)+2.D0*AKOEFF(2,1,1,2)*AKOEFF(1,1,3,1))/(
     & AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1
     & ,1,2)*AKOEFF(1,2,2,1))
      TRACOP(2,4)=(-AKOEFF(2,1,1,2))/(AKOEFF(2,2,1,
     & 1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)*AKOEFF(
     & 1,2,2,1))
      TRACOP(3,1)=(4.D0*AKOEFF(3,1,1,1)*AKOEFF(1,2,2
     & ,1)*AKOEFF(1,1,1,3)-2.D0*AKOEFF(3,1,1,1)*
     & AKOEFF(1,2,1,2)*AKOEFF(1,1,2,2)-2.D0*AKOEFF(2
     & ,2,1,1)*AKOEFF(2,1,2,1)*AKOEFF(1,1,1,3)+
     & AKOEFF(2,2,1,1)*AKOEFF(2,1,1,2)*AKOEFF(1,1
     & ,2,2)+AKOEFF(2,1,2,1)*AKOEFF(2,1,1,2)*
     & AKOEFF(1,2,1,2)-AKOEFF(2,1,1,2)**2.D0*AKOEFF(1,
     & 2,2,1))/(AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)
     & -AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
      TRACOP(3,2)=(-2.D0*AKOEFF(1,2,2,1)*AKOEFF(1,1,
     & 1,3)+AKOEFF(1,2,1,2)*AKOEFF(1,1,2,2))/(
     & AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1
     & ,1,2)*AKOEFF(1,2,2,1))
      TRACOP(3,3)=(-4.D0*AKOEFF(2,2,1,1)*AKOEFF(1,1,
     & 3,1)*AKOEFF(1,1,1,3)+AKOEFF(2,2,1,1)*
     & AKOEFF(1,1,2,2)**2+2.D0*AKOEFF(2,1,2,1)*AKOEFF(
     & 1,2,2,1)*AKOEFF(1,1,1,3)-AKOEFF(2,1,2,
     & 1)*AKOEFF(1,2,1,2)*AKOEFF(1,1,2,2)-AKOEFF(
     & 2,1,1,2)*AKOEFF(1,2,2,1)*AKOEFF(1,1,2,
     & 2)+2.D0*AKOEFF(2,1,1,2)*AKOEFF(1,2,1,2)*
     & AKOEFF(1,1,3,1))/(AKOEFF(2,2,1,1)*AKOEFF(1,
     & 1,2,2)-AKOEFF(2,1,1,2)*AKOEFF(1,2,2,1))
      TRACOP(3,4)=(2.D0*AKOEFF(2,2,1,1)*AKOEFF(1,1,1
     & ,3)-AKOEFF(2,1,1,2)*AKOEFF(1,2,1,2))/(AKOEFF
     & (2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,
     & 2)*AKOEFF(1,2,2,1))
      TRACOP(4,1)=(2.D0*AKOEFF(3,1,1,1)*AKOEFF(1,2,2
     & ,1)-AKOEFF(2,2,1,1)*AKOEFF(2,1,2,1))/(AKOEFF
     & (2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,
     & 2)*AKOEFF(1,2,2,1))
      TRACOP(4,2)=(-AKOEFF(1,2,2,1))/(AKOEFF(2,2,1,
     & 1)*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)*AKOEFF(
     & 1,2,2,1))
      TRACOP(4,3)=(-2.D0*AKOEFF(2,2,1,1)*AKOEFF(1,1,
     & 3,1)+AKOEFF(2,1,2,1)*AKOEFF(1,2,2,1))/(
     & AKOEFF(2,2,1,1)*AKOEFF(1,1,2,2)-AKOEFF(2,1
     & ,1,2)*AKOEFF(1,2,2,1))
      TRACOP(4,4)=AKOEFF(2,2,1,1)/(AKOEFF(2,2,1,1
     & )*AKOEFF(1,1,2,2)-AKOEFF(2,1,1,2)*AKOEFF(1,
     & 2,2,1))

      CALL UTIL_DETERMINANTE(4,TRACOP,DETCOP,IDETFAIL)

C--- LINEAR TRANSFER MATRIX WITHOUT COUPLING

      A1100=AKOEFF(2,2,1,1)
      A2000=AKOEFF(3,1,1,1)
      A0200=AKOEFF(1,3,1,1)

      A0011=AKOEFF(1,1,2,2)
      A0020=AKOEFF(1,1,3,1)
      A0002=AKOEFF(1,1,1,3)

      TRALIN(1,1)=(A1100*A1100-4.D0*A2000*A0200)/A1100
      TRALIN(1,2)=2.D0*A0200/A1100
      TRALIN(2,1)=-2.D0*A2000/A1100
      TRALIN(2,2)=1/A1100

      DETTRA1=TRALIN(1,1)*TRALIN(2,2)-TRALIN(2,1)*TRALIN(1,2)

      IF (I2DIM.EQ.0)  THEN

         TRALIN(3,3)=(A0011*A0011-4.D0*A0020*A0002)/A0011
         TRALIN(3,4)=2.D0*A0002/A0011
         TRALIN(4,3)=-2.D0*A0020/A0011
         TRALIN(4,4)=1/A0011

         DETTRA2=TRALIN(3,3)*TRALIN(4,4)-TRALIN(3,4)*TRALIN(4,3)

      ENDIF

      IF(TRALIN(2,1).NE.0.) FOCEX=-1.D0/TRALIN(2,1)
      IF(TRALIN(4,3).NE.0.) FOCEY=-1.D0/TRALIN(4,3)


C--- CALCULATE TRACKS FROM 2. ORDER TERMS USED FOR LINEAR TRANSFER MATRIX
C    AND COMPARE WITH ORIGINAL TRACKING RESULTS

      DO I=1,NORDNG
      DO J=1,NORDNG
      DO K=1,NORDNG
      DO L=1,NORDNG

      AKOEFFC(I,J,K,L)=0.D0
      IF (I-1+J-1+K-1+L-1 .EQ.2) AKOEFFC(I,J,K,L)=AKOEFF(I,J,K,L)
C     IF (I-1+J-1.EQ.1) AKOEFFC(I,J,K,L)=0.D0 !KILL COUPLING

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL ARESI(      AKOEFFC,MAXTRAP,NORDNGP,
     &                   RSXMT, RSXAVT, RSXT,
     &                   RSYMT, RSYAVT, RSYT,
     &                   RSPXMT,RSPXAVT,RSPXT,
     &                   RSPYMT,RSPYAVT,RSPYT,RESTRA)


C--- CALCULATE EQUIVALENT HALBACH WIGGLER

      IF(IHALBA.NE.0.AND.I2DIM.EQ.0) THEN

C     CALCULATION OF K-VALUES OF MATRIX (REFER BETA MANUAL PAGE 13)
C     ANSATZ: LONG QUADRUPOLE, K-VALUES(FY,FX):

                FY= TRALIN(4,3)/TRALIN(3,4)   ! <0 ,LARGE, FOCUSSING
                FX= TRALIN(2,1)/TRALIN(1,2)   ! >0 ,SMALL, DEFOCUSSING

                IF (FY.GE.0.0.OR.FX+FY.GT.0.0) THEN
                WRITE(6,*)
     &'*** MESSAGE SR TRANPOLY: NO EQUIVALENT HALBACH WIGGLER FOUND ***'
                WRITE(LUNGFO,*)
     &'*** MESSAGE SR TRANPOLY: NO EQUIVALENT HALBACH WIGGLER FOUND ***'
            IHALBA=0
            GOTO 500
         ENDIF

         R2HAL =FX/FY                      ! R2HAL<0 !!
         SINZ2 =-TRALIN(4,3)*TRALIN(3,4)
       SINZ  = DSIGN(DSQRT(SINZ2),TRALIN(3,4))
       COSZ  = (TRALIN(3,3)+TRALIN(4,4))/2.D0
         ZETAZ = ATAN2(SINZ,COSZ)
       IF (ZETAZ.LT.0.D0) ZETAZ=ZETAZ+2.D0*PI1
         ZLENGE= ZETAZ/DSQRT(-FY)
         RHOHAL= 1.D0/DSQRT(-2.D0*(FX+FY))
         ZLHAL = ZLENGE/NPERTRA
         ZKHAL=2.D0*PI1/ZLHAL
       ZKHALR=ZKHAL
       ZKHALI=0.D0

         B0HAL=PEL/CLIGHT1/RHOHAL
       IF (1.D0+R2HAL.GT.0.D0) THEN
          YKHAL=SQRT(ZKHAL*ZKHAL/ABS(1.D0+R2HAL))
        YKHALR=YKHAL
        YKHALI=0.D0
       ELSE IF (1.D0+R2HAL.GT.0.D0) THEN
          YKHAL=SQRT(ZKHAL*ZKHAL/ABS(1.D0+R2HAL))
        YKHALI=YKHAL
        YKHALR=0.D0
       ELSE
        YKHAL=0.D0
        YKHALR=0.D0
        YKHALI=0.D0
       ENDIF

         IF(YKHAL.NE.0.D0) YLHAL=2.D0*PI1/YKHAL
       IF (ZKHAL.GE.YKHAL) THEN
          XKHAL=SQRT((ZKHAL+YKHAL)*(ZKHAL-YKHAL))
        XKHALR=XKHAL
        XKHALI=0.D0
       ELSE
          XKHAL=SQRT((YKHAL+ZKHAL)*(YKHAL-ZKHAL))
        XKHALI=XKHAL
        XKHALR=0.D0
       ENDIF
         IF(XKHAL.NE.0.D0)XLHAL=2.D0*PI1/XKHAL

      ENDIF

C--- CALCULATE EQUIVALENT QUADRUPOLE FROM LINEAR TRANSFER MATRIX

500   CONTINUE

      IF (IQUAD.NE.0) THEN

      T11=0.5D0*(TRALIN(1,1)+TRALIN(2,2))
C     T11=DSQRT(TRALIN(1,1)*TRALIN(2,2))

      IF(T11.GT.1.D0)   THEN

         SIHPHIX=DSQRT((T11+1.D0)*(T11-1.D0))
         PHIX=DLOG(T11+SIHPHIX)

         IF(DABS(TRALIN(2,1)).GT.1.D-6) THEN

                  QUADX=DSQRT(TRALIN(1,2)/TRALIN(2,1)*PHIX*PHIX)
           QUADFX=-QUADX/(PHIX*PHIX)
           QUADR(1,1)=T11
           QUADR(2,2)=QUADR(1,1)
           QUADR(1,2)=QUADX/PHIX*SIHPHIX
           QUADR(2,1)=PHIX/QUADX*SIHPHIX

         ELSE

                  QUADX=0.
           QUADFX=0.
           QUADR(1,1)=T11
           QUADR(2,2)=QUADR(1,1)
           QUADR(1,2)=AKOEFF(1,3,1,1)*2.D0
           QUADR(2,1)=0.

         ENDIF

      ELSEIF (T11.GE.-1.D0) THEN

         SINPHIX=DSQRT(-(T11+1.D0)*(T11-1.D0))
         PHIX=DACOS(T11)

         IF(DABS(TRALIN(2,1)).GT.1.D-6) THEN

             QUADX=DSQRT(-TRALIN(1,2)/TRALIN(2,1)*PHIX*PHIX)
             QUADFX=QUADX/(PHIX*PHIX)
             QUADR(1,1)=T11
             QUADR(2,2)=QUADR(1,1)
             QUADR(1,2)=QUADX/PHIX*SINPHIX
             QUADR(2,1)=-PHIX/QUADX*SINPHIX

         ELSE

             QUADX=0.
             QUADFX=0.
             QUADR(1,1)=T11
             QUADR(2,2)=QUADR(1,1)
             QUADR(1,2)=AKOEFF(1,3,1,1)*2.D0
             QUADR(2,1)=0.

         ENDIF

      ENDIF !T11


      IF (I2DIM.EQ.0) THEN

         T33=0.5D0*(TRALIN(3,3)+TRALIN(4,4))
C        T33=DSQRT(TRALIN(3,3)*TRALIN(4,4))

      IF(T33.GT.1.D0) THEN

         SIHPHIY=DSQRT((T33+1.D0)*(T33-1.D0)) !SINUS-HYP.(PHI)
         PHIY=DLOG(T33+SIHPHIY)  !AREA-COSINUS-HYP.

         IF(DABS(TRALIN(4,3)).GT.1.D-6) THEN

            QUADY=DSQRT(TRALIN(3,4)/TRALIN(4,3)*PHIY*PHIY)
            QUADFY=-QUADY/(PHIY*PHIY)
            QUADR(3,3)=T33
            QUADR(4,4)=QUADR(3,3)
            QUADR(3,4)=QUADY/PHIY*SIHPHIY
            QUADR(4,3)=PHIY/QUADY*SIHPHIY

         ELSE

            QUADY=0.
            QUADFY=0.
            QUADR(3,3)=T33
                 QUADR(4,4)=QUADR(3,3)
                 QUADR(3,4)=AKOEFF(1,1,1,3)*2.D0
                 QUADR(4,3)=0.

         ENDIF

      ELSEIF (T33.GT.-1.D0) THEN

         SINPHIY=DSQRT((T33+1.D0)*(1.D0-T33))
         PHIY=DACOS(T33)

         IF(DABS(TRALIN(4,3)).GT.1.D-6) THEN

             QUADY=DSQRT(-TRALIN(3,4)/TRALIN(4,3)*PHIY*PHIY)
             QUADFY=QUADY/(PHIY*PHIY)
             QUADR(3,3)=T33
             QUADR(4,4)=QUADR(3,3)
             QUADR(3,4)=QUADY/PHIY*SINPHIY
             QUADR(4,3)=-PHIY/QUADY*SINPHIY

         ELSE

             QUADY=0.
             QUADFY=0.
             QUADR(3,3)=T33
             QUADR(4,4)=QUADR(3,3)
             QUADR(3,4)=AKOEFF(1,1,1,3)*2.D0
             QUADR(4,3)=0.

         ENDIF

C        QUADY=DSQRT(-TRALIN(3,4)/TRALIN(4,3)*PHIY*PHIY)
C        QUADFY= QUADY/(PHIY*PHIY)
C        QUADR(3,3)=T33
C        QUADR(4,4)=QUADR(3,3)
C        QUADR(3,4)= QUADY/PHIY*SINPHIY
C        QUADR(4,3)=-PHIY/QUADY*SINPHIY

      ENDIF !T33

      ENDIF !I2DIM

      DO I=1,NORDNG
      DO J=1,NORDNG
      DO K=1,NORDNG
      DO L=1,NORDNG

      AKOEFFC(I,J,K,L)=0.D0

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      IF(QUADR(2,2).NE.0) THEN

         AKOEFFC(2,2,1,1)=1.D0/QUADR(2,2)
         AKOEFFC(1,3,1,1)= QUADR(1,2)/(2.D0*QUADR(2,2))
         AKOEFFC(3,1,1,1)=-QUADR(2,1)/(2.D0*QUADR(2,2))

         IF (I2DIM.EQ.0) THEN

            AKOEFFC(1,1,2,2)=1.D0/QUADR(4,4)
            AKOEFFC(1,1,1,3)= QUADR(3,4)/(2.D0*QUADR(4,4))
            AKOEFFC(1,1,3,1)=-QUADR(4,3)/(2.D0*QUADR(4,4))

         ENDIF

         CALL ARESI(AKOEFFC,MAXTRAP,NORDNGP,
     &                   RSXMTQ, RSXAVTQ, RSXTQ,
     &                   RSYMTQ, RSYAVTQ, RSYTQ,
     &                   RSPXMTQ,RSPXAVTQ,RSPXTQ,
     &                   RSPYMTQ,RSPYAVTQ,RSPYTQ,RESQUAD)

      ENDIF !QUADR(2,2)

      ENDIF !IQUAD


C7.6.93  CALL CLOSEOR(AKOEFF(1,2,1,1),AKOEFF(2,1,1,1),
C7.6.93     &               AKOEFF(2,2,1,1),
C7.6.93     &               AKOEFF(3,1,1,1),AKOEFF(1,3,1,1),
C7.6.93     &               XCLOSE,PCLOSE)


C--- OUTPUT

C030192  WRITE (LUN1,*)'*** WARNUNG: DIE KANONISCHEN VARIABLEN WERDEN OHNE
C030192     & BERUECKSICTIGUNG'
C030192  WRITE (LUN1,*)'             DES VEKTORPOTENTIALES BERECHNET. NUR
C030192     & KORREKT, FALLS'
C030192        WRITE (LUN1,*)'             DAS POTENTIAL NULL IST. ***'

      WRITE (LUN1,*)
      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &  '     Tracks to fit generating function read from file (FILEO):'
      WRITE (LUN1,*)'       ',FILEO
      WRITE (LUN1,*)'     User comment of FILEO:'
      WRITE (LUN1,*)'       ',CODE
      WRITE (LUN1,*)'     Run number of FILEO :  ',kCODE
      WRITE (LUN1,*)'     Run number of TRANPOLY:',ITRANC
      WRITE (LUN1,*)
      WRITE (LUN1,*)'     Number of considered tracks (NTRAJ):',NTRAJ
      WRITE(LUN1,*)

      IF (DRAUSCHX.NE.0.0 .OR. DRAUSCHY.NE.0.0) THEN
      WRITE(LUN1,*)
     &'     Noise amplitude [m] overlaid on Z,ZP (DRAUSCHX):'
     &,SNGL(DRAUSCHX)
      WRITE(LUN1,*)
     &'     Noise amplitude [m] overlaid on Y,YP (DRAUSCHY):'
     &,SNGL(DRAUSCHY)
      WRITE(LUN1,*)
     &'     Noise amplitude on Z,ZP divided by mean length of (XF,PXI):'
     &,SNGL(DRAUSCHX/XFXCAVE)
      IF (I2DIM.EQ.0)
     &   WRITE(LUN1,*)
     &'     Noise amplitude on Y,YP divided by mean length of (YF, PYI):'
     &,SNGL(DRAUSCHY/XFYCAVE)
      ENDIF !DRAUSCH

C     WRITE (6,*)'CODE-WORT AUS WLS.FOR:',CODE
C     WRITE(6,*)'LAUFENDE NUMMER AUS WLS.FOR:',ICODE
C     WRITE(6,*)'LAUFENDE NUMMER AUS TRANPOLY.FOR:',ITRANC
C     WRITE(6,*)'RAUSCHAMPLITUDE IN X,XP:',SNGL(DRAUSCHX)
C     WRITE(6,*)'RAUSCHAMPLITUDE IN Y,YP:',SNGL(DRAUSCHY)

      WRITE (LUN1,*)
      WRITE (LUN1,*)'     Start values of reference orbit (Lab.-System):'
      WRITE (LUN1,*)'     X0,Y0,Z0:   ',SNGL(X0),SNGL(Y0),SNGL(Z0)
      WRITE (LUN1,*)'     ZP0,YP0:    ',SNGL(ZP0),SNGL(YP0)
      WRITE (LUN1,*)'     BX0,BY0,BZ0:',SNGL(BX0),SNGL(BY0),SNGL(BZ0)
      WRITE (LUN1,*)'     AX0,AY0,AZ0:',SNGL(AX0),SNGL(AY0),SNGL(AZ0)
      WRITE (LUN1,*)
      WRITE (LUN1,*)'     Final values of reference orbit (Lab.-System):'
      WRITE (LUN1,*)'     XF0,YF0,ZF0:   ',SNGL(XF0),SNGL(YF0),SNGL(ZF0)
      WRITE (LUN1,*)'     ZPF0,YPF0:    ',SNGL(ZPF0),SNGL(YPF0)
      WRITE (LUN1,*)'     BXF0,BYF0,BZF0:',SNGL(BXF0),SNGL(BYF0),SNGL(BZF0)
      WRITE (LUN1,*)'     AXF0,AYF0,AZF0:',SNGL(AXF0),SNGL(AYF0),SNGL(AZF0)
      WRITE (LUN1,*)

      WRITE (LUN1,*)'     Normal vector of entrance plane:'
      WRITE (LUN1,*)'     ',SNGL(OPNX),SNGL(OPNY),SNGL(OPNZ)
      WRITE (LUN1,*)
      WRITE (LUN1,*)'     Normal vector of exit plane:    '
      WRITE (LUN1,*)'     ',SNGL(OPNFX),SNGL(OPNFY),SNGL(OPNFZ)
      WRITE (LUN1,*)

      WRITE (LUN1,*)
      WRITE (LUN1,*)'     Energy [GeV]:',SNGL(GAMMA*EMASSE1/1.D9)
      WRITE (LUN1,*)'     B*RHO [Tm]:  ',SNGL(BRHOABS)
      WRITE (LUN1,*)

      WRITE (LUN1,*)
      WRITE (LUN1,*)'     Flags:'
      WRITE (LUN1,3300)IA1000,IA0100,IA0010,IA0001
3300  FORMAT('      IA1000,IA0100,IA0010,IA0001:',4I3)
      WRITE (LUN1,3330)I2DIM,ISYM,IA11A20,IWLSHOR
3330  FORMAT('      I2DIM,ISYM,IA11A20,IWLSHOR: ',4I3)
      WRITE (LUN1,*)

      IF (I2DIM.NE.0) THEN
         WRITE (LUN1,*)
         WRITE (LUN1,*)
         WRITE (LUN1,*)
     & '     *** 2-DIMENSIONAL FIT (FLAG I2DIM.NE.0) ***'
         WRITE (LUN1,*)
         WRITE (LUN1,*)
      ENDIF

      WRITE (LUN1,*)
      WRITE (LUN1,*) '     Number of tracks:',MTRAJ
      WRITE (LUN1,*)

      IF (ISELECT.GT.0) THEN

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Selection of initial parameters of tracks:'
      WRITE (LUN1,*)
     &'     ------------------------------------------'
      WRITE (LUN1,*)

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Coordinates and slopes ZI(N), ZPI(N), YI(N), YPI(N):'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
         WRITE (LUN1,1000)I,XI(1,I),XI(2,I),XI(3,I),XI(4,I)
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Canonical variables XI(N), PXF(N), YI(N), PYF(N)):'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
         WRITE (LUN1,1000)I,XIC(I),PXF(I),YIC(I),PYF(I)
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Magnetic field BXI(N), BYI(N), BZI(N)):'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
      WRITE (LUN1,1000)I,BXI(I),BYI(I),BZI(I)
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     XI_lab(N) and vector pot. AXI(N)/BRHO, AYI(N)/BRHO, AZI(N)/BRHO:'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
      WRITE (LUN1,1000)I,SI(I),AXI(I)/BRHOABS,AYI(I)/BRHOABS,AZI(I)/BRHOABS
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Selection of final parameters of tracks:'
      WRITE (LUN1,*)
     &'     ----------------------------------------'
      WRITE (LUN1,*)

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Coordiantes and slopes ZF(N), ZPF(N), YF(N), YPF(N):'
      WRITE (LUN1,*)
      WRITE (LUN1,*)
      DO I=1,MTRAJ,ISELECT
         WRITE (LUN1,1000)I,XF(1,I),XF(2,I),XF(3,I),XF(4,I)
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Canonical variables XF(N), PXI(N), YF(N), PYI(N)):'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
         WRITE (LUN1,1000)I,XFC(I),PXI(I),YFC(I),PYI(I)
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Magnetic field BXI(N), BYI(N), BZI(N)):'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
      WRITE (LUN1,1000)I,BXF(I),BYF(I),BZF(I)
      ENDDO

      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     XF_lab(N) and vector pot. AXF(N)/BRHO, AYF(N)/BRHO, AZF(N)/BRHO:'
      WRITE (LUN1,*)

      DO I=1,MTRAJ,ISELECT
      WRITE (LUN1,1000)I,SF(I),AXF(I)/BRHOABS,AYF(I)/BRHOABS,AZF(I)/BRHOABS
      ENDDO

      ENDIF !ISELECT

      WRITE (LUN1,*)
      WRITE (LUN1,*)
      WRITE (LUN1,*)
     &'     Mean length of canonical vectors (XI, PXF, YI, PYF):',
     &                 SNGL(XICAVE)
      WRITE (LUN1,*)
     &'     Mean length of canonical vectors (XF, PXI, YF, PYI):',
     &                 SNGL(XFCAVE)
      WRITE (LUN1,*)

      WRITE(LUN1,*)
     &'     Scaling factor for system of linear equations:',
     &                 SNGL(DSCALE)

      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Number of fitted coefficients:',NKOEF
      WRITE(LUN1,*)

      WRITE(LUN1,*)
     &'     Check of solution of linear equation system, i.e maximum, mean,'
      WRITE(LUN1,*)'     and rms of deviation:'
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESBM,RESBAV,RESB

      WRITE(LUN1,*)
      IF(IERZ.GT.0) THEN
         WRITE(LUN1,*)
     &'     Run numbers and user comment on file WAVE_ERZFUN.IN:'
         WRITE(LUN1,*)'     ',ICODEERZ,ITRANCE
         WRITE(LUN1,*)'          ',CODEERZ
         WRITE(LUN1,*)
     &'     Number of coefficients AERZ(I,J,K,L):',IERZ
      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean and rms of AKOEFF(I,J,K,L)-AERZ(I,J,K,L):'
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESAM,RESAAV,RESA
      ENDIF !IERZ

      IF(IREF.GT.0) THEN
      WRITE(LUN1,*)
      WRITE(LUN1,*)
         WRITE(LUN1,*)
     &'     Run numbers and user comment on file WAVE_ERZFUN.REF:'
         WRITE(LUN1,*)'     ', ICODEREF,ITRANCR
         WRITE(LUN1,*)'     ',CODEREF
         WRITE(LUN1,*)
     &'     Number of coefficients AREF(I,J,K,L):',IREF
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean and rms of AKOEFF(I,J,K,L)-AREF(I,J,K,L):'
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESAMR,RESAAVR,RESAR
      ENDIF !IREF

      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean, and rms of deviation for XF, PXI, YF, PYI'
      WRITE(LUN1,*)'     recalculated using AKOEFF:'
      WRITE(LUN1,*)
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESXMK, RESXAVK, RESXK
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESPXMK,RESPXAVK,RESPXK
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESYMK, RESYAVK, RESYK
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESPYMK,RESPYAVK,RESPYK

c     WRITE(6,*)'ANZAHL DER GEFITTETEN KOEFF., SQRT(CHI**2/N) DES FITS
c     &:',NKOEF,SNGL(RESAKO)
c     WRITE(6,*)'SQRT(CHI**2/N) DES FITS/MITTLERE LAENGE DER END-VEKTOREN
c     &:',SNGL(RESAKO/XFCAVE)

      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     SQRT(CHI**2/N) corresponding to AKOEFF:         '
     &     ,SNGL(RESAKO)
      WRITE(LUN1,*)
     &'     Dito divided by mean length of (XI,PXF,YI,PYF):'
     &       ,SNGL(RESAKO/XFCAVE)

      IF (IERZ.GT.0) THEN
      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean, and rms of deviation for XF, PXI, YF, PYI'
      WRITE(LUN1,*)'     recalculated using AERZ:'
      WRITE(LUN1,*)
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESXME, RESXAVE, RESXE
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESPXME,RESPXAVE,RESPXE
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESYME, RESYAVE, RESYE
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESPYME,RESPYAVE,RESPYE
      WRITE(LUN1,*)

      WRITE(LUN1,*)
     &'     SQRT(CHI**2/N) corresponding to AERZ:           '
     &      ,SNGL(RESERZ)
      WRITE(LUN1,*)
     &'     Dito divided by mean length of (XI,PXF,YI,PYF):'
     &       ,SNGL(RESERZ/XFCAVE)
      ENDIF !IERZ

      IF (IREF.GT.0) THEN
      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean, and rms of deviation for XF, PXI, YF, PYI'
      WRITE(LUN1,*)'     recalculated using AREF:'
      WRITE(LUN1,*)
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESXMR, RESXAVR, RESXR
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESPXMR,RESPXAVR,RESPXR
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESYMR, RESYAVR, RESYR
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RESPYMR,RESPYAVR,RESPYR
      WRITE(LUN1,*)

      WRITE(LUN1,*)
     &'     SQRT(CHI**2/N) corresponding to AREF:           '
     &     ,SNGL(RESREF)
      WRITE(LUN1,*)
     &'     Dito divided by mean length of (XI,PXF,YI,PYF):'
     &       ,SNGL(RESREF/XFCAVE)
      ENDIF !IREF

C--- WRITE COEFFICIENTS AKOEFF TO OUTPUT FILE WAVE_TRANPOLY.AKO

      OPEN(UNIT=LUNAKO,FILE='wave_tranpoly.ako',STATUS='NEW',FORM='FORMATTED')

C10.6.93 WRITE (LUNAKO,'(I5,4H    ,1A66)')ICODE,CODE
      WRITE (LUNAKO,'(I10,I5,4H    ,1A66)')ICODE,ITRANC,CODE
      WRITE(LUNAKO,*)ICHARGE,NORDNG
      WRITE(LUNAKO,*)ZAPERT,YAPERT,DLAPER,sint,gamma,dgam
      WRITE(LUNAKO,*)X0,Y0,Z0,ZP0,YP0,BX0,BY0,BZ0,AX0,AY0,AZ0
      WRITE(LUNAKO,*)XF0,YF0,ZF0,ZPF0,YPF0,BXF0,BYF0,BZF0,AXF0,AYF0,AZF0
      WRITE(LUNAKO,*)OPNX,OPNY,OPNZ
      WRITE(LUNAKO,*)OPNFX,OPNFY,OPNFZ
      WRITE(LUNAKO,*)A0SCGLOB


      DO N=1,NKOEF
         I=IE(N)
         J=JE(N)
         K=KE(N)
         L=LE(N)
         WRITE(LUNAKO,*) I-1,J-1,K-1,L-1,AKOEFF(I,J,K,L)
      ENDDO

      CLOSE(LUNAKO)

      IF (IERZ.NE.0.OR.IREF.NE.0) THEN

      WRITE (LUN1,*)

      WRITE (LUN1,*)
     &'     Comparison of fitted coefficients AKOEFF with AERZ and AREF;'
      WRITE (LUN1,*)
     &'     AKOEFF, AERZ, 2*(AERZ-AKOEF)/(|AERZ|+|AKOEF|), and also for AREF'
      WRITE (LUN1,*)
     &'     (only printed if deviations larger than 1E-10)'
      WRITE (LUN1,*)

      DO N=2,NORDNG
          WRITE(LUN1,*)
          DO I=1,NORDNG
          DO J=1,NORDNG
          DO K=1,NORDNG
          DO L=1,NORDNG

      IF(I-1+J-1+K-1+L-1.EQ.N-1) THEN
C7.6.93  IF(I-1+J-1+K-1+L-1.EQ.N-1.AND.AKOEFF(I,J,K,L).NE.0.0) THEN

         AK=AKOEFF(I,J,K,L)
         AE=AERZ(I,J,K,L)
         AR=AREF(I,J,K,L)
         AKE=0.D0
         AKR=0.D0

         IF (DABS(AE)+DABS(AK).NE.0.0)
     &             AKE=2.*(AE-AK)/(DABS(AE)+DABS(AK))
         IF (DABS(AR)+DABS(AK).NE.0.0)
     &             AKR=2.*(AR-AK)/(DABS(AR)+DABS(AK))

        IF(DABS(AKE).GT.1.D-10.OR.DABS(AKR).GT.1.D-10) THEN

         WRITE(LUN1,8000) I-1,J-1,K-1,L-1,AK,AE,AKE,AR,AKR
8000    FORMAT('      ','a',4I1,1P,(2D15.7),1P,(D10.2),1P,(D15.7),1P,(D10.2))

C7.6.93     IF (AE*AK.NE.0.D0) AKE=(AK/AE+AE/AK)/2.D0-1.D0
C7.6.93     IF (AR*AK.NE.0.D0) AKR=(AK/AR+AR/AK)/2.D0-1.D0

C7.6.93     IF( AE.NE.0.D0 .AND. AR.NE.0.D0)
C7.6.93     &     WRITE(LUN1,3000) I-1,J-1,K-1,L-1,AK,AE,AKE,AR,AKR
C7.6.93     IF( AE.EQ.0.D0 .AND. AR.NE.0.D0)
C7.6.93     &     WRITE(LUN1,3001) I-1,J-1,K-1,L-1,AK,       AR,AKR
C7.6.93     IF( AE.NE.0.D0 .AND. AR.EQ.0.D0)
C7.6.93     &     WRITE(LUN1,3002) I-1,J-1,K-1,L-1,AK,AE,AKE
C7.6.93     IF( AE.EQ.0.D0 .AND. AR.EQ.0.D0)
C7.6.93     &     WRITE(LUN1,3003) I-1,J-1,K-1,L-1,AK

      ENDIF
      ENDIF

          ENDDO
          ENDDO
          ENDDO
          ENDDO
      ENDDO

      ENDIF

      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean, and rms of deviation for XF, PXI, YF, PYI'
      WRITE(LUN1,*)
     &'     recalculated using only linear terms of AKOEFF, i.e. first'
      WRITE(LUN1,*)'     and second order:'
      WRITE(LUN1,*)
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSXMTK, RSXAVTK, RSXTK
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSPXMTK,RSPXAVTK,RSPXTK
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSYMTK, RSYAVTK, RSYTK
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSPYMTK,RSPYAVTK,RSPYTK
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Corresponding SQRT(CHI**2/N):                   ',SNGL(RESLIN)
      WRITE(LUN1,*)
     &'     Dito divided by mean length of (XI,PXF,YI,PYF):',SNGL(RESLIN/XFCAVE)

      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Linear transfer matrix (for canonical variables):'
      WRITE(LUN1,*)
      DO I=1,4
         WRITE(LUN1,5000) (TRACOP(I,J),J=1,4)
      ENDDO
      WRITE(LUN1,*)
      WRITE(LUN1,*)'      Failure flag and determinant:',IDETFAIL,DETCOP
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Linear transfer matrix (for canonical variables, without coupling):'
      WRITE(LUN1,*)
      DO I=1,4
         WRITE(LUN1,5000) (TRALIN(I,J),J=1,4)
      ENDDO
      WRITE(LUN1,*)
C     WRITE(LUN1,*) '     Deviations of sub-determinants from unity:',
C     &  SNGL(1.D0-DETTRA1),SNGL(1.D0-DETTRA2)
C     WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Horiz. and vert. focal lengths (-1/m21 and -1/m43) [1/m]:'
      WRITE(LUN1,*)'     ', SNGL(FOCEX),SNGL(FOCEY)
      WRITE(LUN1,*)

      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean, and rms of deviation for XF, PXI, YF, PYI recalculated'
      WRITE(LUN1,*)
     &'     using all 2. order terms of AKOEFF (including coupling):'
      WRITE(LUN1,*)
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSXMT, RSXAVT, RSXT
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSPXMT,RSPXAVT,RSPXT
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSYMT, RSYAVT, RSYT
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSPYMT,RSPYAVT,RSPYT
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Corresponding SQRT(CHI**2/N):                   '
     &      ,SNGL(RESTRA)
      WRITE(LUN1,*)
     &'     Dito divided by mean length of (XI,PXF,YI,PYF):'
     &       ,SNGL(RESTRA/XFCAVE)

      WRITE(LUN1,*)
      WRITE(LUN1,*)
      WRITE(LUN1,*)

      IF (IQUAD.NE.0) THEN
      WRITE(LUN1,*) '     Quadrupole matrix (for canonical variables):'
      WRITE(LUN1,*)
      DO I=1,4
         WRITE(LUN1,5000) (QUADR(I,J),J=1,4)
      ENDDO
      WRITE(LUN1,*)
C     WRITE(LUN1,*) '     Deviations of sub-determinants from unity:',
C     &  SNGL(1.D0-(QUADR(1,1)*QUADR(2,2)-QUADR(2,1)*QUADR(1,2))),
C     &  SNGL(1.D0-(QUADR(3,3)*QUADR(4,4)-QUADR(3,4)*QUADR(4,3)))
C     WRITE(LUN1,*)
      WRITE(LUN1,*)'     Hori. and vert. focal length (1/f=k*l) [m]:',
     &  SNGL(QUADFX),SNGL(QUADFY)
C     WRITE(6,*)'FOKAL-LAENGEN IN X UND Y (METER):',QUADFX,QUADFY

C     WRITE(6,*) 'XCLOSE:  ',xclose
C     WRITE(6,*) 'PxfCLOSE:',Pclose
      WRITE(LUN1,*)

      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Maximum, mean, and rms of deviation for XF, PXI, YF, PYI'
      WRITE(LUN1,*)
     &'     recalculated using those 2. order terms of AKOEFF that'
      WRITE(LUN1,*)
     &'     have been used to determine the quadrupole matrix:'
      WRITE(LUN1,*)
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSXMTQ, RSXAVTQ, RSXTQ
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSPXMTQ,RSPXAVTQ,RSPXTQ
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSYMTQ, RSYAVTQ, RSYTQ
      WRITE(LUN1,'(5H     ,1P,(3D15.3))') RSPYMTQ,RSPYAVTQ,RSPYTQ
      WRITE(LUN1,*)
      WRITE(LUN1,*)
     &'     Corresponding SQRT(CHI**2/N):                   ',SNGL(RESQUAD)
      WRITE(LUN1,*)
     &'     Dito divided by mean length of (XI,PXF,YI,PYF):'
     &       ,SNGL(RESQUAD/XFCAVE)
      ENDIF !IQUAD

        IF ((FY.GE.0.0.OR.FX+FY.GT.0.0) .AND. IHALBA.NE.0) THEN
         WRITE (LUN1,*)
         WRITE (LUN1,*)'*** WARNING SR TRANPOLY ***'
         WRITE (LUN1,*)
     & 'EQUIVALENT HALBACH WIGGLER CAN NOT BE DETERMINED'
         WRITE (LUN1,*)
      ENDIF

      IF(IHALBA.NE.0) THEN
         WRITE(LUN1,*)
         WRITE(LUN1,*)
     &'     Parameters of equivalent Halbach wiggler, KX**2+KY**2=KZ**2 (complexe),'
         WRITE(LUN1,*)
     &'     coordinates according to Halbach, z-axis is long. axis:'
         WRITE(LUN1,*)
         WRITE(LUN1,*)'     number of periods: ',NPERTRA
         WRITE(LUN1,*)'     RHO [m], ABS(B0) [T]:'
         WRITE(LUN1,*)'     ',SNGL(RHOHAL),SNGL(B0HAL)
         WRITE(LUN1,*)'     LX, LY, LZ, LZ/2 [m]:     '
         WRITE(LUN1,*)'     ',SNGL(XLHAL),SNGL(YLHAL),SNGL(ZLHAL),SNGL(ZLHAL/2.D0)
         WRITE(LUN1,*)'     KX, KY, KZ (complexe) [m]:'
         WRITE(LUN1,*)'     ',XKHALR,XKHALI
         WRITE(LUN1,*)'     ',YKHALR,YKHALI
         WRITE(LUN1,*)'     ',ZKHALR,ZKHALI
         WRITE(LUN1,*)
         WRITE(LUN1,*)

         OPEN(UNIT=LUNH,FILE='wave_wls.dat',STATUS='NEW',FORM='FORMATTED')
         WRITE(LUNH,*)B0HAL,XLHAL,YLHAL,ZLHAL
         CLOSE(LUNH)

         OPEN(UNIT=LUNSTR,FILE='wave_wls.str',STATUS='NEW',FORM='FORMATTED')
         WRITE(LUNSTR,1234)ZLENGE,RHOHAL,-(YLHAL/XLHAL)**2,ZLENGE
1234  FORMAT(' WLS  ID',E14.6,E14.6,E14.6,E14.6)
         WRITE(LUNSTR,1235)(STRAILEN-ZLENGE)/2.D0
1235  FORMAT(' DWLS SD',E14.6,'  0.000000E+00  0.000000E+00  0.000000E+00')
         CLOSE(LUNSTR)

      ENDIF

      IF(IWARNE.EQ.1) THEN
         WRITE(LUN1,*)
         WRITE(LUN1,*)
         WRITE(LUN1,*)'*** WARNING SR TRANPOLY ***'
         WRITE(LUN1,*)
     &'COEFFICIENTS AERZ EXCEEDING ORDER OF GENERATING FUNCTION IGNORED'
         WRITE(LUN1,*)
         WRITE(LUN1,*)
         WRITE(LUN1,*)
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR TRANPOLY ***'
         WRITE(6,*)
     &'COEFFICIENTS AERZ EXCEEDING ORDER OF GENERATING FUNCTION IGNORED'
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)
      ENDIF

      IF(IWARNR.EQ.1) THEN
         WRITE(LUN1,*)
         WRITE(LUN1,*)
         WRITE(LUN1,*)'*** WARNING SR TRANPOLY ***'
         WRITE(LUN1,*)
     &'COEFFICIENTS AREF EXCEEDING ORDER OF GENERATING FUNCTION IGNORED'
         WRITE(LUN1,*)
         WRITE(LUN1,*)
         WRITE(LUN1,*)
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)'*** WARNING SR TRANPOLY ***'
         WRITE(6,*)
     &'COEFFICIENTS AREF EXCEEDING ORDER OF GENERATING FUNCTION IGNORED'
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)
      ENDIF

C7.6.93  WRITE(LUN1,*) '     X closed orbit, PXF close orbit:'
C7.6.93     &               ,SNGL(xclose),SNGL(Pclose)

1000    FORMAT('      ',I4,' ',1P,(4D15.6))
c2000    FORMAT(' ',I3,'   ',1P,(2D17.8))
C3000      FORMAT(' ',I3,5H     ,'a',4I1,1P,(1D17.5))
c3000    FORMAT('      ','a',4I1,1P,(2D20.10),1P,(D10.2),1P,(D20.10),1P,(D10.2))
c3001    FORMAT(' ','a',4I1,1P,(1D20.10),      1P,(D50.10),1P,(D10.2))
c3002    FORMAT(' ','a',4I1,1P,(2D20.10),1P,(D10.2))
c3003    FORMAT(' ','a',4I1,1P1D20.10)
c4000      FORMAT(4I2,'  ',1P,(D27.20))
5000    FORMAT('     ',4F18.13)

C--- OUTPUT COEFFICIENTS SORTED BY ORDER OF GENERATING FUNCTION

      IF(IKOEFF.NE.0) THEN
          WRITE(LUN1,*)
          WRITE(LUN1,*)'     Coefficients of generating function (AKOEFF):'
          WRITE(LUN1,*)
      DO N=2,NORDNG
          WRITE(LUN1,*)
          WRITE(LUN1,*)'     Order of generating function:',N-1
          WRITE(LUN1,*)
          DO I=1,NORDNG
          DO J=1,NORDNG
          DO K=1,NORDNG
          DO L=1,NORDNG
         IF(I-1+J-1+K-1+L-1.EQ.N-1.AND.AKOEFF(I,J,K,L).NE.0.0)
     &      WRITE (LUN1,1236)I-1,J-1,K-1,L-1,AKOEFF(I,J,K,L)
1236     FORMAT('      a',4I1,1PE30.15)
          ENDDO
          ENDDO
          ENDDO
          ENDDO
      ENDDO
      ENDIF

      IF(LUN1.NE.LUNGFO) CLOSE(LUN1)

C---------------------------------------------------------------

      vx0=v0*v0x
      vy0=v0*v0y
      vz0=v0*v0z
      call idtrmshgf(
     &  x0,y0,z0,vx0,vy0,vz0,
     &  x2,y2,z2,vx2,vy2,vz2,
     &  x0,y0,z0,vx0,vy0,vz0,
     &  xf0,yf0,zf0,vxf0,vyf0,vzf0,"wave_tranpoly.ako")

      write(lungfo,*)
      write(lungfo,*)"Closed orbit test at the end of tranpoly:"
      write(lungfo,*)
      write(lungfo,*)"idtrmshgf tracks from:"
      write(lungfo,*)
      write(lungfo,*)x0,y0,z0
      write(lungfo,*)vx0,vy0,vz0
      write(lungfo,*)
      write(lungfo,*)"to:"
      write(lungfo,*)
      write(lungfo,*)x2,y2,z2
      write(lungfo,*)vx2,vy2,vz2
      write(lungfo,*)

      RETURN
      END
