*CMZ :  4.00/14 22/12/2021  16.40.21  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/01 20/03/2013  10.18.08  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  11.51.36  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.48/04 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.42/04 29/10/2002  10.38.26  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.41/08 12/08/2002  16.38.26  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  18.05.38  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  13.48.28  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  14.54.02  by  Michael Scheer
*CMZ : 00.02/03 31/01/97  15.06.26  by  Michael Scheer
*CMZ :          22/01/97  11.27.20  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97

      SUBROUTINE BPHARMFIT
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

C FIT MAGNETIC FIELD ACCORDING TO SCALAR POTENTIAL:
C     VC:=-B0C/KYC(N,NXY)*COS(NXY*KXC*X)*SINH(KYC(N,NXY)*Y)*COS(N*KZ*Z)$
C     VS:=-B0S/KXS(N,NXY)*COS(NXY*KYS*Y)*SINH(KXS(N,NXY)*X)*SIN(N*KZ*Z)$
C     SEE REDUCE VPOT-HARM, COMMENT FROM 29OCT02

*KEEP,bpharmf90u.
      include 'bpharmf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpharm.
      include 'bpharm.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      INTEGER IPOI,IFAIL,IX,IPAR
      INTEGER MPOI,IREAD
      INTEGER NPHPOIX,NPHPOIY,NPHPOIZ,NPHPOI

      DOUBLE PRECISION
     &  RESBX,RESBY,RESBZ,RESB
     &  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &  ,BERRMX(8),BXERRMX(8),BYERRMX(8),BZERRMX(8)
     &  ,XX,YY,ZZ,BBX,BBY,BBZ,AX,AY,AZ

      DOUBLE PRECISION XOFF

      INTEGER NTUP_P,ICYCLE
      PARAMETER (NTUP_P=9)
      REAL*8 TUP_D(NTUP_P)
      CHARACTER(3) CHTAGS_D(NTUP_P)
      data chtags_d/'x','y','z','bx','by','bz','bxf','byf','bzf'/

      IF (NHARM.GT.NHARMP) THEN
        WRITE(6,*) '*** ERROR IN BPHARMFIT: DIMENSION NHARMP EXCEEDED'
        STOP
      ENDIF
      IF (NTRANS.GT.NTRANSP) THEN
        WRITE(6,*) '*** ERROR IN BPHARMFIT: DIMENSION NTRANSP EXCEEDED'
        STOP
      ENDIF

C--- READ FIELD MAP {

      XOFF=PERLENPH*PHASEPH

      NPHPOIX=0
      NPHPOIY=0
      NPHPOIZ=0
      MPOI=0
      IREAD=0

      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='OLD')

100   CONTINUE
      IF (IWBPHARM.EQ.1) THEN
        READ(LUNBMAP,*,END=90)XX
        READ(LUNBMAP,*)YY
        READ(LUNBMAP,*)ZZ
        READ(LUNBMAP,*)BBX
        READ(LUNBMAP,*)BBY
        READ(LUNBMAP,*)BBZ
        READ(LUNBMAP,*)
      ELSE IF (IWBPHARM.EQ.2) THEN
        READ(LUNBMAP,*,END=90)XX,YY,ZZ,BBX,BBY,BBZ
      ELSE
        WRITE(6,*) '*** ERROR IN BPHARMFIT: FLAG IWBPHARM WRONG, CHECK INPUT'
        STOP
      ENDIF
      IREAD=IREAD+1
      IF (XPHMIN.NE.9999.AND.XX.LT.XPHMIN) GOTO 100
      IF (XPHMAX.NE.9999.AND.XX.GT.XPHMAX) GOTO 100
      IF (YPHMIN.NE.9999.AND.YY.LT.YPHMIN) GOTO 100
      IF (YPHMAX.NE.9999.AND.YY.GT.YPHMAX) GOTO 100
      IF (ZPHMIN.NE.9999.AND.ZZ.LT.ZPHMIN) GOTO 100
      IF (ZPHMAX.NE.9999.AND.ZZ.GT.ZPHMAX) GOTO 100
      IF (BBX.EQ.-9999..AND.BBY.EQ.-9999..AND.BBZ.EQ.-9999.) GOTO 100
      MPOI=MPOI+1

      IF (BBX.NE.-9999.) NPHPOIX=NPHPOIX+1
      IF (BBY.NE.-9999.) NPHPOIY=NPHPOIY+1
      IF (BBZ.NE.-9999.) NPHPOIZ=NPHPOIZ+1
      GOTO 100
90    CLOSE(LUNBMAP)

      ALLOCATE(X(MPOI))
      ALLOCATE(Y(MPOI))
      ALLOCATE(Z(MPOI))
      ALLOCATE(BX(MPOI))
      ALLOCATE(BY(MPOI))
      ALLOCATE(BZ(MPOI))
      ALLOCATE(BXF(MPOI))
      ALLOCATE(BYF(MPOI))
      ALLOCATE(BZF(MPOI))

      NPHPOIX=0
      NPHPOIY=0
      NPHPOIZ=0
      MPOI=0
      IREAD=0
      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='OLD')

101   CONTINUE
      IF (IWBPHARM.EQ.1) THEN
        READ(LUNBMAP,*,END=91)XX
        READ(LUNBMAP,*)YY
        READ(LUNBMAP,*)ZZ
        READ(LUNBMAP,*)BBX
        READ(LUNBMAP,*)BBY
        READ(LUNBMAP,*)BBZ
        READ(LUNBMAP,*)
      ELSE IF (IWBPHARM.EQ.2) THEN
        READ(LUNBMAP,*,END=91)XX,YY,ZZ,BBX,BBY,BBZ
      ELSE
        WRITE(6,*) '*** ERROR IN BPHARMFIT: FLAG IWBPHARM WRONG, CHECK INPUT'
        STOP
      ENDIF
      IREAD=IREAD+1
      IF (XPHMIN.NE.9999.AND.XX.LT.XPHMIN) GOTO 101
      IF (XPHMAX.NE.9999.AND.XX.GT.XPHMAX) GOTO 101
      IF (YPHMIN.NE.9999.AND.YY.LT.YPHMIN) GOTO 101
      IF (YPHMAX.NE.9999.AND.YY.GT.YPHMAX) GOTO 101
      IF (ZPHMIN.NE.9999.AND.ZZ.LT.ZPHMIN) GOTO 101
      IF (ZPHMAX.NE.9999.AND.ZZ.GT.ZPHMAX) GOTO 101
      IF (BBX.EQ.-9999..AND.BBY.EQ.-9999..AND.BBZ.EQ.-9999.) GOTO 101
      MPOI=MPOI+1
      IF (BBX.NE.-9999.) NPHPOIX=NPHPOIX+1
      IF (BBY.NE.-9999.) NPHPOIY=NPHPOIY+1
      IF (BBZ.NE.-9999.) NPHPOIZ=NPHPOIZ+1
      X(MPOI)=XX+XOFF
      Y(MPOI)=YY
      Z(MPOI)=ZZ
      BX(MPOI)=BBX
      BY(MPOI)=BBY
      BZ(MPOI)=BBZ
      GOTO 101
91    CLOSE(LUNBMAP)

      NPHPOI=MPOI

C--- READ FIELD MAP }

      IF (XPHMIN.EQ.9999.) THEN
        XPHMIN=1.D30
        DO IPOI=1,MPOI
          XX=X(IPOI)-XOFF
          IF (XX.LT.XPHMIN) XPHMIN=XX
        ENDDO
        XPHMIN=XPHMIN
      ENDIF

      IF (XPHMAX.EQ.9999.) THEN
        XPHMAX=-1.D30
        DO IPOI=1,MPOI
          XX=X(IPOI)-XOFF
          IF (XX.GT.XPHMAX) XPHMAX=XX
        ENDDO
        XPHMAX=XPHMAX
      ENDIF

      IF (YPHMIN.EQ.9999.) THEN
        YPHMIN=1.D30
        DO IPOI=1,MPOI
          IF (Y(IPOI).LT.YPHMIN) YPHMIN=Y(IPOI)
        ENDDO
        YPHMIN=YPHMIN
      ENDIF

      IF (YPHMAX.EQ.9999.) THEN
        YPHMAX=-1.D30
        DO IPOI=1,MPOI
          IF (Y(IPOI).GT.YPHMAX) YPHMAX=Y(IPOI)
        ENDDO
        YPHMAX=YPHMAX
      ENDIF

      IF (ZPHMIN.EQ.9999.) THEN
        ZPHMIN=1.D30
        DO IPOI=1,MPOI
          IF (Z(IPOI).LT.ZPHMIN) ZPHMIN=Z(IPOI)
        ENDDO
        ZPHMIN=ZPHMIN
      ENDIF

      IF (ZPHMAX.EQ.9999.) THEN
        ZPHMAX=-1.D30
        DO IPOI=1,MPOI
          IF (Z(IPOI).GT.ZPHMAX) ZPHMAX=Z(IPOI)
        ENDDO
        ZPHMAX=ZPHMAX
      ENDIF

      IF (XLENCPH.EQ.9999.) XLENCPH=DABS(ZPHMAX-ZPHMIN)*2.d0
      IF (YLENSPH.EQ.9999.) YLENSPH=DABS(YPHMAX-YPHMIN)*2.d0

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE BPHARMFIT:'
      WRITE(LUNGFO,*)'     ======================'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     period length and long. shift:        '
     &  , SNGL(PERLENPH),SNGL(PHASEPH)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     transversal lengths (hor. and vert.): '
     &  , SNGL(XLENCPH),SNGL(YLENSPH)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     field map read from file:'
      WRITE(LUNGFO,*)'     ',FILEBMAP
      WRITE(LUNGFO,*)

      WRITE(LUNGFO,*)'     number of data points read:      ',IREAD
      WRITE(LUNGFO,*)'     number of data points accepted:  ',MPOI
      WRITE(LUNGFO,*)'     number of Bx, By, Bz .ne. -9999.:'
     &  ,NPHPOIX,NPHPOIY,NPHPOIZ
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     fitted coefficients are written to file:'
      WRITE(LUNGFO,*)'     ',FILEPHFIT
      WRITE(LUNGFO,*)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Considered range in x, y, and z:'
      WRITE(LUNGFO,*)'     ',XPHMIN,XPHMAX
      WRITE(LUNGFO,*)'     ',YPHMIN,YPHMAX
      WRITE(LUNGFO,*)'     ',ZPHMIN,ZPHMAX
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     NTRANS0,NTRANS,NTRANSD:',NTRANS0,NTRANS,NTRANSD
      WRITE(LUNGFO,*)'     NHARM0,NHARM,NHARMD:',NHARM0,NHARM,NHARMD
      WRITE(LUNGFO,*)

      IF (XLENCPH.NE.0.D0) THEN
        XKCPH=2.D0*PI1/XLENCPH
      ELSE
        XKCPH=0.D0
      ENDIF
      IF (YLENSPH.NE.0.D0) THEN
        YKSPH=2.D0*PI1/YLENSPH
      ELSE
        YKSPH=0.D0
      ENDIF

      CALL VPOTPH(MPOI,IFAIL)

      IF (IFAIL.NE.0) THEN
        WRITE(6,*) '*** IFAIL NOT ZERO ***'
        WRITE(6,*) '*** (maybe no transversal field gradient ***'
        WRITE(6,*) '*** (or XLENCPH, YLENSPH zero) ***'
        STOP
      ENDIF

      DO IPOI=1,MPOI
        X(IPOI)=X(IPOI)-XOFF
      ENDDO

      OPEN(UNIT=LUNPHFIT,FILE=FILEPHFIT,STATUS='NEW')
      WRITE(LUNPHFIT,'(I5,1H ,A64)')ICODE,CODE
      WRITE(LUNPHFIT,*)PERLENPH,PHASEPH
      WRITE(LUNPHFIT,*)XLENCPH,YLENSPH
      WRITE(LUNPHFIT,*)XPHMIN,XPHMAX
      WRITE(LUNPHFIT,*)YPHMIN,YPHMAX
      WRITE(LUNPHFIT,*)ZPHMIN,ZPHMAX
      WRITE(LUNPHFIT,*)NTRANS0,NTRANS,NTRANSD
      WRITE(LUNPHFIT,*)NHARM0,NHARM,NHARMD
      WRITE(LUNPHFIT,*)NPARPH

      DO IPAR=1,NPARPH
        WRITE(LUNPHFIT,*) PARPH(IPAR)
      ENDDO
      CLOSE(LUNPHFIT)

      IF (IHBPHARM.NE.0) THEN
        CALL hbookm(NIDBPOLY,'BPHARM FIT',NTUP_P,'//WAVE',mpoi,CHTAGS_D)
      ENDIF

      DO IPOI=1,MPOI

        CALL BPHARM(X(IPOI),Y(IPOI),Z(IPOI)
     &    ,BXF(IPOI),BYF(IPOI),BZF(IPOI),AX,AY,AZ)

        IF (IHBPHARM.NE.0) THEN
          TUP_D(1)=X(IPOI)
          TUP_D(2)=Y(IPOI)
          TUP_D(3)=Z(IPOI)
          TUP_D(4)=BX(IPOI)
          TUP_D(5)=BY(IPOI)
          TUP_D(6)=BZ(IPOI)
          TUP_D(7)=BXF(IPOI)
          TUP_D(8)=BYF(IPOI)
          TUP_D(9)=BZF(IPOI)
          CALL hfm(NIDBPOLY,TUP_D)
        ENDIF

      ENDDO !IPOI

      CALL BRESI      (MPOI,NPHPOIX,NPHPOIY,NPHPOIZ
     &  ,X,Y,Z,BX,BY,BZ,BXF,BYF,BZF
     &  ,RESBX,RESBY,RESBZ,RESB
     &  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &  ,BERRMX,BXERRMX,BYERRMX,BZERRMX)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     RESB:',SNGL(RESB)
      WRITE(LUNGFO,*)'     RESBX, RESBY,RESBZ:',SNGL(RESBX),SNGL(RESBY),SNGL(RESBZ)
      WRITE(LUNGFO,*)
      IF (BAMEAN.NE.0.D0)
     &  WRITE(LUNGFO,*)'     |B|mean,   RESB/|B|mean: ',SNGL(BAMEAN),SNGL(RESB/BAMEAN)
      IF (BXAMEAN.NE.0.D0)
     &  WRITE(LUNGFO,*)'     |BX|mean, RESBX/|BX|mean:',SNGL(BXAMEAN),SNGL(RESBX/BXAMEAN)
      IF (BYAMEAN.NE.0.D0)
     &  WRITE(LUNGFO,*)'     |BY|mean, RESBY/|BY|mean:',SNGL(BYAMEAN),SNGL(RESBY/BYAMEAN)
      IF (BZAMEAN.NE.0.D0)
     &  WRITE(LUNGFO,*)'     |BZ|mean, RESBZ/|BZ|mean:',SNGL(BZAMEAN),SNGL(RESBZ/BZAMEAN)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     MAXIMUM ERRORS OF B, BX, BY AND BZ:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'  X      Y      Z       BX        BY        BZ       B         ERR'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,'(3F8.4,5(1PE10.2))')(SNGL(BERRMX(IX)),IX=1,8)
      WRITE(LUNGFO,'(3F8.4,5(1PE10.2))')(SNGL(BXERRMX(IX)),IX=1,8)
      WRITE(LUNGFO,'(3F8.4,5(1PE10.2))')(SNGL(BYERRMX(IX)),IX=1,8)
      WRITE(LUNGFO,'(3F8.4,5(1PE10.2))')(SNGL(BZERRMX(IX)),IX=1,8)
      WRITE(LUNGFO,*)

      IF (IHBPHARM.NE.0) THEN
        CALL MHROUT(NIDBPOLY,ICYCLE,' ')
        CALL hdeletm(NIDBPOLY)
      ENDIF

      RETURN
      END
