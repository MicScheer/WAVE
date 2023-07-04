*CMZ :  4.00/14 22/12/2021  16.40.21  by  Michael Scheer
*CMZ :  4.00/04 23/08/2019  15.47.38  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  11.54.30  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.42/04 28/10/2002  16.57.50  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.42.22  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  13.48.28  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.02/04 11/02/97  11.55.40  by  Michael Scheer
*CMZ : 00.02/03 31/01/97  15.06.26  by  Michael Scheer
*CMZ :          22/01/97  11.27.20  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97

      SUBROUTINE BPOLY2DHFIT
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

*KEEP,bpoly2dhf90u.
      include 'bpoly2dhf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly2dh.
      include 'bpoly2dh.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      INTEGER IPOI,IFAIL,IX,IPAR
      INTEGER MPOI,IREAD
      INTEGER N2DHPOIX,N2DHPOIY,N2DHPOIZ,N2DHPOI


      DOUBLE PRECISION
     &                   RESBX,RESBY,RESBZ,RESB
     &                  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &                  ,BERRMX(8),BXERRMX(8),BYERRMX(8),BZERRMX(8)
     &                  ,XX,YY,ZZ,BBX,BBY,BBZ,AX,AY,AZ

      DOUBLE PRECISION XOFF

      INTEGER NTUP_P,ICYCLE
        PARAMETER (NTUP_P=9)
        REAL*8 TUP_D(NTUP_P)
        CHARACTER(3) CHTAGS_D(NTUP_P)
      data chtags_d/'x','y','z','bx','by','bz','bxf','byf','bzf'/

      NPAR2DH=NPAR2DHP
      NPARTOT=(NORD2DH+1)/2*NPAR2DHP*2
      IF (NORD2DH.GT.NORD2DHP) THEN
          WRITE(6,*) '*** ERROR IN BPOLY2DHFIT: DIMENSION NORD2DHP EXCEEDED'
          STOP
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE BPOLY2DHFIT:'
      WRITE(LUNGFO,*)'     ======================'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Period length and long. shift:'
     &                 , SNGL(PERLEN2DH),SNGL(PHASE2DH)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     field map read from file:'
      WRITE(LUNGFO,*)'     ',FILEBMAP
      WRITE(LUNGFO,*)

C--- READ FIELD MAP {

      XOFF=PERLEN2DH*PHASE2DH
      PERLEN2DH=PERLEN2DH*XYZ2DH

      N2DHPOIX=0
      N2DHPOIY=0
      N2DHPOIZ=0
      MPOI=0
      IREAD=0

      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='OLD')

100       CONTINUE
          IF (IWBPOLY2DH.EQ.1) THEN
             READ(LUNBMAP,*,END=90)XX
             READ(LUNBMAP,*)YY
             READ(LUNBMAP,*)ZZ
             READ(LUNBMAP,*)BBX
             READ(LUNBMAP,*)BBY
             READ(LUNBMAP,*)BBZ
             READ(LUNBMAP,*)
          ELSEIF (IWBPOLY2DH.EQ.2) THEN
             READ(LUNBMAP,*,END=90)XX,YY,ZZ,BBX,BBY,BBZ
          ELSE
         WRITE(6,*) '*** ERROR IN BPOLY2DHFIT: FLAG IWBPOLY2DH WRONG, CHECK INPUT'
         STOP
          ENDIF
          IREAD=IREAD+1
          IF (X2DHMIN.NE.9999.AND.XX.LT.X2DHMIN) GOTO 100
          IF (X2DHMAX.NE.9999.AND.XX.GT.X2DHMAX) GOTO 100
          IF (Y2DHMIN.NE.9999.AND.YY.LT.Y2DHMIN) GOTO 100
          IF (Y2DHMAX.NE.9999.AND.YY.GT.Y2DHMAX) GOTO 100
          IF (Z2DHMIN.NE.9999.AND.ZZ.LT.Z2DHMIN) GOTO 100
          IF (Z2DHMAX.NE.9999.AND.ZZ.GT.Z2DHMAX) GOTO 100
          IF (BBX.EQ.-9999..AND.BBY.EQ.-9999..AND.BBZ.EQ.-9999.) GOTO 100
          MPOI=MPOI+1
          IF (BBX.NE.-9999.) N2DHPOIX=N2DHPOIX+1
          IF (BBY.NE.-9999.) N2DHPOIY=N2DHPOIY+1
          IF (BBZ.NE.-9999.) N2DHPOIZ=N2DHPOIZ+1
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

      N2DHPOIX=0
      N2DHPOIY=0
      N2DHPOIZ=0
      MPOI=0
      IREAD=0

      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='OLD')

101       CONTINUE
          IF (IWBPOLY2DH.EQ.1) THEN
             READ(LUNBMAP,*,END=91)XX
             READ(LUNBMAP,*)YY
             READ(LUNBMAP,*)ZZ
             READ(LUNBMAP,*)BBX
             READ(LUNBMAP,*)BBY
             READ(LUNBMAP,*)BBZ
             READ(LUNBMAP,*)
          ELSEIF (IWBPOLY2DH.EQ.2) THEN
             READ(LUNBMAP,*,END=91)XX,YY,ZZ,BBX,BBY,BBZ
          ELSE
         WRITE(6,*) '*** ERROR IN BPOLY2DHFIT: FLAG IWBPOLY2DH WRONG, CHECK INPUT'
         STOP
          ENDIF
          IREAD=IREAD+1
          IF (X2DHMIN.NE.9999.AND.XX.LT.X2DHMIN) GOTO 101
          IF (X2DHMAX.NE.9999.AND.XX.GT.X2DHMAX) GOTO 101
          IF (Y2DHMIN.NE.9999.AND.YY.LT.Y2DHMIN) GOTO 101
          IF (Y2DHMAX.NE.9999.AND.YY.GT.Y2DHMAX) GOTO 101
          IF (Z2DHMIN.NE.9999.AND.ZZ.LT.Z2DHMIN) GOTO 101
          IF (Z2DHMAX.NE.9999.AND.ZZ.GT.Z2DHMAX) GOTO 101
          IF (BBX.EQ.-9999..AND.BBY.EQ.-9999..AND.BBZ.EQ.-9999.) GOTO 101
          MPOI=MPOI+1
          IF (BBX.NE.-9999.) N2DHPOIX=N2DHPOIX+1
          IF (BBY.NE.-9999.) N2DHPOIY=N2DHPOIY+1
          IF (BBZ.NE.-9999.) N2DHPOIZ=N2DHPOIZ+1
          X(MPOI)=(XX+XOFF)*XYZ2DH
          Y(MPOI)=YY*XYZ2DH
          Z(MPOI)=ZZ*XYZ2DH
          BX(MPOI)=BBX
          BY(MPOI)=BBY
          BZ(MPOI)=BBZ
          GOTO 101
91    CLOSE(LUNBMAP)

      N2DHPOI=MPOI

      WRITE(LUNGFO,*)'     number of data points read:      ',IREAD
      WRITE(LUNGFO,*)'     number of data points accepted:  ',MPOI
      WRITE(LUNGFO,*)'     number of Bx, By, Bz .ne. -9999.:'
     &                      ,N2DHPOIX,N2DHPOIY,N2DHPOIZ
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     fitted coefficients are written to file:'
      WRITE(LUNGFO,*)'     ',FILE2DHFIT
      WRITE(LUNGFO,*)

C--- READ FIELD MAP }

      CALL VPOT2DH(MPOI,IFAIL)
      IF (IFAIL.NE.0) THEN
          WRITE(6,*) '*** IFAIL NOT ZERO ***'
          STOP
      ENDIF

      DO IPOI=1,MPOI
          X(IPOI)=X(IPOI)-XOFF*XYZ2DH
      ENDDO

      IF (X2DHMIN.EQ.9999.) THEN
          X2DHMIN=1.D30
          DO IPOI=1,MPOI
         IF (X(IPOI).LT.X2DHMIN) X2DHMIN=X(IPOI)
          ENDDO
          X2DHMIN=X2DHMIN/XYZ2DH
      ENDIF

      IF (X2DHMAX.EQ.9999.) THEN
          X2DHMAX=-1.D30
          DO IPOI=1,MPOI
         IF (X(IPOI).GT.X2DHMAX) X2DHMAX=X(IPOI)
          ENDDO
          X2DHMAX=X2DHMAX/XYZ2DH
      ENDIF

      IF (Y2DHMIN.EQ.9999.) THEN
          Y2DHMIN=1.D30
          DO IPOI=1,MPOI
         IF (Y(IPOI).LT.Y2DHMIN) Y2DHMIN=Y(IPOI)
          ENDDO
          Y2DHMIN=Y2DHMIN/XYZ2DH
      ENDIF

      IF (Y2DHMAX.EQ.9999.) THEN
          Y2DHMAX=-1.D30
          DO IPOI=1,MPOI
         IF (Y(IPOI).GT.Y2DHMAX) Y2DHMAX=Y(IPOI)
          ENDDO
          Y2DHMAX=Y2DHMAX/XYZ2DH
      ENDIF

      IF (Z2DHMIN.EQ.9999.) THEN
          Z2DHMIN=1.D30
          DO IPOI=1,MPOI
         IF (Z(IPOI).LT.Z2DHMIN) Z2DHMIN=Z(IPOI)
          ENDDO
          Z2DHMIN=Z2DHMIN/XYZ2DH
      ENDIF

      IF (Z2DHMAX.EQ.9999.) THEN
          Z2DHMAX=-1.D30
          DO IPOI=1,MPOI
         IF (Z(IPOI).GT.Z2DHMAX) Z2DHMAX=Z(IPOI)
          ENDDO
          Z2DHMAX=Z2DHMAX/XYZ2DH
      ENDIF

          PERLEN2DH=PERLEN2DH/XYZ2DH

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Considered range in x, y, and z:'
      WRITE(LUNGFO,*)'     ',X2DHMIN,X2DHMAX
      WRITE(LUNGFO,*)'     ',Y2DHMIN,Y2DHMAX
      WRITE(LUNGFO,*)'     ',Z2DHMIN,Z2DHMAX
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Harmonical order of fit:',NORD2DH
      WRITE(LUNGFO,*)

      OPEN(UNIT=LUN2DHFIT,FILE=FILE2DHFIT,STATUS='NEW')
          WRITE(LUN2DHFIT,'(I5,1H ,A64)')ICODE,CODE
          WRITE(LUN2DHFIT,*)PERLEN2DH,PHASE2DH,XYZ2DH
          WRITE(LUN2DHFIT,*)X2DHMIN,X2DHMAX
          WRITE(LUN2DHFIT,*)Y2DHMIN,Y2DHMAX
          WRITE(LUN2DHFIT,*)Z2DHMIN,Z2DHMAX
          WRITE(LUN2DHFIT,*)NORD2DH,NPAR2DH
          DO IPAR=1,NPARTOT
         WRITE(LUN2DHFIT,*) PAR2DH(IPAR)
          ENDDO
      CLOSE(LUN2DHFIT)

      IF (IHBPOLY2DH.NE.0) THEN
      CALL hbookm(NIDBPOLY,'BPOLY2DH FIT$',NTUP_P,'//WAVE',mpoi,CHTAGS_D)
      ENDIF

      DO IPOI=1,MPOI

          X(IPOI)=X(IPOI)/XYZ2DH
          Y(IPOI)=Y(IPOI)/XYZ2DH
          Z(IPOI)=Z(IPOI)/XYZ2DH

           CALL BPOLY2DH(X(IPOI),Y(IPOI),Z(IPOI)
     &                ,BXF(IPOI),BYF(IPOI),BZF(IPOI),AX,AY,AZ)

      IF (IHBPOLY2DH.NE.0) THEN
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

      CALL BRESI      (MPOI,N2DHPOIX,N2DHPOIY,N2DHPOIZ
     &                  ,X,Y,Z,BX,BY,BZ,BXF,BYF,BZF
     &                  ,RESBX,RESBY,RESBZ,RESB
     &                  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &                  ,BERRMX,BXERRMX,BYERRMX,BZERRMX)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     RESB:',SNGL(RESB)
      WRITE(LUNGFO,*)'     RESBX, RESBY,RESBZ:',SNGL(RESBX),SNGL(RESBY),SNGL(RESBZ)
      WRITE(LUNGFO,*)
      IF (BAMEAN.NE.0.D0)
     &WRITE(LUNGFO,*)'     |B|mean,   RESB/|B|mean: ',SNGL(BAMEAN),SNGL(RESB/BAMEAN)
      IF (BXAMEAN.NE.0.D0)
     &WRITE(LUNGFO,*)'     |BX|mean, RESBX/|BX|mean:',SNGL(BXAMEAN),SNGL(RESBX/BXAMEAN)
      IF (BYAMEAN.NE.0.D0)
     &WRITE(LUNGFO,*)'     |BY|mean, RESBY/|BY|mean:',SNGL(BYAMEAN),SNGL(RESBY/BYAMEAN)
      IF (BZAMEAN.NE.0.D0)
     &WRITE(LUNGFO,*)'     |BZ|mean, RESBZ/|BZ|mean:',SNGL(BZAMEAN),SNGL(RESBZ/BZAMEAN)
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

      IF (IHBPOLY2DH.NE.0) THEN
        CALL MHROUT(NIDBPOLY,ICYCLE,' ')
        CALL hdeletm(NIDBPOLY)
      ENDIF

      RETURN
      END
