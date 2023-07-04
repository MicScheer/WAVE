*CMZ :  4.00/14 22/12/2021  16.40.21  by  Michael Scheer
*CMZ :  4.00/04 23/08/2019  15.47.38  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  11.54.30  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.42/04 28/10/2002  16.57.51  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.30.59  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.39.57  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.27  by  Michael Scheer
*CMZ : 00.02/03 22/01/97  14.23.26  by  Michael Scheer
*CMZ : 00.02/00 25/11/96  10.03.32  by  Michael Scheer
*CMZ : 00.01/09 24/10/95  17.31.19  by  Michael Scheer
*-- Author :    Michael Scheer   22/09/95

      SUBROUTINE BPOLY3DFIT
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

*KEEP,bpoly3df90u.
      include 'bpoly3df90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly3d.
      include 'bpoly3d.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      CHARACTER(60) COMMENT

      INTEGER IPOI,IFAIL,IX,IY,IZ
      INTEGER MPOI,IREAD


      DOUBLE PRECISION
     &                   RESBX,RESBY,RESBZ,RESB
     &                  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &                  ,BERRMX(8),BXERRMX(8),BYERRMX(8),BZERRMX(8)
     &                  ,XX,YY,ZZ,BBX,BBY,BBZ,AX,AY,AZ

C     DOUBLE PRECISION DIV

      INTEGER NTUP_P,ICYCLE
        PARAMETER (NTUP_P=9)
        REAL*8 TUP_D(NTUP_P)
        CHARACTER(3) CHTAGS_D(NTUP_P)
      data chtags_d/'x','y','z','bx','by','bz','bxf','byf','bzf'/

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE BPOLY3DFIT:'
      WRITE(LUNGFO,*)'     ======================'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     order of fit (LORD3D, MORD3D, NDORD3D):'
     &                ,LORD3D, MORD3D, NDORD3D
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     considered range in x, y, and z:'
      WRITE(LUNGFO,*)'     ',X3DMIN,X3DMAX
      WRITE(LUNGFO,*)'     ',Y3DMIN,Y3DMAX
      WRITE(LUNGFO,*)'     ',Z3DMIN,Z3DMAX
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     field map read from file:'
      WRITE(LUNGFO,*)'     ',FILEBMAP
      WRITE(LUNGFO,*)

C--- READ FIELD MAP {

      N3DPOIX=0
      N3DPOIY=0
      N3DPOIZ=0
      MPOI=0
      IREAD=0

      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='OLD')

100       CONTINUE
          IF (IWBPOLY3D.EQ.1) THEN
             READ(LUNBMAP,*,END=90)XX
             READ(LUNBMAP,*)YY
             READ(LUNBMAP,*)ZZ
             READ(LUNBMAP,*)BBX
             READ(LUNBMAP,*)BBY
             READ(LUNBMAP,*)BBZ
             READ(LUNBMAP,*)
          ELSEIF (IWBPOLY3D.EQ.2) THEN
             READ(LUNBMAP,*,END=90)XX,YY,ZZ,BBX,BBY,BBZ
          ELSE
         WRITE(6,*) '*** ERROR IN BPOLY3DFIT: FLAG IWBPOLY3D WRONG, CHECK INPUT'
         STOP
          ENDIF
          IREAD=IREAD+1
          IF (X3DMIN.NE.9999.AND.XX.LT.X3DMIN) GOTO 100
          IF (X3DMAX.NE.9999.AND.XX.GT.X3DMAX) GOTO 100
          IF (Y3DMIN.NE.9999.AND.YY.LT.Y3DMIN) GOTO 100
          IF (Y3DMAX.NE.9999.AND.YY.GT.Y3DMAX) GOTO 100
          IF (Z3DMIN.NE.9999.AND.ZZ.LT.Z3DMIN) GOTO 100
          IF (Z3DMAX.NE.9999.AND.ZZ.GT.Z3DMAX) GOTO 100
          IF (BBX.EQ.-9999..AND.BBY.EQ.-9999..AND.BBZ.EQ.-9999.) GOTO 100
          MPOI=MPOI+1
          IF (BBX.NE.-9999.) N3DPOIX=N3DPOIX+1
          IF (BBY.NE.-9999.) N3DPOIY=N3DPOIY+1
          IF (BBZ.NE.-9999.) N3DPOIZ=N3DPOIZ+1
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

      N3DPOIX=0
      N3DPOIY=0
      N3DPOIZ=0
      MPOI=0
      IREAD=0


      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='OLD')

101       CONTINUE
          IF (IWBPOLY3D.EQ.1) THEN
             READ(LUNBMAP,*,END=91)XX
             READ(LUNBMAP,*)YY
             READ(LUNBMAP,*)ZZ
             READ(LUNBMAP,*)BBX
             READ(LUNBMAP,*)BBY
             READ(LUNBMAP,*)BBZ
             READ(LUNBMAP,*)
          ELSEIF (IWBPOLY3D.EQ.2) THEN
             READ(LUNBMAP,*,END=91)XX,YY,ZZ,BBX,BBY,BBZ
          ELSE
         WRITE(6,*) '*** ERROR IN BPOLY3DFIT: FLAG IWBPOLY3D WRONG, CHECK INPUT'
         STOP
          ENDIF
          IREAD=IREAD+1
          IF (X3DMIN.NE.9999.AND.XX.LT.X3DMIN) GOTO 101
          IF (X3DMAX.NE.9999.AND.XX.GT.X3DMAX) GOTO 101
          IF (Y3DMIN.NE.9999.AND.YY.LT.Y3DMIN) GOTO 101
          IF (Y3DMAX.NE.9999.AND.YY.GT.Y3DMAX) GOTO 101
          IF (Z3DMIN.NE.9999.AND.ZZ.LT.Z3DMIN) GOTO 101
          IF (Z3DMAX.NE.9999.AND.ZZ.GT.Z3DMAX) GOTO 101
          IF (BBX.EQ.-9999..AND.BBY.EQ.-9999..AND.BBZ.EQ.-9999.) GOTO 101
          MPOI=MPOI+1
          IF (BBX.NE.-9999.) N3DPOIX=N3DPOIX+1
          IF (BBY.NE.-9999.) N3DPOIY=N3DPOIY+1
          IF (BBZ.NE.-9999.) N3DPOIZ=N3DPOIZ+1
          X(MPOI)=XX*XYZ3DSC
          Y(MPOI)=YY*XYZ3DSC
          Z(MPOI)=ZZ*XYZ3DSC
          BX(MPOI)=BBX
          BY(MPOI)=BBY
          BZ(MPOI)=BBZ
          GOTO 101
91    CLOSE(LUNBMAP)

      N3DPOI=MPOI

      WRITE(LUNGFO,*)'     number of data points read:    ',IREAD
      WRITE(LUNGFO,*)'     number of data points accepted:',MPOI
      WRITE(LUNGFO,*)'     number of Bx, By, Bz .ne. -9999.:'
     &                      ,N3DPOIX,N3DPOIY,N3DPOIZ
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     fitted coefficients are written to file:'
      WRITE(LUNGFO,*)'     ',FILE3DFIT
      WRITE(LUNGFO,*)

C--- READ FIELD MAP }

      CALL VPOT3D(MPOI,X,Y,Z,BX,BY,BZ,IFAIL,COMMENT)

      IF (IFAIL.NE.0) THEN
          WRITE(6,*) '*** IFAIL NOT ZERO ***'
          STOP
      ENDIF

      IF (X3DMIN.EQ.9999.) THEN
          X3DMIN=1.D30
          DO IPOI=1,MPOI
         IF (X(IPOI).LT.X3DMIN) X3DMIN=X(IPOI)
          ENDDO
          X3DMIN=X3DMIN/XYZ3DSC
      ENDIF

      IF (X3DMAX.EQ.9999.) THEN
          X3DMAX=-1.D30
          DO IPOI=1,MPOI
         IF (X(IPOI).GT.X3DMAX) X3DMAX=X(IPOI)
          ENDDO
          X3DMAX=X3DMAX/XYZ3DSC
      ENDIF

      IF (Y3DMIN.EQ.9999.) THEN
          Y3DMIN=1.D30
          DO IPOI=1,MPOI
         IF (Y(IPOI).LT.Y3DMIN) Y3DMIN=Y(IPOI)
          ENDDO
          Y3DMIN=Y3DMIN/XYZ3DSC
      ENDIF

      IF (Y3DMAX.EQ.9999.) THEN
          Y3DMAX=-1.D30
          DO IPOI=1,MPOI
         IF (Y(IPOI).GT.Y3DMAX) Y3DMAX=Y(IPOI)
          ENDDO
          Y3DMAX=Y3DMAX/XYZ3DSC
      ENDIF

      IF (Z3DMIN.EQ.9999.) THEN
          Z3DMIN=1.D30
          DO IPOI=1,MPOI
         IF (Z(IPOI).LT.Z3DMIN) Z3DMIN=Z(IPOI)
          ENDDO
          Z3DMIN=Z3DMIN/XYZ3DSC
      ENDIF

      IF (Z3DMAX.EQ.9999.) THEN
          Z3DMAX=-1.D30
          DO IPOI=1,MPOI
         IF (Z(IPOI).GT.Z3DMAX) Z3DMAX=Z(IPOI)
          ENDDO
          Z3DMAX=Z3DMAX/XYZ3DSC
      ENDIF

      OPEN(UNIT=LUN3DFIT,FILE=FILE3DFIT,STATUS='NEW')
          WRITE(LUN3DFIT,'(I5,1H ,A64)')ICODE,CODE
          WRITE(LUN3DFIT,'(A60)')COMMENT
          WRITE(LUN3DFIT,*)XYZ3DSC
          WRITE(LUN3DFIT,*)X3DMIN,X3DMAX
          WRITE(LUN3DFIT,*)Y3DMIN,Y3DMAX
          WRITE(LUN3DFIT,*)Z3DMIN,Z3DMAX
          DO IZ=1,MORD3D+1
          DO IY=1,MORD3D+1
          DO IX=1,MORD3D+1
             IF (C(IX,IY,IZ).NE.0.0D0) WRITE(LUN3DFIT,*)IX,IY,IZ,C(IX,IY,IZ)
          ENDDO
          ENDDO
          ENDDO
      CLOSE(LUN3DFIT)

      CALL hbookm(NIDBPOLY,'BPOLY3D FIT',NTUP_P,'//WAVE',mpoi,CHTAGS_D)

      DO IPOI=1,MPOI

          X(IPOI)=X(IPOI)/XYZ3DSC
          Y(IPOI)=Y(IPOI)/XYZ3DSC
          Z(IPOI)=Z(IPOI)/XYZ3DSC

           CALL BPOLY3D(X(IPOI),Y(IPOI),Z(IPOI)
     &                ,BXF(IPOI),BYF(IPOI),BZF(IPOI),AX,AY,AZ)

C       CALL DIVB(X(IPOI),Y(IPOI),Z(IPOI),NDIMC+1,C
C     &             ,LORD3D,MORD3D,NDORD3D,DIV)

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

      ENDDO !IPOI

      CALL BRESI      (MPOI,N3DPOIX,N3DPOIY,N3DPOIZ
     &                  ,X,Y,Z,BX,BY,BZ,BXF,BYF,BZF
     &                  ,RESBX,RESBY,RESBZ,RESB
     &                  ,BXAMEAN,BYAMEAN,BZAMEAN,BAMEAN
     &                  ,BERRMX,BXERRMX,BYERRMX,BZERRMX)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     RESB:',SNGL(RESB)
      WRITE(LUNGFO,*)'     RESBX, RESBY,RESBZ:',SNGL(RESBX),SNGL(RESBY),SNGL(RESBZ)
      WRITE(LUNGFO,*)
      IF (BAMEAN.NE.0.D0)
     &WRITE(LUNGFO,*)'     |B|mean, RESB/|B|mean:',SNGL(BAMEAN),SNGL(RESB/BAMEAN)
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

      CALL MHROUT(NIDBPOLY,ICYCLE,' ')
      CALL hdeletm(NIDBPOLY)

      RETURN
      END
