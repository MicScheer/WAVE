*CMZ :  2.44/01 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  22.20.23  by  Michael Scheer
*CMZ :  1.00/00 24/07/97  16.38.25  by  Michael Scheer
*CMZ : 00.02/00 21/11/96  14.40.10  by  Michael Scheer
*CMZ :  1.00/03 27/09/95  17.15.24  by  Michael Scheer
*CMZ :  1.00/01 25/09/95  18.20.26  by  Michael Scheer
*CMZ :  1.00/00 21/09/95  17.13.57  by  Michael Scheer
*-- Author :    Michael Scheer   21/09/95
      SUBROUTINE BMPOT3D(NPOI,X,Y,Z,BX,BY,BZ,IFAIL,COMMENT)
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

*KEEP,bmpot3du.
      include 'bmpot3du.cmn'
*KEND.

C---  TO FIT 3D POTENTIAL V=SUM( CG(I,J,K) * X**(I-1) * Y**(J-1) * Z**(K-1))
C     OF A MAGNETIC FIELD B=(BX,BY,BZ)=-GRAD(V)
C

C--- INPUT:

C     NPOI  : NUMBER OF DATA POINTS X,Y,Z,BX,BY,BZ

C--- OUTPUT:

C     IFAIL : FAILURE FLAG
C     COMMENT : COMMENT ON COEF-FILE (SR BMPOT3DINIT)

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEND.

      INTEGER NPOI,IFAIL,I,LUN,MTOT,MFIT

      DOUBLE PRECISION X(NPOI),Y(NPOI),Z(NPOI),BX(NPOI),BY(NPOI),BZ(NPOI)

      INTEGER NSTAKP
      PARAMETER (NSTAKP=4)

      CHARACTER(60) COMMENT,FILE

      INTEGER ICAL,NCOEF,ICOEF,IX,IY,IZ,IS,NS
     &         ,JX,JY,JZ


      DATA ICAL/0/
      DATA LUN/99/
      DATA FILE/'BMESS.COEF'/

      IF (ICAL.EQ.0) THEN

         ALLOCATE(ICG(NDIMCG,NDIMCG,NDIMCG))
         ALLOCATE(XPOW(NDIMCIPG))
         ALLOCATE(YPOW(NDIMCIPG))
         ALLOCATE(ZPOW(NDIMCIPG))
         ALLOCATE(VC(NDIMCIPG))
         ALLOCATE(A(NDIMCIPG,NDIMCIPG))
         ALLOCATE(WA(NDIMCIPG,NDIMCIPG))
         ALLOCATE(B(NDIMCIPG))
         ALLOCATE(WS(2*NDIMCIPG))
         ALLOCATE(WB(NDIMCIPG))
         ALLOCATE(INDEX(4,NDIMCIPG))
         ALLOCATE(FSTAK(4,NSTAKP,NDIMCIPG))

         MTOT=(MORD3DG*(MORD3DG+1)*(2*MORD3DG+1)/6
     &              +MORD3DG*(MORD3DG+1)*3/2+2*MORD3DG)/2

         MFIT=(MORD3DG+1)**2-1

      IF (MTOT.GT.NDIMCIPG) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'ERROR IN BMPOT3D: DIMENSION NDIMCIPG EXCEEDED ***'
          WRITE(LUNGFO,*)'CHECK INPUT OR BMPOLY3D.CMN'
          WRITE(6,*)
          WRITE(6,*)'ERROR IN BMPOT3D: DIMENSION NDIMPCIPG EXCEEDED ***'
          WRITE(6,*)'CHECK INPUT OR BMPOLY3D.CMN'
          WRITE(6,*)
          STOP
      ENDIF

      IF (MFIT.GT.NDIMCIPG) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'ERROR IN BMPOT3D: DIMENSION NDIMCIPG EXCEEDED ***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'ERROR IN BMPOT3D: DIMENSION NDIMCIPG EXCEEDED ***'
          WRITE(6,*)
          STOP
      ENDIF

          CALL BMPOT3DINIT(NDIMCIPG,LORD3DG,MORD3DG,NDORD3DG,
     &                       INDEX,NCOEF
     &                      ,NSTAKP,FSTAK,LUN,FILE,COMMENT)

      ENDIF !ICAL

      DO IZ=1,NDIMCG
      DO IY=1,NDIMCG
      DO IX=1,NDIMCG
          CG(IX,IY,IZ)=0.D0
          ICG(IX,IY,IZ)=0
      ENDDO
      ENDDO
      ENDDO

      CALL BMPOT3DFIT(NPOI,X,Y,Z,BX,BY,BZ,NDIMCIPG,NCOEF,MORD3DG
     &                       ,INDEX,VC,NSTAKP,FSTAK,A,B,WS,WA,WB
     &                       ,XPOW,YPOW,ZPOW,IFAIL)

C--- GET FITTED COEFFICIENTS

      DO ICOEF=1,NCOEF
          IX=INDEX(1,ICOEF)
          IY=INDEX(2,ICOEF)
          IZ=INDEX(3,ICOEF)
          CG(IX+1,IY+1,IZ+1)=VC(ICOEF)
          ICG(IX+1,IY+1,IZ+1)=1
      ENDDO

C--- CALCULATE OTHER COEFFICIENTS

      DO ICOEF=1,NCOEF
            NS=INDEX(4,ICOEF)
          IF (NS.GT.1) THEN
             IX=INDEX(1,ICOEF)
             IY=INDEX(2,ICOEF)
             IZ=INDEX(3,ICOEF)
             DO IS=1,NS
            JX=NINT(FSTAK(1,IS,ICOEF))
            JY=NINT(FSTAK(2,IS,ICOEF))
            JZ=NINT(FSTAK(3,IS,ICOEF))
            IF (IX.NE.JX.OR.IY.NE.JY.OR.IZ.NE.JZ) THEN
            IF (IX+1.GT.NDIMCG.OR.IY+1.GT.NDIMCG.OR.IZ+1.GT.NDIMCG) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'ERROR IN BMPOT3D: DIMENSION NDIMCG EXCEEDED ***'
                WRITE(LUNGFO,*)
                WRITE(6,*)
                WRITE(6,*)'ERROR IN BMPOT3D: DIMENSION NDIMCG EXCEEDED ***'
                WRITE(6,*)
                STOP
            ENDIF
            ICG(JX+1,JY+1,JZ+1)=2
            CG(JX+1,JY+1,JZ+1)=CG(JX+1,JY+1,JZ+1)
     &         +FSTAK(4,IS,ICOEF)*CG(IX+1,IY+1,IZ+1)
                 ENDIF
             ENDDO   !IS
         ENDIF
      ENDDO

      IF (ICAL.EQ.0) THEN

         NCINDG=0
              DO IZ=1,NDIMCG
              DO IY=1,NDIMCG
              DO IX=1,NDIMCG
                  IF (ICG(IX,IY,IZ).NE.0) THEN
            NCINDG=NCINDG+1
            ICINDG(1,NCINDG)=IX
            ICINDG(2,NCINDG)=IY
            ICINDG(3,NCINDG)=IZ
             ENDIF
              ENDDO
              ENDDO
              ENDDO

         WRITE(LUNGFO,*)'     SR BMPOT3D:'
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'     NUMBER OF FITTED COEFFICIENTS:    ',NCOEF
         WRITE(LUNGFO,*)'     TOTAL NUMBER OF COEFFICIENTS:',NCINDG
         WRITE(LUNGFO,*)

         IF (NCINDG.NE.MTOT.OR.NCOEF.NE.MFIT) THEN
             WRITE (6,*)
     &              '*** ERROR IN BMPOT3D: WRONG NUMBER OF COEFFICIENTS'
             WRITE (6,*)
     &              '*** ERROR IN BMPOT3D: CHECK FILE BMESS.COEF'
             WRITE (LUNGFO,*)
     &              '*** ERROR IN BMPOT3D: WRONG NUMBER OF COEFFICIENTS'
             WRITE (LUNGFO,*)
     &              '*** ERROR IN BMPOT3D: CHECK FILE BMESS.COEF'
             STOP
         ENDIF


         ICAL=1

      ENDIF !ICAL

      DO I=1,NCINDG
          CINDG(I)=CG(ICINDG(1,I),ICINDG(2,I),ICINDG(3,I))
      ENDDO

      RETURN
      END
