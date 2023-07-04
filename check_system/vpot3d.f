*CMZ :  2.69/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.41/10 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  22.20.56  by  Michael Scheer
*CMZ : 00.02/00 13/11/96  16.46.41  by  Michael Scheer
*CMZ :  1.00/03 27/09/95  17.15.24  by  Michael Scheer
*CMZ :  1.00/01 25/09/95  18.20.26  by  Michael Scheer
*CMZ :  1.00/00 21/09/95  17.13.57  by  Michael Scheer
*-- Author :    Michael Scheer   21/09/95
      SUBROUTINE VPOT3D(NPOI,X,Y,Z,BX,BY,BZ,IFAIL,COMMENT)
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

C---  TO FIT 3D POTENTIAL V=SUM( C(I,J,K) * X**(I-1) * Y**(J-1) * Z**(K-1))
C     OF A MAGNETIC FIELD B=(BX,BY,BZ)=-GRAD(V)
C

C--- INPUT:

C     NPOI  : NUMBER OF DATA POINTS X,Y,Z,BX,BY,BZ

C--- OUTPUT:

C     IFAIL : FAILURE FLAG
C     COMMENT : COMMENT ON COEF-FILE (SR VPOT3DINIT)

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly3d.
      include 'bpoly3d.cmn'
*KEND.

      INTEGER NPOI,IFAIL,MCOEF

      DOUBLE PRECISION X(NPOI),Y(NPOI),Z(NPOI),BX(NPOI),BY(NPOI),BZ(NPOI)


      INTEGER NDIMP,NSTAKP
      PARAMETER (NDIMP=256,NSTAKP=16)

        CHARACTER(60) COMMENT

      INTEGER JDIMP
      PARAMETER (JDIMP=NDIMC+1)
      INTEGER ICAL,NCOEF,INDEX(4,NDIMP),ICOEF,IX,IY,IZ,IS,NS
     &         ,JX,JY,JZ,JNDEX(JDIMP,JDIMP,JDIMP)

      DOUBLE PRECISION VC(NDIMP),A(NDIMP,NDIMP),B(NDIMP),WS(2*NDIMP)
     &                  ,WA(NDIMP,NDIMP),WB(NDIMP)
     &                  ,XPOW(NDIMP),YPOW(NDIMP),ZPOW(NDIMP)
     &                  ,FSTAK(4,NSTAKP,NDIMP)


      DATA ICAL/0/

      DO IZ=1,MORD3D+1
      DO IY=1,MORD3D+1
      DO IX=1,MORD3D+1
          JNDEX(IX,IY,IZ)=-1
      ENDDO
      ENDDO
      ENDDO

      IF (ICAL.EQ.0) THEN
          CALL VPOT3DINIT(NDIMP,LORD3D,MORD3D,NDORD3D,INDEX,NCOEF,NSTAKP,FSTAK
     &                        ,LUN3DCOE,FILE3DCOE,COMMENT)
      ENDIF

      IF (NCOEF.GT.NDIMP) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'ERROR IN VPOT3D: DIMENSION NDIMP EXCEEDED ***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'ERROR IN VPOT3D: DIMENSION NDIMP EXCEEDED ***'
          WRITE(6,*)
          STOP
      ENDIF

      CALL VPOT3DFIT(NPOI,X,Y,Z,BX,BY,BZ,NDIMP,NCOEF,MORD3D
     &                       ,INDEX,VC,NSTAKP,FSTAK,A,B,WS,WA,WB
     &                       ,XPOW,YPOW,ZPOW,IFAIL)

      DO IZ=1,NDIMC
      DO IY=1,NDIMC
      DO IX=1,NDIMC
          C(IX,IY,IZ)=0.D0
      ENDDO
      ENDDO
      ENDDO

C--- GET FITTED COEFFICIENTS

      DO ICOEF=1,NCOEF
          IX=INDEX(1,ICOEF)
          IY=INDEX(2,ICOEF)
          IZ=INDEX(3,ICOEF)
          IF (IX+1.GT.NDIMC.OR.IY+1.GT.NDIMC.OR.IZ+1.GT.NDIMC) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'ERROR IN VPOT3D: DIMENSION NDIMC EXCEEDED ***'
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)'ERROR IN VPOT3D: DIMENSION NDIMC EXCEEDED ***'
             WRITE(6,*)
             STOP
          ENDIF
          C(IX+1,IY+1,IZ+1)=VC(ICOEF)
          JNDEX(IX+1,IY+1,IZ+1)=1
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
            IF (IX+1.GT.NDIMC.OR.IY+1.GT.NDIMC.OR.IZ+1.GT.NDIMC) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'ERROR IN VPOT3D: DIMENSION NDIMC EXCEEDED ***'
                WRITE(LUNGFO,*)
                WRITE(6,*)
                WRITE(6,*)'ERROR IN VPOT3D: DIMENSION NDIMC EXCEEDED ***'
                WRITE(6,*)
                STOP
            ENDIF
            C(JX+1,JY+1,JZ+1)=C(JX+1,JY+1,JZ+1)
     &                 +FSTAK(4,IS,ICOEF)*C(IX+1,IY+1,IZ+1)
            JNDEX(JX+1,JY+1,JZ+1)=1
                 ENDIF
             ENDDO   !IS
         ENDIF
      ENDDO

      MCOEF=0
      DO IZ=1,MORD3D+1
      DO IY=1,MORD3D+1
      DO IX=1,MORD3D+1
          IF (JNDEX(IX,IY,IZ).EQ.1) MCOEF=MCOEF+1
      ENDDO
      ENDDO
      ENDDO

      IF (ICAL.EQ.0) THEN
      WRITE(LUNGFO,*)'     SR VPOT3D:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     number of fitted coefficients:    ',NCOEF
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     number of fitted and calc. coeff.:',MCOEF
      WRITE(LUNGFO,*)
          ICAL=1
      ENDIF

      RETURN
      END
