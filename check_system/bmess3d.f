*CMZ :  2.44/01 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  1.00/00 28/07/97  11.04.15  by  Michael Scheer
*CMZ : 00.02/00 15/11/96  11.30.41  by  Michael Scheer
*CMZ : 00.01/09 25/10/95  17.37.17  by  Michael Scheer
*-- Author :    Michael Scheer   29/09/95

      SUBROUTINE BMESS3D(X,Y,Z,BX,BY,BZ)

C--- ALL INDICES ACCORDING TO FORTRAN, BUT LORD3DG AND MORD3DG REFER TO MATH. INDICES

      IMPLICIT NONE

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
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEND.

      INTEGER IX,IY,IZ,NPOWP
     &         ,IWARNXMN,IWARNXMX,IWARNYMN,IWARNYMX,IWARNZMN,IWARNZMX
     &         ,IND,ICAL

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,XPOW1,YPOW1,ZPOW1,AX,AY,AZ

      PARAMETER (NPOWP=NDIMCG+1)
      DOUBLE PRECISION XPOW(NPOWP),YPOW(NPOWP),ZPOW(NPOWP)

      DATA IWARNXMN,IWARNXMX,IWARNYMN,IWARNYMX,IWARNZMN,IWARNZMX/6*0/
      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

      IF (MORD3DG.GT.NDIMCG-1) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BMESS3D: DIMENSION NDIMCG EXCEEDED ***'
          WRITE(6,*)
          STOP
      ENDIF

      IF (MORD3DG.GT.NPOWP) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BMESS3D: DIMENSION NPOWP EXCEEDED ***'
          WRITE(6,*)
          STOP
      ENDIF

      AX=0.D0
      AY=0.D0
      AZ=0.D0

      ICAL=1

      ENDIF

      IF (IWARNXMN.EQ.0.AND.X.LT.X3DMING) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BMESS3D: X OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'X .LT. X3DMING'
          WRITE(LUNGFO,*)'X, X3DMING:',X,X3DMING
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BMESS3D: X OUT OF FIT RANGE'
          WRITE(6,*)'X .LT. X3DMING'
          WRITE(6,*)'X, X3DMING:',X,X3DMING
          WRITE(6,*)
          IWARNXMN=1
      ENDIF

      IF (IWARNXMX.EQ.0.AND.X.GT.X3DMAXG) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BMESS3D: X OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'X .GT. X3DMAXG'
          WRITE(LUNGFO,*)'X, X3DMAXG:',X,X3DMAXG
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BMESS3D: X OUT OF FIT RANGE'
          WRITE(6,*)'X .GT. X3DMAXG'
          WRITE(6,*)'X, X3DMAXG:',X,X3DMAXG
          WRITE(6,*)
          IWARNXMX=1
      ENDIF

      IF (IWARNYMN.EQ.0.AND.Y.LT.Y3DMING) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BMESS3D: Y OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Y .LT. Y3DMING'
          WRITE(LUNGFO,*)'Y, Y3DMING:',Y,Y3DMING
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BMESS3D: Y OUT OF FIT RANGE'
          WRITE(6,*)'Y .LT. Y3DMING'
          WRITE(6,*)'Y, Y3DMING:',Y,Y3DMING
          WRITE(6,*)
          IWARNYMN=1
      ENDIF

      IF (IWARNYMX.EQ.0.AND.Y.GT.Y3DMAXG) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BMESS3D: Y OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Y .GT. Y3DMAXG'
          WRITE(LUNGFO,*)'Y, Y3DMAXG:',Y,Y3DMAXG
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BMESS3D: Y OUT OF FIT RANGE'
          WRITE(6,*)'Y .GT. Y3DMAXG'
          WRITE(6,*)'Y, Y3DMAXG:',Y,Y3DMAXG
          WRITE(6,*)
          IWARNYMX=1
      ENDIF

      IF (IWARNZMN.EQ.0.AND.Z.LT.Z3DMING) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BMESS3D: Z OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Z .LT. Z3DMING'
          WRITE(LUNGFO,*)'Z, Z3DMING:',Z,Z3DMING
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BMESS3D: Z OUT OF FIT RANGE'
          WRITE(6,*)'Z .LT. Z3DMING'
          WRITE(6,*)'Z, Z3DMING:',Z,Z3DMING
          WRITE(6,*)
          IWARNZMN=1
      ENDIF

      IF (IWARNZMX.EQ.0.AND.Z.GT.Z3DMAXG) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING SR BMESS3D: Z OUT OF FIT RANGE'
          WRITE(LUNGFO,*)'Z .GT. Z3DMAXG'
          WRITE(LUNGFO,*)'Z, Z3DMAXG:',Z,Z3DMAXG
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING SR BMESS3D: Z OUT OF FIT RANGE'
          WRITE(6,*)'Z .GT. Z3DMAXG'
          WRITE(6,*)'Z, Z3DMAXG:',Z,Z3DMAXG
          WRITE(6,*)
          IWARNZMX=1
      ENDIF

      XPOW(1)=1.D0
      DO IX=2,MORD3DG
          XPOW(IX)=XPOW(IX-1)*X
      ENDDO

      YPOW(1)=1.D0
      DO IY=2,MORD3DG
          YPOW(IY)=YPOW(IY-1)*Y
      ENDDO

      ZPOW(1)=1.D0
      DO IZ=2,MORD3DG
          ZPOW(IZ)=ZPOW(IZ-1)*Z
      ENDDO

      BX=0.D0
      BY=0.D0
      BZ=0.D0

      DO IND=1,NCINDG

             IX=ICINDG(1,IND)
             IY=ICINDG(2,IND)
             IZ=ICINDG(3,IND)

             IF (IX.GT.1) THEN
            XPOW1=XPOW(IX-1)
             ELSE
            XPOW1=1.D0
             ENDIF
             BX=BX-(IX-1)*CINDG(IND)*XPOW1*YPOW(IY)*ZPOW(IZ)

             IF (IY.GT.1) THEN
            YPOW1=YPOW(IY-1)
             ELSE
            YPOW1=1.D0
             ENDIF

             BY=BY-(IY-1)*CINDG(IND)*XPOW(IX)*YPOW1*ZPOW(IZ)

             IF (IZ.GT.1) THEN
            ZPOW1=ZPOW(IZ-1)
             ELSE
            ZPOW1=1.D0
             ENDIF

             BZ=BZ-(IZ-1)*CINDG(IND)*XPOW(IX)*YPOW(IY)*ZPOW1

      ENDDO

      RETURN
      END
