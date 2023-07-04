*CMZ :  2.45/03 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.45/01 13/12/2002  17.24.23  by  Michael Scheer
*CMZ :  2.44/02 12/12/2002  16.04.18  by  Michael Scheer
*CMZ :  2.44/01 12/12/2002  11.07.45  by  Michael Scheer
*CMZ :  2.44/00 07/11/2002  15.12.01  by  Michael Scheer
*CMZ :  2.42/04 14/09/2002  07.16.48  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.39/01 15/01/2002  16.47.24  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.11  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.02/03 23/01/97  17.29.29  by  Michael Scheer
*CMZ : 00.02/00 25/11/96  09.51.41  by  Michael Scheer
*-- Author :    Michael Scheer   28/09/95
      SUBROUTINE WBMAP4
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

C--- TO WRITE 3D FIELD MAP IN TERMS OF POLYNOMIAL COEFFICIENTS TO FILE


      IMPLICIT NONE

      EXTERNAL DCOSD,DSIND
      DOUBLE PRECISION DCOSD,DSIND

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEEP,bmap.
      include 'bmap.cmn'
*KEND.

      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: X,Y,Z,BX,BY,BZ
     &                                                ,DXX,DYY,DZZ

      DOUBLE PRECISION DX,DY,DZ,AX,AY,AZ,SCALXYZ

      INTEGER IX,IY,IZ,NPOI,IX1,IY1,IZ1,IX2,IY2,IZ2,JX,JY,JZ,IND,IFAIL
      INTEGER LX1,LY1,LZ1,LX2,LY2,LZ2

      CHARACTER(60) COMMENT

      DATA DX,DY,DZ/0.D0,0.D0,0.D0/
      DATA SCALXYZ/100.D0/

      IF (IBMRADIAL.NE.0) THEN
          WRITE(6,*) '*** ERROR IN WBMAP4: IBMRADIAL.NE.0'
          WRITE(LUNGFO,*) '*** ERROR IN WBMAP4: IBMRADIAL.NE.0'
          STOP
      ENDIF

      LORD3DG=1
      NDORD3DG=1

        IF (NMAPX.EQ.-9999) NMAPX=NINT((XSTOP-XSTART)*MYINUM)+1

      ALLOCATE(X(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(Y(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(Z(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(DXX(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(DYY(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(DZZ(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(BX(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(BY(NBMDATX*NBMDATY*NBMDATZ))
      ALLOCATE(BZ(NBMDATX*NBMDATY*NBMDATZ))

        IF (XMAPMN.EQ.9999.) XMAPMN=XSTART
        IF (XMAPMX.EQ.9999.) XMAPMX=XSTOP

        IF (IWBMAPEXT.EQ.0.AND.NMAPX.LT.NBMDATX) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: NMAPX.LT.NBMDATX'
            WRITE(LUNGFO,*)'CHECK NAMELISTS BMAP AND BGRIDN IN WAVE.IN'
            WRITE(6,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WBMAP4: NMAPX.LT.NBMDATX'
            WRITE(6,*)'CHECK NAMELISTS BMAP AND BGRIDN IN WAVE.IN'
            WRITE(6,*)
            STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

        IF (IWBMAPEXT.EQ.0.AND.NMAPY.LT.NBMDATY) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: NMAPY.LT.NBMDATY'
            WRITE(LUNGFO,*)'CHECK NAMELISTS BMAP AND BGRIDN IN WAVE.IN'
            WRITE(6,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WBMAP4: NMAPY.LT.NBMDATY'
            WRITE(6,*)'CHECK NAMELISTS BMAP AND BGRIDN IN WAVE.IN'
            WRITE(6,*)
            STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

        IF (IWBMAPEXT.EQ.0.AND.NMAPZ.LT.NBMDATZ) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: NMAPZ.LT.NBMDATZ'
            WRITE(LUNGFO,*)'CHECK NAMELISTS BMAP AND BGRIDN IN WAVE.IN'
            WRITE(6,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WBMAP4: NMAPZ.LT.NBMDATZ'
            WRITE(6,*)'CHECK NAMELISTS BMAP AND BGRIDN IN WAVE.IN'
            WRITE(6,*)
            STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

        IF (XMAPMX.LT.XMAPMN) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: XMAPMX.LT.XMAPMN'
            WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WBMAP4: XMAPMX.LT.XMAPMN'
            WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
            WRITE(6,*)
            STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

        IF (YMAPMX.LT.YMAPMN) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: YMAPMX.LT.YMAPMN'
            WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WBMAP4: YMAPMX.LT.YMAPMN'
            WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
            WRITE(6,*)
            STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

        IF (ZMAPMX.LT.ZMAPMN) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: ZMAPMX.LT.ZMAPMN'
            WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN WBMAP4: ZMAPMX.LT.ZMAPMN'
            WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
            WRITE(6,*)
            STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF

      IF (NMAPX.EQ.-9999) NMAPX=NINT((XSTOP-XSTART)*MYINUM)+1
      IF (NMAPX.GT.1)   DX=(XMAPMX-XMAPMN)/(NMAPX-1)*SCALXYZ
      IF (NMAPY.GT.1) DY=(YMAPMX-YMAPMN)/(NMAPY-1)*SCALXYZ
      IF (NMAPZ.GT.1) DZ=(ZMAPMX-ZMAPMN)/(NMAPZ-1)*SCALXYZ

      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='NEW'
     & ,FORM='UNFORMATTED')

      WRITE(LUNBMAP)'MAP OF COEFFS'
      WRITE(LUNBMAP)ICODE,CODE
      WRITE(LUNBMAP)NMAPX,XMAPMN,XMAPMX
      WRITE(LUNBMAP)NMAPY,YMAPMN,YMAPMX
      WRITE(LUNBMAP)NMAPZ,ZMAPMN,ZMAPMX

C DUMMY CALL TO INITIALIZE ICINDG

      NPOI=0
        DO IX=1,NBMDATX
        DO IY=1,NBMDATY
        DO IZ=1,NBMDATZ

         NPOI=NPOI+1

         X(NPOI)=IX
         Y(NPOI)=IY
         Z(NPOI)=IZ
         BX(NPOI)=0.D0
         BY(NPOI)=0.D0
         BZ(NPOI)=0.D0
      ENDDO
      ENDDO
      ENDDO

      CALL BMPOT3D4(NPOI,X,Y,Z,BX,BY,BZ,IFAIL,COMMENT)

      WRITE(LUNBMAP)NCINDG,MORD3DG
      DO IND=1,NCINDG
          WRITE(LUNBMAP)ICINDG(1,IND),ICINDG(2,IND),ICINDG(3,IND)
      ENDDO

      DO JX=1,NMAPX-1
      DO JY=1,NMAPY-1
      DO JZ=1,NMAPZ-1

          IF (IWBMAPEXT.NE.0) THEN
             IX1=-(NBMDATX-1)/2-1
             IX2=IX1+NBMDATX-1
             LX1=IX1
             LX2=IX2
            ELSE
             IX1=-(NBMDATX-1)/2-1
             LX1=IX1
             LX2=LX1+NBMDATX-1
             IF (IX1.LT.-JX) IX1=-JX
             IX2=IX1+NBMDATX-1
             IF (IX2+JX.GT.NMAPX-1) THEN
            IX2=NMAPX-1-JX
            IX1=IX2-NBMDATX+1
             ENDIF
          ENDIF

          IF (IWBMAPEXT.NE.0) THEN
             IY1=-(NBMDATY-1)/2-1
             IY2=IY1+NBMDATY-1
             LY1=IY1
             LY2=IY2
            ELSE
             IY1=-(NBMDATY-1)/2-1
             LY1=IY1
             LY2=LY1+NBMDATY-1
             IF (IY1.LT.-JY) IY1=-JY
             IY2=IY1+NBMDATY-1
             IF (IY2+JY.GT.NMAPY-1) THEN
            IY2=NMAPY-1-JY
            IY1=IY2-NBMDATY+1
             ENDIF
          ENDIF

          IF (IWBMAPEXT.NE.0) THEN
             IZ1=-(NBMDATZ-1)/2-1
             IZ2=IZ1+NBMDATZ-1
             LZ1=IZ1
             LZ2=IZ2
            ELSE
             IZ1=-(NBMDATZ-1)/2-1
             LZ1=IZ1
             LZ2=LZ1+NBMDATZ-1
             IF (IZ1.LT.-JZ) IZ1=-JZ
             IZ2=IZ1+NBMDATZ-1
             IF (IZ2+JZ.GT.NMAPZ-1) THEN
            IZ2=NMAPZ-1-JZ
            IZ1=IZ2-NBMDATZ+1
             ENDIF
          ENDIF

          NPOI=0
          DO IX=IX1,IX2
          DO IY=IY1,IY2
          DO IZ=IZ1,IZ2

         NPOI=NPOI+1

         X(NPOI)=XMAPMN+(JX+IX)*DX/SCALXYZ
         Y(NPOI)=YMAPMN+(JY+IY)*DY/SCALXYZ
         Z(NPOI)=ZMAPMN+(JZ+IZ)*DZ/SCALXYZ

         CALL MYBFELD(X(NPOI),Y(NPOI),Z(NPOI)
     &         ,BX(NPOI),BY(NPOI),BZ(NPOI),AX,AY,AZ)

            ENDDO
          ENDDO
          ENDDO

          NPOI=0
          DO IX=LX1,LX2
          DO IY=LY1,LY2
          DO IZ=LZ1,LZ2

         NPOI=NPOI+1

         DXX(NPOI)=(IX-LX1-1)*DX
         DYY(NPOI)=(IY-LY1-1)*DY
         DZZ(NPOI)=(IZ-LZ1-1)*DZ

            ENDDO
          ENDDO
          ENDDO

          CALL BMPOT3D4(NPOI,DXX,DYY,DZZ,BX,BY,BZ,IFAIL,COMMENT)

          IF (IFAIL.NE.0) THEN
             WRITE(6,*)'*** ERROR IN WBMAP4: FIT FAILED'
             WRITE(6,*)'X,Y,Z:'
             WRITE(6,*)X,Y,Z
             WRITE(6,*)'XPOI(NPOI),YPOI(NPOI),ZPOI(NPOI)'
             WRITE(6,*)'BXPOI(NPOI),BYPOI(NPOI),BZPOI(NPOI)'
             WRITE(LUNGFO,*)'*** ERROR IN WBMAP4: FIT FAILED'
             WRITE(LUNGFO,*)'X,Y,Z:'
             WRITE(LUNGFO,*)X,Y,Z
             WRITE(LUNGFO,*)'XPOI(NPOI),YPOI(NPOI),ZPOI(NPOI)'
             WRITE(LUNGFO,*)'BXPOI(NPOI),BYPOI(NPOI),BZPOI(NPOI)'
             STOP
          ENDIF

          WRITE(LUNBMAP)(CINDG4(IND),IND=1,NCINDG)

      ENDDO
      ENDDO
      ENDDO

      WRITE(LUNBMAP)SCALXYZ

      CLOSE (LUNBMAP)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR WBMAP4: FIELD MAP WRITTEN TO FILE'
      WRITE(LUNGFO,*)'     ',FILEBMAP
      WRITE(LUNGFO,*)
      WRITE(6,*)
      WRITE(6,*)'     SR WBMAP4: FIELD MAP WRITTEN TO FILE'
      WRITE(6,*)'     ',FILEBMAP
      WRITE(6,*)

      DEALLOCATE (X)
      DEALLOCATE (Y)
      DEALLOCATE (Z)
      DEALLOCATE (BX)
      DEALLOCATE (BY)
      DEALLOCATE (BZ)

      RETURN
      END
