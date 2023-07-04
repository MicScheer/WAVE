*CMZ :  4.00/04 17/05/2019  14.17.20  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.19.07  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  08.28.26  by  Michael Scheer
*CMZ :  2.54/07 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.47/03 12/03/2003  16.01.17  by  Michael Scheer
*CMZ :  2.45/03 16/12/2002  17.52.31  by  Michael Scheer
*CMZ :  2.44/01 10/12/2002  17.49.10  by  Michael Scheer
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
      SUBROUTINE WBMAP5
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

*KEEP,bmessf90u.
      include 'bmessf90u.cmn'
*KEND.

C--- CONVERT COLUMN FORMAT OF FILEB0 TO BMESS-FORMAT FOR BMESS


      IMPLICIT NONE

      EXTERNAL DCOSD,DSIND
      DOUBLE PRECISION DCOSD,DSIND

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bmap.
      include 'bmap.cmn'
*KEEP,bmessf90.
      include 'bmessf90.cmn'
*KEND.

      DOUBLE PRECISION DX,DY,DZ,BX,BY,BZ,X,Y,Z,EPSXYZ,DXX,DYY,DZZ

      INTEGER I,IX,IY,IZ

      DATA DX,DY,DZ/0.D0,0.D0,0.D0/
      DATA EPSXYZ/1.D-6/

      IF (IBMRADIAL.NE.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: IBMRADIAL.NE.0'
          WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WBMAP5: IBMRADIAL.NE.0'
          WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(6,*)
          STOP
      ENDIF !IBMRADIAL

      IF (NMAPX.LT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: NMAPX.LT.2'
          WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WBMAP5: NMAPX.LT.2'
          WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(6,*)
          STOP
      ENDIF !NMAPX

      IF (NMAPY.LT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: NMAPY.LT.2'
          WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WBMAP5: NMAPY.LT.2'
          WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(6,*)
          STOP
      ENDIF !NMAPY

      IF (NMAPZ.LT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: NMAPZ.LT.2'
          WRITE(LUNGFO,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WBMAP5: NMAPZ.LT.2'
          WRITE(6,*)'CHECK NAMELIST BMAP IN WAVE.IN'
          WRITE(6,*)
          STOP
      ENDIF !NMAPZ

      NBMESSX=NMAPX
      NBMESSY=NMAPY
      NBMESSZ=NMAPZ

      DX=(XMAPMX-XMAPMN)/(NMAPX-1)
      DY=(YMAPMX-YMAPMN)/(NMAPY-1)
      DZ=(ZMAPMX-ZMAPMN)/(NMAPZ-1)

        ALLOCATE(BDATA(3,NBMESSZ,NBMESSY,NBMESSX))

      DO IX=1,NBMESSX
        DO IY=1,NBMESSY
        DO IZ=1,NBMESSZ
        DO I=1,3
            BDATA(I,IZ,IY,IX)=-9999.
        ENDDO
        ENDDO
        ENDDO
        ENDDO

      OPEN(UNIT=LUNB0,FILE=FILEB0,STATUS='OLD',FORM='FORMATTED')

      I=0
1     READ(LUNB0,*,END=9)X,Y,Z,BX,BY,BZ

         I=I+1

         IX=NINT((X-XMAPMN)/DX)+1
         IY=NINT((Y-YMAPMN)/DY)+1
         IZ=NINT((Z-ZMAPMN)/DZ)+1

         DXX=X-(XMAPMN+(IX-1)*DX)
         DYY=Y-(YMAPMN+(IY-1)*DY)
         DZZ=Z-(ZMAPMN+(IZ-1)*DZ)

         IF (ABS(DXX).LT.EPSXYZ) THEN
             X=XMAPMN+(IX-1)*DX
         ENDIF

         IF (ABS(DYY).LT.EPSXYZ) THEN
             Y=YMAPMN+(IY-1)*DY
         ENDIF

         IF (ABS(DZZ).LT.EPSXYZ) THEN
             Z=ZMAPMN+(IZ-1)*DZ
         ENDIF

         IF (IX.LT.1.OR.IX.GT.NMAPX
     &         .OR.X.LT.XMAPMN-EPSXYZ.OR.X.GT.XMAPMX+EPSXYZ) THEN
             WRITE(6,*)'*** ERROR IN WBMAP5: BAD X'
             WRITE(6,*)'X: ',X
             WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: BAD X'
             WRITE(LUNGFO,*)'IX,X: ',IX,X
             STOP
         ENDIF

         IF (IY.LT.1.OR.IY.GT.NMAPY
     &         .OR.Y.LT.YMAPMN-EPSXYZ.OR.Y.GT.YMAPMX+EPSXYZ) THEN
             WRITE(6,*)'*** ERROR IN WBMAP5: BAD Y'
             WRITE(6,*)'Y: ',Y
             WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: BAD Y'
             WRITE(LUNGFO,*)'IY,Y: ',IY,Y
             STOP
         ENDIF

         IF (IZ.LT.1.OR.IZ.GT.NMAPZ
     &         .OR.Z.LT.ZMAPMN-EPSXYZ.OR.Z.GT.ZMAPMX+EPSXYZ) THEN
             WRITE(6,*)'*** ERROR IN WBMAP5: BAD Z'
             WRITE(6,*)'Z: ',Z
             WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: BAD Z'
             WRITE(LUNGFO,*)'IZ,Z: ',IZ,Z
             STOP
         ENDIF

         BDATA(1,IZ,IY,IX)=BX
         BDATA(2,IZ,IY,IX)=BY
         BDATA(3,IZ,IY,IX)=BZ

      GOTO 1

9     CLOSE (LUNBMAP)

      IF (I.NE.NMAPX*NMAPY*NMAPZ) THEN
             WRITE(6,*)'*** ERROR IN WBMAP5: BAD NUMBER OF DATA'
             WRITE(6,*)'EXPECTED: ',NMAPX*NMAPY*NMAPZ
             WRITE(6,*)'FOUND: ',I
             WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: BAD NUMBER OF DATA'
             WRITE(LUNGFO,*)'EXPECTED: ',NMAPX*NMAPY*NMAPZ
             WRITE(LUNGFO,*)'FOUND: ',I
             STOP
          ENDIF

      OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='NEW'
     &      ,FORM='UNFORMATTED')

                WRITE(LUNBMAP)ICODE,CODE
                WRITE(LUNBMAP)NMAPX,XMAPMN,XMAPMX
                WRITE(LUNBMAP)NMAPY,YMAPMN,YMAPMX
                WRITE(LUNBMAP)NMAPZ,ZMAPMN,ZMAPMX

      DO IX=1,NBMESSX
        DO IY=1,NBMESSY
        DO IZ=1,NBMESSZ
        DO I=1,3
            IF (BDATA(I,IZ,IY,IX).EQ.-9999.) THEN
             WRITE(6,*)'*** ERROR IN WBMAP5: FIELD DATA MISSING'
             WRITE(6,*)'IX,IY,IZ:'
             WRITE(6,*)IX,IY,IZ
             WRITE(6,*)'DX,DY,DZ:'
             WRITE(6,*)DX
             WRITE(6,*)DY
             WRITE(6,*)DZ
             WRITE(6,*)'X,Y,Z:'
             WRITE(6,*)XMAPMN+(IX-1)*DX
             WRITE(6,*)YMAPMN+(IY-1)*DY
             WRITE(6,*)ZMAPMN+(IZ-1)*DZ
             WRITE(LUNGFO,*)'*** ERROR IN WBMAP5: FIELD DATA MISSING'
             WRITE(LUNGFO,*)'X,Y,Z:'
             WRITE(LUNGFO,*)XMAPMN+(IX-1)*DX
             WRITE(LUNGFO,*)YMAPMN+(IY-1)*DY
             WRITE(LUNGFO,*)ZMAPMN+(IZ-1)*DZ
             STOP
          ENDIF
        ENDDO
        WRITE(LUNBMAP)(BDATA(I,IZ,IY,IX),I=1,3)
        ENDDO
        ENDDO
        ENDDO

      CLOSE (LUNB0)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'       WBMAP5: FIELD MAP CONVERTED'
      WRITE(LUNGFO,*)'     INPUT FILE: ',FILEB0
      WRITE(LUNGFO,*)'     OUTPUT FILE: ',FILEBMAP
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     NX, XMIN, XMAX:  ',NMAPX,SNGL(XMAPMN),SNGL(XMAPMX)
      WRITE(LUNGFO,*)'     NY, YMIN, YMAX:  ',NMAPY,SNGL(YMAPMN),SNGL(YMAPMX)
      WRITE(LUNGFO,*)'     NZ, ZMIN, ZMAX:  ',NMAPZ,SNGL(ZMAPMN),SNGL(ZMAPMX)
      WRITE(LUNGFO,*)

      WRITE(6,*)
     &'       WBMAP5: FIELD MAP CONVERTED'
      WRITE(6,*)'     INPUT FILE: ',FILEB0
      WRITE(6,*)'     OUTPUT FILE: ',FILEBMAP
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'     NX, XMIN, XMAX:  ',NMAPX,SNGL(XMAPMN),SNGL(XMAPMX)
      WRITE(6,*)'     NY, YMIN, YMAX:  ',NMAPY,SNGL(YMAPMN),SNGL(YMAPMX)
      WRITE(6,*)'     NZ, ZMIN, ZMAX:  ',NMAPZ,SNGL(ZMAPMN),SNGL(ZMAPMX)
      WRITE(6,*)


      RETURN
      END
