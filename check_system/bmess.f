*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/02 02/07/2012  11.43.27  by  Michael Scheer
*CMZ :  2.57/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.54/05 17/05/2005  15.01.37  by  Michael Scheer
*CMZ :  2.45/01 13/12/2002  14.40.32  by  Michael Scheer
*CMZ :  2.44/01 10/12/2002  17.52.22  by  Michael Scheer
*CMZ :  2.44/00 07/11/2002  17.13.48  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.34/09 20/09/2001  12.05.02  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.13.01  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.58.52  by  Michael Scheer
*CMZ :  1.00/00 06/08/97  17.14.18  by  Michael Scheer
*CMZ : 00.02/05 16/04/97  17.05.38  by  Michael Scheer
*CMZ : 00.02/00 22/11/96  15.24.30  by  Michael Scheer
*CMZ : 00.01/09 05/10/95  13.50.07  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.31.15  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.54.08  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.33.07  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.00  by  Michael Scheer
*-- Author : Michael Scheer
C**********************************************************************
      SUBROUTINE BMESS(X,Y,Z,BX,BY,BZ)

C     INTERPOLATION OF 3D EQUALLY SPACED FIELD MAP
C**********************************************************************

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
*KEEP,bmessf90u.
      include 'bmessf90u.cmn'
*KEND.

      IMPLICIT NONE

      CHARACTER(64) COMMENT

      INTEGER ICAL,IX,IX1,IX2,IY,IY1,IY2,IZ,IZ1,IZ2,NPOIP,NPOI,IFAIL
      INTEGER IX1O,IY1O,IZ1O,I,NX,NY,NZ,IWARNX
      PARAMETER (NPOIP=4*4*4)

      DOUBLE PRECISION BXO,BYO,BZO
      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,STEP,DXX,DYY,DZZ
      DOUBLE PRECISION DX,DY,DZ,SCALXYZ
      DOUBLE PRECISION XPOI(NPOIP),YPOI(NPOIP),ZPOI(NPOIP)
      DOUBLE PRECISION BXPOI(NPOIP),BYPOI(NPOIP),BZPOI(NPOIP)
      DOUBLE PRECISION B111(3),B211(3), B121(3),B221(3)
      DOUBLE PRECISION B112(3),B212(3), B122(3),B222(3)
      DOUBLE PRECISION B112111(3),B212211(3), B122121(3),B222221(3)
      DOUBLE PRECISION BLOW(3),BHIG(3),B(3)

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,bmessf90.
      include 'bmessf90.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEND.

      DATA ICAL/0/
      DATA SCALXYZ/100.D0/
      DATA IWARNX/0/

      IF (IRFILB0.EQ.-4) THEN
        CALL BMESS4(X,Y,Z,BX,BY,BZ)
        RETURN
      else IF (IRFILB0.EQ.6.or.irfilb0.eq.-6) THEN
        CALL bmap(X,Y,Z,BX,BY,BZ)
        RETURN
      ENDIF

      IF (ICAL.EQ.0) THEN

        STEP=1.D0/MYINUM
        LORD3DG=1
        NDORD3DG=1

        CALL BMESSINI

        B3DSCALEG=SCALXYZ

        X3DMING=BMXMIN*SCALXYZ
        X3DMAXG=BMXMAX*SCALXYZ
        Y3DMING=BMYMIN*SCALXYZ
        Y3DMAXG=BMYMAX*SCALXYZ
        Z3DMING=BMZMIN*SCALXYZ
        Z3DMAXG=BMZMAX*SCALXYZ

        NX=NBMDATX
        NY=NBMDATY
        NZ=NBMDATZ

        IF (IRFILB0.EQ.-1.OR.IRFILB0.EQ.-2) THEN
          NX=2
          NY=2
          NZ=2
        ENDIF

        ICAL=1

      ENDIF

      IF (IWARNX.EQ.0.AND.(X.LT.BMXMIN-STEP.OR.X.GT.BMXMAX+STEP)) THEN
        WRITE(6,*)'*** WARNING: IN BMESS: X OUT OF RANGE'
        WRITE(6,*)'X:',X
        WRITE(6,*)'Y:',Y
        WRITE(6,*)'Z:',Z
        WRITE(LUNGFO,*)'*** WARNING IN BMESS: X OUT OF RANGE'
        WRITE(LUNGFO,*)'X:',X
        WRITE(LUNGFO,*)'Y:',Y
        WRITE(LUNGFO,*)'Z:',Z
        IWARNX=1
      ENDIF

      IF (IWARNX.NE.0.AND.(X.LT.BMXMIN-STEP.OR.X.GT.BMXMAX+STEP)) THEN
        BX=0.0D0
        BY=0.0D0
        BZ=0.0D0
        RETURN
      ENDIF

      IF (Y.LT.BMYMIN-STEP.OR.Y.GT.BMYMAX+STEP) THEN
        WRITE(6,*)'*** ERROR IN BMESS: Y OUT OF RANGE'
        WRITE(6,*)'X:',X
        WRITE(6,*)'Y:',Y
        WRITE(6,*)'Z:',Z
        WRITE(LUNGFO,*)'*** ERROR IN BMESS: Y OUT OF RANGE'
        WRITE(LUNGFO,*)'X:',X
        WRITE(LUNGFO,*)'Y:',Y
        WRITE(LUNGFO,*)'Z:',Z
        STOP
      ENDIF

      IF (Z.LT.BMZMIN-STEP.OR.Z.GT.BMZMAX+STEP) THEN
        WRITE(6,*)'*** ERROR IN BMESS: Z OUT OF RANGE'
        WRITE(6,*)'X:',X
        WRITE(6,*)'Y:',Y
        WRITE(6,*)'Z:',Z
        WRITE(LUNGFO,*)'*** ERROR IN BMESS: Z OUT OF RANGE'
        WRITE(LUNGFO,*)'X:',X
        WRITE(LUNGFO,*)'Y:',Y
        WRITE(LUNGFO,*)'Z:',Z
        STOP
      ENDIF

      IF (X.LT.BMXMIN) THEN
        IX1=1
      ELSE
C        IX1=NINT((X-BMXMIN+BMESSEPS)/BMESSDX)+1
        IX1=(X-BMXMIN)/BMESSDX+1-(NX-1)/2
        IX1=MAX(1,IX1)
      ENDIF
      IX2=IX1+NX-1
      IF (IX2.GT.NBMESSX) THEN
        IX2=NBMESSX
        IX1=IX2-NX+1
      ENDIF

      IF (Y.LT.BMYMIN) THEN
        IY1=1
      ELSE
C        IY1=NINT((Y-BMYMIN+BMESSEPS)/BMESSDY)+1
        IY1=(Y-BMYMIN)/BMESSDY+1-(NY-1)/2
        IY1=MAX(1,IY1)
      ENDIF
      IY2=IY1+NY-1
      IF (IY2.GT.NBMESSY) THEN
        IY2=NBMESSY
        IY1=IY2-NY+1
      ENDIF

      IF (Z.LT.BMZMIN) THEN
        IZ1=1
      ELSE
C        IZ1=NINT((Z-BMZMIN+BMESSEPS)/BMESSDZ)+1
        IZ1=(Z-BMZMIN)/BMESSDZ+1-(NZ-1)/2
        IZ1=MAX(1,IZ1)
      ENDIF
      IZ2=IZ1+NZ-1
      IF (IZ2.GT.NBMESSZ) THEN
        IZ2=NBMESSZ
        IZ1=IZ2-NZ+1
      ENDIF

      IF (IRFILB0.EQ.-1) THEN

        IF (IX1.EQ.IX1O.AND.IY1.EQ.IY1O.AND.IZ1.EQ.IZ1O) THEN
          BX=BXO
          BY=BYO
          BZ=BZO
        ELSE
          BX=BDATA(1,IZ1,IY1,IX1)
          BY=BDATA(2,IZ1,IY1,IX1)
          BZ=BDATA(3,IZ1,IY1,IX1)
          BXO=BX
          BYO=BY
          BZO=BZ
        ENDIF   !OLD VOXL

      ELSE IF (IRFILB0.EQ.-2) THEN

        DXX=(X-(BMXMIN+(IX1-1)*BMESSDX))/BMESSDX
        DYY=(Y-(BMYMIN+(IY1-1)*BMESSDY))/BMESSDY
        DZZ=(Z-(BMZMIN+(IZ1-1)*BMESSDZ))/BMESSDZ

        DO I=1,3

          B111(I)=BDATA(I,IZ1,IY1,IX1)
          B211(I)=BDATA(I,IZ2,IY1,IX1)
          B121(I)=BDATA(I,IZ1,IY2,IX1)
          B221(I)=BDATA(I,IZ2,IY2,IX1)
          B112(I)=BDATA(I,IZ1,IY1,IX2)
          B212(I)=BDATA(I,IZ2,IY1,IX2)
          B122(I)=BDATA(I,IZ1,IY2,IX2)
          B222(I)=BDATA(I,IZ2,IY2,IX2)

          B112111(I)=B111(I)+(B112(I)-B111(I))*DXX
          B122121(I)=B121(I)+(B122(I)-B121(I))*DXX
          B212211(I)=B211(I)+(B212(I)-B211(I))*DXX
          B222221(I)=B221(I)+(B222(I)-B221(I))*DXX

          BLOW(I)=B112111(I)+(B212211(I)-B112111(I))*DZZ
          BHIG(I)=B122121(I)+(B222221(I)-B122121(I))*DZZ

          B(I)=BLOW(I)+(BHIG(I)-BLOW(I))*DYY

        ENDDO

        BX=B(1)
        BY=B(2)
        BZ=B(3)

      ELSE !(IRFILB0.EQ.-1)

        DX=(X-(BMXMIN+(IX1-1)*BMESSDX))*SCALXYZ
        DY=(Y-(BMYMIN+(IY1-1)*BMESSDY))*SCALXYZ
        DZ=(Z-(BMZMIN+(IZ1-1)*BMESSDZ))*SCALXYZ

C SKIP FIT, IF STILL IN OLD VOXEL
        IF (IX1.EQ.IX1O.AND.IY1.EQ.IY1O.AND.IZ1.EQ.IZ1O) GOTO 500

        NPOI=0
        DO IX=IX1,IX2
          DO IY=IY1,IY2
            DO IZ=IZ1,IZ2

              NPOI=NPOI+1
              IF (NPOI.GT.NPOIP) THEN
                WRITE(6,*)'*** ERROR IN BMESS: DIMENSION EXCEEDED ***'
                STOP
              ENDIF


              XPOI(NPOI)=BMESSDX*(IX-IX1)*SCALXYZ
              YPOI(NPOI)=BMESSDY*(IY-IY1)*SCALXYZ
              ZPOI(NPOI)=BMESSDZ*(IZ-IZ1)*SCALXYZ

C ATTENTION: IZ,IY,IX !! BE AWARE
              BXPOI(NPOI)=BDATA(1,IZ,IY,IX)
              BYPOI(NPOI)=BDATA(2,IZ,IY,IX)
              BZPOI(NPOI)=BDATA(3,IZ,IY,IX)

            ENDDO
          ENDDO
        ENDDO

        CALL BMPOT3D(NPOI,XPOI,YPOI,ZPOI,BXPOI,BYPOI,BZPOI,IFAIL,COMMENT)

        IF (IFAIL.NE.0) THEN
          WRITE(6,*)'*** ERROR IN BMESS: FIT FAILED'
          WRITE(6,*)'X,Y,Z:'
          WRITE(6,*)X,Y,Z
          WRITE(6,*)'XPOI(NPOI),YPOI(NPOI),ZPOI(NPOI)'
          WRITE(6,*)'BXPOI(NPOI),BYPOI(NPOI),BZPOI(NPOI)'
          WRITE(LUNGFO,*)'*** ERROR IN BMESS: FIT FAILED'
          WRITE(LUNGFO,*)'X,Y,Z:'
          WRITE(LUNGFO,*)X,Y,Z
          WRITE(LUNGFO,*)'XPOI(NPOI),YPOI(NPOI),ZPOI(NPOI)'
          WRITE(LUNGFO,*)'BXPOI(NPOI),BYPOI(NPOI),BZPOI(NPOI)'
          NPOI=0
          DO IX=IX1,IX2
            DO IY=IY1,IY2
              DO IZ=IZ1,IZ2
                NPOI=NPOI+1
                WRITE(6,*)XPOI(NPOI),YPOI(NPOI),ZPOI(NPOI)
                WRITE(6,*)BXPOI(NPOI),BYPOI(NPOI),BZPOI(NPOI)
                WRITE(LUNGFO,*)XPOI(NPOI),YPOI(NPOI),ZPOI(NPOI)
                WRITE(LUNGFO,*)BXPOI(NPOI),BYPOI(NPOI),BZPOI(NPOI)
              ENDDO
            ENDDO
          ENDDO
          STOP
        ENDIF

500   CONTINUE

      CALL BMESS3D(DX,DY,DZ,BX,BY,BZ)

      ENDIF !(IRFILB0.EQ.-1)

      IX1O=IX1
      IY1O=IY1
      IZ1O=IZ1

      RETURN
      END
