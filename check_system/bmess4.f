*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.46/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.44/02 12/12/2002  14.09.23  by  Michael Scheer
*CMZ :  2.44/01 06/12/2002  11.10.22  by  Michael Scheer
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
      SUBROUTINE BMESS4(X,Y,Z,BX,BY,BZ)

C     INTERPOLATION OF 3D EQUALLY SPACED COEFFICIENT MAP
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
*KEND.

*KEEP,bmessf90u.
      include 'bmessf90u.cmn'
*KEND.

      IMPLICIT NONE

      INTEGER ICAL,IX1,IY1,IZ1
      INTEGER IX1O,IY1O,IZ1O

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,STEP
      DOUBLE PRECISION DX,DY,DZ

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

      IF (ICAL.EQ.0) THEN

            STEP=1.D0/MYINUM

          CALL BMESSINI4

          X3DMING=BMXMIN*B3DSCALEG
          X3DMAXG=BMXMAX*B3DSCALEG
          Y3DMING=BMYMIN*B3DSCALEG
          Y3DMAXG=BMYMAX*B3DSCALEG
          Z3DMING=BMZMIN*B3DSCALEG
          Z3DMAXG=BMZMAX*B3DSCALEG

          ICAL=1

      ENDIF

      IF (X.LT.BMXMIN-STEP.OR.X.GT.BMXMAX+STEP) THEN
          WRITE(6,*)'*** ERROR IN BMESS4: X OUT OF RANGE'
          WRITE(6,*)'X:',X
          WRITE(LUNGFO,*)'*** ERROR IN BMESS4: X OUT OF RANGE'
          WRITE(LUNGFO,*)'X:',X
          STOP
      ENDIF

      IF (Y.LT.BMYMIN-STEP.OR.Y.GT.BMYMAX+STEP) THEN
          WRITE(6,*)'*** ERROR IN BMESS4: Y OUT OF RANGE'
          WRITE(6,*)'Y:',Y
          WRITE(LUNGFO,*)'*** ERROR IN BMESS4: Y OUT OF RANGE'
          WRITE(LUNGFO,*)'Y:',Y
          STOP
      ENDIF

      IF (Z.LT.BMZMIN-STEP.OR.Z.GT.BMZMAX+STEP) THEN
          WRITE(6,*)'*** ERROR IN BMESS4: Z OUT OF RANGE'
          WRITE(6,*)'Z:',Z
          WRITE(LUNGFO,*)'*** ERROR IN BMESS4: Z OUT OF RANGE'
          WRITE(LUNGFO,*)'Z:',Z
          STOP
      ENDIF

      IX1=(X-BMXMIN)/BMESSDX+1
      IX1=MAX(1,IX1)
      IX1=MIN(NBMESSX-1,IX1)

      IY1=(Y-BMYMIN)/BMESSDY+1
      IY1=MAX(1,IY1)
      IY1=MIN(NBMESSY-1,IY1)

      IZ1=(Z-BMZMIN)/BMESSDZ+1
      IZ1=MAX(1,IZ1)
      IZ1=MIN(NBMESSZ-1,IZ1)

500   CONTINUE

        DX=(X-(BMXMIN+(IX1-1)*BMESSDX))*B3DSCALEG
        DY=(Y-(BMYMIN+(IY1-1)*BMESSDY))*B3DSCALEG
        DZ=(Z-(BMZMIN+(IZ1-1)*BMESSDZ))*B3DSCALEG

      CALL BMESS3D4(DX,DY,DZ,BX,BY,BZ,IX1,IY1,IZ1)

      IX1O=IX1
      IY1O=IY1
      IZ1O=IZ1

      RETURN
      END
