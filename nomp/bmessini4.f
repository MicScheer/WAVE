*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.45/01 13/12/2002  17.24.57  by  Michael Scheer
*CMZ :  2.44/02 12/12/2002  14.09.47  by  Michael Scheer
*CMZ :  2.44/01 06/12/2002  14.56.41  by  Michael Scheer
*CMZ :  2.44/00 08/11/2002  12.15.59  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.34/08 17/09/2001  19.11.02  by  Michael Scheer
*CMZ :  2.34/05 23/08/2001  17.35.07  by  Michael Scheer
*CMZ :  2.34/01 02/07/2001  17.22.17  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.25.47  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.02.52  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  18.22.03  by  Michael Scheer
*CMZ :  1.00/00 29/07/97  16.33.37  by  Michael Scheer
*CMZ : 00.02/05 17/04/97  17.59.48  by  Michael Scheer
*CMZ : 00.02/04 12/02/97  14.52.46  by  Michael Scheer
*CMZ : 00.02/00 28/11/96  17.04.35  by  Michael Scheer
*-- Author :    Michael Scheer   11/11/96

      SUBROUTINE BMESSINI4
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

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,bmessf90.
      include 'bmessf90.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.


      INTEGER IBMESSX,IBMESSY,IBMESSZ
      INTEGER ICODEBM,IND,ICAL

      REAL*4 X,Y,Z

        CHARACTER(13) COMMENT4
        CHARACTER(65) COMMENT

      DATA ICAL/0/

      IF (ICAL.NE.0) RETURN

      ICAL=1

      BMESSEPS=0.00001

      IF (NTUPGRID.NE.0) THEN
          WRITE(6,*)'*** ERROR IN BMESSINI4: NTUPGRID.NE.0'
          WRITE(LUNGFO,*)'*** ERROR IN BMESSINI4: NTUPGRID.NE.0'
          STOP
      ENDIF

      OPEN(UNIT=LUNB0,FILE=FILEB0,STATUS='OLD'
     &    ,FORM='UNFORMATTED')
        REWIND(LUNB0)

        READ(LUNB0) COMMENT4
      IF (COMMENT4.NE.'MAP OF COEFFS') THEN
          WRITE(6,*)'*** ERROR IN BMESSINI4: NOT A COEFFICIENT MAP ON FILE'
          WRITE(6,*)FILEBMAP
          WRITE(LUNGFO,*)'*** ERROR IN BMESSINI4: NOT A COEFFICIENT MAP ON FILE'
          WRITE(LUNGFO,*)FILEBMAP
          STOP
      ENDIF

      READ(LUNB0)ICODEBM,COMMENT
      READ(LUNB0)NBMESSX,BMXMIN,BMXMAX
      READ(LUNB0)NBMESSY,BMYMIN,BMYMAX
      READ(LUNB0)NBMESSZ,BMZMIN,BMZMAX

      READ(LUNB0)NCINDG,MORD3DG
      IF (NCINDG.GT.NDIMCIPG) THEN
          WRITE(6,*)'*** ERROR IN BMESSINI4: DIMENSION NDIMCIPG EXCEEDED'
          WRITE(LUNGFO,*)'*** ERROR IN BMESSINI4: DIMENSION NDIMCIPG EXCEEDED'
          STOP
      ENDIF
      DO IND=1,NCINDG
        READ(LUNB0) ICINDG(1,IND),ICINDG(2,IND),ICINDG(3,IND)
      ENDDO

      ALLOCATE(BMCOEF(NCINDG,NBMESSZ,NBMESSY,NBMESSX))

          DO IBMESSX=1,NBMESSX-1
          DO IBMESSY=1,NBMESSY-1
          DO IBMESSZ=1,NBMESSZ-1
          READ(LUNB0)(BMCOEF(IND,IBMESSZ,IBMESSY,IBMESSX),IND=1,NCINDG)
          ENDDO
          ENDDO
          ENDDO

          READ(LUNB0)B3DSCALEG

      CLOSE(LUNB0)

          BMESSDX=(BMXMAX-BMXMIN)/(NBMESSX-1)
          BMESSDY=(BMYMAX-BMYMIN)/(NBMESSY-1)
          BMESSDZ=(BMZMAX-BMZMIN)/(NBMESSZ-1)

      X=BMXMAX
          IBMESSX=NINT((X-BMXMIN)/BMESSDX)
          IF (ABS(X-(BMXMIN+IBMESSX*BMESSDX)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI4:'
            WRITE(6,*)'BAD X-VALUE:',X
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI4:'
            WRITE(LUNGFO,*)'BAD X-VALUE:',X
            STOP
          ENDIF
      IF (NBMESSX.EQ.9999) THEN
          NBMESSX=IBMESSX+1
      ENDIF

      Y=BMYMAX
          IBMESSY=NINT((Y-BMYMIN)/BMESSDY)
          IF (ABS(Y-(BMYMIN+IBMESSY*BMESSDY)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI4:'
            WRITE(6,*)'BAD Y-VALUE:',Y
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI4:'
            WRITE(LUNGFO,*)'BAD Y-VALUE:',Y
            STOP
          ENDIF
      IF (NBMESSY.EQ.9999) THEN
          NBMESSY=IBMESSY+1
      ENDIF

      Z=BMZMAX
          IBMESSZ=NINT((Z-BMZMIN)/BMESSDZ)
          IF (ABS(Z-(BMZMIN+IBMESSZ*BMESSDZ)).GT.BMESSEPS) THEN
            WRITE(6,*)'*** ERROR IN BMESSINI4:'
            WRITE(6,*)'BAD Z-VALUE:',Z
            WRITE(LUNGFO,*)'*** ERROR IN BMESSINI4:'
            WRITE(LUNGFO,*)'BAD Z-VALUE:',Z
            STOP
          ENDIF
      IF (NBMESSZ.EQ.9999) THEN
          NBMESSZ=IBMESSZ+1
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE BMESSINI4:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     DATA FILE:'
      WRITE(LUNGFO,*)'     ',FILEB0
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     NUMBER OF COEFFICIENTS READ:',NCINDG
      WRITE(LUNGFO,*)'     XMIN, XMAX:  ',BMXMIN,BMXMAX
      WRITE(LUNGFO,*)'     YMIN, YMAX:  ',BMYMIN,BMYMAX
      WRITE(LUNGFO,*)'     ZMIN, ZMAX:  ',BMZMIN,BMZMAX
      WRITE(LUNGFO,*)'     DX, DY, DZ:  '
      WRITE(LUNGFO,*)'     ',BMESSDX,BMESSDY,BMESSDZ
      WRITE(LUNGFO,*)'     EPSILON:     ',BMESSEPS

      IF (XSTART.EQ.9999.) THEN
          IF (IPERIODG.NE.0.AND.SIGNG.LT.0.D0) THEN
         XSTART=2.D0*BMXMIN
          ELSE
         XSTART=BMXMIN
          ENDIF
      ENDIF

      IF (XSTOP.EQ.9999.) THEN
          IF (IPERIODG.NE.0.AND.SIGNG.LT.0.D0) THEN
         XSTOP=2.D0*BMXMAX
          ELSE
         XSTOP=BMXMAX
          ENDIF
      ENDIF

      RETURN
      END
