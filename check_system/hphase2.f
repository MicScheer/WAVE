*CMZ :  3.00/00 11/03/2013  15.13.36  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/00 11/12/2012  12.00.31  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  09.13.47  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.43  by  Michael Scheer
*CMZ :  2.49/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.16/00 08/06/2000  12.11.00  by  Michael Scheer
*CMZ :  1.03/06 23/09/98  17.08.00  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.01/06 13/02/95  13.14.40  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  10.41.22  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.53.40  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.27  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.49  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HPHASE2(ID,TIT,FILL)
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

C--- STORE RESULTS OF SPECTRUM CALCULATION ON HISTOGRAM FILE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      INTEGER IOB,IOBZ,IOBY,ID
      INTEGER ICYCLE,ICYCLE1

      REAL*4 ZMIN,ZMAX,YMIN,YMAX,ZFILL,YFILL,DZ,DY
      REAL*8 FILL(MPHASEZ*MPHASEY)

      CHARACTER(80) TIT,TIT1

      IF (mphaseZ.GT.1.AND.mphaseY.GT.1) THEN

        DZ=PHWID/(nphaseZ-1)
        DY=PHHIG/(nphaseY-1)
        ZMIN=PHCENZ-(mphasez-1)*dz/2.0d0
        ZMax=PHCENZ+(mphasez-1)*dz/2.0d0
        yMIN=PHCENy-(mphasey-1)*dy/2.0d0
        yMax=PHCENy+(mphasey-1)*dy/2.0d0

        call hbook2m(ID,TIT,mphaseZ,ZMIN,ZMAX,mphaseY,YMIN,YMAX,VMX)
        TIT1=TIT(1:40)//'(HORIZONTAL CUT)'
        call hbook1m(ID-1,TIT1,mphaseZ,ZMIN,ZMAX,VMX)
        TIT1=TIT(1:40)//'(VERTICAL CUT)'
        call hbook1m(ID-2,TIT1,mphaseY,YMIN,YMAX,VMX)

        IOB=0
        DO IOBY=1,mphaseY
          DO IOBZ=1,mphaseZ

            ZFILL=PHCENZ-PHWID/2.+(IOBZ-1)*DZ
            YFILL=PHCENY-PHHIG/2.+(IOBY-1)*DY
            IOB=IOB+1

            CALL hfillm(ID,ZFILL,YFILL,FILL(IOB))
            IF (IOBZ.EQ.mphaseZ/2+1) CALL hfillm(ID-2,YFILL,0.,FILL(IOB))
            IF (IOBY.EQ.mphaseY/2+1) CALL hfillm(ID-1,ZFILL,0.,FILL(IOB))

          ENDDO !IOBZ
        ENDDO !IOBY

        ICYCLE1=ICYCLE
        CALL MHROUT(ID,ICYCLE,' ')
        CALL hdeletm(ID)
        CALL MHROUT(ID-1,ICYCLE1,' ')
        ICYCLE1=ICYCLE
        CALL hdeletm(ID-1)
        CALL MHROUT(ID-2,ICYCLE1,' ')
        CALL hdeletm(ID-2)

      ELSE IF (mphaseZ.GT.1) THEN

        DZ=PHWID/(mphaseZ-1)
        ZMIN=PHCENZ-PHWID/2.-DZ/2.
        ZMAX=PHCENZ+PHWID/2.+DZ/2.

        TIT1=TIT(1:40)//'(HORIZONTAL CUT)'
        call hbook1m(ID-1,TIT1,mphaseZ,ZMIN,ZMAX,VMX)

        DO IOBZ=1,mphaseZ
          ZFILL=PHCENZ-PHWID/2.+(IOBZ-1)*DZ
          CALL hfillm(ID-1,ZFILL,0.,FILL(IOB))
        ENDDO !IOBZ

        CALL MHROUT(ID-1,ICYCLE,' ')
        CALL hdeletm(ID-1)

      ELSE IF (mphaseY.GT.1) THEN

        DY=PHHIG/(mphaseY-1)
        YMIN=PHCENY-PHHIG/2.-DY/2.
        YMAX=PHCENY+PHHIG/2.+DY/2.

        TIT1=TIT(1:40)//'(VERTICAL CUT)'
        call hbook1m(ID-2,TIT1,mphaseY,YMIN,YMAX,VMX)

        DO IOBY=1,mphaseY
            YFILL=PHCENY-PHHIG/2.+(IOBY-1)*DY
            CALL hfillm(ID-2,YFILL,0.,FILL(IOB))
        ENDDO !IOBY

        CALL MHROUT(ID-2,ICYCLE,' ')
        CALL hdeletm(ID-2)

      ENDIF !mphaseZ, mphaseY

      RETURN
      END
