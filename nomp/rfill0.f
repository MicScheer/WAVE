*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.52/13 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.38/00 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.33/01 03/05/2001  12.28.19  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  12.21.53  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.35.07  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.16.26  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.41  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.38  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE RFILL0
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.

C--- SOURCE POINTS ARE READ FROM FILE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,specdip.
      include 'specdip.cmn'
*KEND.

      INTEGER ICODL0,ISOUR,I,J

C--- READ FILE FILEL0

      OPEN (UNIT=LUNL0,FILE = FILEL0,FORM = 'FORMATTED',STATUS = 'OLD')

      READ(LUNL0,*)ICODL0
      READ(LUNL0,*)NSOURCE
      READ(LUNL0,*)WGWINFC,WBL0CUT,WBL0HYS
      READ(LUNL0,*)CX1,CY1,CZ1,WID1,HIG1
      READ(LUNL0,*)CX2,CY2,CZ2,WID2,HIG2

        IF (ISPECDIP.GT.0.AND.NDIP.GT.NSOURCE) NSOURCE=NDIP
      ALLOCATE(SOURCEA(3,4,NSOURCE))
      ALLOCATE(SOURCEE(3,4,NSOURCE))
      ALLOCATE(SOURCEAO(3,4,NSOURCE))
      ALLOCATE(SOURCEEO(3,4,NSOURCE))
      ALLOCATE(SOURCEN(3,4,NSOURCE))
      ALLOCATE(SOURCET(3,NSOURCE))

      DO ISOUR=1,NSOURCE
        READ(LUNL0,*)   ((SOURCEA(I,J,ISOUR),I=1,3),J=1,4)
        READ(LUNL0,*)   ((SOURCEE(I,J,ISOUR),I=1,3),J=1,4)
      ENDDO

      CLOSE(LUNL0)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     *** SR RFILL0: SOURCE READ FROM FILE ***'
      WRITE(LUNGFO,*)'     FILE:,',FILEL0
      WRITE(LUNGFO,*)'     RUN NUMBER:',ICODL0
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     ATTENTION: INPUT PARAMETERS OVERWRITTEN'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     PARAMETERS:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     WGWINFC:',WGWINFC
      WRITE(LUNGFO,*)'     WBL0CUT:',WBL0CUT
      WRITE(LUNGFO,*)'     WBL0HYS:',WBL0HYS
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     PINHOLES OF COLLIMATOR:'
      WRITE(LUNGFO,*)'     (X,Y,Z), WIDTH, HIGHT:'
      WRITE(LUNGFO,*)'     ',SNGL(CX1),SNGL(CY1),SNGL(CZ1)
     &                        ,SNGL(WID1),SNGL(HIG1)
      WRITE(LUNGFO,*)'     ',SNGL(CX2),SNGL(CY2),SNGL(CZ2)
     &                       ,SNGL(WID2),SNGL(HIG2)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     NUMBER OF SOURCES:',NSOURCE
      WRITE(LUNGFO,*)'     BEGIN AND END OF SOURCES (XYZ):'
      DO ISOUR=1,NSOURCE
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     ',(SOURCEA(I,1,ISOUR),I=1,3)
        WRITE(LUNGFO,*)'     ',(SOURCEE(I,1,ISOUR),I=1,3)
      ENDDO
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     INTEGRATION IS LIMITED TO X BETWEEN XIANF AND XIEND'
      WRITE(LUNGFO,*)'     XIANF, XIEND:',SNGL(XIANF),SNGL(XIEND)
      WRITE(LUNGFO,*)
      DO ISOUR=1,NSOURCE
      WRITE(LUNGFO,*)'     SLOPES AT BEGIN AND END OF SOURCES (HOR.,VER.):'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &     '     ',SOURCEA(3,2,ISOUR)/SOURCEA(1,2,ISOUR)
     &    ,SOURCEA(2,2,ISOUR)/SOURCEA(1,2,ISOUR)
        WRITE(LUNGFO,*)
     &     '     ',SOURCEE(3,2,ISOUR)/SOURCEE(1,2,ISOUR)
     &    ,SOURCEE(2,2,ISOUR)/SOURCEE(1,2,ISOUR)
      ENDDO
      WRITE(LUNGFO,*)

      RETURN
      END
