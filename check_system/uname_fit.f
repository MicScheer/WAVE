*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/07 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.16/08 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.12  by  Michael Scheer
*CMZ :  1.00/00 08/07/97  10.51.15  by  Michael Scheer
*CMZ : 00.01/10 02/06/96  12.03.03  by  Michael Scheer
*CMZ : 00.01/06 20/02/95  16.06.18  by  Michael Scheer
*CMZ : 00.01/05 01/02/95  14.05.08  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  17.53.53  by  Michael Scheer
*CMZ : 00.00/07 24/05/94  09.48.27  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UNAME_FIT
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

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,modulator.
      include 'modulator.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      INTEGER I
      DOUBLE PRECISION VARFIT(100)

      OPEN(UNIT=99,FILE='UNAME.FIT',STATUS='OLD')
      DO I=1,100
          READ(99,*,END=99)VARFIT(I)
            VARFIT(I)=ASIN(VARFIT(I))*radgra1
      ENDDO
99    CLOSE(99)

      DO I=1,NMAGMOD

          IF (THEROT(I).EQ.11111.) THEN
         THEROT(I)=VARFIT(1)
          else if (THEROT(I).EQ.-11111.) THEN
         THEROT(I)=-VARFIT(1)

          else if (THEROT(I).EQ.22222.) THEN
         THEROT(I)=VARFIT(2)
          else if (THEROT(I).EQ.-22222.) THEN
         THEROT(I)=-VARFIT(2)

          else if (THEROT(I).EQ.33333.) THEN
         THEROT(I)=VARFIT(3)
          else if (THEROT(I).EQ.-33333.) THEN
         THEROT(I)=-VARFIT(3)

          ENDIF

      ENDDO

      RETURN
      END
