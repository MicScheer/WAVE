*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 12/11/2009  16.27.11  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.08.53  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.02.58  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.27  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************************
      SUBROUTINE POWINT(X,Y,IWALL,IMODE,IPOL,ISTART,N)
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER IWALL,IPOL,IMODE,IWALLO,IMODEO,IPOLO,N,ISTART

      DATA IMODEO/0/
      DATA IPOLO/0/
      DATA IWALLO/0/

      REAL*4 X,Y

C--- SPLINE COEFFICIENTS

      IF (
     &      IWALLO.NE.IWALL
     &  .OR.
     &      IPOLO.NE.IPOL
     &  .OR.
     &      IMODEO.NE.IMODE
     &  ) THEN

          CALL POWSPL(IWALL,IMODE,IPOL,ISTART,N,1.D30,1.D30)
          IWALLO=IWALL
          IPOLO=IPOL
          IMODEO=IMODE

      ENDIF !ICAL

          CALL POWSPI(X,Y,IWALL,IMODE,IPOL,ISTART,N)

      RETURN
      END
