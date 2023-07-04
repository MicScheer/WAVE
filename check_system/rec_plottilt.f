*CMZ :  3.01/09 12/08/2014  15.10.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.05  by  Michael Scheer
*CMZ :  1.00/00 23/05/97  12.32.14  by  Michael Scheer
*CMZ : 00.01/08 21/06/95  10.07.32  by  Michael Scheer
*CMZ : 00.01/07 09/03/95  14.40.13  by  Michael Scheer
*CMZ : 00.00/01 03/03/95  15.54.17  by  Michael Scheer
*-- Author : Michael Scheer
C*****************************************************************
      SUBROUTINE REC_PLOTTILT(IMAG,MODE)
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

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEND.

      INTEGER IMAG,MODE
      REAL*4 X(8),Y(8),Z(8),XP(2),YP(2),X0,Y0,Z0
      REAL*4 XXXLEN,YYYLEN,XXLEN,YYLEN
        REAL*4 SINTHE,COSTHE
      IF (PHI(IMAG).NE.0.0) THEN
        WRITE(6,*) '*** ERROR IN REC_PLOTTILT: PHI NOT ZERO ***'
        STOP
      ENDIF

      X0=DX(IMAG)
      Y0=DY(IMAG)
      Z0=DZ(IMAG)

      XXLEN=XLEN(IMAG)
      YYLEN=YLEN(IMAG)

      SINTHE=DSIN(THETA(IMAG))
      COSTHE=DCOS(THETA(IMAG))
      XXXLEN=ABS(XXLEN*costhe-YYLEN*sinthe)
      YYYLEN=ABS(YYLEN*costhe+XXLEN*sinthe)
      X(1)=X0-XXXLEN/2.
      Y(1)=Y0-YYYLEN/2.
      X(6)=X0+XXXLEN/2.
      Y(6)=Y0+YYYLEN/2.

      SINTHE=DSIN(-THETA(IMAG))
      COSTHE=DCOS(-THETA(IMAG))
      XXXLEN=ABS(XXLEN*costhe-YYLEN*sinthe)
      YYYLEN=ABS(YYLEN*costhe+XXLEN*sinthe)
      X(2)=X0+XXXLEN/2.
      Y(2)=Y0-YYYLEN/2.
      X(5)=X0-XXXLEN/2.
      Y(5)=Y0+YYYLEN/2.

      X(4)=X(1)
      Y(4)=Y(1)
      X(8)=X(5)
      Y(8)=Y(5)
      X(3)=X(2)
      Y(3)=Y(2)
      X(7)=X(6)
      Y(7)=Y(6)

      Z(1)=Z0+ZLEN(IMAG)/2.
      Z(4)=Z0-ZLEN(IMAG)/2.
      Z(2)=Z(1)
      Z(5)=Z(1)
      Z(6)=Z(1)
      Z(3)=Z(4)
      Z(7)=Z(4)
      Z(8)=Z(4)

      IF (MODE.EQ.0) THEN
          XP(1)=X(1)
          YP(1)=Z(1)
          XP(2)=X(4)
          YP(2)=Z(4)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(2)
          YP(1)=Z(2)
          XP(2)=X(3)
          YP(2)=Z(3)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(5)
          YP(1)=Z(5)
          XP(2)=X(8)
          YP(2)=Z(8)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(6)
          YP(1)=Z(6)
          XP(2)=X(7)
          YP(2)=Z(7)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(1)
          YP(1)=Z(1)
          XP(2)=X(2)
          YP(2)=Z(2)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(3)
          YP(1)=Z(3)
          XP(2)=X(4)
          YP(2)=Z(4)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(5)
          YP(1)=Z(5)
          XP(2)=X(6)
          YP(2)=Z(6)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(7)
          YP(1)=Z(7)
          XP(2)=X(8)
          YP(2)=Z(8)
          call mshplt_pline(2,xp,yp)
        ELSE   !MODE
          XP(1)=X(1)
          YP(1)=Y(1)
          XP(2)=X(2)
          YP(2)=Y(2)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(6)
          YP(1)=Y(6)
          XP(2)=X(2)
          YP(2)=Y(2)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(5)
          YP(1)=Y(5)
          XP(2)=X(6)
          YP(2)=Y(6)
          call mshplt_pline(2,xp,yp)
          XP(1)=X(1)
          YP(1)=Y(1)
          XP(2)=X(5)
          YP(2)=Y(5)
          call mshplt_pline(2,xp,yp)
      ENDIF !MODE

      RETURN
      END
