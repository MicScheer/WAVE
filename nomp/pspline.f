*CMZ :  3.07/00 11/03/2019  13.47.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.13/03 15/12/99  16.42.06  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.15  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.52  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PSPLINE(ISOUR,IFREQ)
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
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- CALCULATES COEFFICIENTS OF CUBIC SPLINES THAT INTERPOLATE THE
C    INTENSITY INSIDE THE PINHOLE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER ISOUR,IFREQ,IOBSV,IZ,IY

      DOUBLE PRECISION S(NSPLINEP),S2(NSPLINEP)
      DOUBLE PRECISION YPP0,YPPN

C- SPLINES IN Z

      IOBSV=0
      DO IY=1,NOBSVY

        DO IZ=1,NOBSVZ
          IOBSV=IOBSV+1
          S(IZ)=SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))
        ENDDO      !IZ

c        YPP0=1.D30
c        YPPN=1.D30

C       CALL FSPLINEZ(OBSVZ,S,NOBSVZ,YPP0,YPPN,S2)
C060793    CALL FSPLINDX(OBSVZ, S,NOBSVZ,YPP0,YPPN,S2)
        CALL FSPLINDX(OBSVDZ,S,NOBSVZ,0.D0,0.D0,S2)

        DO IZ=1,NOBSVZ
          SPCOEF(IZ+(IY-1)*NOBSVZ)=S2(IZ)
        ENDDO      !IZ
      ENDDO      !IY

      RETURN
      END
