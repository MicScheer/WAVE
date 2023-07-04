*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/09 08/03/2000  17.52.23  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.04  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BFORM

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

C--- BERECHNET HARD-EDGE B-FELDKONFIGURATION


      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bfeld.
      include 'bfeld.cmn'
*KEND.

      DOUBLE PRECISION B0LM

      IF (IBSYM.EQ.0) STOP '*** S/R BFORM: IBSYM.EQ.0 ***'

CC    FB0M=FB0MFUN(FB0N)

      BBY1=B0FORM
      BBY2=-B0FORM/FB0M
      BBY3=B0FORM/FB0M/FB0N
      BBY4=-B0FORM/FB0N
      BBY5=BBY3

      BBY6=0.0
        BBY7=0.0

      B0LM=FB0N/2.*B0LP

      XM1=0.0
      XP1=B0LP/8.
      XM2=XP1
      XP2=XM2+B0LP/8.
      XM3=XP2
      XP3=XM3+B0LM/8.
      XM4=XP3
      XP4=XM4+B0LM/4.
      XM5=XP4
      XP5=XM5+B0LM/8.

C250991  XM6=XP5+1.
C250991  XP6=XM6+1.
      XM6=XP5  !C250991
      XP6=XM6  !C250991
        XM7=XP7
        XP7=XM7

      RETURN
      END
