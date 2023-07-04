*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.23.06  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.13  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.12  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ERZANA(OPSTARTX,OPENDX,
     &                  OPNX,OPNY,OPNZ,
     &                  OPNFX,OPNFY,OPNFZ,
     &                  ZI,ZPI,YI,YPI,ZF,ZPF,YF,YPF,
     &                  BXI,BYI,BZI,BXF,BYF,BZF,
     &                  AXI,AYI,AZI,AXF,AYF,AZF)
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

      INTEGER ICAL

      DOUBLE PRECISION ZI,ZPI,ZF,ZPF,YI,YPI,YF,YPF
     &        ,OPSTARTX,OPENDX,OPNX,OPNY,OPNZ,OPNFX,OPNFY,OPNFZ
     &        ,BXI,BYI,BZI,BXF,BYF,BZF
     &        ,AXI,AYI,AZI,AXF,AYF,AZF

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

          OPSTARTX=0.
          OPENDX=1.

          OPNX=1.
          OPNY=0.
          OPNZ=0.

          OPNFX=1.
          OPNFY=0.
          OPNFZ=0.

          ICAL=1

      ENDIF

      BXI=0.
      BYI=0.
      BZI=0.

      BXF=0.
      BYF=0.
      BZF=0.

      AXI=0.
      AYI=0.
      AZI=0.

      AXF=0.
      AYF=0.
      AZF=0.

      ZF=ZI+ZPI
      ZPF=ZPI

      YF=YI+YPI
      YPF=YPI

      RETURN
      END
