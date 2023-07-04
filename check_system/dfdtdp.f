*CMZ :  4.00/13 01/12/2021  09.06.02  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  12.53.41  by  Michael Scheer
*CMZ :  2.33/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.00/00 15/12/98  10.35.46  by  Michael Scheer
*CMZ : 00.02/02 15/01/97  15.18.41  by  Michael Scheer
*CMZ : 00.00/07 18/05/94  14.54.15  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.12  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.24  by  Michael Scheer
*-- Author : Michael Scheer
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
      DOUBLE PRECISION FUNCTION DFDTDP
     &  (Y,PSI,DMYGAMMA,DMYCUR,BANWID,PAR,PER,POWR)

C--- CALCULATE DIPOL SPECTRUM

      IMPLICIT NONE

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      integer ical

      DOUBLE PRECISION Y,PSI,DMYGAMMA,DMYCUR,BANWID,PAR,PER,POWR
      DOUBLE PRECISION X,XX,XX1,XI,BK13,BK23,DBSKR3,CONST

      data ical/0/

      if (ical.eq.0) then
        CONST=3.0d0/4.0d0*ALPHA1/PI1**2*DMYGAMMA**2*BANWID*DMYCUR/ECHARGE1
        ical=1
      endif

      IF (Y.EQ.0.0d0) THEN
        WRITE(6,*)'     E/Ec is zero in function DFDTDP'
        STOP '--- PROGRAM WAVE ABORTED ---'
      ENDIF

      X=DMYGAMMA*PSI
      XX=X*X
      XX1=XX+1.0d0
      XI=Y*DSQRT(XX1)**3/2.0d0
      BK13=DBSKR3(XI,1)
      BK23=DBSKR3(XI,2)
      PAR=CONST*(Y*XX1)**2*BK23**2
      PER=CONST*Y*Y*XX*XX1*BK13**2
C      DFDTDP=CONST*(Y*XX1)**2*(BK23**2+XX/XX1*BK13**2)
      DFDTDP=PAR+PER
C      POWR=7.0d0/16.0d0*(ECHARGE1*(CLIGHT1*10.0d0))**2
C     &       /(SQRT(PSI*PSI+1.0d0/DMYGAMMA/DMYGAMMA))**5
C     &       *(1.0d0+5.0d0/7.0d0*XX/XX1)*1.0d-9*DMYCUR/ECHARGE1
      POWR=7.0d0/16.0d0*ECHARGE1*ECHARGE1/4.0d0/PI1/EPS01
     &  /(SQRT(PSI*PSI+1.0d0/DMYGAMMA/DMYGAMMA))**5
     &  *(1.0d0+5.0d0/7.0d0*XX/XX1)*DMYCUR/ECHARGE1

      RETURN
      END
