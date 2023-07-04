*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  1.01/01 10/12/97  13.24.52  by  Michael Scheer
*CMZ : 00.02/05 24/03/97  11.32.00  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.22.12  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.03  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPBUFINT(X,Y,N,RESULT,SIMPLE)
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

C--- INTEGRATION USING SPLINES

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.

      INTEGER N,IFREQ
      DOUBLE PRECISION X(N),Y(N),RESULT,S2(NDFREQP),SIMPLE,DFREQ2,YSUM

C--- SPLINE COEFFICIENTS

      CALL FSPLINEF(X,Y,N,0.D0,0.D0,S2)

C--- INTEGRATION

      RESULT=0.0D0
      SIMPLE=0.D0

      DO IFREQ=1,N-1

      DFREQ2=(X(IFREQ+1)-X(IFREQ))/2.D0
      YSUM=Y(IFREQ)+Y(IFREQ+1)
      SIMPLE=SIMPLE+YSUM*DFREQ2

      RESULT=RESULT
     &          +DFREQ2
     &          *YSUM
     &          -DFREQ2**3/3.D0
     &          *(S2(IFREQ)+S2(IFREQ+1))

      ENDDO !IFREQ


      RETURN
      END
