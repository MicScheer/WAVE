*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.47/04 12/03/2003  15.52.58  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.47  by  Michael Scheer
*-- Author : Michael Scheer
C**********************************************************************
      DOUBLE PRECISION FUNCTION FBETFUN(NN,D0M,PHI,R1)
C**********************************************************************
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

C     Formfactor Fbeta, calculated with REDUCE (26.2.91)
c     revised 31.03.92

      IMPLICIT NONE

      DOUBLE PRECISION NN,D0,D0M,PHI,R1,ANS

      D0=-D0M  !GODE RECHNET MIT NEGATIVER ABLAGE UND NEGATIVER DISPERSION

        ANS=(0.00066314559621623d0*(1795444365.9387d0*D0**2*NN**2    !310392
     . +1795444365.9387d0*D0**2+71438461.471411d0*D0*R1*NN**3
     . *PHI**2+100460336.44417d0*D0*R1*NN**2*PHI**2+46323377.360368d0
     &  *D0*R1*NN*PHI**2+71438461.471411d0*D0*R1*PHI**
     . 2+710611.51687843d0*R1**2*NN**4*PHI**4+1998594.8912206d0
     . *R1**2*NN**3*PHI**4+1755729.1592639d0*R1**2*NN**2*
     . PHI**4+1063696.6143274d0*R1**2*NN*PHI**4+852733.82025412d0
     &*R1**2*PHI**4)**0.5)/((4.*NN**2+6.)**0.5*R1*PHI
     . **2)

      FBETFUN=ANS

      RETURN
      END
