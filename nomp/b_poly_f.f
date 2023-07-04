*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/03 15/11/94  19.17.57  by  Michael Scheer
*CMZ :  0.00/03 15/11/94  18.15.46  by  Michael Scheer
*-- Author :Johannes Bahrdt
c*******************************************************************
      Subroutine B_FIELD(X,Y,Z,Bx,By,Bz,Qa0,Qa,XKX,YKY,width,
     &  NFIRSTX,NORDX,NSTEPX,NFIRSTY,NORDY,NSTEPY,NORDP)
c*******************************************************************
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

      INTEGER NFIRSTX,NORDX,NSTEPX,NFIRSTY,NORDY,NSTEPY,NORDP
      INTEGER IORDX,IORDY


      DOUBLE PRECISION X,Y,Z,Bx,By,Bz,Qa0(NORDP),Qa(NORDP,NORDP)
      DOUBLE PRECISION XKX,YKY,width

      DOUBLE PRECISION c0,s0,xn,xm,xknm

      Bx=0.
      By=0.
      Bz=0.

      DO IORDY=NFIRSTY,NORDY,NSTEPY    ! step size = 2
         xn=iordy
           c0=dcos(xn*yky*z)
           s0=-xn*yky*dsin(xn*yky*z)
         bx=bx
         by=by+c0*qa0(iordy)*xn*yky*dcosh(xn*yky*y)
           bz=bz+s0*qa0(iordy)*dsinh(xn*yky*y)
      DO IORDX=NFIRSTX,NORDX,NSTEPX    ! step size = 1
         xm=iordx
         xknm=dsqrt(xn**2*yky**2+xm**2*xkx**2)
         bx=bx+c0*dsin(xm*xkx*y)*xknm*dsinh(xknm*x)*
     &                     (qa(iordy,iordx)/dcosh(xknm*width))
         by=by+c0*xm*xkx*dcos(xm*xkx*y)*dcosh(xknm*x)*
     &                     (qa(iordy,iordx)/dcosh(xknm*width))
         bz=bz+s0*dsin(xm*xkx*y)*dcosh(xknm*x)*
     &                     (qa(iordy,iordx)/dcosh(xknm*width))
      ENDDO
      ENDDO

      RETURN
      END
