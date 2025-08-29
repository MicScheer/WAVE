*CMZ :  3.03/02 19/11/2015  13.37.52  by  Michael Scheer
*CMZ :  3.02/04 11/03/2015  13.19.06  by  Michael Scheer
*CMZ :  2.70/00 12/11/2012  11.53.09  by  Michael Scheer
*CMZ :  2.68/04 03/09/2012  13.24.02  by  Michael Scheer
*CMZ :  2.68/03 25/08/2012  14.52.34  by  Michael Scheer
*CMZ :  2.66/21 22/11/2011  13.51.05  by  Michael Scheer
*-- Author :    Michael Scheer
      subroutine uradestep(icharge,x,y,z,vx,vy,vz,
     &  efieldx,efieldy,efieldz,
     &  dt,gamma,dgamma)

c Author: Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de

c NO WARRANTY

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

c Simple approach to electrical fields into account
c The energy and gamma of the particle is not changed, but returned in dgamma

      double precision x,y,z,vx,vy,vz,dt,gamma,dgamma
      double precision px,py,pz,qdt,eg,vn,pn,
     &  px0,py0,pz0,pn0,pn2,pn02,de,dpx,dpy,dpz,dtm,
     &  efieldx,efieldy,efieldz,
     &  echarge1,emasskg1,emassg1

      integer icharge

      data echarge1/1.602176462D-19/
      data emasskg1/9.10938188D-31/
      data emassg1/0.510998902D-3/

      if (efieldx.eq.0.0d0.and.efieldy.eq.0.0d0.and.efieldz.eq.0.0d0)
     &    then
        dgamma=0.0d0
        return
      endif

      qdt=icharge*echarge1*dt
      eg=emasskg1*gamma

      px0=vx*eg !SI-units
      py0=vy*eg
      pz0=vz*eg
      pn02=px0*px0+py0*py0+pz0*pz0
      pn0=sqrt(pn02)

      dpx=efieldx*qdt !SI-units
      dpy=efieldy*qdt
      dpz=efieldz*qdt

      px=px0+dpx !SI-units
      py=py0+dpy
      pz=pz0+dpz

      vn=sqrt(vx*vx+vy*vy+vz*vz)
      pn2=px*px+py*py+pz*pz
      pn=sqrt(pn2)

c total momentum and energy are kept!!
      vx=px/pn*vn
      vy=py/pn*vn
      vz=pz/pn*vn

      dtm=0.5d0*dt/(emasskg1*gamma) ! F=dp/dt, m=m_e*gamma, a=F/m*dp/dt/m_e/gamma

      x=x+dtm*dpx
      y=y+dtm*dpy
      z=z+dtm*dpz

      de=(vx*dpx+dpy*vy+dpz*vz)/echarge1/1.0d9 !GeV
      dgamma=de/emassg1

      return
      end
