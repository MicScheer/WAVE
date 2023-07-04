*CMZ :  2.67/00 16/02/2012  13.16.00  by  Michael Scheer
*CMZ :  2.53/05 24/02/2005  14.03.18  by  Michael Scheer
*CMZ :  0.99/07 16/02/2004  17.22.16  by  Michael Scheer
*CMZ :  0.99/03 12/02/2004  10.50.28  by  Michael Scheer
*CMZ :  0.99/02 12/02/2004  10.20.35  by  Michael Scheer
*CMZ :  0.99/01 11/02/2004  13.54.20  by  Michael Scheer
*CMZ :  0.99/00 29/01/2004  14.14.29  by  Michael Scheer
*CMZ :  0.00/08 23/01/2004  15.31.25  by  Michael Scheer
*CMZ :  0.00/06 14/01/2004  16.32.20  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  16.08.27  by  Michael Scheer
*CMZ :  0.00/04 19/12/2003  18.32.08  by  Michael Scheer
*-- Author :    Michael Scheer   19/12/2003
      subroutine bpeq(x1,x2,a,b,zin,qx,qy,qz,tiny,
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
     &  reverse,iwarn)

      implicit none

      double precision x1,x2,a,b,zin,z,tiny,twopi,pi,f1,f2,dum1,dum2
      double precision qx,qy,qz,xi1,xi2,abz2,a2z2b2,
     &  are1,aim1,are2,aim2,arei1,aimi1,arei2,aimi2,dphi,dphi1,dphi2,phi1,phi2,
     &  sdum,dum,xpi(5),reverse,x1r,x2r,y1,y2,a21,ab,z2,a2,b2,az,a2z2,
     &  arem,aimm,arep,aimp,phim,phip,xm1,xp1,xm2,xp2

      double precision ra21,rho1,rho2,x2rxz,x1rxz,
     &  arg1,arg2,arg31,arg32,arg4

      integer iwarn
      parameter (pi=3.1415926535897932385d0)
      parameter (twopi=6.2831853071795864769d0)
c      if (zin.lt.0.d0) print *,'neg. z: ',zin

      xm1=1.d30
      xp1=1.d30
      xm2=1.d30
      xp2=1.d30
      dphi1=0.d0
      dphi2=0.d0
      iwarn=0

      if (abs(zin).lt.tiny) then
        z=sign(tiny,zin)
        iwarn=1
      else
        z=zin
      endif

      if (x1.gt.x2) then
        x1r=x2
        x2r=x1
        reverse=-1.d0
      else
        x1r=x1
        x2r=x2
        reverse=1.d0
      endif

c-------------------------------------------------------------------------

      y1=a*x1+b
      y2=a*x2+b
      z2=z*z

      ab=a*b
      a2=a*a
      b2=b*b
      az=a*z
      a2z2=a2*z2
      abz2=ab*z2
      a2z2b2=a2z2+b2

      a21=1.d0+a2
      ra21=Sqrt(a21)

      rho1=Sqrt(x1**2+y1**2+z2)
      rho2=Sqrt(x2**2+y2**2+z2)

      arg1=(ab+x1*a21)/ra21+rho1
      arg2=(ab+x2*a21)/ra21+rho2

      arg31=((b*(y1+rho1)+z2)**2 + (z*(a21*x1+a*(b+rho1)))**2)/(x1**2+z**2)
      arg32=((b*(y2+rho2)+z2)**2 + (z*(a21*x2+a*(b+rho2)))**2)/(x2**2+z**2)

      if (
     &    ra21.eq.0.d0.or.arg2.eq.0.d0.or.arg32.eq.0.d0
     &    .or.sign(1.d0,arg1).ne.sign(1.d0,arg2)
     &    .or.sign(1.d0,arg31).ne.sign(1.d0,arg32)
     &   ) then
        iwarn=10
        return
      endif

      arg4=Log(arg1/arg2)/ra21

      qx=Log(arg31/arg32)/2.d0-a*arg4

c-------------------------------------------------------------------------

      x2rxz=x2+Sqrt(x2**2+z2)
      x1rxz=x1+Sqrt(x1**2+z2)

      if (x1rxz.eq.0.d0.or.sign(1.d0,x2rxz).ne.sign(1.d0,x1rxz)) then
        iwarn=11
        return
      endif

      qy=Log(x2rxz/x1rxz)+arg4

c-------------------------------------------------------------------------

      f1=((az+b)*(az-b))**2*a2z2b2-4.0d0*abz2**2
      f2=a21*z2+b2;

      sdum=f1*f2

      dum=(((1.0d0+3.0d0*a21)*z2+b2)*a2z2+(a2z2-b2)*b2)*b2

      if (sdum.gt.0.d0.and.dum.ne.0.d0) then

        dum1=2.0d0*sqrt(sdum)*abs(ab)*z2
        dum2=abz2*a2z2b2**2

        xi1=(+dum1-dum2)/dum
        xi2=(-dum1-dum2)/dum

      else

        xi1=-1.d30
        xi2=-1.d30

      endif

c--------------

      if (xi1.gt.xi2) then
        xpi(1)=xi1
        xi1=xi2
        xi2=xpi(1)
      else if (xi1.eq.xi2) then
        xi2=-1.d30
      endif

c On the way from x1 to x2, the atan2 may be not continuous. So we check,
c xi1, and xi2 and add +/- pi, respectivly

      if (xi1.lt.x1r.or.xi1.gt.x2r) then
        xi1=-1.d30
      endif

      if (xi2.lt.x1r.or.xi2.gt.x2r) then
        xi2=-1.d30
      endif

c----------

      if (xi1.ne.-1.d30) then

c is it real null or pi?

        xm1=xi1-tiny
        xp1=xi1+tiny

        call areim(xm1,a,b,z,arem,aimm)
        call areim(xi1,a,b,z,arei1,aimi1)
        call areim(xp1,a,b,z,arep,aimp)

        if (arem.eq.0.d0.and.aimm.eq.0.d0) then
          stop '*** Error in BPEQ: 0/0'
        endif

        if (arep.eq.0.d0.and.aimp.eq.0.d0) then
          stop '*** Error in BPEQ: 0/0'
        endif

        phim=atan2(aimm,arem)
        phip=atan2(aimp,arep)

        if (sign(1.d0,phim).eq.sign(1.d0,phip)) then
          xi1=-1.d30
        else
          dphi1=pi*nint((phim-phip)/pi)
        endif
        if (arei1.gt.tiny) then
          xi1=-1.d30
        endif

      endif !xi1.ne.-1.d30

c--------------

      if (xi2.ne.-1.d30) then

        xm2=xi2-tiny
        xp2=xi2+tiny

        call areim(xm2,a,b,z,arem,aimm)
        call areim(xi2,a,b,z,arei2,aimi2)
        call areim(xp2,a,b,z,arep,aimp)

        if (arem.eq.0.d0.and.aimm.eq.0.d0) then
          stop '*** Error in BPEQ: 0/0'
        endif

        if (arep.eq.0.d0.and.aimp.eq.0.d0) then
          stop '*** Error in BPEQ: 0/0'
        endif

        phim=atan2(aimm,arem)
        phip=atan2(aimp,arep)

        if (sign(1.d0,phim).eq.sign(1.d0,phip)) then
          xi2=-1.d30
        else
          dphi2=pi*nint((phim-phip)/pi)
        endif

        if (arei2.eq.0.d0.and.aimi2.eq.0.d0) then
          stop '*** Error in BPEQ: 0/0'
        endif


        if (arei2.gt.tiny) then
          xi2=-1.d30
        endif

      endif !xi2.ne.-1.d30

      call areim(x1r,a,b,z,are1,aim1)

      if (are1.eq.0.d0.and.aim1.eq.0.d0) then
        iwarn=4
        phi1=0.d0
      else
        phi1=atan2(aim1,are1)
      endif


      call areim(x2r,a,b,z,are2,aim2)

      if (are2.eq.0.d0.and.aim2.eq.0.d0) then
        iwarn=5
        phi2=0.d0
      else
        phi2=atan2(aim2,are2)
      endif


      dphi=dphi1+dphi2

      if (abs(dphi).gt.twopi) then

c probably two nulls detected due to numerical problems, but actually only one


        if (dphi.gt.twopi) dphi=dphi-twopi
        if (dphi.lt.-twopi) dphi=dphi+twopi

        iwarn=6

      endif

      dphi=phi2-phi1+dphi

      qz=-reverse*dphi

      return
      end
