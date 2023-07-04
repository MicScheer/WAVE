*CMZ :  4.00/13 01/09/2021  14.42.02  by  Michael Scheer
*CMZ :  3.04/00 19/01/2018  12.10.01  by  Michael Scheer
*CMZ :  3.01/00 16/07/2013  09.32.23  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/02 27/06/2012  16.34.34  by  Michael Scheer
*CMZ :  2.54/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.53/05 11/02/2005  09.55.20  by  Michael Scheer
*CMZ :  2.52/14 20/12/2004  17.10.56  by  Michael Scheer
*CMZ :  2.52/09 21/10/2004  15.47.48  by  Michael Scheer
*CMZ :  2.52/06 14/10/2004  09.16.20  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.34/07 04/09/2001  16.15.01  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.26.51  by  Michael Scheer
*CMZ : 00.01/11 11/09/96  17.24.24  by  Michael Scheer
*CMZ : 00.01/10 11/09/96  12.42.14  by  Michael Scheer
*CMZ : 00.01/07 16/03/95  14.21.07  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.45.58  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.05.04  by  Michael Scheer
*CMZ : 00.00/03 29/04/94  10.18.17  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.37  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BFOUmag(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,axout,ayout,azout,im)

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

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.
c+SEQ,CMPARA.
*KEEP,fourier.
      include 'fourier.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION DNULL
      COMPLEX*16 CDEXPOMX,CEXPOMZ,CDEXPOMZ

      double precision xin,yin,zin,xcen,bxout,byout,bzout
      double precision
     &  xk0four,xl0four,yk0four,zk0four,zl0four,
     &  xkfour,xlfour,ykfour,zkfour,zlfour,
     &  DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ,axout,ayout,azout
     &  ,BXH,BYH,BZH,AXH,AYH,AZH,AN,AM,X,ak,expomy,expomy1,dexpomy

      integer imag,im,i,kmag,k

      DATA DNULL/0.0D0/

      if (im.lt.0) then
        imag=-im
      else
        imag=im
      endif

      kmag=0
      do i=1,nfoumags
        if (nint(xfoubounds(1,i)).eq.imag) then
          kmag=i
          exit
        endif
      enddo

      if (kmag.eq.0) then
        stop "*** Error in BFOUMAG: Magnet not found ***"
      endif

      xcen=(xfoubounds(3,kmag)+xfoubounds(2,kmag))/2.0d0
      zl0four=xfoubounds(3,kmag)-xfoubounds(2,kmag)
      if (zl0four.ne.0.0d0) then
        zk0four=twopi1/zl0four
      else
        zk0four=0.0d0
      endif

      xl0four=xfoubounds(4,kmag)
      if (xl0four.ne.0.0d0) then
        xk0four=twopi1/xl0four
      else
        xk0four=0.0d0
      endif
      yk0four=sqrt(xk0four**2+zk0four**2)
      x=dmod(xin-xcen,zl0four)

      if (x.gt.zl0four/2.0d0) then
        x=x-zl0four
      else if (x.lt.-zl0four/2.0d0) then
        x=x+zl0four
      endif

      cdexpomx=cdexp(dcmplx(dnull,xk0four*(-zin)))
      dcsxkx=dreal(cdexpomx)
      dsnxkx=dimag(cdexpomx)

      dexpomy=dexp(yk0four*yin)
      expomy=1.0D0

      cdexpomz=cdexp(dcmplx(dnull,zk0four*x))
      cexpomz=dcmplx(1.0d0,dnull)

      bxh=0.0d0
      byh=foumags(1,kmag)/2.0d0
      bzh=0.0d0

      axh=0.0d0
      ayh=0.0d0
      azh=0.0d0

      do k=2,nint(xfoubounds(5,kmag))

        zkfour=zk0four*k
        xkfour=xk0four
        ykfour=sqrt(zkfour**2+xkfour**2)

        if (xk0four.ne.0.0d0) then
          expomy=dexp(ykfour*yin)
        else
          expomy=expomy*dexpomy
        endif
        expomy1=1.0d0/expomy
        dchyky=(expomy+expomy1)*0.5d0
        dshyky=(expomy-expomy1)*0.5d0

        cexpomz=cexpomz*cdexpomz
        dcszkz=dreal(cexpomz)
        dsnzkz=dimag(cexpomz)

        ak=foumags(k,kmag)

        bxh=bxh-ak*xkfour/ykfour*dsnxkx*dshyky*dcszkz
        byh=byh+ak*                    dcsxkx*dchyky*dcszkz
        bzh=bzh-ak*zkfour/ykfour*dcsxkx*dshyky*dsnzkz

        AXH=AXH+ak/ZKFOUR*DCSXKX*DCHYKY*DSNZKZ
        AZH=AZH+0.0
        AYH=AYH+ak/ZKFOUR*XKFOUR/YKFOUR*DSNXKX*DSHYKY*DSNZKZ

      enddo

      if (im.gt.0) then
        bzout=-bxh
        byout= byh
        bxout= bzh
        azout=-axh
        ayout= ayh
        axout= azh
      else
        bzout= byh
        byout= bxh
        bxout= bzh
        azout= ayh
        ayout= axh
        axout= azh
      endif

      return
      end
