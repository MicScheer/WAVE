*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  3.03/02 14/01/2016  16.34.05  by  Michael Scheer
*CMZ :  3.00/01 02/04/2013  14.03.46  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.09.17  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.68/02 03/07/2012  15.35.20  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.43  by  Michael Scheer
*CMZ :  2.66/07 25/06/2010  12.15.46  by  Michael Scheer
*CMZ :  2.16/08 23/10/2009  09.19.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.05.51  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.33  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HBETA
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

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,wbetaf90u.
      include 'wbetaf90u.cmn'
*KEND.

C--- HISTOGRAMS FOR OPTICAL FUNCTIONS

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

*KEEP,track.
      include 'track.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      REAL*8 BETH,BETPH,BETV,BETPV,ETA
      REAL*8 TANPHI,BANA,BANAC,BPANA,BPANAC,RI,BY

      real x,XI,XE
      real*8 xid,xed,xd,dxd
      INTEGER I,nbin,icycle,k

      XId=WBETA(1,1)-DS0/2.
      XEd=WBETA(1,NCO)+DS0/2.
      XI=xid
      XE=xed

      nbin=nco/ihtrsmp
      call hbook1m(IDBETAH-1 ,'hits in beta histos',nbin,XI,XE,VMX)

      call hbook1m(IDBETAH ,'HOR. BETA-FUNCTION',nbin,XI,XE,VMX)
      call hbook1m(IDBETHPC,'HOR. BETA-FUNCTION, ANALYTICALLY',nbin,XI,XE,VMX)
      call hbook1m(IDBETHP, 'HOR. BETA-FUNCTION, PARABOLIC ANSATZ',nbin,XI,XE,VMX)

      call hbook1m(IDBETAPH,'DERIV. OF HOR. BETA-FUNCTION',nbin,XI,XE,VMX)
      call hbook1m(IDBPHPC, 'DERIV. OF HOR. BETA-FUNCTION, ANA.',nbin,XI,XE,VMX)
      call hbook1m(IDBPHP,  'DERIV. OF HOR. BETA-FUNCTION, PARAB.',nbin,XI,XE,VMX)

      call hbook1m(IDBETAV ,'VERT. BETA-FUNCTION',nbin,XI,XE,VMX)
      call hbook1m(IDBETAPV,'DERIV. OF VERT. BETA-FUNCTION',nbin,XI,XE,VMX)
      call hbook1m(IDETA   ,'DISPERSION',nbin,XI,XE,VMX)

      dxd=(xed-xid)/nbin

      xd=xid+ds0/2.-dxd
      do i=1,nbin
        xd=xd+dxd
        do k=1,nco-1
          if (xd.ge.wbeta(1,k).and.xd.lt.wbeta(1,k+1)) then
            goto 9
          endif
        enddo
9       continue

        BETH =wbeta(2,k)
        BETPH=wbeta(3,k)
        BETV =wbeta(4,k)
        BETPV=wbeta(5,k)
        ETA  =wbeta(6,k)

        BY    =wtra(2,3,k)
        RI    =CLIGHT1*BY/EMOM
        TANPHI=WTRA(3,2,k)/WTRA(1,2,k)
        BANA  =BETFUN+X**2./BETFUN
        BPANA =2.*X/BETFUN
        BANAC =BANA/(1.+TANPHI**2)
        BPANAC=BPANA/(1+TANPHI**2)+2.*BANA*TANPHI/(1.+TANPHI**2)
     &    *RI/(1+TANPHI**2)

        x=xd

        CALL hfillm(IDBETAH-1, X,0.,1.0d0)
        CALL hfillm(IDBETAH, X,0.,BETH)
        CALL hfillm(IDBETHPC,X,0.,BANAC)
        CALL hfillm(IDBETHP ,X,0.,BANA)

        CALL hfillm(IDBETAPH,X,0.,BETPH)
        CALL hfillm(IDBPHPC ,X,0.,BPANAC)
        CALL hfillm(IDBPHP  ,X,0.,BPANA)

        CALL hfillm(IDBETAV,X,0.,BETV)
        CALL hfillm(IDBETAPV,X,0.,BETPV)
        CALL hfillm(IDETA,X,0.,ETA)

      enddo

      call hoperam(idbetah,'/',idbetah-1,idbetah,1.,1.)
      call hoperam(idbethpc,'/',idbetah-1,idbethpc,1.,1.)
      call hoperam(idbethp,'/',idbetah-1,idbethp,1.,1.)
      call hoperam(idbetav,'/',idbetah-1,idbetav,1.,1.)
      call hoperam(idbetapv,'/',idbetah-1,idbetapv,1.,1.)
      call hoperam(ideta,'/',idbetah-1,ideta,1.,1.)

      call mhrout(idbetah-1,icycle,' ')
      call mhrout(idbetah,icycle,' ')
      call mhrout(idbethpc,icycle,' ')
      call mhrout(idbethp,icycle,' ')
      call mhrout(idbetav,icycle,' ')
      call mhrout(idbetapv,icycle,' ')
      call mhrout(ideta,icycle,' ')

      call hdeletm(idbetah-1)
      call hdeletm(idbetah)
      call hdeletm(idbethpc)
      call hdeletm(idbethp)
      call hdeletm(idbetav)
      call hdeletm(idbetapv)
      call hdeletm(ideta)

      RETURN
      END
