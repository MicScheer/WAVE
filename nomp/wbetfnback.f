*CMZ :  4.00/15 05/04/2022  11.53.19  by  Michael Scheer
*CMZ :  3.04/00 11/01/2018  11.45.58  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.09.17  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.68/02 08/06/2012  09.23.03  by  Michael Scheer
*CMZ :  2.67/02 26/04/2012  14.49.41  by  Michael Scheer
*CMZ :  2.66/09 25/06/2010  12.15.46  by  Michael Scheer
*CMZ :  2.66/07 20/01/2010  16.25.35  by  Michael Scheer
*CMZ :  2.63/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.55.32  by  Michael Scheer
*CMZ :  2.47/12 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  16.15.48  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.36  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.03  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  10.46.38  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.09  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.10  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WBETFNback
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

      implicit none

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,tralin.
      include 'tralin.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEND.

      INTEGER IT,IP,IPP,IPP2,jfail,ical,ifail

      DOUBLE PRECISION X0,Y0,Z0,YP0,ZP0
      DOUBLE PRECISION X1,Y1,Z1,YP1,ZP1
      DOUBLE PRECISION X2,Y2,Z2,t2
      DOUBLE PRECISION XF0,YF0,ZF0
      DOUBLE PRECISION X1T,Y1T,Z1T

      DOUBLE PRECISION V0,VX0,VY0,VZ0
      DOUBLE PRECISION VX1,VY1,VZ1
      DOUBLE PRECISION VX2,VY2,VZ2
      DOUBLE PRECISION VXF0,VYF0,VZF0
      DOUBLE PRECISION VX1T,VY1T,VZ1T

      DOUBLE PRECISION UNX,UNY,UNZ
      DOUBLE PRECISION VNX,VNY,VNZ,VN
      DOUBLE PRECISION WNX,WNY,WNZ
      DOUBLE PRECISION EWSFX,EWSFY,EWSFZ,BSHIFT

      DOUBLE PRECISION EPS0HO,EPS0VO,DUMP,DUMM,SPLYP0,SPLYPN,gammal

      DOUBLE PRECISION DS2,ALPHA,ALPHAP,BETA,GAMA,
     &  betmxh,betmnh,betcenh,
     &  betmxv,betmnv,betcenv

      double precision X2E,Y2E,Z2E,VX2E,VY2E,VZ2E,GAMMAE,
     &  xpar(3),ypar(3),a(3),yp(3),xopt,yopt

      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::
     &  XBUFF,YBUFF,Y2BUFF,AABUFF,BBBUFF,CCBUFF,CBUFF

      DATA BSHIFT/0.5D0/,ical/0/

      if (ical.ne.0) return
      ical=1

      IF (IENELOSS.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    ' *** WARNING IN WBETFNBACK: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
        WRITE(LUNGFO,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
     &    ' *** WARNING IN WBETFNBACK: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
        WRITE(6,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(6,*)
      ENDIF

      X0=WSXYZ(1,nco)
      Y0=WSXYZ(2,nco)
      Z0=WSXYZ(3,nco)
      XF0=WSXYZ(1,1)

      VX0=-WVXYZ(1,nco)
      VY0=-WVXYZ(2,nco)
      VZ0=-WVXYZ(3,nco)
      V0=SQRT(VX0**2+VY0**2+VZ0**2)

      YP0=VY0/VX0
      ZP0=VZ0/VX0

      IF (YP0.NE.0.0 .OR. Y0.NE.0
     &   .OR.ZP0.NE.0.0 .OR. Z0.NE.0) THEN

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** WARNING IN WBETFNBACK ***'
          WRITE(LUNGFO,*) 'START VALUES OF REFERENCE ORBIT NOT ZERO !!??'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*) '*** WARNING IN WBETFNBACK ***'
          WRITE(6,*) 'START VALUES OF REFERENCE ORBIT NOT ZERO !!??'
          WRITE(6,*)
C         STOP '*** ERROR IN WBETFNBACK ***'

      ENDIF

      EPS0HO=EPS0H
      EPS0VO=EPS0V

      IF (EPS0H.LT.5.D-9) THEN

          EPS0H=5.D-9

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** WARNING SR WBETFNBACK ***'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) 'HORIZ. EMITTANCE LOWER THAN 5.E-9'
          WRITE(LUNGFO,*) 'FOR CALCULATION OF BETA-FUNCTION SET 5.E-9'
          WRITE(LUNGFO,*) 'OLD VALUE RESTORED AFTERWARD'
          WRITE(LUNGFO,*)

      ENDIF

      IF (EPS0V.LT.5.D-9) THEN

          EPS0V=5.D-9

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** WARNING SR WBETFNBACK ***'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) 'VERTICAL. EMITTANCE LOWER THAN 5.E-9'
          WRITE(LUNGFO,*) 'FOR CALCULATION OF BETA-FUNCTION SET 5.E-9'
          WRITE(LUNGFO,*) 'OLD VALUE RESTORED AFTERWARD'
          WRITE(LUNGFO,*)

      ENDIF

1     continue

C--- CALCULATE PRINCIPAL TRAJECTORIES

      DO IT=1,4

C--- REFERENCE ORBIT

        X0=WSXYZ(1,nco)
        Y0=WSXYZ(2,nco)
        Z0=WSXYZ(3,nco)

        VX0=-WVXYZ(1,nco)
        VY0=-WVXYZ(2,nco)
        VZ0=-WVXYZ(3,nco)

        YP0=VY0/VX0
        ZP0=VZ0/VX0

C--- NORMAL VECTOR OF ENTRANCE PLANE I.E Ex-VECTOR OF REFERENCE SYSTEM

        UNX=VX0/V0
        UNY=VY0/V0
        UNZ=VZ0/V0

C--- VECTOR Ey OF REFERENCE SYSTEM

        VN=SQRT(UNX**2+UNY**2)
        VNX=-UNY/VN
        VNY= UNX/VN
        VNZ=0.0d0

C--- VECTOR Ez OF REFERENCE SYSTEM

        WNX=UNY*VNZ-UNZ*VNY
        WNY=UNZ*VNX-UNX*VNZ
        WNZ=UNX*VNY-UNY*VNX

C COORDINATES IN REFERENCE SYSTEM

        IF (IT.EQ.1) THEN

          ALPHA=-.5D0*BETAPH

          X1=0.0d0
          Y1=0.0d0
          Z1=SQRT(EPS0H*BETAH)

          YP1=0.0d0
          ZP1=-SQRT(EPS0H/BETAH)*ALPHA

        ELSE IF (IT.EQ.2) THEN

          X1=0.0d0
          Y1=0.0d0
          Z1=0.0d0

          YP1=0.0d0
          ZP1=-SQRT(EPS0H/BETAH)

        ELSE IF (IT.EQ.3) THEN

          ALPHA=-.5D0*BETAPV

          X1=0.0d0
          Y1=SQRT(EPS0V*BETAV)
          Z1=0.0d0

          YP1=-SQRT(EPS0V/BETAV)*ALPHA
          ZP1=0.0d0

        ELSE IF (IT.EQ.4) THEN

          X1=0.0d0
          Y1=0.0d0
          Z1=0.0d0

          YP1=-SQRT(EPS0V/BETAV)
          ZP1=0.0d0

        ENDIF

C ABSOLUTE COORDINATES (LAB.-SYSTEM)

        X1T=X0+X1*UNX+Y1*VNX+Z1*WNX
        Y1T=Y0+X1*UNY+Y1*VNY+Z1*WNY
        Z1T=Z0+X1*UNZ+Y1*VNZ+Z1*WNZ

        X1=X1T
        Y1=Y1T
        Z1=Z1T

        VX1=V0/SQRT(1.D0+ZP1**2+YP1**2)
        VY1=VX1*YP1
        VZ1=VX1*ZP1

        VX1T=VX1*UNX+VY1*VNX+VZ1*WNX
        VY1T=VX1*UNY+VY1*VNY+VZ1*WNY
        VZ1T=VX1*UNZ+VY1*VNZ+VZ1*WNZ

        VX1=VX1T
        VY1=VY1T
        VZ1=VZ1T

        ZP1=VZ1/VX1
        YP1=VY1/VX1

        gammae=dmygamma*(1.0d0+deltae)

        IF (IT.LE.2) THEN
          WBETA(IT+1,1)=(X1-X0)**2+(Z1-Z0)**2
          WBETA(9+IT+1,1)=(X1-X0)**2+(Z1-Z0)**2
        ELSE
          WBETA(IT+1,1)=(X1-X0)**2+(Y1-Y0)**2
          WBETA(9+IT+1,1)=(X1-X0)**2+(Y1-Y0)**2
        ENDIF

C--- LOOP OVER POINTS

        WBETA(1,1)=x0

        DO IP=1,NCO-1

C--- MOVE PARTICLE TO PLANE THAT CORRESPOND TO NEXT POINT OF REF.ORBIT

          IPP=IP+1

          XF0=WSXYZ(1,IPP)
          YF0=WSXYZ(2,IPP)
          ZF0=WSXYZ(3,IPP)

          VXF0=WVXYZ(1,IPP)
          VYF0=WVXYZ(2,IPP)
          VZF0=WVXYZ(3,IPP)

          EWSFX=VXF0/V0
          EWSFY=VYF0/V0
          EWSFZ=VZF0/V0

C For the particle with energy deviation, we also take
c the normal reference orbit!?

          nutrack=-3

          CALL TRACKSHORT(ISNORDER,X1,Y1,Z1,VX1,VY1,VZ1,
     &      XF0,YF0,ZF0,EWSFX,EWSFY,EWSFZ,
     &      X2E,Y2E,Z2E,T2,VX2E,VY2E,VZ2E,DTIM0,BSHIFT,GAMMAE,BMOVECUT,IUSTEP,
     &      IENELOSS,GAMMAL)

          nutrack=-3

          CALL TRACKSHORT(ISNORDER,X1,Y1,Z1,VX1,VY1,VZ1,
     &      XF0,YF0,ZF0,EWSFX,EWSFY,EWSFZ,
     &      X2,Y2,Z2,T2,VX2,VY2,VZ2,DTIM0,BSHIFT,DMYGAMMA,BMOVECUT,IUSTEP,
     &      IENELOSS,GAMMAL)

          IF (IT.LE.2) THEN
            WBETA(IT+1,IPP)=(X2-XF0)**2+(Z2-ZF0)**2
            WBETA(9+IT+1,IPP)=(X2e-XF0)**2+(Z2e-ZF0)**2
          ELSE
            WBETA(IT+1,IPP)=(X2-XF0)**2+(Y2-YF0)**2
            WBETA(9+IT+1,IPP)=(X2e-XF0)**2+(Y2e-YF0)**2
          ENDIF

          if(it.eq.1) WBETA(1,IPP)=wbeta(1,ip)+(x2-x1)

          X1=X2
          Y1=Y2
          Z1=Z2

          VX1=VX2
          VY1=VY2
          VZ1=VZ2

        ENDDO !IP

      ENDDO !IT

c{30.4.2010: Correct for different step size of last step

      if (nco.gt.2) then
        x2=wbeta(1,nco-1)+ds0
        xpar(1)=wbeta(1,nco-2)
        xpar(2)=wbeta(1,nco-1)
        xpar(3)=wbeta(1,nco)
        do it=1,4
          ypar(1)=wbeta(it+1,nco-2)
          ypar(2)=wbeta(it+1,nco-1)
          ypar(3)=wbeta(it+1,nco)
          call UTIL_PARABEL(xpar,ypar,A,YP,XOPT,yopt,IFAIL)
          if (ifail.ne.0) then
            wbeta(it+1,nco)=wbeta(it+1,nco-1)+
     &        (wbeta(it+1,nco)-wbeta(it+1,nco-1))/
     &        (wbeta(1,nco)-wbeta(1,nco-1))*ds0
          else
            wbeta(it+1,nco)=a(1)+a(2)*x2+a(3)*x2**2
          endif
        enddo
        do it=10,13
          ypar(1)=wbeta(it+1,nco-2)
          ypar(2)=wbeta(it+1,nco-1)
          ypar(3)=wbeta(it+1,nco)
          call UTIL_PARABEL(xpar,ypar,A,YP,XOPT,yopt,IFAIL)
          if (ifail.ne.0) then
            wbeta(it+1,nco)=wbeta(it+1,nco-1)+
     &        (wbeta(it+1,nco)-wbeta(it+1,nco-1))/
     &        (wbeta(1,nco)-wbeta(1,nco-1))*ds0
          else
            wbeta(it+1,nco)=a(1)+a(2)*x2+a(3)*x2**2
          endif
        enddo
        wbeta(1,nco)=x2
      endif

c}30.4.2010: Correct for different step size of last step

C--- CALCULATE BETA-FUNCTIONS

      DO IP=1,NCO

C21.7.92     WBETA(1,IP)=WSXYZ(1,IP)
c16.12.2009          WBETA(1,IP)=-DFLOAT(NCO-1)*DS0/2.+DFLOAT(IP-1)*DS0
c30.4.2010        WBETA(1,IP)=x0+DFLOAT(IP-1)*DS0

        WBETA(2,IP)=(WBETA(2,IP)+WBETA(3,IP))/EPS0H
        WBETA(4,IP)=(WBETA(4,IP)+WBETA(5,IP))/EPS0V
        WBETA(11,IP)=(WBETA(11,IP)+WBETA(12,IP))/EPS0H
        WBETA(13,IP)=(WBETA(13,IP)+WBETA(14,IP))/EPS0V

      ENDDO !IP

C{ CALCULATE DERIVATIVES AND INTEGRALS

      ALLOCATE(XBUFF(NCO))
      ALLOCATE(YBUFF(NCO))
      ALLOCATE(Y2BUFF(NCO))
      ALLOCATE(AABUFF(NCO))
      ALLOCATE(BBBUFF(NCO))
      ALLOCATE(CCBUFF(NCO))
      ALLOCATE(CBUFF(NCO))

c{ e-dev

      DO IPP=0,1

        IPP2=IPP*2

        DO IP=1,NCO
          XBUFF(IP)=WBETA(1,IP)
          YBUFF(IP)=WBETA(11+IPP2,IP)
        ENDDO

        if (ipp.eq.0) then
          call util_min_parabel(nco,xbuff,ybuff,xmx,betmnh,aabuff,bbbuff,jfail)
          call util_max_parabel(nco,xbuff,ybuff,xmx,betmxh,aabuff,bbbuff,jfail)
          betcenh=ybuff(nco/2+1)
        endif

        if (ipp.eq.1) then
          call util_min_parabel(nco,xbuff,ybuff,xmx,betmnv,aabuff,bbbuff,jfail)
          call util_max_parabel(nco,xbuff,ybuff,xmx,betmxv,aabuff,bbbuff,jfail)
          betcenv=ybuff(nco/2+1)
        endif

        DUMP=(YBUFF(3)-YBUFF(2))/DS0
        DUMM=(YBUFF(2)-YBUFF(1))/DS0
        SPLYP0=(DUMP-DUMM)/DS0

        DUMP=(YBUFF(NCO)-YBUFF(NCO-1))/DS0
        DUMM=(YBUFF(NCO-1)-YBUFF(NCO-2))/DS0
        SPLYPN=(DUMP-DUMM)/DS0

        CALL UTIL_SPLINE_COEF(XBUFF,YBUFF,NCO,SPLYP0,SPLYPN,Y2BUFF,
     &    AABUFF,BBBUFF,CCBUFF,CBUFF)

        CALL UTIL_SPLINE_INTER_DERIV(XBUFF,YBUFF,Y2BUFF,NCO,
     &    XBUFF(1),WBETA(11+IPP2,1),WBETA(12+IPP2,1),-1)

        DO IP=2,NCO

          CALL UTIL_SPLINE_INTER_DERIV(XBUFF,YBUFF,Y2BUFF,NCO,
     &      XBUFF(IP),WBETA(11+IPP2,IP),WBETA(12+IPP2,IP),0)

        ENDDO

      ENDDO ! HORI and VERT

      DO IPP=0,1

        IPP2=IPP*2

        DO IP=1,NCO
          XBUFF(IP)=WBETA(1,IP)
          YBUFF(IP)=1.D0/WBETA(11+IPP2,IP)
        ENDDO

        DUMP=(YBUFF(3)-YBUFF(2))/DS0
        DUMM=(YBUFF(2)-YBUFF(1))/DS0
        SPLYP0=(DUMP-DUMM)/DS0

        DUMP=(YBUFF(NCO)-YBUFF(NCO-1))/DS0
        DUMM=(YBUFF(NCO-1)-YBUFF(NCO-2))/DS0
        SPLYPN=(DUMP-DUMM)/DS0

      ENDDO ! HORI and VERT

c} e-dev

      DO IPP=0,1

        IPP2=IPP*2

        DO IP=1,NCO
          XBUFF(IP)=WBETA(1,IP)
          YBUFF(IP)=WBETA(2+IPP2,IP)
        ENDDO

        if (ipp.eq.0) then
          call util_min_parabel(nco,xbuff,ybuff,xmx,betmnh,aabuff,bbbuff,jfail)
          call util_max_parabel(nco,xbuff,ybuff,xmx,betmxh,aabuff,bbbuff,jfail)
          betcenh=ybuff(nco/2+1)
        endif

        if (ipp.eq.1) then
          call util_min_parabel(nco,xbuff,ybuff,xmx,betmnv,aabuff,bbbuff,jfail)
          call util_max_parabel(nco,xbuff,ybuff,xmx,betmxv,aabuff,bbbuff,jfail)
          betcenv=ybuff(nco/2+1)
        endif

        DUMP=(YBUFF(3)-YBUFF(2))/DS0
        DUMM=(YBUFF(2)-YBUFF(1))/DS0
        SPLYP0=(DUMP-DUMM)/DS0

        DUMP=(YBUFF(NCO)-YBUFF(NCO-1))/DS0
        DUMM=(YBUFF(NCO-1)-YBUFF(NCO-2))/DS0
        SPLYPN=(DUMP-DUMM)/DS0

        CALL UTIL_SPLINE_COEF(XBUFF,YBUFF,NCO,SPLYP0,SPLYPN,Y2BUFF,
     &    AABUFF,BBBUFF,CCBUFF,CBUFF)

        CALL UTIL_SPLINE_INTER_DERIV(XBUFF,YBUFF,Y2BUFF,NCO,
     &    XBUFF(1),WBETA(2+IPP2,1),WBETA(3+IPP2,1),-1)

        DO IP=2,NCO

          CALL UTIL_SPLINE_INTER_DERIV(XBUFF,YBUFF,Y2BUFF,NCO,
     &      XBUFF(IP),WBETA(2+IPP2,IP),WBETA(3+IPP2,IP),0)

        ENDDO

      ENDDO ! HORI and VERT

c{e-dev
      DO IPP=0,1

        IPP2=IPP*2

        DO IP=1,NCO
          XBUFF(IP)=WBETA(1,IP)
          YBUFF(IP)=1.D0/WBETA(11+IPP2,IP)
        ENDDO

        DUMP=(YBUFF(12)-YBUFF(2))/DS0
        DUMM=(YBUFF(11)-YBUFF(1))/DS0
        SPLYP0=(DUMP-DUMM)/DS0

        DUMP=(YBUFF(NCO)-YBUFF(NCO-1))/DS0
        DUMM=(YBUFF(NCO-1)-YBUFF(NCO-2))/DS0
        SPLYPN=(DUMP-DUMM)/DS0

        CALL UTIL_SPLINE_COEF(XBUFF,YBUFF,NCO,SPLYP0,SPLYPN,Y2BUFF,
     &    AABUFF,BBBUFF,CCBUFF,CBUFF)

        WBETA(14+IPP,1)=0.D0

        DO IP=2,NCO

          WBETA(14+IPP,IP) = WBETA(14+IPP,IP-1)
     &      +(XBUFF(IP)-XBUFF(IP-1))*0.5D0
     &      *(YBUFF(IP-1)+YBUFF(IP))
     &      -(XBUFF(IP)-XBUFF(IP-1))**3/24.D0
     &      *(Y2BUFF(IP-1)+Y2BUFF(IP))

        ENDDO

      ENDDO ! HORI and VERT

c}e-dev

      DO IPP=0,1

        IPP2=IPP*2

        DO IP=1,NCO
          XBUFF(IP)=WBETA(1,IP)
          YBUFF(IP)=1.D0/WBETA(2+IPP2,IP)
        ENDDO

        DUMP=(YBUFF(3)-YBUFF(2))/DS0
        DUMM=(YBUFF(2)-YBUFF(1))/DS0
        SPLYP0=(DUMP-DUMM)/DS0

        DUMP=(YBUFF(NCO)-YBUFF(NCO-1))/DS0
        DUMM=(YBUFF(NCO-1)-YBUFF(NCO-2))/DS0
        SPLYPN=(DUMP-DUMM)/DS0

        CALL UTIL_SPLINE_COEF(XBUFF,YBUFF,NCO,SPLYP0,SPLYPN,Y2BUFF,
     &    AABUFF,BBBUFF,CCBUFF,CBUFF)

        WBETA(8+IPP,1)=0.D0

        DO IP=2,NCO

          WBETA(8+IPP,IP) = WBETA(8+IPP,IP-1)
     &      +(XBUFF(IP)-XBUFF(IP-1))*0.5D0
     &      *(YBUFF(IP-1)+YBUFF(IP))
     &      -(XBUFF(IP)-XBUFF(IP-1))**3/24.D0
     &      *(Y2BUFF(IP-1)+Y2BUFF(IP))

        ENDDO

      ENDDO ! HORI and VERT

      DEALLOCATE(XBUFF)
      DEALLOCATE(YBUFF)
      DEALLOCATE(Y2BUFF)
      DEALLOCATE(AABUFF)
      DEALLOCATE(BBBUFF)
      DEALLOCATE(CCBUFF)
      DEALLOCATE(CBUFF)

C} CALCULATE DERIVATIVES AND INTEGRALS

      DS2=DS0*2.D0

      DO IP=3,NCO-2
        WBETAK(1,IP)=(WBETA(3,IP+1)-WBETA(3,IP-1))/DS2
        WBETAK(2,IP)=(WBETA(5,IP+1)-WBETA(5,IP-1))/DS2
      ENDDO

      WBETAK(1,2)=WBETAK(1,3)-(WBETAK(1,4)-WBETAK(1,3))
      WBETAK(1,NCO-1)=WBETAK(1,NCO-2)+(WBETAK(1,NCO-2)-WBETAK(1,NCO-3))
      WBETAK(2,2)=WBETAK(2,3)-(WBETAK(2,4)-WBETAK(2,3))
      WBETAK(2,NCO-1)=WBETAK(2,NCO-2)+(WBETAK(2,NCO-2)-WBETAK(2,NCO-3))

      WBETAK(1,1)=WBETAK(1,2)-(WBETAK(1,3)-WBETAK(1,2))
      WBETAK(1,NCO)=WBETAK(1,NCO-1)+(WBETAK(1,NCO-1)-WBETAK(1,NCO-2))
      WBETAK(2,1)=WBETAK(2,2)-(WBETAK(2,3)-WBETAK(2,2))
      WBETAK(2,NCO)=WBETAK(2,NCO-1)+(WBETAK(2,NCO-1)-WBETAK(2,NCO-2))

C--- DETERMINATION OF GRADIENT Kz=-(kz-1/rho**2)
C    AND Ky=ky=-kz

      DO IP=1,NCO
        DO IT=1,2

          ALPHA=-WBETA(1+IT*2,IP)/2.D0
          ALPHAP=-WBETAK(IT,IP)/2.D0
          BETA=WBETA(2*IT,IP)
          GAMA=(1+ALPHA**2)/BETA
          WBETAK(IT,IP)=(ALPHAP+GAMA)/BETA

        ENDDO !IT
      ENDDO !IP

C--- RECALCULATE MAGNETIC FIELD FROM Ky+Kz=1/rho**2

      DO IP=1,NCO

        IF(WBETAK(2,IP)+WBETAK(1,IP).GE.0.0d0) THEN
          WBETAK(3,IP)=SQRT(WBETAK(2,IP)+WBETAK(1,IP))*EMOM/CLIGHT1
        ELSE
          WBETAK(3,IP)=0.0d0
        ENDIF

      ENDDO !IP

      EPS0H=EPS0HO
      EPS0V=EPS0VO

      write(lungfo,*)
      write(lungfo,*)'      Subroutine WBETFNBACK:'
      write(lungfo,*)' '

      write(lungfo,*)'      Max., min., and value of horizontal beta-function in the center:'
      write(lungfo,*)' '
      write(lungfo,*)'      ',sngl(betmxh),sngl(betmnh),sngl(betcenh)
      write(lungfo,*)' '

      write(lungfo,*)'      Max., min., and value of vertical beta-function in the center:'
      write(lungfo,*)' '
      write(lungfo,*)'      ',sngl(betmxv),sngl(betmnv),sngl(betcenv)
      write(lungfo,*)' '

      RETURN
      END
