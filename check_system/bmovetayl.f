*CMZ :  3.06/00 14/02/2019  21.59.37  by  Michael Scheer
*CMZ :  3.05/20 31/10/2018  10.50.14  by  Michael Scheer
*CMZ :  3.05/13 19/09/2018  08.53.32  by  Michael Scheer
*CMZ :  3.05/04 04/07/2018  13.24.53  by  Michael Scheer
*CMZ :  3.05/00 02/05/2018  11.05.41  by  Michael Scheer
*CMZ :  3.04/01 08/03/2018  16.10.29  by  Michael Scheer
*CMZ :  3.04/00 01/03/2018  15.48.20  by  Michael Scheer
*CMZ :  3.02/04 02/12/2014  16.21.05  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  16.41.15  by  Michael Scheer
*CMZ :  3.01/04 20/05/2014  12.53.10  by  Michael Scheer
*CMZ :  3.01/00 18/07/2013  13.23.00  by  Michael Scheer
*CMZ :  2.68/03 22/08/2012  12.47.44  by  Michael Scheer
*CMZ :  2.68/02 06/06/2012  09.45.08  by  Michael Scheer
*CMZ :  2.68/00 24/05/2012  15.54.07  by  Michael Scheer
*CMZ :  2.67/06 24/05/2012  14.24.46  by  Michael Scheer
*CMZ :  2.67/02 09/05/2012  09.06.48  by  Michael Scheer
*CMZ :  2.66/20 22/11/2011  10.34.00  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/02 14/03/2007  16.37.46  by  Michael Scheer
*CMZ :  2.53/02 24/01/2005  13.48.52  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  12.59.53  by  Michael Scheer
*CMZ :  2.52/16 21/01/2005  17.08.32  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  11.03.35  by  Michael Scheer
*CMZ :  2.20/11 11/04/2001  11.16.32  by  Michael Scheer
*CMZ :  2.20/01 05/12/2000  14.03.36  by  Michael Scheer
*CMZ :  2.16/08 20/10/2000  17.40.03  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  18.01.33  by  Michael Scheer
*CMZ :  2.11/00 11/05/99  14.23.02  by  Michael Scheer
*-- Author :    Michael Scheer   11/05/99
C*******************************************************************************
      Subroutine BMOVETAYL(Xin,yin,zin,VX1,VY1,VZ1,BXIN,BYIN,BZIN,DTIM,
     &  X2,Y2,Z2,VX2,VY2,VZ2,VXP,VYP,VZP,GAMMA,ICHARGE,BMOVECUT,
     &  IUSTEP,IENELOSS,DGAMMA)
C*******************************************************************************
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
C
C
C BMOVE OPTIMIERT, ABER NACH WIE VOR WIRD (VXP,VYP,VZP) AM EINGANG DES ZEIT-
C INTERVALLES BERECHNET
C

C BERECHNET DIE 3-DIM Trajektorie EINES TEILCHENS IN EINEM
C 3-DIM HOMOGENEN Magnetfeld
C
C LAENGEN IN METERN
C GeschwindigkeitEN IN M/SEC
C ZEIT IN SEKUNDEN
C B-FELDER IN TESLA V SEC/M**2
C
C*******************************************************************

      IMPLICIT NONE
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,ustep.
      include 'ustep.cmn'
*KEND.

      DOUBLE PRECISION BBET1,VSENK1,VSENK2,TZ,TZ2,TZ3,TZ4,TZ5,VSENKZYK
      !31.10.2018
      double precision xin,yin,zin

      DOUBLE PRECISION BX2,BY2,BZ2,BSQ,BBET,BUX,BUY,BUZ,V0SQ,VX1,VY1,VZ1,
     &  V0BET,VPAR,VPARX,VPARY,VPARZ,VPARSQ,VSENK,DTIM,BXIN,BYIN,BZIN
      DOUBLE PRECISION X1N,Y1N,Z1N,X2N,Y2N,Z2N,
     &  vx12,vy12,vz12

      DOUBLE PRECISION X2,Y2,Z2,VX2,VY2,VZ2,X1,Y1,Z1,VXP,VYP,VZP,ZYK,GAMMA,SZ,CZ
     &  ,DGAMMA,PDUM

      DOUBLE PRECISION BMOVECUT,dgammao,v0,v12n(3),emom,pel(3),
     &  zp,yp,brho,dz,dy,acc,y12,z12,dk

      INTEGER ICHARGE,IUSTEP,IENELOSS,icorry,icorrz,isbig

*KEEP,photon.
      include 'photon.cmn'
*KEEP,efield.
      include 'efield.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      !31.10.2018
      x1=0.0d0
      y1=0.0d0
      z1=0.0d0

      isbig=0

      PDUM=CGAM1/2.0D0/PI1*CLIGHT1*(CLIGHT1/1.0D9)**2*EMASSG1

      DGAMMA=0.0D0

      IF (ICHARGE.GT.0) THEN
        BX2=-BXIN
        BY2=-BYIN
        BZ2=-BZIN
      ELSE
        BX2=BXIN
        BY2=BYIN
        BZ2=BZIN
      ENDIF

      IF(DABS(BX2).LT.BMOVECUT.AND.DABS(BY2).LT.BMOVECUT
     &    .AND.DABS(BZ2).LT.BMOVECUT) THEN

        VXP=0.0D0
        VYP=0.0D0
        VZP=0.0D0

        X2=X1+VX1*DTIM
        Y2=Y1+VY1*DTIM
        Z2=Z1+VZ1*DTIM

        VX2=VX1
        VY2=VY1
        VZ2=VZ1

        GOTO 999

      ENDIF !B-CUT

      BSQ=BX2*BX2+BY2*BY2+BZ2*BZ2
      BBET=DSQRT(BSQ)
      BBET1=1.0D0/BBET

      BUX=BX2*BBET1
      BUY=BY2*BBET1
      BUZ=BZ2*BBET1

C
C   BETRAG VON V0 PARALELL UND SENKRECHT
C
      V0SQ=VX1*VX1+VY1*VY1+VZ1*VZ1
      v0=sqrt(v0sq)
      V0BET=DSQRT(V0SQ)

C
C  VPAR
C
      VPAR=VX1*BUX+VY1*BUY+VZ1*BUZ
      VPARSQ=VPAR*VPAR
      VPARX=VPAR*BUX
      VPARY=VPAR*BUY
      VPARZ=VPAR*BUZ

c      VSENK2=(V0SQ-VPARSQ)
      vsenk2=(v0bet-vpar)*(v0bet+vpar)

      IF (VSENK2.LE.0.0D0) THEN

C
C  ZEITABLEITUNG DER Geschwindigkeit
C
        VXP=0.D0
        VYP=0.D0
        VZP=0.D0

C
C V(DTIM),X(DTIM) BERECHNEN
C

        X2=X1+VX1*DTIM
        Y2=Y1+VY1*DTIM
        Z2=Z1+VZ1*DTIM

        VX2=VX1
        VY2=VY1
        VZ2=VZ1

        GOTO 999

      ELSE    !(VSENK2.LT.0.0)

        VSENK=DSQRT(VSENK2)
        VSENK1=1.0D0/VSENK

C
C   VEKTOR N1 BERECHNEN
C

        X1N=(VX1-VPAR*BUX)*VSENK1
        Y1N=(VY1-VPAR*BUY)*VSENK1
        Z1N=(VZ1-VPAR*BUZ)*VSENK1

C
C  VEKTOR N2=(BUX,BUY,BUZ) KREUZ N1
C

        X2N = BUY*Z1N - BUZ*Y1N
        Y2N = BUZ*X1N - BUX*Z1N
        Z2N = BUX*Y1N - BUY*X1N

C
C ZYKLOTRONFREQUENZ
C
        ZYK=(ECHARGE1/(GAMMA*EMASSKG1))*BBET
C
C

        TZ=ZYK*DTIM

        IF (TZ.LE.0.03D0) THEN
          TZ2=TZ*TZ
          TZ3=TZ2*TZ
          TZ4=TZ3*TZ
          TZ5=TZ4*TZ
          CZ=1.0D0-TZ2/2.0D0+TZ4/24.0D0
          SZ=TZ-TZ3/6.0D0+TZ5/120.0D0
        ELSE
          CZ=COS(TZ)
          SZ=SIN(TZ)
          isbig=1
        ENDIF

C
C  ZEITABLEITUNG DER Geschwindigkeit
C
cerror 4.7.2018, this is at the beginning of the step
c        VXP=VSENK*ZYK*X2N
c        VYP=VSENK*ZYK*Y2N
c        VZP=VSENK*ZYK*Z2N

C
C V(DTIM),X(DTIM) BERECHNEN
C

        VX2=VPARX + VSENK*(X1N*CZ+X2N*SZ)
        VY2=VPARY + VSENK*(Y1N*CZ+Y2N*SZ)
        VZ2=VPARZ + VSENK*(Z1N*CZ+Z2N*SZ)

C
C X(DTIM) BERECHNEN
C

        VSENKZYK=VSENK/ZYK

        X2=X1+VSENKZYK*(X2N+X1N*SZ-X2N*CZ)+VPARX*DTIM
        Y2=Y1+VSENKZYK*(Y2N+Y1N*SZ-Y2N*CZ)+VPARY*DTIM
        Z2=Z1+VSENKZYK*(Z2N+Z1N*SZ-Z2N*CZ)+VPARZ*DTIM

        emom=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
        brho=emom*1.0d9/clight1

        if (ieneloss.ne.0) then
          dgamma=-pdum*bsq*vsenk2/v0sq*gamma**2*dtim
          if (ieneloss.eq.-1) then
c            v12n(1)=(vx1+vx2)/2.0d0/v0
c            v12n(2)=(vy1+vy2)/2.0d0/v0
c            v12n(3)=(vz1+vz2)/2.0d0/v0
            v12n(1)=vx2/v0
            v12n(2)=vy2/v0
            v12n(3)=vz2/v0
            call photon(x2+xin,y2+yin,z2+zin,
     &        v12n,gamma,bx2,by2,bz2,dgamma,dtim,-1) !dgamma will be overwritten
            if (dgamma.ne.0.0d0) then
c              dgamma=0.0d0
c              dpphoton=0.0d0
              pel=emom*v12n+dpphoton
              dgamma=
     &          sqrt(1.0d0+(pel(1)**2+pel(2)**2+pel(3)**2)/emassg1**2)-
     &          gamma
c              emom=emassg1*dsqrt((gamma+dgamma-1.0d0)*(gamma+dgamma+1.0d0)) !GeV
              emom=sqrt(pel(1)**2+pel(2)**2+pel(3)**2)
              vx2=pel(1)/((gamma+dgamma)*emassg1)*clight1
              vy2=pel(2)/((gamma+dgamma)*emassg1)*clight1
              vz2=pel(3)/((gamma+dgamma)*emassg1)*clight1
            endif !dgamma
          endif !ieneloss .eq. -1
        endif !ieneloss

      ENDIF   !(VSENK.LT.0.0)


cccccccccccccccccccccccccccccccccccccccccc

      acc=icharge*echarge1/(gamma*emasskg1)
cerror 4.7.2018, due to error above, now here
      vx12=(vx1+vx2)/2.0d0
      vy12=(vy1+vy2)/2.0d0
      vz12=(vz1+vz2)/2.0d0

      vxp=acc*(vy12*bzin-vz12*byin)
      vyp=acc*(vz12*bxin-vx12*bzin)
      vzp=acc*(vx12*byin-vy12*bxin)

      if (isbig.eq.0) then

        dk=(clight1*dtim)**2/brho

        dy=abs(bz2*dk)
        y12=abs(y1-y2)
        icorry=0
        if (dy.ge.2.0d0*y12
     &      .or.gamma.gt.1.0d10
     &      ) then
          y2=y1+vy12*dtim
          icorry=1
        endif

        dz=abs(by2*dk)
        z12=abs(z1-z2)
        icorrz=0
        if (dz.ge.2.0d0*z12
     &      .or.gamma.gt.1.0d10
     &      ) then
          z2=z1+vz12*dtim
          icorrz=1
        endif

        if (icorry.eq.1.and.icorrz.eq.1) then
          x2=x1+vx1*dtim
        endif

      endif !(isbig.eq.0) then

cccccccccccccccccccccccccccccccccccccccccc

999   CONTINUE

      IF (kefield.NE.0) THEN
        dgammao=dgamma
        CALL estep(X2,Y2,Z2,VX2,VY2,VZ2,DTIM,gamma,dgamma)
        dgamma=dgammao+dgamma
        VXP=(VX2-VX1)/DTIM
        VYP=(VY2-VY1)/DTIM
        VZP=(VZ2-VZ1)/DTIM
      ENDIF

      gammaustep=gamma
      dgammaustep=dgamma
      bxustep=bx2
      byustep=by2
      bzustep=bz2
      IF (IUSTEP.NE.0) THEN
        dgammao=dgamma
        CALL USTEP(X2,Y2,Z2,VX2,VY2,VZ2,vxp,vyp,vzp,DTIM,gamma,dgamma)
        dgamma=dgammao+dgamma
      ENDIF

      !31.10.2018

      x2=x2+xin
      y2=y2+yin
      z2=z2+zin

      RETURN
      End
