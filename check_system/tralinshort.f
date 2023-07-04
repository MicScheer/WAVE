*CMZ :  3.03/04 11/10/2017  11.28.21  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.61/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.61/00 30/01/2007  18.56.07  by  Michael Scheer
*CMZ :  2.59/02 24/01/2007  14.57.52  by  Michael Scheer
*CMZ :  2.59/01 24/01/2007  14.30.15  by  Michael Scheer
*CMZ :  2.58/00 16/01/2007  16.51.31  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.47/14 01/08/2003  13.37.51  by  Michael Scheer
*CMZ :  2.47/12 01/07/2003  14.03.43  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.34/07 06/09/2001  15.28.35  by  Michael Scheer
*CMZ :  2.34/05 23/08/2001  17.35.09  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ : 00.01/12 10/10/96  15.17.26  by  Michael Scheer
*CMZ : 00.01/10 04/06/96  10.40.08  by  Michael Scheer
*-- Author :    Michael Scheer   01/06/96
      subroutine tralinshort(
     &  xi0,yi0,zi0,ypi0,zpi0,xe0,gamma,tfmh,tfmv,tfmdeh,tfmdev)
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

C SIMPLE ESTIMATE OF LINEAR TRANSFER MATRIX, light version of tralin

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,tralin.
      include 'tralin.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,phasetrack.
      include 'phasetrack.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

        double precision xe0,ydum,zdum,gamma,
     &    tfmh(2,2),tfmv(2,2),
     &    tfmdeh(2,2),tfmdev(2,2)

        DOUBLE PRECISION DTIM,BSHIFT,BETA,V0
      DOUBLE PRECISION XI,YI,ZI,XF,YF,ZF,YPF,ZPF
      DOUBLE PRECISION XI0,YI0,ZI0,YPI0,ZPI0
      DOUBLE PRECISION VXI,VYI,VZI,VXF,VYF,VZF,XF0,YF0,ZF0,ZPF0,YPF0
      DOUBLE PRECISION VXF0,VYF0,VZF0
        DOUBLE PRECISION WXI0,WYI0,WZI0
      DOUBLE PRECISION EWS(3),EWY(3),EWZ(3),EWSF(3),EWYF(3),EWZF(3)
        DOUBLE PRECISION DX,DY,DZ,DYP,DZP,WX,WY,WZ
     &    ,GAMMAL
     &    ,tfmo(4,4),gammao

      DATA BSHIFT/0.5D0/

        tfmo=tfm !if tralin has calculated tfm already
        gammao=gamma

        IF (IENELOSS.NE.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      ' *** WARNING IN TRALINSHORT: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
          WRITE(LUNGFO,*)
     &      ' *** BE CAREFUL!! ***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)
     &      ' *** WARNING IN TRALINSHORT: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
          WRITE(6,*)
     &      ' *** BE CAREFUL!! ***'
          WRITE(6,*)
        ENDIF

        IF (DELTAZ.EQ.0.D0) DELTAZ=0.00001
        IF (DELTAY.EQ.0.D0) DELTAY=0.00001
        IF (DELTAZP.EQ.0.D0) DELTAZP=0.00001
        IF (DELTAYP.EQ.0.D0) DELTAYP=0.00001
        IF (DELTAE.EQ.0.D0) DELTAE=0.01

1     continue

      BETA=DSQRT((1.D0-1.D0/gamma)*(1.D0+1.D0/gamma))
        V0=CLIGHT1*BETA
        DTIM=1.D0/(v0*myinum)

        WXI0=V0/sqrt(1.0d0+(zpi0**2+ypi0**2))
        WYI0=wxi0*ypi0
        WZI0=wxi0*zpi0

C--- CLOSED ORBIT

        EWS(1)=wXI0/v0
        EWS(2)=wYI0/v0
        EWS(3)=wZI0/v0

C EWS = [EWS,(0,1,0)]

      EWZ(1)=-EWS(3)
      EWZ(2)=0.D0
      EWZ(3)=+EWS(1)

C EWY = [EWZ,EWS]

      EWY(1)=EWZ(2)*EWS(3)-EWZ(3)*EWS(2)
      EWY(2)=EWZ(3)*EWS(1)-EWZ(1)*EWS(3)
      EWY(3)=EWZ(1)*EWS(2)-EWZ(2)*EWS(1)

        CALL track(XI0,YI0,ZI0,WXI0,WYI0,WZI0
     &    ,xe0,ydum,zdum,1.0d0,0.0d0,0.0d0
     &    ,XF0,YF0,ZF0,VXF0,VYF0,VZF0,DTIM,BSHIFT,gamma,GAMMAL)

      ZPF0=VZF0/VXF0
      YPF0=VYF0/VXF0

      EWSF(1)=VXF0/V0
      EWSF(2)=VYF0/V0
      EWSF(3)=VZF0/V0

C EWSF = [EWSF,(0,1,0)]

      EWZF(1)=-EWSF(3)
      EWZF(2)=0.D0
      EWZF(3)=+EWSF(1)

C EWYF = [EWZF,EWSF]

      EWYF(1)=EWZF(2)*EWSF(3)-EWZF(3)*EWSF(2)
      EWYF(2)=EWZF(3)*EWSF(1)-EWZF(1)*EWSF(3)
      EWYF(3)=EWZF(1)*EWSF(2)-EWZF(2)*EWSF(1)

C--T(1:4,1):

      DZ=DELTAZ
      DY=0.0D0
      DZP=0.D0
      DYP=0.D0

      XI=XI0+DY*EWY(1)+DZ*EWZ(1)
      YI=YI0+DY*EWY(2)+DZ*EWZ(2)
      ZI=ZI0+DY*EWY(3)+DZ*EWZ(3)

      WX=V0/SQRT(1.D0+DZP**2+DYP**2)
      WY=DYP*WX
      WZ=DZP*WX

      VXI=WX*EWS(1)+WY*EWY(1)+WZ*EWZ(1)
      VYI=WX*EWS(2)+WY*EWY(2)+WZ*EWZ(2)
      VZI=WX*EWS(3)+WY*EWY(3)+WZ*EWZ(3)

      CALL TRACK(XI,YI,ZI,VXI,VYI,VZI
     &            ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &            ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,gamma,GAMMAL)

      DX=XF-XF0
      DY=YF-YF0
      DZ=ZF-ZF0

      XF=DX*EWSF(1)+DY*EWSF(2)+DZ*EWSF(3)
      YF=DX*EWYF(1)+DY*EWYF(2)+DZ*EWYF(3)
      ZF=DX*EWZF(1)+DY*EWZF(2)+DZ*EWZF(3)

      WX=VXF*EWSF(1)+VYF*EWSF(2)+VZF*EWSF(3)
      WY=VXF*EWYF(1)+VYF*EWYF(2)+VZF*EWYF(3)
      WZ=VXF*EWZF(1)+VYF*EWZF(2)+VZF*EWZF(3)

      ZPF=WZ/WX
      YPF=WY/WX

      TFM(1,1)=ZF/DELTAZ
      TFM(2,1)=ZPF/DELTAZ
      TFM(3,1)=YF/DELTAZ
      TFM(4,1)=YPF/DELTAZ

C--T(1:4,2):

      DZ=0.D0
      DY=0.0D0
      DZP=DELTAZP
      DYP=0.D0

      XI=XI0+DY*EWY(1)+DZ*EWZ(1)
      YI=YI0+DY*EWY(2)+DZ*EWZ(2)
      ZI=ZI0+DY*EWY(3)+DZ*EWZ(3)

      WX=V0/SQRT(1.D0+DZP**2+DYP**2)
      WY=DYP*WX
      WZ=DZP*WX

      VXI=WX*EWS(1)+WY*EWY(1)+WZ*EWZ(1)
      VYI=WX*EWS(2)+WY*EWY(2)+WZ*EWZ(2)
      VZI=WX*EWS(3)+WY*EWY(3)+WZ*EWZ(3)

      CALL TRACK(XI,YI,ZI,VXI,VYI,VZI
     &            ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &            ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,gamma,GAMMAL)


      DX=XF-XF0
      DY=YF-YF0
      DZ=ZF-ZF0

      XF=DX*EWSF(1)+DY*EWSF(2)+DZ*EWSF(3)
      YF=DX*EWYF(1)+DY*EWYF(2)+DZ*EWYF(3)
      ZF=DX*EWZF(1)+DY*EWZF(2)+DZ*EWZF(3)

      WX=VXF*EWSF(1)+VYF*EWSF(2)+VZF*EWSF(3)
      WY=VXF*EWYF(1)+VYF*EWYF(2)+VZF*EWYF(3)
      WZ=VXF*EWZF(1)+VYF*EWZF(2)+VZF*EWZF(3)

      ZPF=WZ/WX
      YPF=WY/WX

      TFM(1,2)=ZF/DELTAZP
      TFM(2,2)=ZPF/DELTAZP
      TFM(3,2)=YF/DELTAZP
      TFM(4,2)=YPF/DELTAZP

C--T(1:4,3):

      DZ=0.D0
      DY=DELTAY
      DZP=0.D0
      DYP=0.D0

      XI=XI0+DY*EWY(1)+DZ*EWZ(1)
      YI=YI0+DY*EWY(2)+DZ*EWZ(2)
      ZI=ZI0+DY*EWY(3)+DZ*EWZ(3)

      WX=V0/SQRT(1.D0+DZP**2+DYP**2)
      WY=DYP*WX
      WZ=DZP*WX

      VXI=WX*EWS(1)+WY*EWY(1)+WZ*EWZ(1)
      VYI=WX*EWS(2)+WY*EWY(2)+WZ*EWZ(2)
      VZI=WX*EWS(3)+WY*EWY(3)+WZ*EWZ(3)

      CALL TRACK(XI,YI,ZI,VXI,VYI,VZI
     &            ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &            ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,gamma,GAMMAL)


      DX=XF-XF0
      DY=YF-YF0
      DZ=ZF-ZF0

      XF=DX*EWSF(1)+DY*EWSF(2)+DZ*EWSF(3)
      YF=DX*EWYF(1)+DY*EWYF(2)+DZ*EWYF(3)
      ZF=DX*EWZF(1)+DY*EWZF(2)+DZ*EWZF(3)

      WX=VXF*EWSF(1)+VYF*EWSF(2)+VZF*EWSF(3)
      WY=VXF*EWYF(1)+VYF*EWYF(2)+VZF*EWYF(3)
      WZ=VXF*EWZF(1)+VYF*EWZF(2)+VZF*EWZF(3)

      ZPF=WZ/WX
      YPF=WY/WX

      TFM(1,3)=ZF/DELTAY
      TFM(2,3)=ZPF/DELTAY
      TFM(3,3)=YF/DELTAY
      TFM(4,3)=YPF/DELTAY

C--T(1:4,4):

      DZ=0.D0
      DY=0.0D0
      DZP=0.D0
      DYP=DELTAYP

      XI=XI0+DY*EWY(1)+DZ*EWZ(1)
      YI=YI0+DY*EWY(2)+DZ*EWZ(2)
      ZI=ZI0+DY*EWY(3)+DZ*EWZ(3)

      WX=V0/SQRT(1.D0+DZP**2+DYP**2)
      WY=DYP*WX
      WZ=DZP*WX

      VXI=WX*EWS(1)+WY*EWY(1)+WZ*EWZ(1)
      VYI=WX*EWS(2)+WY*EWY(2)+WZ*EWZ(2)
      VZI=WX*EWS(3)+WY*EWY(3)+WZ*EWZ(3)

      CALL TRACK(XI,YI,ZI,VXI,VYI,VZI
     &            ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &            ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,gamma,GAMMAL)


      DX=XF-XF0
      DY=YF-YF0
      DZ=ZF-ZF0

      XF=DX*EWSF(1)+DY*EWSF(2)+DZ*EWSF(3)
      YF=DX*EWYF(1)+DY*EWYF(2)+DZ*EWYF(3)
      ZF=DX*EWZF(1)+DY*EWZF(2)+DZ*EWZF(3)

      WX=VXF*EWSF(1)+VYF*EWSF(2)+VZF*EWSF(3)
      WY=VXF*EWYF(1)+VYF*EWYF(2)+VZF*EWYF(3)
      WZ=VXF*EWZF(1)+VYF*EWZF(2)+VZF*EWZF(3)

      ZPF=WZ/WX
      YPF=WY/WX

      TFM(1,4)=ZF/DELTAYP
      TFM(2,4)=ZPF/DELTAYP
      TFM(3,4)=YF/DELTAYP
      TFM(4,4)=YPF/DELTAYP

        if (gamma.eq.gammao) then
          tfmh(1:2,1:2)=tfm(1:2,1:2)
          tfmv(1:2,1:2)=tfm(3:4,3:4)
          gamma=gamma*(1.0d0+deltae)
          goto 1
        else
          tfmdeh(1:2,1:2)=tfm(1:2,1:2)
          tfmdev(1:2,1:2)=tfm(3:4,3:4)
c          tfmdeh=tfmh
c          tfmdev=tfmv
        endif

        tfm=tfmo !if tralin has calculated tfm already
        gamma=gammao

        RETURN
      END
