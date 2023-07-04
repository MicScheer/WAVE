*CMZ :  4.00/07 07/06/2020  14.37.38  by  Michael Scheer
*CMZ :  3.05/16 09/10/2018  15.04.43  by  Michael Scheer
*CMZ :  3.05/10 08/08/2018  14.39.43  by  Michael Scheer
*CMZ :  3.04/00 12/01/2018  12.51.17  by  Michael Scheer
*CMZ :  3.03/04 04/12/2017  12.15.53  by  Michael Scheer
*CMZ :  3.00/01 02/04/2013  16.13.11  by  Michael Scheer
*CMZ :  2.68/02 01/08/2012  09.56.07  by  Michael Scheer
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
      SUBROUTINE TRALIN
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

C SIMPLE ESTIMATE OF LINEAR TRANSFER MATRIX

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bmessf90.
      include 'bmessf90.cmn'
*KEEP,tralin.
      include 'tralin.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,phasetrack.
      include 'phasetrack.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER I,IFAIL,J,IPHI,iw2(2)

      DOUBLE PRECISION DTIM,BSHIFT,DET,BETA,V0,D10
      DOUBLE PRECISION XI,YI,ZI,XF,YF,ZF,YPF,ZPF
      DOUBLE PRECISION XI0,YI0,ZI0,YPI0,ZPI0
      DOUBLE PRECISION VXI,VYI,VZI,VXF,VYF,VZF,XF0,YF0,ZF0,ZPF0,YPF0
      DOUBLE PRECISION VXI0,VYI0,VZI0,VXF0,VYF0,VZF0
      DOUBLE PRECISION WXI0,WYI0,WZI0
      DOUBLE PRECISION ERRZ,ERRZP,ERRY,ERRYP,ERRZ10,ERRZP10,ERRY10,ERRYP10
      DOUBLE PRECISION EWS(3),EWY(3),EWZ(3),EWSF(3),EWYF(3),EWZF(3)
      DOUBLE PRECISION DX,DY,DZ,DYP,DZP,WX,WY,WZ
      DOUBLE PRECISION RKXKY
      DOUBLE PRECISION XTM(4,2),eps,PHI,DPHI,CS,SN,dz0,dzp0,dy0,dyp0
     &  ,ZIPH,ZPIPH,YIPH,YPIPH,GAMMAL,a,b,g,beta0m(2,2),
     &  betahm(2,2),betavm(2,2),
     &  betahmi(2,2),betavmi(2,2),betahb,betahpb,betavb,betavpb,
     &  tfm22(2,2),tfmtr(2,2),dum22(2,2),bt(2,2),xrefold

      DATA BSHIFT/0.5D0/
      DATA D10/10.D0/

      write(lungfo,'(/a)')"      --- Entered TRALIN ---"
      write(lungfo,*)

      if (irfilb0.ne.0) then
        if (deltay*d10.gt.bmymax) then
          write(6,*)'*** Warning in TRALIN: DELTAY*10 out of range of field map ***'
          write(6,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
          write(lungfo,*)'*** Warning in TRALIN: DELTAY*10 out of range of field map ***'
          write(lungfo,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
        endif
        if (deltayp*(xstop-xstart)*d10.gt.bmymax) then
          write(6,*)'*** Warning in TRALIN: DELTAYP*(XSTOP-XSTART)*10 out of range of field map ***'
          write(6,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
          write(lungfo,*)'*** Warning in TRALIN: DELTAYP*(XSTOP-XSTART)*10 out of range of field map ***'
          write(lungfo,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
        endif
        if (deltaz*d10.gt.bmzmax) then
          write(6,*)'*** Warning in TRALIN: DELTAZ*10 out of range of field map ***'
          write(6,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
          write(lungfo,*)'*** Warning in TRALIN: DELTAZ*10 out of range of field map ***'
          write(lungfo,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
        endif
        if (deltazp*(xstop-xstart)*d10.gt.bmzmax) then
          write(6,*)'*** Warning in TRALIN: DELTAZP*(XSTOP-XSTART)*10 out of range of field map ***'
          write(6,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
          write(lungfo,*)'*** Warning in TRALIN: DELTAZP*(XSTOP-XSTART)*10 out of range of field map ***'
          write(lungfo,*)'*** Estimation of uncertainty of the transfer matrix will fail ***'
        endif
        if (deltay.gt.bmymax) then
          write(6,*)'*** Error in TRALIN: DELTAY out of range of field map ***'
          write(6,*)'*** Program WAVE aborted ***'
          write(lungfo,*)'*** Error in TRALIN: DELTAY out of range of field map ***'
          write(lungfo,*)'*** Program WAVE aborted ***'
          stop
        endif
        if (deltayp*(xstop-xstart).gt.bmymax) then
          write(6,*)'*** Error in TRALIN: DELTAYP*(XSTOP-XSTART) out of range of field map ***'
          write(6,*)'*** Program WAVE aborted ***'
          write(lungfo,*)'*** Error in TRALIN: DELTAYP*(XSTOP-XSTART) out of range of field map ***'
          write(lungfo,*)'*** Program WAVE aborted ***'
          stop
        endif
        if (deltaz.gt.bmzmax) then
          write(lungfo,*)'*** Error in TRALIN: DELTAZ out of range of field map ***'
          write(lungfo,*)'*** Program WAVE aborted ***'
          write(6,*)'*** Error in TRALIN: DELTAZ out of range of field map ***'
          write(6,*)'*** Program WAVE aborted ***'
          stop
        endif
        if (deltazp*(xstop-xstart).gt.bmzmax) then
          write(lungfo,*)'*** Error in TRALIN: DELTAZP*(XSTOP-XSTART) out of range of field map ***'
          write(lungfo,*)'*** Program WAVE aborted ***'
          write(6,*)'*** Error in TRALIN: DELTAZP*(XSTOP-XSTART) out of range of field map ***'
          write(6,*)'*** Program WAVE aborted ***'
          stop
        endif
      endif !(irfilb0.ne.0) then

      IF (IENELOSS.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    ' *** WARNING IN TRALIN: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
        WRITE(LUNGFO,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
     &    ' *** WARNING IN TRALIN: IENELOSS .NE. 0, NOT YET IMPLEMENTED ***'
        WRITE(6,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(6,*)
      ENDIF

      IF (DELTAZ.EQ.0.D0) DELTAZ=0.00001
      IF (DELTAY.EQ.0.D0) DELTAY=0.00001
      IF (DELTAZP.EQ.0.D0) DELTAZP=0.00001
      IF (DELTAYP.EQ.0.D0) DELTAYP=0.00001

      BETA=DSQRT((1.D0-1.D0/DMYGAMMA)*(1.D0+1.D0/DMYGAMMA))
      V0=CLIGHT1*BETA
      WXI0=VXIN*V0
      WYI0=VYIN*V0
      WZI0=VZIN*V0
      DTIM=1.D0/(v0*myinum)

      XI=XSTART
      XI0=XI

C--- CLOSED ORBIT

      ZI0=ZSTART
      YI0=YSTART
      VXI0=VXIN
      VYI0=VYIN
      VZI0=VZIN
      ZPI0=VZIN/VXIN
      YPI0=VYIN/VXIN

      EWS(1)=VXI0
      EWS(2)=VYI0
      EWS(3)=VZI0

C EWS = [EWS,(0,1,0)]

      EWZ(1)=-EWS(3)
      EWZ(2)=0.D0
      EWZ(3)=+EWS(1)

C EWY = [EWZ,EWS]

      EWY(1)=EWZ(2)*EWS(3)-EWZ(3)*EWS(2)
      EWY(2)=EWZ(3)*EWS(1)-EWZ(1)*EWS(3)
      EWY(3)=EWZ(1)*EWS(2)-EWZ(2)*EWS(1)

      CALL TRACK(XI0,YI0,ZI0,WXI0,WYI0,WZI0
     &  ,XSTOP,0.D0,0.D0,1.D0,0.D0,0.D0
     &  ,XF0,YF0,ZF0,VXF0,VYF0,VZF0,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

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
     &  ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &  ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

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
     &  ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &  ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

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
     &  ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &  ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

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
     &  ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &  ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

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

C--- CHECK

      DZ=DELTAZ
      DY=DELTAY
      DZP=DELTAZP
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
     &  ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &  ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

      ZPF=VZF/VXF
      YPF=VYF/VXF

      ERRZ=( ZF-ZF0
     &-(TFM(1,1)*DELTAZ+TFM(1,2)*DELTAZP
     &+ TFM(1,3)*DELTAY+TFM(1,4)*DELTAYP))
     &/DELTAZ
      ERRZP=( ZPF-ZPF0
     &-(TFM(2,1)*DELTAZ+TFM(2,2)*DELTAZP
     &+ TFM(2,3)*DELTAY+TFM(2,4)*DELTAYP))
     &/DELTAZP
      ERRY=( YF-YF0
     &-(TFM(3,1)*DELTAZ+TFM(3,2)*DELTAZP
     &+ TFM(3,3)*DELTAY+TFM(3,4)*DELTAYP))
     &/DELTAY
      ERRYP=( YPF-YPF0
     &-(TFM(4,1)*DELTAZ+TFM(4,2)*DELTAZP
     &+ TFM(4,3)*DELTAY+TFM(4,4)*DELTAYP))
     &/DELTAYP

      DELTAZ=DELTAZ*D10
      DELTAZP=DELTAZP*D10
      DELTAY=DELTAY*D10
      DELTAYP=DELTAYP*D10

      DZ=DELTAZ
      DY=DELTAY
      DZP=DELTAZP
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
     &  ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &  ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

      ZPF=VZF/VXF
      YPF=VYF/VXF

      ERRZ10=( ZF-ZF0
     &-(TFM(1,1)*DELTAZ+TFM(1,2)*DELTAZP
     &+ TFM(1,3)*DELTAY+TFM(1,4)*DELTAYP))
     &/DELTAZ
      ERRZP10=( ZPF-ZPF0
     &-(TFM(2,1)*DELTAZ+TFM(2,2)*DELTAZP
     &+ TFM(2,3)*DELTAY+TFM(2,4)*DELTAYP))
     &/DELTAZP
      ERRY10=( YF-YF0
     &-(TFM(3,1)*DELTAZ+TFM(3,2)*DELTAZP
     &+ TFM(3,3)*DELTAY+TFM(3,4)*DELTAYP))
     &/DELTAY
      ERRYP10=( YPF-YPF0
     &-(TFM(4,1)*DELTAZ+TFM(4,2)*DELTAZP
     &+ TFM(4,3)*DELTAY+TFM(4,4)*DELTAYP))
     &/DELTAYP

      DELTAZ=DELTAZ/D10
      DELTAZP=DELTAZP/D10
      DELTAY=DELTAY/D10
      DELTAYP=DELTAYP/D10

      CALL UTIL_DETERMINANTE(4,TFM,DET,IFAIL)

      IF (TFM(1,2).NE.0.D0) THEN
        TRAXKX=TFM(2,1)/TFM(1,2)
      ELSE
        TRAXKX=0.D0
      ENDIF

      IF (TFM(3,4).NE.0.D0) THEN
        TRAXKY=TFM(4,3)/TFM(3,4)
      ELSE
        TRAXKY=0.D0
      ENDIF

      RKXKY=TRAXKX+TRAXKY
      IF (RKXKY.NE.0.D0) THEN
        TRAXKX=TRAXKX/RKXKY
        TRAXKY=TRAXKY/RKXKY
      ENDIF

C APPLY TRANSFER MATRIX TO YSTART,ZSTART,VYIN/VXIN,VZIN/VXIN

      XTM(1,1)=PHTRZ0 ! SIGN NOT CHECKED
      XTM(2,1)=PHTRZP0 ! SIGN NOT CHECKED
      XTM(3,1)=PHTRY0
      XTM(4,1)=PHTRYP0

      DO I=1,4
        XTM(I,2)=0.D0
        DO J=1,4
          XTM(I,2)=XTM(I,2)+TFM(I,J)*XTM(J,1)
        ENDDO
      ENDDO

c{ invert transfer matrixces

      tfminv=0.0d0

      tfm22(1,1)=tfm(1,1)
      tfm22(1,2)=tfm(1,2)
      tfm22(2,1)=tfm(2,1)
      tfm22(2,2)=tfm(2,2)

      dum22(1,1)=1.0d0
      dum22(1,2)=0.0d0
      dum22(2,1)=0.0d0
      dum22(2,2)=1.0d0

      call deqinv(2,tfm22,2,iw2,ifail,2,dum22)

      tfminv(1,1)=tfm22(1,1)
      tfminv(1,2)=tfm22(1,2)
      tfminv(2,1)=tfm22(2,1)
      tfminv(2,2)=tfm22(2,2)

      tfm22(1,1)=tfm(3,3)
      tfm22(1,2)=tfm(3,4)
      tfm22(2,1)=tfm(4,3)
      tfm22(2,2)=tfm(4,4)

      dum22(1,1)=1.0d0
      dum22(1,2)=0.0d0
      dum22(2,1)=0.0d0
      dum22(2,2)=1.0d0

      call deqinv(2,tfm22,2,iw2,ifail,2,dum22)

      tfminv(3,3)=tfm22(1,1)
      tfminv(3,4)=tfm22(1,2)
      tfminv(4,3)=tfm22(2,1)
      tfminv(4,4)=tfm22(2,2)

c} invert transfer matrixces
      if (ibetback.ne.0) then

        write(lungfo,*)
        write(lungfo,*)'       IBETBACK:',ibetback
        write(lungfo,*)
C{ APPLY TRANSFER MATRIX TO BETA-MATRICES, BUT BACKWARD

        a=-betaph/2.0d0
        g=(1.0d0+a**2)/betah
        b=betah

        beta0m(1,1)=b
        beta0m(1,2)=-a
        beta0m(2,1)=-a
        beta0m(2,2)=g

        tfm22(1:2,1:2)=tfminv(1:2,1:2)

        tfmtr(1,1)=tfm22(1,1)
        tfmtr(1,2)=tfm22(2,1)
        tfmtr(2,1)=tfm22(1,2)
        tfmtr(2,2)=tfm22(2,2)

        call util_matrix_multiplication(2,2,2,tfm22,beta0m,bt,dum22)
        call util_matrix_multiplication(2,2,2,bt,tfmtr,betahmi,dum22)

c          WRITE(LUNGFO,*)
c          WRITE(LUNGFO,*)'       BETAH and BETAPH transformed backward to XSTART:'
c          WRITE(LUNGFO,*)'       ',sngl(betahmi(1,1)),sngl(betahmi(1,2)*2.0d0)
c          WRITE(LUNGFO,*)

        betah=betahmi(1,1)
        betaph=2.0d0*betahmi(1,2)

        a=-betapv/2.0d0
        g=(1.0d0+a**2)/betav
        b=betav

        beta0m(1,1)=b
        beta0m(1,2)=-a
        beta0m(2,1)=-a
        beta0m(2,2)=g

        tfm22(1:2,1:2)=tfminv(3:4,3:4)

        tfmtr(1,1)=tfm22(1,1)
        tfmtr(1,2)=tfm22(2,1)
        tfmtr(2,1)=tfm22(1,2)
        tfmtr(2,2)=tfm22(2,2)

        call util_matrix_multiplication(2,2,2,tfm22,beta0m,bt,dum22)
        call util_matrix_multiplication(2,2,2,bt,tfmtr,betavmi,dum22)

c          WRITE(LUNGFO,*)
c          WRITE(LUNGFO,*)'       BETAV and BETAPV transformed backward to XSTART:'
c          WRITE(LUNGFO,*)'       ',sngl(betavmi(1,1)),sngl(betavmi(1,2)*2.0d0)
c          WRITE(LUNGFO,*)

        betav=betavmi(1,1)
        betapv=2.0d0*betavmi(1,2)

C} APPLY TRANSFER MATRIX TO BETA-MATRICES, BUT BACKWARD

      endif

      if (betfun.eq.-9999.) then

C Calculate periodic solution if required:

        tfm22(1,1)=tfm(1,1)
        tfm22(1,2)=tfm(1,2)
        tfm22(2,1)=tfm(2,1)
        tfm22(2,2)=tfm(2,2)

        a=
     &    tfm22(1,1)**2*tfm22(2,2)**2-tfm22(1,1)**2-
     &    2.0d0*tfm22(1,1)*tfm22(1,2)*tfm22(2,1)*tfm22(2,2)+
     &    tfm22(1,2)**2*tfm22(2,1)**2-
     &    2.0d0*tfm22(1,2)*tfm22(2,1)-tfm22(2,2)**2+1.0d0

        if (a.gt.0.0d0) then

          betah=abs((tfm22(1,1)*tfm22(1,2)*tfm22(2,2)-tfm22(1,2)**2*tfm22(2,1)+
     &      tfm22(1,2))/sqrt(a))
          a=(tfm22(1,1)**2*tfm22(2,2)-tfm22(1,1)*tfm22(1,2)*tfm22(2,1)-
     &      tfm22(2,2))/sqrt(a)
          betaph=-2.0d0*a

        else

          write(6,*)
          write(6,*)
          write(6,*)'*** Warning in TRALIN: No periodic solution found for '
          write(6,*)'*** horizontal beta-function, old values of BETAH and BETAHP kept'
          write(lungfo,*)
          write(lungfo,*)'*** Warning in TRALIN: No periodic solution found for '
          write(lungfo,*)'*** horizontal beta-function, old values of BETAH and BETAHP kept'
          write(lungfo,*)

        endif
      endif

      if (betfunv.eq.-9999.) then

        tfm22(1,1)=tfm(3,3)
        tfm22(1,2)=tfm(3,4)
        tfm22(2,1)=tfm(4,3)
        tfm22(2,2)=tfm(4,4)

        a=
     &    tfm22(1,1)**2*tfm22(2,2)**2-tfm22(1,1)**2-
     &    2.0d0*tfm22(1,1)*tfm22(1,2)*tfm22(2,1)*tfm22(2,2)+
     &    tfm22(1,2)**2*tfm22(2,1)**2-
     &    2.0d0*tfm22(1,2)*tfm22(2,1)-tfm22(2,2)**2+1.0d0

        if (a.gt.0.0d0) then
          betav=abs((tfm22(1,1)*tfm22(1,2)*tfm22(2,2)-tfm22(1,2)**2*tfm22(2,1)+
     &      tfm22(1,2))/sqrt(a))
          a=(tfm22(1,1)**2*tfm22(2,2)-tfm22(1,1)*tfm22(1,2)*tfm22(2,1)-
     &      tfm22(2,2))/sqrt(a)
          betapv=-2.0d0*a
        else
          write(6,*)
          write(6,*)'*** Warning in TRALIN: No periodic solution found for '
          write(6,*)'*** vertical beta-function, old values of BETAV and BETAVP kept'
          write(6,*)
          write(lungfo,*)
          write(lungfo,*)'*** Warning in TRALIN: No periodic solution found for '
          write(lungfo,*)'*** vertical beta-function, old values of BETAV and BETAVP kept'
          write(lungfo,*)
        endif
      endif



C{ APPLY TRANSFER MATRIX TO BETA-MATRICES

      a=-betaph/2.0d0
      g=(1.0d0+a**2)/betah
      b=betah

      beta0m(1,1)=b
      beta0m(1,2)=-a
      beta0m(2,1)=-a
      beta0m(2,2)=g

      tfm22(1,1)=tfm(1,1)
      tfm22(1,2)=tfm(1,2)
      tfm22(2,1)=tfm(2,1)
      tfm22(2,2)=tfm(2,2)

      tfmtr(1,1)=tfm22(1,1)
      tfmtr(1,2)=tfm22(2,1)
      tfmtr(2,1)=tfm22(1,2)
      tfmtr(2,2)=tfm22(2,2)

      call util_matrix_multiplication(2,2,2,tfm22,beta0m,bt,dum22)
      call util_matrix_multiplication(2,2,2,bt,tfmtr,betahm,dum22)

      a=-betapv/2.0d0
      g=(1.0d0+a**2)/betav
      b=betav

      beta0m(1,1)=b
      beta0m(1,2)=-a
      beta0m(2,1)=-a
      beta0m(2,2)=g

      tfm22(1,1)=tfm(3,3)
      tfm22(1,2)=tfm(3,4)
      tfm22(2,1)=tfm(4,3)
      tfm22(2,2)=tfm(4,4)

      tfmtr(1,1)=tfm22(1,1)
      tfmtr(1,2)=tfm22(2,1)
      tfmtr(2,1)=tfm22(1,2)
      tfmtr(2,2)=tfm22(2,2)

      call util_matrix_multiplication(2,2,2,tfm22,beta0m,bt,dum22)
      call util_matrix_multiplication(2,2,2,bt,tfmtr,betavm,dum22)

C} APPLY TRANSFER MATRIX TO BETA-MATRICES

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Linear transfer matrix (raw estimate):'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       deltaz, deltazp, deltay, deltayp:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',SNGL(DELTAZ),SNGL(DELTAZP),SNGL(DELTAY),SNGL(DELTAYP)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',(SNGL(TFM(1,I)),I=1,4)
      WRITE(LUNGFO,*)'       ',(SNGL(TFM(2,I)),I=1,4)
      WRITE(LUNGFO,*)'       ',(SNGL(TFM(3,I)),I=1,4)
      WRITE(LUNGFO,*)'       ',(SNGL(TFM(4,I)),I=1,4)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Relative error for track (deltaz,deltazp,deltay,deltayp)'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',SNGL(ERRZ),SNGL(ERRZP),SNGL(ERRY),SNGL(ERRYP)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Relative error for track 10 X (deltaz,deltazp,deltay,deltayp)'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',SNGL(ERRZ10),SNGL(ERRZP10),SNGL(ERRY10),SNGL(ERRYP10)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Failure flag and determinant:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',IFAIL,DET
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      normalized kx, ky:',SNGL(TRAXKX),SNGL(TRAXKY)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '      Transfer Matrix applied to (PHTRZ0, PHTRZP0, PHTRY0, PHTRYP0):'
      WRITE(LUNGFO,*)'        ',(SNGL(XTM(I,2)),I=1,4)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Inverse linear transfer matrix (raw estimate):'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',(SNGL(tfminv(1,I)),I=1,4)
      WRITE(LUNGFO,*)'       ',(SNGL(tfminv(2,I)),I=1,4)
      WRITE(LUNGFO,*)'       ',(SNGL(tfminv(3,I)),I=1,4)
      WRITE(LUNGFO,*)'       ',(SNGL(tfminv(4,I)),I=1,4)
      WRITE(LUNGFO,*)

      DO I=1,2
        XTM(I,2)=0.D0
        DO J=1,2
          XTM(I,2)=XTM(I,2)+TFM(I,J)*XTM(J,1)
        ENDDO
      ENDDO

      DO I=3,4
        XTM(I,2)=0.D0
        DO J=3,4
          XTM(I,2)=XTM(I,2)+TFM(I,J)*XTM(J,1)
        ENDDO
      ENDDO
      WRITE(LUNGFO,*)
     &  '      dito but without coupling of planes:'
      WRITE(LUNGFO,*)'        ',(SNGL(XTM(I,2)),I=1,4)
      WRITE(LUNGFO,*)

      b=betah
      a=-betaph/2.0d0
      g=(1.0d0+a**2)/betah

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'         Horizontal beta matrix, i.e.:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'              ( beta  -alpha )'
      WRITE(LUNGFO,*)'              ( -alpha gamma )'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',SNGL(b),sngl(-a)
      WRITE(LUNGFO,*)'       ',sngl(-a),sngl(g)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       Beta, BetaP:',sngl(betah),sngl(betaph)
      WRITE(LUNGFO,*)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Transformed hor. beta matrix,'
      WRITE(LUNGFO,*)'      (must be the same for the periodic solution):'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',(SNGL(betahm(1,I)),I=1,2)
      WRITE(LUNGFO,*)'       ',(SNGL(betahm(2,I)),I=1,2)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       Beta, BetaP:',
     &  sngl(betahm(1,1)),sngl(betahm(1,2)*2.0d0)
      WRITE(LUNGFO,*)

      b=betav
      a=-betapv/2.0d0
      g=(1.0d0+a**2)/betav

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'         Vertical beta matrix, i.e.:'
      WRITE(LUNGFO,*)'              ( beta  -alpha )'
      WRITE(LUNGFO,*)'              ( -alpha gamma )'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',SNGL(b),sngl(-a)
      WRITE(LUNGFO,*)'       ',sngl(-a),sngl(g)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       Beta, BetaP:',sngl(betav),sngl(betapv)
      WRITE(LUNGFO,*)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      Transformed ver. beta matrix,'
      WRITE(LUNGFO,*)'      (must be the same for the periodic solution):'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ',(SNGL(betavm(1,I)),I=1,2)
      WRITE(LUNGFO,*)'       ',(SNGL(betavm(2,I)),I=1,2)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       Beta, BetaP:',
     &  sngl(betavm(1,1)),sngl(betavm(1,2)*2.0d0)
      WRITE(LUNGFO,*)

      OPEN(UNIT=99,FILE='tralin.wav',STATUS='unknown')
      REWIND(99)
      WRITE(99,'(4(1PE18.10))')((TFM(I,J),J=1,4),I=1,4)
      CLOSE(99)

C-- phasespace ellipes

      if (iphellip.ne.0.0d0) then

        iphellip=abs(iphellip)

        if (phbetah.ne.0.0d0) then
          eps=phdisph**2/phbetah
          dz0=  sqrt(eps*phbetah)
          dzp0=-sqrt(eps/phbetah)
        else
          eps=0.0d0
          dz0=0.0d0
          dzp0=0.0d0
        endif

        dphi=2.0d0*pi1/iphellip

        DY=PHTRY0
        DYP=PHTRYP0

        YIPH=DY
        YPIPH=DYP

        open(unit=99,file='wave_phase_ellipse_hori.wva',
     &    status='unknown',recl=256)

        write(99,'(a,i10,a,a)')'* ',icode,' ',code
        write(99,'(a,e15.5,a,e15.5)')'* y =',PHTRY0,' yp =',PHTRYP0
        write(99,'(a)')'* Meaning of columns phi, zi, zpi, zf, zpf'

        do iphi=1,iphellip

          phi=(iphi-1)*dphi
          cs=cos(phi)
          sn=sin(phi)

          DZ=  dz0*cs
          DZP=dzp0*sn
          DY=PHTRY0
          DYP=PHTRYP0

          ZIPH=DZ
          ZPIPH=DZP

          XI=XI0+DY*EWY(1)+DZ*EWZ(1)
          YI=YI0+DY*EWY(2)+DZ*EWZ(2)
          ZI=ZI0+DY*EWY(3)+DZ*EWZ(3)

          WX=V0/SQRT(1.D0+DZP**2+DYP**2)
          WY=DYP*WX
          WZ=DZP*WX

          VXI=WX*EWS(1)+WY*EWY(1)+WZ*EWZ(1)
          VYI=WX*EWS(2)+WY*EWY(2)+WZ*EWZ(2)
          VZI=WX*EWS(3)+WY*EWY(3)+WZ*EWZ(3)

          IF (IERZFUN.EQ.0) THEN

            CALL TRACK(XI,YI,ZI,VXI,VYI,VZI
     &        ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &        ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

          ELSE !IERZFUN

            CALL IDTRMSHGF(XI,YI,ZI,VXI,VYI,VZI,XF,YF,ZF,VXF,VYF,VZF,
     &        XI0,YI0,ZI0,WXI0,WYI0,WZI0,
     &        XF0,YF0,ZF0,VXF0,VYF0,VZF0,"wave_erzfun.in")

          ENDIF

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

          WRITE(99,*)PHI,ZIPH,ZPIPH,ZF,ZPF

        enddo !iphi=0,iphellip

        close(99)

c vertically

        if (phbetav.ne.0.0d0) then
          eps=phdispv**2/phbetav
          dy0=  sqrt(eps*phbetav)
          dyp0=-sqrt(eps/phbetav)
        else
          eps=0.0d0
          dy0=0.0d0
          dyp0=0.0d0
        endif

        DZ=PHTRZ0
        DZP=PHTRZP0
        ZIPH=DZ
        ZPIPH=DZP

        open(unit=99,file='wave_phase_ellipse_vert.wva',
     &    status='unknown',recl=256)
        write(99,'(a,i10,a,a)')'* ',icode,' ',code
        write(99,'(a,e15.5,a,e15.5)')'* z =',PHTRZ0,' zp =',PHTRZP0
        write(99,'(a)')'* Meaning of columns phi, yi, ypi, yf, ypf'

        do iphi=1,iphellip

          phi=(iphi-1)*dphi
          cs=cos(phi)
          sn=sin(phi)

          DZ=PHTRZ0
          DZP=PHTRZP0
          DY=  dy0*cs
          DYP=dyp0*sn

          YIPH=DY
          YPIPH=DYP

          XI=XI0+DY*EWY(1)+DZ*EWZ(1)
          YI=YI0+DY*EWY(2)+DZ*EWZ(2)
          ZI=ZI0+DY*EWY(3)+DZ*EWZ(3)

          WX=V0/SQRT(1.D0+DZP**2+DYP**2)
          WY=DYP*WX
          WZ=DZP*WX

          VXI=WX*EWS(1)+WY*EWY(1)+WZ*EWZ(1)
          VYI=WX*EWS(2)+WY*EWY(2)+WZ*EWZ(2)
          VZI=WX*EWS(3)+WY*EWY(3)+WZ*EWZ(3)

          IF (IERZFUN.EQ.0) THEN

            CALL TRACK(XI,YI,ZI,VXI,VYI,VZI
     &        ,XF0,YF0,ZF0,VXF0/V0,VYF0/V0,VZF0/V0
     &        ,XF,YF,ZF,VXF,VYF,VZF,DTIM,BSHIFT,DMYGAMMA,GAMMAL)

          ELSE !IERZFUN

            CALL IDTRMSHGF(XI,YI,ZI,VXI,VYI,VZI,XF,YF,ZF,VXF,VYF,VZF,
     &        XI0,YI0,ZI0,WXI0,WYI0,WZI0,
     &        XF0,YF0,ZF0,VXF0,VYF0,VZF0,"wave_erzfun.in")

          ENDIF

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

          WRITE(99,*)PHI,YIPH,YPIPH,YF,YPF

        enddo !iphi=0,iphellip

        close(99)

        write(lungfo,*)
        write(lungfo,*)'      Tralin: Phasespace ellipes written to files:'

        write(lungfo,*)'      wave_phase_ellipse_hori.wva'
        write(lungfo,*)'      wave_phase_ellipse_vert.wva'

        write(lungfo,*)

        write(lungfo,*)'      PHTRZ0,PHTRZP0,PHTRY0,PHTRYP0:'
        write(lungfo,*)'      ', PHTRZ0,PHTRZP0
        write(lungfo,*)'      ', PHTRY0,PHTRYP0

        write(lungfo,*)

        write(lungfo,*)'      PHTRZ0,PHTRZP0,PHTRY0,PHTRYP0:'
        write(lungfo,*)'      ', PHTRZ0,PHTRZP0
        write(lungfo,*)'      ', PHTRY0,PHTRYP0

        write(lungfo,*)

        write(lungfo,*)'      PHDISPH,PBETAH,PHDISPV,PBETAV:'
        write(lungfo,*)'      ',PHDISPH,PHBETAH
        write(lungfo,*)'      ',PHDISPV,PHBETAV

        write(lungfo,*)

      endif !(iphellip.ne.0.0d0)

      write(lungfo,'(/a)')"      --- Leaving TRALIN ---"
      write(lungfo,*)" "

      RETURN
      END
