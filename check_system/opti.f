*CMZ :  4.00/14 22/12/2021  16.40.21  by  Michael Scheer
*CMZ :  4.00/11 22/11/2020  13.38.35  by  Michael Scheer
*CMZ :  4.00/07 09/06/2020  09.56.54  by  Michael Scheer
*CMZ :  3.04/00 19/01/2018  15.53.17  by  Michael Scheer
*CMZ :  3.03/04 19/12/2017  10.39.56  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  10.45.00  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.43.52  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.29.06  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/00 05/10/2012  08.15.04  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.61/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.61/00 29/01/2007  17.00.55  by  Michael Scheer
*CMZ :  2.54/07 16/06/2005  12.05.48  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.39/02 22/01/2002  11.29.06  by  Michael Scheer
*CMZ :  2.34/07 05/09/2001  12.06.28  by  Michael Scheer
*CMZ :  2.34/05 23/08/2001  17.35.48  by  Michael Scheer
*CMZ :  2.16/04 23/08/2001  15.27.45  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.01  by  Michael Scheer
*CMZ : 00.02/00 22/11/96  16.22.27  by  Michael Scheer
*CMZ : 00.01/10 02/09/96  16.02.38  by  Michael Scheer
*CMZ : 00.01/08 03/04/95  10.04.39  by  Michael Scheer
*CMZ : 00.01/07 23/03/95  12.23.36  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.02.55  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.48  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.48  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE OPTI(X0,Y0,Z0,BETX0,BETY0,BETZ0,
     &  XF0,YF0,ZF0,BETXF0,BETYF0,BETZF0,
     &  DTIM,BSHIFT,GAMMA)
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

C     TRACKS A SET OF PARTICLES AND WRITES INITIAL
C     AND FINAL COORDINATES TO FILE "FILEO"
C     THE FILE IS READ BY THE PROGRAM TRANPOLY TO CALCULATE THE
C     GENERATING FUNCTION OF THE DEVICE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,trackf90.
      include 'trackf90.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      INTEGER ICOUNT,ICOUNT10,KCOUNT,K10,IZFAIL,IYFAIL,IZ,IY,IZP,IYP,I,J,NTOT
      INTEGER IOPEND,IOPNF,JERZFUN,nkoef,luni,istat,ntread

*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,genfun.
      include 'genfun.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION OPTBUF(23,MAXTRA),BSHIFT,DTIM,GAMMA
     &  ,AX1,AY1,AZ1,AX2,AY2,AZ2,dgam
     &  ,AXIO,AYIO,AZIO,AXFO,AYFO,AZFO
     &  ,BX1,BY1,BZ1,BX2,BY2,BZ2
     &  ,BXIO,BYIO,BZIO,BXFO,BYFO,BZFO
     &  ,AXFIRST,AYFIRST,AZFIRST,BXFIRST,BYFIRST,BZFIRST
     &  ,XFIRST,YFIRST,ZFIRST,YPFIRST,ZPFIRST
     &  ,X1,Y1,Z1,ZP1,YP1,X2,Y2,Z2,ZP2,YP2
     &  ,XIO,YIO,ZIO,ZPIO,YPIO,XFO,YFO,ZFO,ZPFO,YPFO
     &  ,X0,Y0,Z0,ZP0,YP0,XF0,YF0,ZF0,ZPF0,YPF0,
     &  bxf0,byf0,bzf0,axf0,ayf0,azf0
     &  ,VX1,VY1,VZ1,VX2,VY2,VZ2,V0,VF
     &  ,VX0,VY0,VZ0,VXF0,VYF0,VZF0
     &  ,VXIO,VYIO,VZIO,VXFO,VYFO,VZFO
     &  ,BETX0,BETY0,BETZ0,BETXF0,BETYF0,BETZF0
     &  ,W0,WS1,WZ1,WY1,EWS(3),EWY(3),EWZ(3),S1,S1S,Y1S,Z1S
     &  ,EWSF(3),EWYF(3),EWZF(3),EN
     &  ,GAMMALOSS,phi,dum

      DOUBLE PRECISION BETA0H,betahh,ALPHAH,alphav,GAMMAH,EPSAH,EPSBH,EPSHZ,EPSHY
      DOUBLE PRECISION BETA0L,BETAL,ALPHAL,GAMMAL,EPSAL,EPSBL,EPSLZ,EPSLY
      DOUBLE PRECISION ZAPERTP,YAPERTP
      DOUBLE PRECISION XBMAXI,YBMAXI,ZBMAXI,VXBMAXI,VYBMAXI,VZBMAXI,BOLD
      DOUBLE PRECISION BSTORE,XSTORE,YSTORE,ZSTORE,VN

      DOUBLE PRECISION EPS
      DATA EPS/1.D-5/

      REAL*4 RNDM,rr

      INTEGER NTUP
      PARAMETER (NTUP=23)
      CHARACTER(4) CHTAGS(NTUP)
      REAL*8 TUP(NTUP)
      DATA CHTAGS/
     &  'xi',
     &  'zi','zpi',
     &  'yi','ypi',
     &  'bxi','byi','bzi',
     &  'axi','ayi','azi',
     &  'xf',
     &  'zf','zpf',
     &  'yf','ypf',
     &  'bxf','byf','bzf',
     &  'axf','ayf','azf',
     &  'dl'
     &  /

      INTEGER NTUP1
      PARAMETER (NTUP1=4)
      CHARACTER(4) CHTAGS1(NTUP1)
      DATA CHTAGS1/
     &  'zi','zpi',
     &  'yi','ypi'
     &  /

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SUBROUTINE OPTI'
      WRITE(LUNGFO,*)'     ==============='

      IF (IENELOSS.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    ' *** WARNING IN OPTI: IENELOSS .NE. 0, NOT YET FULLY IMPLEMENTED ***'
        WRITE(LUNGFO,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
     &    ' *** WARNING IN OPTI: IENELOSS .NE. 0, NOT YET FULLY IMPLEMENTED ***'
        WRITE(6,*)
     &    ' *** BE CAREFUL!! ***'
        WRITE(6,*)
      ENDIF

      if (ioptic.eq.-2) then
        optbuf=0.0d0
        icount=0
        open(newunit=luni,file="wave_opti_in.dat",status='old')
        call util_skip_comment(luni)
        do ntread=1,maxtra
          read(luni,*,end=123)i,
     &      x1,y1,z1,vx1,vy1,vz1,
     &      x2,y2,z2,vx2,vy2,vz2,
     &      istat
          if (istat.ne.0) cycle
          icount=icount+1
          optbuf(1,icount)=x1
          optbuf(2,icount)=z1
          optbuf(3,icount)=vz1/vx1
          optbuf(4,icount)=y1
          optbuf(5,icount)=vy1/vx1
          optbuf(12,icount)=x2
          optbuf(13,icount)=z2
          optbuf(14,icount)=vz2/vx2
          optbuf(15,icount)=y2
          optbuf(16,icount)=vy2/vx2
        enddo
123     close(luni)
      endif

      BSTORE=BMAXGL2
      XSTORE=XBMAXGL
      YSTORE=YBMAXGL
      ZSTORE=ZBMAXGL

      BMAXGL2=-1.D30
      BOLD=BMAXGL2

      IF (IHPHSPAC.NE.0) THEN
        CALL hbookm(NIDTRAC,'ZI,ZPI,YI,YPI ON GRID$',NTUP1,'//WAVE',
     &    (nzopt*2+1)*(nyopt*2+1)*(nzpopt*2+1)*(nypopt*2+1),
     &    CHTAGS1)
        CALL hbookm(NIDTRAC+1,'TRANS. TRACKS$',NTUP,'//WAVE',
     &    (nzopt*2+1)*(nyopt*2+1)*(nzpopt*2+1)*(nypopt*2+1),
     &    CHTAGS)
      ENDIF

C--- DEFAULTS

      IF (DLAPER.EQ.9999.) DLAPER=XF0-X0

      IF (DZOPT.EQ.9999.) DZOPT=ZAPERT/DFLOAT(NZOPT)
      IF (DZPOPT.EQ.9999.) DZPOPT=ZAPERT/DLAPER/DFLOAT(NZPOPT)
      IF (DYOPT.EQ.9999.) DYOPT=YAPERT/DFLOAT(NYOPT)
      IF (DYPOPT.EQ.9999.) DYPOPT=YAPERT/DLAPER/DFLOAT(NYPOPT)

      IF(KHALBA**2+IBHELM**2
     &  +KBFELD**2+IRFILB**2+IRFILP**2 .GT.1) STOP
     &  '*** ERROR IN OPTI: IERZFUN**2+IBHTRACK**2+... .GT.1 ***'

      JERZFUN=IERZFUN
      IF (IERZFUN.EQ.100) IERZFUN=1
      IF(IERZANA**2+IERZFUN**2+IBHTRACK**2+IBHARD**2 .GT.1) STOP
     &  '*** ERROR IN OPTI: IERZANA+IERZFUN+IBHTRACK+IBHARD .GT.1 ***'
      IERZFUN=JERZFUN

      ZP0=BETZ0/BETX0
      YP0=BETY0/BETX0
      V0=CLIGHT1*DSQRT( (1.D0-1.D0/GAMMA)*(1.D0+1.D0/GAMMA) )
      VX0=BETX0*CLIGHT1
      VY0=BETY0*CLIGHT1
      VZ0=BETZ0*CLIGHT1
      ZPF0=BETZF0/BETXF0
      YPF0=BETYF0/BETXF0
      VXF0=BETXF0*CLIGHT1
      VYF0=BETYF0*CLIGHT1
      VZF0=BETZF0*CLIGHT1
      NTOT=(NZOPT+1)*(NYOPT+1)*(NZPOPT+1)*(NYPOPT+1)

C--- FIND STARTING POINT

C- DEFAULTS

      IF (OPSTARTX.EQ.9999.) THEN
        OPSTARTX=X0
        OPSTARTY=Y0
        OPSTARTZ=Z0
      ENDIF

      IF (OPENDX.EQ.9999.)  THEN
        IOPEND=9999
        OPENDX=XF0
        OPENDY=YF0
        OPENDZ=ZF0
      ELSE
        IOPEND=0
      ENDIF

      IF(OPNX.EQ.9999.) THEN
        OPNX=VX0
        OPNY=VY0
        OPNZ=VZ0
      ENDIF

      EN=DSQRT(OPNX**2+OPNY**2+OPNZ**2)
      OPNX=OPNX/EN
      OPNY=OPNY/EN
      OPNZ=OPNZ/EN

      IF (OPNFX.EQ.9999) THEN
        IOPNF=9999
        OPNFX=VXF0
        OPNFY=VYF0
        OPNFZ=VZF0
      ELSE
        IOPNF=0
      ENDIF

      EN=DSQRT(OPNFX**2+OPNFY**2+OPNFZ**2)
      OPNFX=OPNFX/EN
      OPNFY=OPNFY/EN
      OPNFZ=OPNFZ/EN

      IF (OPENDX.LE.OPSTARTX) THEN
        WRITE(LUNGFO,*)'*** ERROR IN OPTI ***'
        WRITE(LUNGFO,*)'OPENDX .LE. OPSTARTX'
        WRITE(LUNGFO,*)'CHECK INPUT FILE'
        WRITE(6,*)'*** ERROR IN OPTI ***'
        WRITE(6,*)'OPENDX .LE. OPSTARTX'
        WRITE(6,*)'CHECK INPUT FILE'
        STOP
      ENDIF

      XIO=OPSTARTX
      YIO=OPSTARTY
      ZIO=OPSTARTZ
      VXIO=V0*OPNX
      VYIO=V0*OPNY
      VZIO=V0*OPNZ
      ZPIO=VZIO/VXIO
      YPIO=VYIO/VXIO

      CALL TRACK(XIO,YIO,ZIO,VXIO,VYIO,VZIO,
     &  OPENDX,OPENDY,OPENDZ,OPNFX,OPNFY,OPNFZ,
     &  XFO,YFO,ZFO,VXFO,VYFO,VZFO,DTIM,BSHIFT,GAMMA,GAMMALOSS)

      dgam=gammaloss

      IF (IOPEND.EQ.0) THEN
        XFO=OPENDX
        YFO=OPENDY
        ZFO=OPENDZ
      ENDIF   !(IOPEND.EQ.9999)

      IF (IOPNF.EQ.0) THEN
        VXFO=V0*OPNFX
        VYFO=V0*OPNFY
        VZFO=V0*OPNFZ
      ENDIF   !(IOPNF.EQ.9999)

      ZPFO=VZFO/VXFO
      YPFO=VYFO/VXFO

C- VECTOR-POTENTIALS

      CALL MYBFELD(XIO,YIO,ZIO,
     &  BXIO,BYIO,BZIO,AXIO,AYIO,AZIO)
c18.1.2019{
      CALL MYBFELD(XF0,YF0,ZF0,
     &  BXF0,BYF0,BZF0,AXF0,AYF0,AZF0)
c18.1.2019}
      CALL MYBFELD(XFO,YFO,ZFO,
     &  BXFO,BYFO,BZFO,AXFO,AYFO,AZFO)

C--- COORDINATE SYSTEM OF THE REFERENCE ORBIT (STARTING POINT)

      EWS(1)=VXIO/V0 !UNIT-VECTOR EWS
      EWS(2)=VYIO/V0
      EWS(3)=VZIO/V0

C        UNIT-VECTOR EWZ~[EWS,(0,1,0)] (CROSS-PRODUCT)

      EWZ(1)=-EWS(3)/DSQRT(EWS(3)*EWS(3)+EWS(1)*EWS(1))
      EWZ(2)= 0.
      EWZ(3)= EWS(1)/DSQRT(EWS(3)*EWS(3)+EWS(1)*EWS(1))

      EWY(1)= EWZ(2)*EWS(3) - EWZ(3)*EWS(2)  !UNIT-VECTOR EWY=[EWZ,EWS]
      EWY(2)= EWZ(3)*EWS(1) - EWZ(1)*EWS(3)
      EWY(3)= EWZ(1)*EWS(2) - EWZ(2)*EWS(1)


C--- CHECK ENTRANCE PLANE: ENTRANCE PLANE SHOULD BE PERPENDICULAR
C                          TO REFERENCE ORBIT

      IF (
     &    ABS(EWS(1)-OPNX).GT.1.D-10
     &    .OR.
     &    ABS(EWS(2)-OPNY).GT.1.D-10
     &    .OR.
     &    ABS(EWS(3)-OPNZ).GT.1.D-10
     &    ) THEN
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN OPTI ***'
        WRITE(6,*)'START VECTOR OF REFERENCE ORBIT AND NORMAL VECTOR'
        WRITE(6,*)'OF ENTRANCE PLANE OF DEVICE ARE DIFFERENT'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN OPTI ***'
        WRITE(LUNGFO,*)'START VECTOR OF REFERENCE ORBIT AND NORMAL VECTOR'
        WRITE(LUNGFO,*)'OF ENTRANCE PLANE OF DEVICE ARE DIFFERENT'
        WRITE(LUNGFO,*)
      ENDIF

C--- COORDINATE SYSTEM OF THE REFERENCE ORBIT (EXIT PLAIN)

      VF=DSQRT(VXFO*VXFO+VYFO*VYFO+VZFO*VZFO)
      EWSF(1)=VXFO/VF
      EWSF(2)=VYFO/VF
      EWSF(3)=VZFO/VF

      EWZF(1)=-EWSF(3)/DSQRT(EWSF(3)*EWSF(3)+EWSF(1)*EWSF(1))
      EWZF(2)= 0.
      EWZF(3)= EWSF(1)/DSQRT(EWSF(3)*EWSF(3)+EWSF(1)*EWSF(1))

      EWYF(1)= EWZF(2)*EWSF(3) - EWZF(3)*EWSF(2)
      EWYF(2)= EWZF(3)*EWSF(1) - EWZF(1)*EWSF(3)
      EWYF(3)= EWZF(1)*EWSF(2) - EWZF(2)*EWSF(1)

      IF (I2DIM.NE.0 .AND. (DABS( EWY(2)-1.D0).GT.1.D-15
     &                   .OR. DABS(EWYF(2)-1.D0).GT.1.D-15)) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN OPTI ***'
          WRITE(LUNGFO,*)
     &      'FLAG I2DIM SET, BUT REFERENCE ORBIT NOT IN X,Z - PLANE'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN OPTI ***'
          WRITE(6,*)
     &      'FLAG I2DIM SET, BUT REFERENCE ORBIT NOT IN X,Z - PLANE'
          STOP
      ENDIF

      IF (
     &  ABS(EWSF(1)-OPNFX).GT.1.D-10
     &  .OR.
     &  ABS(EWSF(2)-OPNFY).GT.1.D-10
     &  .OR.
     &  ABS(EWSF(3)-OPNFZ).GT.1.D-10
     &    ) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN OPTI ***'
          WRITE(LUNGFO,*)'START VECTOR OF REFERENCE ORBIT AND NORMAL VECTOR '
          WRITE(LUNGFO,*)'ENTRANCE PLANE OF DEVICE ARE DIFFERENT'
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN OPTI ***'
          WRITE(6,*)'START VECTOR OF REFERENCE ORBIT AND NORMAL VECTOR '
          WRITE(6,*)'ENTRANCE PLANE OF DEVICE ARE DIFFERENT'
      ENDIF

      if (ioptic.eq.-2) goto 1301

C- CHECK NUMBER OF TRAJECTORIES THAT PASS APERTURE
C  (ROUGH CHECK, TRAJECTORIES ASSUMED TO BE STRAIGHT LINES)

      ICOUNT=0
      IZFAIL=0
      IYFAIL=0

      DO IZ=-NZOPT,NZOPT,2
        DO IZP=-NZPOPT,NZPOPT,2
          DO IY=-NYOPT,NYOPT,2
            DO IYP=-NYPOPT,NYPOPT,2

              RR=real(IZ+IZP+IY+IYP)
              S1 =0.

              if (ioptic.gt.0) then

                Z1 =DZOPT *(DFLOAT (IZ) -DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))
                Y1 =DYOPT *(DFLOAT (IY) -DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))
                ZP1=DZPOPT*(DFLOAT (IZP)-DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))
                YP1=DYPOPT*(DFLOAT (IYP)-DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))

              else if (ioptic.eq.-1) then

c phasespace ellipse (alpha=0): gammah*z**2+2*alphah*z*zp+betah*zp**2=eps0h

                phi=DBLE(RNDM(RR))*twopi1
                alphah=-betaph/2.0d0
                z1= sqrt(eps0h*betah)*cos(phi)
                zp1=-sqrt(eps0h/betah)*(alphah*cos(phi)+sin(phi))

                phi=DBLE(RNDM(RR))*twopi1
                alphav=-betapv/2.0d0
                y1= sqrt(eps0v*betav)*cos(phi)
                yp1=-sqrt(eps0v/betav)*(alphav*cos(phi)+sin(phi))

              endif

              IF (I2DIM.NE.0) THEN
                Y1=0.
                YP1=0.
              ENDIF

              IF (IPHASPAC.NE.0) THEN

C--- HORIZONTAL HIGH BETA  (BETA0ZH)

                BETA0H=BETA0ZH
                betahh=BETA0H+DLAPER/BETA0H**2
                ALPHAH=DLAPER/BETA0H
                GAMMAH=(1.D0+ALPHAH**2)/betahh
                EPSBH=ZAPERT**2/betahh
                ZAPERTP=ZAPERT/2.D0/DLAPER
                EPSAH=(2.D0*ZAPERTP)**2/GAMMAH
                EPSHZ=DMIN1(EPSAH,EPSBH)

C--- HORIZONTAL LOW BETA  (BETA0ZL)

                BETA0L=BETA0ZL
                BETAL=BETA0L+DLAPER/BETA0L**2
                ALPHAL=DLAPER/BETA0L
                GAMMAL=(1.D0+ALPHAL**2)/BETAL
                EPSBL=ZAPERT**2/BETAL
                ZAPERTP=ZAPERT/2.D0/DLAPER
                EPSAL=(2.D0*ZAPERTP)**2/GAMMAL
                EPSLZ=DMIN1(EPSAL,EPSBL)

                IF(
     &              Z1**2*GAMMAH+2.D0*ALPHAH*Z1*ZP1+betahh*ZP1**2.GT.EPSHZ
     &              .OR.
     &              Z1**2*GAMMAL+2.D0*ALPHAL*Z1*ZP1+BETAL*ZP1**2.GT.EPSLZ
     &              ) THEN
                  IZFAIL=IZFAIL+1
                  GOTO 131 !I.E. SKIP
                ENDIF

C--- VERTICAL HIGH BETA  (BETA0YH)

                BETA0H=BETA0YH
                betahh=BETA0H+DLAPER/BETA0H**2
                ALPHAH=DLAPER/BETA0H
                GAMMAH=(1.D0+ALPHAH**2)/betahh
                EPSBH=YAPERT**2/betahh
                YAPERTP=YAPERT/2.D0/DLAPER
                EPSAH=(2.D0*YAPERTP)**2/GAMMAH
                EPSHY=DMIN1(EPSAH,EPSBH)

C--- VERTICAL LOW BETA  (BETA0YL)

                BETA0L=BETA0YL
                BETAL=BETA0L+DLAPER/BETA0L**2
                ALPHAL=DLAPER/BETA0L
                GAMMAL=(1.D0+ALPHAL**2)/BETAL
                EPSBL=YAPERT**2/BETAL
                YAPERTP=YAPERT/2.D0/DLAPER
                EPSAL=(2.D0*YAPERTP)**2/GAMMAL
                EPSLY=DMIN1(EPSAL,EPSBL)

                IF(
     &              Y1**2*GAMMAH+2.D0*ALPHAH*Y1*YP1+betahh*YP1**2.GT.EPSHY
     &              .OR.
     &              Y1**2*GAMMAL+2.D0*ALPHAL*Y1*YP1+BETAL*YP1**2.GT.EPSLY
     &              ) THEN
                  IYFAIL=IYFAIL+1
                  GOTO 131 !I.E. SKIP
                ENDIF

              ELSE   !IPHASPAC

                IF       (DABS(Z1)           .GT.DABS(ZAPERT+EPS)
     &              .OR.  DABS(Z1+ZP1*DLAPER).GT.DABS(ZAPERT+EPS)) THEN
                  IZFAIL=IZFAIL+1
                  GOTO 131 !I.E. SKIP
                ENDIF

                IF       (DABS(Y1)           .GT.DABS(YAPERT+EPS)
     &              .OR.  DABS(Y1+YP1*DLAPER).GT.DABS(YAPERT+EPS)) THEN
                  IYFAIL=IYFAIL+1
                  GOTO 131 !I.E. SKIP
                ENDIF

              ENDIF  !IPHASPAC

              ICOUNT=ICOUNT+1

              IF (icount.GT.MAXTRA) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** Error in OPTI ***'
                WRITE(LUNGFO,*)'TOO MANY TRACKS'
                WRITE(LUNGFO,*)'INCREASE PARAMETER MAXTRA IN FILE GENFUN.CMN'
                WRITE(LUNGFO,*)
                WRITE(6,*)
                WRITE(6,*)'*** ERROR  in OPTI ***'
                WRITE(6,*)'TOO MANY TRACKS'
                WRITE(6,*)'INCREASE PARAMETER MAXTRA IN FILE GENFUN.CMN'
                WRITE(6,*)
                STOP
              ENDIF

              if (ioptic.lt.0) then
                optbuf(2,icount)=z1
                optbuf(3,icount)=zp1
                optbuf(4,icount)=y1
                optbuf(5,icount)=yp1
              endif

131           continue
            ENDDO
          ENDDO
        ENDDO
      ENDDO

132   WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'*** MESSAGE SR OPTI ***'
      WRITE(6,*)
     &  'NUMBER OF PARTICLES TO BE TRACKED:',ICOUNT
      WRITE(6,*)
     &  'LOSSES DUE TO HORI. AND VERT. APERTURE:'
     &  ,IZFAIL,IYFAIL
      WRITE(6,*)
      WRITE(6,*)

      ICOUNT10=ICOUNT/10
      KCOUNT=1

      IF(ICOUNT.GT.MAXTRA) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** SR OPTI ***'
        WRITE(LUNGFO,*)'TOO MANY TRAJECTORIES'
        WRITE(LUNGFO,*)'INCREASE MAXTRA IN FILE GENFUN.CMN'
        WRITE(6,*)
        WRITE(6,*)'*** SR OPTI ***'
        WRITE(6,*)'TOO MANY TRAJECTORIES'
        WRITE(6,*)'INCREASE MAXTRA IN FILE GENFUN.CMN'
        STOP
      ENDIF

C---LOOP OVER RASTER

      ICOUNT=0
      IZFAIL=0
      IYFAIL=0
      DO 100 IZ=-NZOPT,NZOPT,2
        DO 110 IZP=-NZPOPT,NZPOPT,2
          DO 120 IY=-NYOPT,NYOPT,2
            DO 130 IYP=-NYPOPT,NYPOPT,2

C--- START OF TRAJEKTORY

              RR=real(IZ+IZP+IY+IYP)
              S1 =0.

              if (ioptic.ge.0) then

                Z1 =DZOPT *(DFLOAT(IZ) -DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))
                Y1 =DYOPT *(DFLOAT(IY) -DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))
                ZP1=DZPOPT*(DFLOAT(IZP)-DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))
                YP1=DYPOPT*(DFLOAT(IYP)-DRANDO+2.D0*DRANDO*DBLE(RNDM(RR)))

              else if (ioptic.eq.-1) then

                icount=icount+1
                z1=optbuf(2,icount)
                zp1=optbuf(3,icount)
                y1=optbuf(4,icount)
                yp1=optbuf(5,icount)
                icount=icount-1

              endif

              IF (I2DIM.NE.0) THEN
                Y1=0.
                YP1=0.
              ENDIF

              IF (IPHASPAC.NE.0) THEN

C--- HORIZONTAL HIGH BETA  (BETA0ZH)

                BETA0H=BETA0ZH
                betahh=BETA0H+DLAPER/BETA0H**2
                ALPHAH=DLAPER/BETA0H
                GAMMAH=(1.D0+ALPHAH**2)/betahh
                EPSBH=ZAPERT**2/betahh
                ZAPERTP=ZAPERT/2.D0/DLAPER
                EPSAH=(2.D0*ZAPERTP)**2/GAMMAH
                EPSHZ=DMIN1(EPSAH,EPSBH)

C--- HORIZONTAL LOW BETA  (BETA0ZL)

                BETA0L=BETA0ZL
                BETAL=BETA0L+DLAPER/BETA0L**2
                ALPHAL=DLAPER/BETA0L
                GAMMAL=(1.D0+ALPHAL**2)/BETAL
                EPSBL=ZAPERT**2/BETAL
                ZAPERTP=ZAPERT/2.D0/DLAPER
                EPSAL=(2.D0*ZAPERTP)**2/GAMMAL
                EPSLZ=DMIN1(EPSAL,EPSBL)

                IF(
     &              Z1**2*GAMMAH+2.D0*ALPHAH*Z1*ZP1+betahh*ZP1**2.GT.EPSHZ
     &              .OR.
     &              Z1**2*GAMMAL+2.D0*ALPHAL*Z1*ZP1+BETAL*ZP1**2.GT.EPSLZ
     &              ) THEN
                  IZFAIL=IZFAIL+1
                  GOTO 130 !I.E. SKIP
                ENDIF

C--- VERTICAL HIGH BETA  (BETA0YH)

                BETA0H=BETA0YH
                betahh=BETA0H+DLAPER/BETA0H**2
                ALPHAH=DLAPER/BETA0H
                GAMMAH=(1.D0+ALPHAH**2)/betahh
                EPSBH=YAPERT**2/betahh
                YAPERTP=YAPERT/2.D0/DLAPER
                EPSAH=(2.D0*YAPERTP)**2/GAMMAH
                EPSHY=DMIN1(EPSAH,EPSBH)

C--- VERTICAL LOW BETA  (BETA0YL)

                BETA0L=BETA0YL
                BETAL=BETA0L+DLAPER/BETA0L**2
                ALPHAL=DLAPER/BETA0L
                GAMMAL=(1.D0+ALPHAL**2)/BETAL
                EPSBL=YAPERT**2/BETAL
                YAPERTP=YAPERT/2.D0/DLAPER
                EPSAL=(2.D0*YAPERTP)**2/GAMMAL
                EPSLY=DMIN1(EPSAL,EPSBL)

                IF(
     &              Y1**2*GAMMAH+2.D0*ALPHAH*Y1*YP1+betahh*YP1**2.GT.EPSHY
     &              .OR.
     &              Y1**2*GAMMAL+2.D0*ALPHAL*Y1*YP1+BETAL*YP1**2.GT.EPSLY
     &              ) THEN
                  IYFAIL=IYFAIL+1
                  GOTO 130 !I.E. SKIP
                ENDIF

              ELSE   !IPHASPAC

                IF       (DABS(Z1)           .GT.DABS(ZAPERT+EPS)
     &              .OR.  DABS(Z1+ZP1*DLAPER).GT.DABS(ZAPERT+EPS)) THEN
                  IZFAIL=IZFAIL+1
                  GOTO 130 !I.E. SKIP
                ENDIF

                IF       (DABS(Y1)           .GT.DABS(YAPERT+EPS)
     &              .OR.  DABS(Y1+YP1*DLAPER).GT.DABS(YAPERT+EPS)) THEN
                  IYFAIL=IYFAIL+1
                  GOTO 130 !I.E. SKIP
                ENDIF

              ENDIF  !IPHASPAC

              IF (IHPHSPAC.NE.0) THEN

                TUP(1)=Z1
                TUP(2)=ZP1
                TUP(3)=Y1
                TUP(4)=YP1
                CALL hfm(NIDTRAC,TUP)

              ENDIF

              ICOUNT=ICOUNT+1

              IF(ICOUNT.EQ.KCOUNT) THEN

                CALL ZEIT(6)
                WRITE(6,*)'OPTI: NUMBER OF TRACKS CALCULATED SO FAR:',ICOUNT

                IF (KCOUNT.LT.ICOUNT10) THEN
                  KCOUNT=KCOUNT*10
                  K10=KCOUNT
                ELSE
                  KCOUNT=KCOUNT+K10
                ENDIF

              ENDIF

              IF(ICOUNT.GT.MAXTRA) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** SUBROUTINE OPTI ***'
                WRITE(LUNGFO,*)'TOO MANY TRAJECTORIES'
                WRITE(LUNGFO,*)'INCREASE MAXTRA IN FILE GENFUN.CMN'
                WRITE(6,*)
                WRITE(6,*)'*** SUBROUTINE OPTI ***'
                WRITE(6,*)'TOO MANY TRAJECTORIES'
                WRITE(6,*)'INCREASE MAXTRA IN FILE GENFUN.CMN'
                STOP
              ENDIF

              IF ((IERZFUN.EQ.0 .or. IERZFUN.EQ.100) .AND. IERZANA.EQ.0) THEN

C--- KINEMATICS IN THE TRAJEKTORY SYSTEM

                W0=V0

C        VELOCITY IN THE REFERENCE ORBIT SYSTEM

                WS1=W0/DSQRT(1.D0+ZP1*ZP1+YP1*YP1)
                WZ1=WS1*ZP1
                WY1=WS1*YP1

C--- TRANSFORM VELOCITY-VECTOR W INTO LAB.SYSTEM

                VX1=EWS(1)*WS1+EWY(1)*WY1+EWZ(1)*WZ1
                VY1=EWS(2)*WS1+EWY(2)*WY1+EWZ(2)*WZ1
                VZ1=EWS(3)*WS1+EWY(3)*WY1+EWZ(3)*WZ1

C--- REDEFINE ZP1,YP1

                ZP1=VZ1/VX1
                YP1=VY1/VX1

C--- TRANSFORM COORDINATES INTO LAB.SYSTEM

                S1S=S1
                Y1S=Y1
                Z1S=Z1

                X1=XIO+EWS(1)*S1S+EWY(1)*Y1S+EWZ(1)*Z1S
                Y1=YIO+EWS(2)*S1S+EWY(2)*Y1S+EWZ(2)*Z1S
                Z1=ZIO+EWS(3)*S1S+EWY(3)*Y1S+EWZ(3)*Z1S

C--- NOW WE ARE IN THE LAB

                XFIRST=X1
                YFIRST=Y1
                ZFIRST=Z1
                ZPFIRST=ZP1
                YPFIRST=YP1

C--- CALCULATE VECTOR POTENTIAL AT THE BEGINNING

                CALL MYBFELD(X1,Y1,Z1,BX1,BY1,BZ1,AX1,AY1,AZ1)

              ELSE

                CALL MYBFELD(OPENDX,Y1,Z1,BX1,BY1,BZ1,AX1,AY1,AZ1)

              ENDIF !(IERZFUN)

              AXFIRST=AX1
              AYFIRST=AY1
              AZFIRST=AZ1

              BXFIRST=BX1
              BYFIRST=BY1
              BZFIRST=BZ1

              IF (IERZFUN.NE.0.AND.IERZFUN.NE.100) THEN

                X1=OPSTARTX
                CALL ERZFUN(GAMMA,X1,BX1,BY1,BZ1,BX2,BY2,BZ2,
     &            AX1,AY1,AZ1,AX2,AY2,AZ2,Z1,ZP1,Y1,YP1,
     &            Z2,ZP2,Y2,YP2,OPSTARTX,OPENDX)
                X2=OPENDX

              ELSE IF (ABS(IERZFUN).EQ.100) THEN

                CALL IDTRMSHGF(
     &            X1,Y1,Z1,VX1,VY1,VZ1,
     &            X2,Y2,Z2,VX2,VY2,VZ2,
     &            X0,Y0,Z0,VX0,VY0,VZ0,
     &            XF0,YF0,ZF0,VXF0,VYF0,VZF0,"wave_erzfun.in")

                ZP2=VZ2/VX2
                YP2=VY2/VX2

                CALL MYBFELD(X2,Y2,Z2,BX2,BY2,BZ2,AX2,AY2,AZ2)

              ELSE IF (KBFELD.NE.0.AND.IBHARD.NE.0) THEN

C        CALL BHARD(X0,XSTOP,GAMMA,Z1,ZP1,Y1,YP1,Z2,ZP2,Y2,YP2)

                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** ERROR IN OPTI ***'
                WRITE(LUNGFO,*)'SR BHARD HERE NOT AVAILABLE'
                WRITE(LUNGFO,*)'CHECK INPUT FILE (FLAG IBHARD)'
                WRITE(6,*)
                WRITE(6,*)'*** ERROR IN OPTI ***'
                WRITE(6,*)'SR BHARD HERE NOT AVAILABLE'
                WRITE(6,*)'CHECK INPUT FILE (FLAG IBHARD)'
                STOP

              ELSE IF (IERZANA.NE.0) THEN

C START- AND END PLANE OF MAGNET ARE OVERWRITTEN IN SR ERZANA

                CALL ERZANA(OPSTARTX,OPENDX,
     &            OPNX,OPNY,OPNZ,
     &            OPNFX,OPNFY,OPNFZ,
     &            Z1,ZP1,Y1,YP1,Z2,ZP2,Y2,YP2,
     &            BX1,BY1,BZ1,BX2,BY2,BZ2,
     &            AX1,AY1,AZ1,AX2,AY2,AZ2)

              ELSE IF (KBFELD.NE.0.AND.IBHTRACK.NE.0) THEN

C CHECK OPSTART, OPEND
C     CALL BHTRACK (GAMMA,X1,Y1,Z1,V0,VX1,VY1,VZ1,
C     &                       X2,Y2,Z2,ZP2,YP2)

                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** ERROR IN OPTI ***'
                WRITE(LUNGFO,*)'SR BHTRACK HERE NOT AVAILABLE'
                WRITE(LUNGFO,*)'CHECK INPUT FILE (FLAG IBHTRACK)'
                WRITE(6,*)
                WRITE(6,*)'*** ERROR IN OPTI ***'
                WRITE(6,*)'SR BHTRACK HERE NOT AVAILABLE'
                WRITE(6,*)'CHECK INPUT FILE (FLAG IBHTRACK)'
                STOP

              ELSE !IERZFUN

C NORMAL CASE

                CALL TRACK(X1,Y1,Z1,VX1,VY1,VZ1,
     &            XFO,YFO,ZFO,OPNFX,OPNFY,OPNFZ,
     &            X2,Y2,Z2,VX2,VY2,VZ2,DTIM,BSHIFT,GAMMA,GAMMALOSS)

                if (iwarnmyb.ne.0) then
                  icount=icount-1
                  goto 130
                endif

                IF (BMAXGL2.GT.BOLD) THEN
                  BOLD=BMAXGL2
                  XBMAXI=X1
                  YBMAXI=Y1
                  ZBMAXI=Z1
                  VXBMAXI=VX1
                  VYBMAXI=VY1
                  VZBMAXI=VZ1
                ENDIF

                ZP2=VZ2/VX2
                YP2=VY2/VX2

                CALL MYBFELD(X2,Y2,Z2,BX2,BY2,BZ2,AX2,AY2,AZ2)

              ENDIF  !IERZFUN

              IF (IERZFUN.EQ.0 .AND. IERZANA.EQ.0) THEN

                OPTBUF(1,ICOUNT)=XFIRST
                OPTBUF(2,ICOUNT)=ZFIRST
                OPTBUF(3,ICOUNT)=ZPFIRST
                OPTBUF(4,ICOUNT)=YFIRST
                OPTBUF(5,ICOUNT)=YPFIRST
                OPTBUF(6,ICOUNT)=BXFIRST
                OPTBUF(7,ICOUNT)=BYFIRST
                OPTBUF(8,ICOUNT)=BZFIRST
                OPTBUF(9,ICOUNT)=AXFIRST
                OPTBUF(10,ICOUNT)=AYFIRST
                OPTBUF(11,ICOUNT)=AZFIRST

              ELSE   !IERZFUN

                OPTBUF(1,ICOUNT)=X1
                OPTBUF(2,ICOUNT)=Z1
                OPTBUF(3,ICOUNT)=ZP1
                OPTBUF(4,ICOUNT)=Y1
                OPTBUF(5,ICOUNT)=YP1
                OPTBUF(6,ICOUNT)=BX1
                OPTBUF(7,ICOUNT)=BY1
                OPTBUF(8,ICOUNT)=BZ1
                OPTBUF(9,ICOUNT)=AX1
                OPTBUF(10,ICOUNT)=AY1
                OPTBUF(11,ICOUNT)=AZ1

              ENDIF !IERZFUN

              OPTBUF(12,ICOUNT)=X2
              OPTBUF(13,ICOUNT)=Z2
              OPTBUF(14,ICOUNT)=ZP2
              OPTBUF(15,ICOUNT)=Y2
              OPTBUF(16,ICOUNT)=YP2
              OPTBUF(17,ICOUNT)=BX2
              OPTBUF(18,ICOUNT)=BY2
              OPTBUF(19,ICOUNT)=BZ2
              OPTBUF(20,ICOUNT)=AX2
              OPTBUF(21,ICOUNT)=AY2
              OPTBUF(22,ICOUNT)=AZ2
              OPTBUF(23,ICOUNT)=WTRA2IC
c              OPTBUF(24,ICOUNT)=tint

130         CONTINUE
120       CONTINUE
110     CONTINUE
100   CONTINUE

1301  NTOT=ICOUNT

      OPEN(UNIT=LUNO,FILE=FILEO,STATUS='unknown',recl=1024) !,FORM='UNFORMATTED')

      write(luno,*)CODE
      write(luno,*)ICODE
      write(luno,*)GAMMA,dgam
      write(luno,*)NTOT

      IF (IERZFUN.NE.0) THEN

          XIO=OPSTARTX
          YIO=0.
          ZIO=0.
          ZPIO=0.
          YPIO=0.
          BXIO=0.
          BYIO=0.
          BZIO=0.
          AXIO=0.
          AYIO=0.
          AZIO=0.

          XFO=OPENDX
          YFO=0.
          ZFO=0.
          ZPFO=0.
          YPFO=0.
          BXFO=0.
          BYFO=0.
          BZFO=0.
          AXFO=0.
          AYFO=0.
          AZFO=0.

          OPNX=1.
          OPNY=0.
          OPNZ=0.

          OPNFX=1.
          OPNFY=0.
          OPNFZ=0.

      ENDIF

      IF (IERZANA.NE.0) THEN

          XIO=OPSTARTX
          YIO=0.
          ZIO=0.
          ZPIO=0.
          YPIO=0.
          XFO=OPENDX

          CALL ERZANA(OPSTARTX,OPENDX,OPNX,OPNY,OPNZ, !REFERENCE ORBIT
     &                  OPNFX,OPNFY,OPNFZ,
     &                  ZIO,ZPIO,YIO,YPIO,ZFO,ZPFO,YFO,YPFO,
     &                  BXIO,BYIO,BZIO,BXFO,BYFO,BZFO,
     &                  AXIO,AYIO,AZIO,AXFO,AYFO,AZFO)

      ENDIF

      write(luno,*)XIO,YIO,ZIO,ZPIO,YPIO,BXIO,BYIO,BZIO,AXIO,AYIO,AZIO
      write(luno,*)XFO,YFO,ZFO,ZPFO,YPFO,BXFO,BYFO,BZFO,AXFO,AYFO,AZFO
      write(luno,*)xf0,yf0,zf0,vzf0/vxf0,vyf0/vxf0,axf0,ayf0,azf0 !18.1.2018
      write(luno,*)OPNX,OPNY,OPNZ
      write(luno,*)OPNFX,OPNFY,OPNFZ
      write(luno,*)ZAPERT,YAPERT,DLAPER,bint0

      DO I=1,NTOT
        write(luno,*) (OPTBUF(J,I),J=1,23)
c        write(luno,*) (OPTBUF(J,I),J=1,11)
c        write(luno,*) (OPTBUF(J,I),J=12,23)
        IF (IHPHSPAC.NE.0) THEN
          DO J=1,23
            TUP(J)=OPTBUF(J,I)
          ENDDO
          CALL hfm(NIDTRAC+1,TUP)
        ENDIF
      ENDDO

      CLOSE(LUNO)

C     WRITE(6,*)'*** SR OPTI:',ICOUNT,' ANFANGS- UND ENDZUSTAENDE AUF',
C     &            FILEO,' GESCHRIEBEN ***'
C     WRITE(6,*)
C     WRITE(6,*) 'VERLUSTE DURCH HORZ. UND VERT. APERTUR:',IZFAIL,IYFAIL
C     WRITE(6,*)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     Number of tracks written to file:',ICOUNT
      WRITE(LUNGFO,*)'     Filename: ',FILEO

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     Reference point of entrance plane (OPSTARTX,OPSTARTY,OPSTARTZ):'
      WRITE(LUNGFO,*)'     ',SNGL(OPSTARTX),SNGL(OPSTARTY),SNGL(OPSTARTZ)
      WRITE(LUNGFO,*)
     &  '     Normal vector of entrance plane (OPNX,OPNY,OPNZ):'
      WRITE(LUNGFO,*)'     ',SNGL(OPNX),SNGL(OPNY),SNGL(OPNZ)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     Reference point of exit plane (OPENDX,OPENDY,OPENDZ):'
      WRITE(LUNGFO,*)'     ',SNGL(OPENDX),SNGL(OPENDY),SNGL(OPENDZ)
      WRITE(LUNGFO,*)
     &  '     Normal vector of exit plane (OPNFX,OPNFY,OPNFZ):'
      WRITE(LUNGFO,*)'     ',SNGL(OPNFX),SNGL(OPNFY),SNGL(OPNFZ)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &'     Parameter NZOPT, NYOPT, NZPOPT, NYPOPT to define number of grid points:'
      WRITE(LUNGFO,*)'     ',NZOPT, NYOPT, NZPOPT, NYPOPT
      WRITE(LUNGFO,*)
     &'     Mash sizes 2.*DZOPT, 2.*DYOPT, 2.*DZPOPT, 2.*DYPOPT of grid:'
      WRITE(LUNGFO,*)'     ',SNGL(2.*DZOPT), SNGL(2.*DYOPT),
     &                         SNGL(2.*DZPOPT),SNGL(2.*DYPOPT)

      WRITE(LUNGFO,*)
      IF (IPHASPAC.NE.0) THEN
      WRITE(LUNGFO,*)
     &'     Horizontal low and high beta functions:',SNGL(BETA0ZL),SNGL(BETA0ZH)
      WRITE(LUNGFO,*)
     &'     Vertical low and high beta functions:',SNGL(BETA0YL),SNGL(BETA0YH)
      WRITE(LUNGFO,*)
      ENDIF
      WRITE(LUNGFO,*)
     &'     Losses due to hori. and vert. aperture cut:',IZFAIL,IYFAIL
      WRITE(LUNGFO,*)
     &'     Hori. and vert. aperture and collimator length (ZAPERT,YAPERT,DLAPER):'
      WRITE(LUNGFO,*)'     ',SNGL(ZAPERT),SNGL(YAPERT),SNGL(DLAPER)
      WRITE(LUNGFO,*)

      IF (IERZANA.NE.0.AND.IERZFUN.EQ.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     BMAXGL found while in SR OPTI:'
     &    ,DSQRT(DABS(BOLD))
        WRITE(LUNGFO,*)
     &    'start of corresponding track:'
        WRITE(LUNGFO,*)'     XBMAXI, YBMAXI, ZBMAXI:'
        WRITE(LUNGFO,*)'     ',XBMAXI, YBMAXI, ZBMAXI
        WRITE(LUNGFO,*)'     VXBMAXI, VYBMAXI, VZBMAXI:'
        VN=DSQRT(VXBMAXI**2+VYBMAXI**2+VZBMAXI**2)
        WRITE(LUNGFO,*)'     ',VXBMAXI/VN, VYBMAXI/VN, VZBMAXI/VN
        WRITE(LUNGFO,*)

        BMAXGL2=BSTORE
        XBMAXGL=XSTORE
        YBMAXGL=YSTORE
        ZBMAXGL=ZSTORE

      ENDIF

      RETURN
      END
