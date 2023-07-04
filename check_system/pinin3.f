*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.68/05 17/09/2012  12.26.08  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  15.19.07  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PININ3
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- INITIALIZE GRID OF OBERSERVATION POINTS OF PINHOLE for ipin.eq.3

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      if (ipin.ne.3) stop '*** ERROR IN PININ3: IPIN.NE.3'

      IF (IPINCIRC.NE.0) THEN
        PINW=2.D0*PINR
        PINH=2.D0*PINR
      ENDIF !IPINCIRC

      IF (OBSVDZ.NE.0.D0) THEN
        MPINZ=NINT(PINW/OBSVDZ)+1
        PINW=OBSVDZ*(MPINZ-1)
      ENDIF !(OBSVDZ.NE.0.D0)

      IF (OBSVDY.NE.0.D0) THEN
        MPINY=NINT(PINH/OBSVDY)+1
        PINH=OBSVDY*(MPINY-1)
      ENDIF !(OBSVDZ.NE.0.D0)

      mpinzorig=mpinz
      mpinyorig=mpiny

      IF (PINCEN(2).EQ.9999.) THEN
        IF     (IPBRILL.EQ.0) THEN
          PINCEN(2)=0.0
        ELSE IF (IPBRILL.EQ.1) THEN
          PINCEN(2)=PINH/2.D0
        ELSE IF (IPBRILL.EQ.2) THEN
          PINCEN(2)=PINH/2.D0
        ELSE IF (IPBRILL.EQ.3) THEN
          PINCEN(2)=-PINH/2.D0
        ELSE IF (IPBRILL.EQ.4) THEN
          PINCEN(2)=-PINH/2.D0
        ENDIF !IPBRILL
      ELSE IF (PINCEN(2).EQ.-9999.) THEN
        PINCEN(2)=YSTART+VYIN/VXIN*(PINCEN(1)-XSTART)
      ENDIF !PINCEN(2)

      IF (PINCEN(3).EQ.9999.) THEN
        IF     (IPBRILL.EQ.0) THEN
          PINCEN(3)=0.0
        ELSE IF (IPBRILL.EQ.1) THEN
          PINCEN(3)=PINW/2.D0
        ELSE IF (IPBRILL.EQ.2) THEN
          PINCEN(3)=-PINW/2.D0
        ELSE IF (IPBRILL.EQ.3) THEN
          PINCEN(3)=-PINW/2.D0
        ELSE IF (IPBRILL.EQ.4) THEN
          PINCEN(3)=PINW/2.D0
        ENDIF !IPBRILL
      ELSE IF (PINCEN(3).EQ.-9999.) THEN
        PINCEN(3)=ZSTART+VZIN/VXIN*(PINCEN(1)-XSTART)
      ENDIF !PINCEN(3)

      IF(NDOBSVZ*NDOBSVY.NE.NDOBSV) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*) '*** ERROR IN PININ ***'
        WRITE(LUNGFO,*) 'DIMENSION DECLARATIONS NOT CONSISTENT'
        WRITE(LUNGFO,*) 'NDOBSV MUST BE EQUAL TO NDOBSVZ*NDOBSVY'
        WRITE(LUNGFO,*) 'CHANGE PARAMETER IN CMPARA.CMN'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
        WRITE(6,*) '*** ERROR IN PININ ***'
        WRITE(6,*) 'DIMENSION DECLARATIONS NOT CONSISTENT'
        WRITE(6,*) 'NDOBSV MUST BE EQUAL TO NDOBSVZ*NDOBSVY'
        WRITE(6,*) 'CHANGE PARAMETER IN CMPARA.CMN'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF

      IF (IF1DIM.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININ3 ***'
        WRITE(LUNGFO,*)'FLAG IF1DIM SET BUT IPIN =3'
        WRITE(LUNGFO,*)'IF1DIM SET TO ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN PININ3 ***'
        WRITE(6,*)'FLAG IF1DIM SET BUT IPIN =3'
        WRITE(6,*)'IF1DIM SET TO ZERO'
        WRITE(6,*)
        IF1DIM=0
      ENDIF

      IF (mpinr.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININ3 ***'
        WRITE(LUNGFO,*)'FLAG MPINR SET BUT IPIN=3'
        WRITE(LUNGFO,*)'MPINR SET TO ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN PININ3 ***'
        WRITE(6,*)'FLAG MPINR SET BUT IPIN=3'
        WRITE(6,*)'MPINR SET TO ZERO'
        WRITE(6,*)
        MPINR=0
      ENDIF

C- DETERMINE RMS VALUES FOR FOLDING

      if (iemit.ne.0) CALL WSIGFOL

      IF (IFOLD.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININ3 ***'
        WRITE(LUNGFO,*)'FLAG IFOLD SET BUT IPIN=3'
        WRITE(LUNGFO,*)'IFOLD SET TO ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN PININ3 ***'
        WRITE(6,*)'FLAG IFOLD SET BUT IPIN=3'
        WRITE(6,*)'IFOLD SET TO ZERO'
        WRITE(6,*)
      ENDIF

      IF (IEFOLD.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININ3 ***'
        WRITE(LUNGFO,*)'FLAG IEFOLD SET BUT IPIN=3'
        WRITE(LUNGFO,*)'IEFOLD SET TO ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN PININ3 ***'
        WRITE(6,*)'FLAG IEFOLD SET BUT IPIN=3'
        WRITE(6,*)'IEFOLD SET TO ZERO'
        WRITE(6,*)
      ENDIF

      ifold=0
      iefold=0

      mpiny=1
      mpinz=1

      ihpin=0
      ihfold=0

      nobsv=1
      mobsv=1
      nobsvz=1
      nobsvy=1
      mobsvy=1
      mobsvz=1

      ipbrill=1
      icbrill=1

      obsv(1,1)=obs1x
      obsv(2,1)=obs1y
      obsv(3,1)=obs1z

      obsvy(1)=obs1y
      obsvz(1)=obs1z

      if (ihisini_c.ne.0) call hisini3

      RETURN
      END
