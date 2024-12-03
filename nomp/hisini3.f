*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  3.03/00 18/08/2015  15.44.00  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.13.36  by  Michael Scheer
*CMZ :  2.68/05 17/09/2012  12.42.05  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HISINI3
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
*KEND.

      use bunchmod

C--- INITIALIZES HBOOK

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,photon.
      include 'photon.cmn'
*KEEP,ntuple2.
      include 'ntuple2.cmn'
*KEND.

      real zmin,zmax,ymin,ymax
      CHARACTER(4) CHTAGS(22)
      CHARACTER(96) tit
      CHARACTER(4) CHSPEC(21)
      CHARACTER(4) CHPOW(8)
      CHARACTER(4) CHSTOK(12)

      data chspec/'isou','iobs','x','y','z','ener','spec'
     &  ,'iz','iy','iene'
     &  ,'re_x','im_x','re_y','im_y','re_z','im_z','rf_y','if_y','rf_z','if_z',
     &  'phi0'
     &  /

      data chpow/'x','y','z','pow','iz','iy','iobs','isou'/

      data chstok/'iobs','x','y','z','ener','s0','s1','s2','s3'
     &  ,'iz','iy','iene'
     &  /

      if (ispec.ne.0.and.ipin.eq.3) then

        if (ihisini_c.eq.0) call hisini

        write(6,*)' '
        write(6,*)'      HISINI3:'
        write(6,*)'      '
        write(6,*)'      For IPIN=3, not all histograms and ntuples are available'
        write(6,*)'      '

        write(lungfo,*)' '
        write(lungfo,*)'      HISINI3:'
        write(lungfo,*)'      '
        write(lungfo,*)'      For IPIN=3, not all histograms and ntuples are available'
        write(lungfo,*)'      '

        IF (IPINCIRC.NE.0) THEN
          PINW=2.D0*PINR
          PINH=2.D0*PINR
        ENDIF !IPINCIRC

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

        ymin=pincen(2)-pinh/2.
        ymax=pincen(2)+pinh/2.
        zmin=pincen(3)-pinw/2.
        zmax=pincen(3)+pinw/2.
        tit='dist. in pinhole (ipin=3)'
        call hbook2m(idspec,tit,mpinzorig,zmin,zmax,mpinyorig,ymin,ymax,0.0)
        tit='hori. dist. in pinhole (ipin=3)'
        call hbook1m(idspec-1,tit,mpinzorig,zmin,zmax,0.0)
        tit='vert. dist. in pinhole (ipin=3)'
        call hbook1m(idspec-2,tit,mpinyorig,ymin,ymax,0.0)
        call hbookm(nidspec,'arrays spect, reaima (ipin=3)',21
     &    ,'//WAVE',nsource*nobsv*nfreq,chspec)
        call hbookm(nidpow,'power density from spectrum (ipin=3)',8
     &    ,'//WAVE',nsource*nobsv,chpow)
        chtags(1)='egam'
        chtags(2)='flux'
        call hbookm(nidfreqp,'photon flux in pinhole (ipin=3)',2,
     &    '//WAVE',nfreq,chtags)
        if (istokes.ne.0) then
          call hbookm(nidstok,'stokes arrays (ipin=3)',12
     &      ,'//WAVE',1000,chstok)
          chtags(1)='egam'
          chtags(2)='s0'
          chtags(3)='s1'
          chtags(4)='s2'
          chtags(5)='s3'
          call hbookm(4600,'stokes flux in pinhole (ipin=3)',5,
     &      '//WAVE',nfreq,chtags)
        endif
      endif !ipin.eq.3

      RETURN
      END
