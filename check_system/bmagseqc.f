*CMZ :  4.00/11 26/07/2021  09.08.58  by  Michael Scheer
*CMZ :  3.06/00 11/02/2019  12.49.34  by  Michael Scheer
*CMZ :  3.04/00 19/01/2018  16.33.13  by  Michael Scheer
*CMZ :  3.03/02 16/02/2016  12.18.47  by  Michael Scheer
*CMZ :  3.01/00 15/07/2013  08.04.32  by  Michael Scheer
*CMZ :  2.68/03 07/08/2012  13.09.30  by  Michael Scheer
*CMZ :  2.66/07 04/12/2009  16.11.19  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.02/00 19/12/97  16.15.53  by  Michael Scheer
*CMZ :  1.01/00 28/10/97  12.14.09  by  Michael Scheer
*CMZ : 00.01/08 01/04/95  16.54.24  by  Michael Scheer
*CMZ : 00.01/07 10/03/95  11.22.55  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.21.20  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.03  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.42  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BMAGSEQC(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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

      INTEGER IM,ISTORE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEEP,fourier.
      include 'fourier.cmn'
*KEND.

      DOUBLE PRECISION BX,BY,BZ,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION XIN,YIN,ZIN,xlen2,xr,yr,zr,ex,ey,ez,dist,fint,gap,
     &  bxr,byr,bzr,axr,ayr,azr,xbend,xshift,xsym,ybend,zbend,fringeout,
     &  strength,pin(3),center(3),pout(3),vnin(3),vnout(3),b(3),edge(2)

      integer ical,istatus,modus
      data ical/0/

      if (ical.eq.0) then

        bmsqbounds(1)=1.0d30
        bmsqbounds(2)=-1.0d30

        do im=1,mmag

          if (ctyp(im).eq.'QP') then
            xlen2=pmag(1,im)/2.0d0
            qfbounds(1,im)=pmag(3,im)-xlen2
            qfbounds(2,im)=pmag(3,im)+xlen2
            if (qfbounds(2,im).lt.qfbounds(1,im)) then
              write(lungfo,*)'*** Error in BMAGSEQC: Bounderies of QP ',im,'bad! ***'
              write(6,*)'*** Error in BMAGSEQC: Bounderies of QP ',im,'bad! ***'
              stop '*** Program WAVE aborted ***'
            endif
            if (qfbounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=qfbounds(1,im)
            if (qfbounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=qfbounds(2,im)
          else if (ctyp(im).eq.'QF') then
            xlen2=pmag(1,im)/2.0d0
            qfbounds(1,im)=pmag(3,im)-70.0d0/1000.0d0-xlen2
            qfbounds(2,im)=pmag(3,im)+70.0d0/1000.0d0+xlen2
            if (qfbounds(2,im).lt.qfbounds(1,im)) then
              write(lungfo,*)'*** Error in BMAGSEQC: Bounderies of QF ',im,'bad! ***'
              write(6,*)'*** Error in BMAGSEQC: Bounderies of QF ',im,'bad! ***'
              stop '*** Program WAVE aborted ***'
            endif
            if (qfbounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=qfbounds(1,im)
            if (qfbounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=qfbounds(2,im)
          else if (ctyp(im).eq.'SX') then
            xlen2=pmag(1,im)/2.0d0
            sxbounds(1,im)=pmag(3,im)-70.0d0/1000.0d0-xlen2
            sxbounds(2,im)=pmag(3,im)+70.0d0/1000.0d0+xlen2
            if (sxbounds(2,im).lt.sxbounds(1,im)) then
              write(lungfo,*)'*** Error in BMAGSEQC: Bounderies of SX ',im,'bad! ***'
              write(6,*)'*** Error in BMAGSEQC: Bounderies of SX ',im,'bad! ***'
              stop '*** Program WAVE aborted ***'
            endif
            if (sxbounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=sxbounds(1,im)
            if (sxbounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=sxbounds(2,im)
          else if (ctyp(im).eq.'DI'.or.ctyp(im).eq.'DIF') then
            xlen2=dabs(pmag(2,im)*sin(pmag(1,im)/2.0d0))
            dibounds(1,im)=pmag(3,im)-70.0d0/pmag(4,im)-xlen2
            dibounds(2,im)=pmag(3,im)+70.0d0/pmag(4,im)+xlen2
            if (dibounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=dibounds(1,im)
            if (dibounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=dibounds(2,im)
          else if (ctyp(im).eq.'DIL') then
            dibounds(1,im)=pmag(3,im)-pmag(11,im)/2.0d0
            dibounds(2,im)=pmag(3,im)+pmag(11,im)/2.0d0
            if (dibounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=dibounds(1,im)
            if (dibounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=dibounds(2,im)
          else if (ctyp(im).eq.'BEND') then
            if (dibounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=dibounds(1,im)
            if (dibounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=dibounds(2,im)
          else if (ctyp(im).eq.'DCS') then
            dibounds(1,im)=pmag(3,im)-pmag(11,im)/2.0d0
            dibounds(2,im)=pmag(3,im)+pmag(11,im)/2.0d0
            if (dibounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=dibounds(1,im)
            if (dibounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=dibounds(2,im)
          else if (ctyp(im).eq.'DQS') then
            dibounds(1,im)=pmag(3,im)-pmag(11,im)/2.0d0
            dibounds(2,im)=pmag(3,im)+pmag(11,im)/2.0d0
            if (dibounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=dibounds(1,im)
            if (dibounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=dibounds(2,im)
          else if (ctyp(im).eq.'DH'.or.ctyp(im).eq.'DHF') then
            xlen2=dabs(pmag(2,im)*sin(pmag(1,im)/2.0d0))
            dhbounds(1,im)=pmag(3,im)-70.0d0/pmag(4,im)-xlen2
            dhbounds(2,im)=pmag(3,im)+70.0d0/pmag(4,im)+xlen2
            if (dhbounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=dhbounds(1,im)
            if (dhbounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=dhbounds(2,im)
          else if (ctyp(im).eq.'UE') then
            if (uebounds(2,im).lt.uebounds(1,im)) then
              write(lungfo,*)'*** Error in BMAGSEQC: Bounderies of UE ',im,'bad! ***'
              write(6,*)'*** Error in BMAGSEQC: Bounderies of UE ',im,'bad! ***'
              stop '*** Program WAVE aborted ***'
            endif
            if (uebounds(1,im).lt.bmsqbounds(1)) bmsqbounds(1)=uebounds(1,im)
            if (uebounds(2,im).gt.bmsqbounds(2)) bmsqbounds(2)=uebounds(2,im)
          endif
        enddo

        ical=1

      endif

      BXOUT=0.0d0
      BYOUT=0.0d0
      BZOUT=0.0d0

      AXOUT=0.0d0
      AYOUT=0.0d0
      AZOUT=0.0d0

      BX=0.0d0
      BY=0.0d0
      BZ=0.0d0

C- FIELD

      IF(IWFILF.NE.99.AND.IMGSQF.NE.0) THEN
        ISTORE=IRFILF
        IRFILF=99
        CALL BFOUR(XIN,YIN,ZIN,BX,BY,BZ,AXOUT,AYOUT,AZOUT)
        IRFILF=ISTORE
      ENDIF !IWFILF

      DO IM=1,mmag

        IF (ctyp(IM).EQ.'QP') THEN
          CALL BQP(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,IM)
        ELSE IF (ctyp(IM).EQ.'QF') THEN
          if (xin+1.0d-10.ge.qfbounds(1,im).and.xin-1.0d-10.le.qfbounds(2,im)) then
            CALL BQF(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif
        ELSE IF (ctyp(IM).EQ.'SX') THEN
          if (xin+1.0d-10.ge.sxbounds(1,im).and.xin-1.0d-10.le.sxbounds(2,im)) then
            CALL bsx(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif

        ELSE IF (ctyp(IM).EQ.'BEND') THEN

          bxout=0.0d0
          byout=0.0d0
          bzout=0.0d0

          if (xin+1.0d-10.ge.dibounds(1,im).and.xin-1.0d-10.le.dibounds(2,im)) then

            pin=[pmag(1,im),0.0d0, pmag(2,im)]
            center=[pmag(3,im), 0.0d0, pmag(4,im)]
            pout=[pmag(5,im), 0.0d0,pmag(6,im)]

            strength=pmag(7,im)
            edge=pmag(8:9,im)
            fint=pmag(10,im)
            gap=pmag(11,im)

            b=[0.0d0,strength,0.0d0]

            vnin=[pmag(14,im), 0.0d0, pmag(15,im)]
            vnout=[pmag(16,im), 0.0d0, pmag(17,im)]

            modus=int(pmag(12,im))

            call bbend(xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &        fint,gap,center,b,
     &        pin,vnin,pout,vnout,
     &        modus,istatus,0)
          endif

        ELSE IF (ctyp(IM).EQ.'DI') THEN
          if (xin+1.0d-10.ge.dibounds(1,im).and.xin-1.0d-10.le.dibounds(2,im)) then
            CALL BDI(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif
        ELSE IF (ctyp(IM).EQ.'DCS') THEN
          bxout=0.0d0
          byout=0.0d0
          bzout=0.0d0
          if (xin+1.0d-10.ge.dibounds(1,im).and.xin-1.0d-10.le.dibounds(2,im)) then
            xshift=pmag(11,im)/2.0d0
            xr=xin-pmag(3,im)+xshift
            yr=yin-pmag(4,im)
            zr=zin-pmag(5,im)
            ex=pmag(8,im)
            ey=pmag(9,im)
            ez=pmag(10,im)
            dist=xr*ex+yr*ey+zr*ez
            if (dist.lt.0.0d0) return
            xbend=dist
            ybend=yr
            zbend=0.0d0 ! Or change mrad_csbend etc.
            xsym=xbend
            gap=pmag(7,im)
            fint=pmag(6,im)
            strength=pmag(12,im)
            if (xshift.gt.0.0d0) then
              if (xbend.gt.xshift) then
                xsym=xshift-(xbend-xshift)
              endif
              call mrad_fringe_cubic_spline(xsym,ybend,zbend,
     &          bxr,byr,bzr,axr,ayr,azr,fint,gap,fringeout,
     &          istatus)
            endif
            if (xbend.gt.xshift) then
              bxr=-bxr
            endif
            bxout=bxr*strength
            byout=byr*strength
            bzout=bzr*strength
          endif
        ELSE IF (ctyp(IM).EQ.'DIL') THEN
          bxout=0.0d0
          byout=0.0d0
          bzout=0.0d0
          if (xin+1.0d-10.ge.dibounds(1,im).and.xin-1.0d-10.le.dibounds(2,im)) then
            xshift=pmag(11,im)/2.0d0
            xr=xin-pmag(3,im)+xshift
            yr=yin-pmag(4,im)
            zr=zin-pmag(5,im)
            ex=pmag(8,im)
            ey=pmag(9,im)
            ez=pmag(10,im)
            dist=xr*ex+yr*ey+zr*ez
            if (dist.lt.0.0d0) return
            xbend=dist
            ybend=yr
            zbend=0.0d0 ! Or change mrad_csbend etc.
            xsym=xbend
            gap=pmag(7,im)
            fint=pmag(6,im)
            strength=pmag(12,im)
            if (xshift.gt.0.0d0) then
              if (xbend.gt.xshift) then
                xsym=xshift-(xbend-xshift)
              endif
              call mrad_fringe_linear(xsym,ybend,zbend,
     &          bxr,byr,bzr,axr,ayr,azr,fint,gap,fringeout,
     &          istatus)
            endif
            if (xbend.gt.xshift) then
              bxr=-bxr
            endif
            bxout=bxr*strength
            byout=byr*strength
            bzout=bzr*strength
          endif
        ELSE IF (ctyp(IM).EQ.'DQS') THEN
          bxout=0.0d0
          byout=0.0d0
          bzout=0.0d0
          if (xin+1.0d-10.ge.dibounds(1,im).and.xin-1.0d-10.le.dibounds(2,im)) then
            xshift=pmag(11,im)/2.0d0
            xr=xin-pmag(3,im)+xshift
            yr=yin-pmag(4,im)
            zr=zin-pmag(5,im)
            ex=pmag(8,im)
            ey=pmag(9,im)
            ez=pmag(10,im)
            dist=xr*ex+yr*ey+zr*ez
            if (dist.lt.0.0d0) return
            xbend=dist
            ybend=yr
            zbend=0.0d0 ! Or change mrad_csbend etc.
            xsym=xbend
            gap=pmag(7,im)
            fint=pmag(6,im)
            strength=pmag(12,im)
            if (xshift.gt.0.0d0) then
              if (xbend.gt.xshift) then
                xsym=xshift-(xbend-xshift)
              endif
              call mrad_fringe_quintic_spline(xsym,ybend,zbend,
     &          bxr,byr,bzr,axr,ayr,azr,fint,gap,fringeout,
     &          istatus)
            endif
            if (xbend.gt.xshift) then
              bxr=-bxr
            endif
            bxout=bxr*strength
            byout=byr*strength
            bzout=bzr*strength
          endif
        ELSE IF (ctyp(IM).EQ.'DIF') THEN
          if (xin+1.0d-10.ge.dibounds(1,im).and.xin-1.0d-10.le.dibounds(2,im)) then
            CALL bfoumag(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,axout,ayout,azout,IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif
        ELSE IF (ctyp(IM).EQ.'DHF') THEN
          if (xin+1.0d-10.ge.dhbounds(1,im).and.xin-1.0d-10.le.dhbounds(2,im)) then
            CALL bfoumag(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,axout,ayout,azout,-IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif
        ELSE IF (ctyp(IM).EQ.'DH') THEN
          if (xin+1.0d-10.ge.dhbounds(1,im).and.xin-1.0d-10.le.dhbounds(2,im)) then
            CALL BDH(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif
        ELSE IF (ctyp(IM).EQ.'UE') THEN
          if (xin+1.0d-10.ge.uebounds(1,im).and.xin-1.0d-10.le.uebounds(2,im)) then
            CALL BUE(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,axout,ayout,azout,IM)
          else
            bxout=0.0d0
            byout=0.0d0
            bzout=0.0d0
          endif
        ENDIF !CTYP

        BX=BX+BXOUT
        BY=BY+BYOUT
        BZ=BZ+BZOUT

      ENDDO !IM

      BXOUT=BX
      BYOUT=BY
      BZOUT=BZ

      RETURN
      END
