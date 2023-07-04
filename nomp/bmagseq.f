*CMZ :  4.00/11 26/07/2021  08.38.41  by  Michael Scheer
*CMZ :  4.00/07 09/07/2020  12.27.02  by  Michael Scheer
*CMZ :  3.06/00 11/02/2019  12.49.00  by  Michael Scheer
*CMZ :  3.04/00 16/01/2018  17.30.02  by  Michael Scheer
*CMZ :  3.03/04 11/10/2017  11.26.09  by  Michael Scheer
*CMZ :  3.03/02 18/03/2016  14.37.19  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  08.47.45  by  Michael Scheer
*CMZ :  3.01/01 29/07/2013  14.03.24  by  Michael Scheer
*CMZ :  3.01/00 17/07/2013  16.10.24  by  Michael Scheer
*CMZ :  3.00/01 28/03/2013  12.44.42  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/03 07/08/2012  11.36.24  by  Michael Scheer
*CMZ :  2.68/02 06/07/2012  13.23.04  by  Michael Scheer
*CMZ :  2.67/06 24/05/2012  12.31.33  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  09.41.56  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.61/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.53/01 24/01/2005  10.56.03  by  Michael Scheer
*CMZ :  2.52/11 06/12/2004  15.54.00  by  Michael Scheer
*CMZ :  2.47/10 30/05/2003  11.44.20  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.20/03 23/02/2001  11.01.50  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  1.02/03 07/01/98  11.52.22  by  Michael Scheer
*CMZ :  1.02/00 06/01/98  15.08.07  by  Michael Scheer
*CMZ :  1.01/00 28/10/97  12.14.09  by  Michael Scheer
*CMZ : 00.01/08 01/04/95  16.54.24  by  Michael Scheer
*CMZ : 00.01/07 10/03/95  11.22.55  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.21.20  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.03  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.42  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BMAGSEQ(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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

C--- READ DATA-FILE BMAGSEQ.IN, THAT CONTAINS MAGNET CONFIGURATION
C    STRUCTURE IS CENTERED AROUND ORIGIN

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL,IM,ieof,imag

*KEEP,fourier.
      include 'fourier.cmn'
*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEND.

      DOUBLE PRECISION BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION XIN,YIN,ZIN
      double precision shift,xcen,perlen,pern,xlamb,ahwpol,totlen

      DOUBLE PRECISION VN,BETA,V0,X1,Y1,Z1,X2,Y2,Z2
     &  ,VX1,VY1,VZ1,VX2,VY2,VZ2,ANG1z,ANG2z,DANGz,ang1y,ang2y,dangy
     &  ,DTIM,BSHIFT,xlen2,dint,bx,by,bz,
     &  xfour(nfoumagcp+2+2),dxfour,
     &  posi(4,3),edge(2),strength,angle,dlength,seclen,
     &  fint,hgap,de,ds,dum

      COMPLEX CKOEF(nfoumagcp/2+1+2)
      REAL*4  YFOUR(nfoumagcp+2+2)
      EQUIVALENCE (CKOEF,YFOUR)

      integer ip,i1,i,mfour,k,istatus

      character(32) cbmodel
      CHARACTER(3) CDUM2
      CHARACTER(5) CDUM1
      CHARACTER(256) cnam(nmgsqp)

      save ical, cnam

      DATA ICAL/0/

C--- INITIALISATION

      IF (ICAL.EQ.0) THEN

C- OPEN FILE, READ FIRST TIME IN ORDER TO DECODE MAGNET-TYPES

        OPEN(UNIT=LUNMG,FILE=FILEMG,FORM='FORMATTED',STATUS='OLD')

        DO IM=1,NMGSQP
          call util_skip_comment_end(lunmg,ieof)
          if (ieof.ne.0) goto 99
          READ(LUNMG,*,END=99) cnam(im),CTYP(IM)
          if (ctyp(im).eq.'BEND') then
            READ(LUNMG,*,END=99) DUM
          endif
        ENDDO !IM

99      CONTINUE
        mmag=IM-1

C- REWIND FILE AND READ AGAIN TO GET PARAMETERS

        REWIND(LUNMG)

        nfoumags=0

        DO IM=1,mmag

          call util_skip_comment_end(lunmg,ieof)
          CORR(IM)=1.0D0

          IF     (CTYP(IM).EQ.'DI') THEN
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1,IM),PMAG(2,IM),PMAG(3,IM),PMAG(4,IM)

          else IF (CTYP(IM).EQ.'BEND') THEN

            ! 1-2: x,z of point in entrance plane
            ! 3-4: x,z of point in orbit plane
            ! 5-6: x,z of point in exit plane
            ! 7-9: field strength, angles of entrance and exit planes
            ! 10-11: fint, gap
            ! 12: mode

            read(lunmg,*)cdum1,cdum2,pmag(7:12,im)
            read(lunmg,*)pmag(1:6,im)

            if (pmag(12,im).ne.3.0d0.and.pmag(12,im).ne.5.0d0)
     &        pmag(10:11,im)=0.0d0

c            if (strength.ne.0.0d) then
c              r=dbrho/strength
c            else
c              r=1.0d30
c            endif

            dibounds(1,im)=pmag(1,im)-abs(sin(pmag(8,im)))
            dibounds(2,im)=pmag(5,im)+abs(sin(pmag(9,im)))

            pmag(14,im)=cos(pmag(8,im))
            pmag(15,im)=sin(pmag(8,im))
            pmag(16,im)=cos(pmag(9,im))
            pmag(17,im)=sin(pmag(9,im))

          else IF (CTYP(IM).EQ.'DIL') THEN
            READ(LUNMG,*)CDUM1,CDUM2,PMAG(1:10,IM)
            ds=sqrt(pmag(8,im)**2+pmag(9,im)**2+pmag(10,im)**2)
            if (ds.eq.0.0d0) then
              print*," "
              print*,"*** Error in bmagseq: Bad normal vector for entrance plane of DCS element ***"
              print*," "
              stop
            else
              pmag(8:10,im)=pmag(8:10,im)/ds
            endif
            ds=1.0D0/dble(myinum)
            cbmodel="linear"
            angle=pmag(1,im)*radgra1
            dlength=pmag(1,im)*pmag(2,im)
            edge=angle/2.0d0
            fint=pmag(6,im)
            hgap=pmag(7,im)/2.0d0
            call csbend(cbmodel,strength,angle,dlength,edge,seclen,
     &        posi,fint,hgap,de,dmyenergy,bmovecut,ds,istatus)
            pmag(11,im)=seclen
            pmag(12,im)=strength
          else IF (CTYP(IM).EQ.'DQS') THEN
            READ(LUNMG,*)CDUM1,CDUM2,PMAG(1:10,IM)
            ds=sqrt(pmag(8,im)**2+pmag(9,im)**2+pmag(10,im)**2)
            if (ds.eq.0.0d0) then
              print*," "
              print*,"*** Error in bmagseq: Bad normal vector for entrance plane of DQS element ***"
              print*," "
              stop
            else
              pmag(8:10,im)=pmag(8:10,im)/ds
            endif
            ds=1.0D0/dble(myinum)
            cbmodel="quintic-spline"
            angle=pmag(1,im)*radgra1
            dlength=pmag(1,im)*pmag(2,im)
            edge=angle/2.0d0
            fint=pmag(6,im)
            hgap=pmag(7,im)/2.0d0
            call csbend(cbmodel,strength,angle,dlength,edge,seclen,
     &        posi,fint,hgap,de,dmyenergy,bmovecut,ds,istatus)
            pmag(11,im)=seclen
            pmag(12,im)=strength
          else IF (CTYP(IM).EQ.'DH') THEN
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1,IM),PMAG(2,IM),PMAG(3,IM),PMAG(4,IM)
          else IF (CTYP(IM).EQ.'DIF') THEN
            nfoumags=nfoumags+1
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1,IM),PMAG(2,IM),PMAG(3,IM),PMAG(4,IM),
     &        pmag(5,im),pmag(6,im)
            xfoubounds(4,nfoumags)=pmag(5,im)
            xfoubounds(5,nfoumags)=pmag(6,im)
            if (pmag(6,im).gt.nfoumagcp) then
              print*,"*** Error in BMAGSEQ: Number of Fourier coefficients exceeds dimension nfoumagcp =",nfoumagcp,"  ***"
              stop "*** Program WAVE aborted ***"
            endif
          else IF (CTYP(IM).EQ.'DHF') THEN
            nfoumags=nfoumags+1
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1,IM),PMAG(2,IM),PMAG(3,IM),PMAG(4,IM),
     &        pmag(5,im),pmag(6,im)
            xfoubounds(4,nfoumags)=pmag(5,im)
            xfoubounds(5,nfoumags)=pmag(6,im)
            if (pmag(6,im).gt.nfoumagcp) then
              print*,"*** Error in BMAGSEQ: Number of Fourier coefficients exceeds dimension nfoumagcp =",nfoumagcp,"  ***"
              stop "*** Program WAVE aborted ***"
            endif
          ELSE IF (CTYP(IM).EQ.'QP'.OR.CTYP(IM).EQ.'QF') THEN
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1,IM),PMAG(2,IM),PMAG(3,IM),PMAG(4,IM),PMAG(5,IM)
          ELSE IF (CTYP(IM).EQ.'SX') THEN
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1,IM),PMAG(2,IM),PMAG(3,IM),PMAG(4,IM),PMAG(5,IM)
          ELSE IF (CTYP(IM).EQ.'UE') THEN
            READ(LUNMG,*)CDUM1,CDUM2,
     &        PMAG(1:13,IM)
            ! 1. 2.   3.     4.    5.      6.     7.     8.       9.   10.
            ! K  B0V, B0H, Shift, XCen, PerLen, NPer, Lambda_X, Nharm Eharm
            ! 11.     12.    13.
            ! ctaper z-shift ang
            shift=pmag(4,im)
            xcen=pmag(5,im)
            perlen=pmag(6,im)
            pern=pmag(7,im)
            xlamb=pmag(8,im)
            ahwpol=((pern-1.0d0)*2+1.0d0)
            totlen=perlen*((ahwpol-1.0d0)/2.0d0+1.0d0)+shift
            uebounds(1,im)=xcen-0.5d0*totlen
            uebounds(2,im)=xcen+0.5d0*totlen
          ELSE

            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BMAGSEQ ***'
            WRITE(LUNGFO,*)'ILLEGAL MAGNET TYP ',ctyp(im), ' ON FILE FILEMG'
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BMAGSEQ ***'
            WRITE(6,*)'ILLEGAL MAGNET TYP ',ctyp(im), ' ON FILE FILEMG'
            WRITE(6,*)

            STOP

          ENDIF !CTYP

        ENDDO   !IM

        CLOSE(LUNMG)

        if (nfoumags.gt.maxfoumagp) then
          print*
          print*,"*** Error in BMAGSEQ: Too many Fourier magnets ***"
          print*,"*** Please, increase MAXFOUMAGSP and recompile WAVE"
          print*,"*** Program WAVE aborted ***"
          print*
          stop
        endif

C---{ CORRECT FOR FRINGE-FIELD-EFFECTS

        IF (KMAGCOR.NE.0) THEN

          BSHIFT=0.5D0          !DONT WORRY
          BETA=DSQRT((1.D0-1.D0/DMYGAMMA)*(1.D0+1.D0/DMYGAMMA))
          V0=CLIGHT1*BETA
          DTIM=1.0D0/(v0*myinum)   !TIME INTERVALLS FOR TRACKING

          X1=XSTART
          Y1=YSTART
          Z1=ZSTART
          VX1=VXIN
          VY1=VYIN
          VZ1=VZIN

          DO IM=1,mmag
            IF (CTYP(IM).EQ.'DI'.or.ctyp(im).eq.'DIF') THEN
c              CALL BDI(XSTART,Y1,Z1,BXOUT,BYOUT,BZOUT,IM)
              ! The correction has only an effect, if the plateau is not reached
              xlen2=dabs(pmag(2,im)*sin(pmag(1,im)/2.0d0))
              dint=exp(2.0d0*pmag(4,im)*xlen2)/
     &          (pmag(4,im)*(exp(2.0d0*pmag(4,im)*xlen2)-1.0d0))*
     &          2.0d0*pmag(4,im)*xlen2
              corr(im)=corr(im)/dint*(2.0d0*xlen2)
            else IF (CTYP(IM).EQ.'DH'.or.ctyp(im).eq.'DHF') THEN
c              CALL BDH(XSTART,Y1,Z1,BXOUT,BYOUT,BZOUT,IM)
              ! The correction has only an effect, if the plateau is not reached
              xlen2=dabs(pmag(2,im)*sin(pmag(1,im)/2.0d0))
              dint=exp(2.0d0*pmag(4,im)*xlen2)/
     &          (pmag(4,im)*(exp(2.0d0*pmag(4,im)*xlen2)-1.0d0))*
     &          2.0d0*pmag(4,im)*xlen2
              corr(im)=corr(im)/dint*(2.0d0*xlen2)
c            ELSE !DI
c              BXOUT=0.0D0
c              BYOUT=0.0D0
c              BZOUT=0.0D0
            ENDIF !DI
c            IF (BYOUT.GT.1.0D-6) THEN
c              WRITE(LUNGFO,*)
c              WRITE(LUNGFO,*)
c     &          '*** WARNING IN BMAGSEQ: MAGNETIC FIELD FOR XSTART NOT ZERO'
c              WRITE(6,*)
c              WRITE(6,*)
c     &          '*** WARNING IN BMAGSEQ: MAGNETIC FIELD FOR XSTART NOT ZERO'
c            STOP 'PROGRAM ABORTED'
c            ENDIF
          ENDDO   !mmag

          if (kmagcor.gt.0) then

            VN=V0/DSQRT(VX1*VX1+VY1*VY1+VZ1*VZ1)

            VX1=VX1*VN
            VY1=VY1*VN
            VZ1=VZ1*VN

            ANG2y=DATAN(Vy1/VX1)
            ANG2z=DATAN(VZ1/VX1)

            DO IMag=1,mmag

              IF (CTYP(IM).EQ.'DI'.or.ctyp(im).eq.'DIF') THEN
                im=imag
              else if (CTYP(IM).EQ.'DH'.or.ctyp(im).eq.'DHF') THEN
                im=-imag
              else
                cycle
              endif

              ANG1y=ANG2y
              ANG1z=ANG2z

              CALL TRACKBMAG(1,X1,Y1,Z1,VX1,VY1,VZ1,
     &          xstop,0.D0,0.D0,1.D0,0.D0,0.D0,
     &          X2,Y2,Z2,VX2,VY2,VZ2,DTIM,BSHIFT,DMYGAMMA,IM,BMOVECUT,
     &          IUSTEP,IENELOSS)

              ANG2y=DATAN(Vy2/VX2)
              ANG2z=DATAN(VZ2/VX2)
              DANGy=ANG2y-ANG1y
              DANGz=ANG2z-ANG1z

              IF (im.lt.0.and.DANGy.NE.0.0) THEN
                CORR(IMag)=corr(imag)*DABS(PMAG(1,IMag)/DANGy)
              else IF (im.gt.0.and.DANGz.NE.0.0) THEN
                CORR(IMag)=corr(imag)*DABS(PMAG(1,IMag)/DANGz)
              ENDIF

              CALL TRACKBMAG(1,X1,Y1,Z1,VX1,VY1,VZ1,
     &          xstop,0.D0,0.D0,1.D0,0.D0,0.D0,
     &          X2,Y2,Z2,VX2,VY2,VZ2,DTIM,BSHIFT,DMYGAMMA,IM,BMOVECUT,
     &          IUSTEP,IENELOSS)

              ANG2y=DATAN(Vy2/VX2)
              ANG2z=DATAN(VZ2/VX2)
              DANGy=ANG2y-ANG1y
              DANGz=ANG2z-ANG1z

              IF (im.lt.0.and.DANGy.NE.0.0) THEN
                CORR(IMag)=corr(imag)*DABS(PMAG(1,IMag)/DANGy)
              else IF (im.gt.0.and.DANGz.NE.0.0) THEN
                CORR(IMag)=corr(imag)*DABS(PMAG(1,IMag)/DANGz)
              ENDIF

            ENDDO !mmag

          ENDIF   !(KMAGCOR.NE.0)

          WRITE(LUNGFO,*)'     after corrections:'
          WRITE(LUNGFO,*)
          DO IM=1,mmag
            WRITE(LUNGFO,1200) CTYP(IM),
     &        PMAG(1,IM),PMAG(2,IM)/CORR(IM),PMAG(3:13,IM)
1200        FORMAT('      ',A,13E14.6)
          ENDDO !IM
          WRITE(LUNGFO,*)

        ENDIF  !(KMAGCOR.NE.0)

C---} CORRECT FOR FRINGE-FIELD-EFFECTS

        if (nfoumags.gt.0) then

          mfour=nint(alog(float(nfoumagcp))/alog(2.0E0))
          nfoumags=0

          DO imag=1,mmag

            IF (CTYP(IMag).EQ.'DIF') THEN
              nfoumags=nfoumags+1
              xlen2=dabs(pmag(2,imag)*sin(pmag(1,imag)/2.0d0))
     &          +70.0d0/pmag(4,imag)
              dxfour=4.0d0*xlen2/nfoumagcp
              do i=1,nfoumagcp/2+1
                xfour(i)=-dxfour*(nfoumagcp/2+1-i)
                xfour(nfoumagcp+1-i+1)=-xfour(i)
              enddo
              do i=1,nfoumagcp/2+1
                i1=i-1
                ip=nfoumagcp/2+1+i1
                im=nfoumagcp/2+1-i1
                call bdi(pmag(3,imag)-xfour(ip),0.0d0,0.0d0,bx,by,bz,imag)
                yfour(ip)=by
                yfour(im)=by
              enddo
              xfoubounds(1,nfoumags)=imag
              xfoubounds(2,nfoumags)=xfour(1)+pmag(3,imag)
              xfoubounds(3,nfoumags)=-xfour(1)+pmag(3,imag)
              call rfft(ckoef,-mfour) !fft mit cern-routine d703
              do k=1,nint(xfoubounds(5,nfoumags))  !reelle koeffizienten
                foumags(k,nfoumags)=(-1.)**(k-1)*2.0*real(ckoef(k))
              enddo
            else IF (CTYP(IMag).EQ.'DHF') THEN
              nfoumags=nfoumags+1
              xlen2=dabs(pmag(2,imag)*sin(pmag(1,imag)/2.0d0))
     &          +70.0d0/pmag(4,imag)
              dxfour=4.0d0*xlen2/nfoumagcp
              do i=1,nfoumagcp/2+1
                xfour(i)=-dxfour*(nfoumagcp/2+1-i)
                xfour(nfoumagcp+1-i+1)=-xfour(i)
              enddo
              do i=1,nfoumagcp/2+1
                i1=i-1
                ip=nfoumagcp/2+1+i1
                im=nfoumagcp/2+1-i1
                call bdh(pmag(3,imag)-xfour(ip),0.0d0,0.0d0,bx,by,bz,imag)
                yfour(ip)=bz
                yfour(im)=bz
              enddo
              xfoubounds(1,nfoumags)=imag
              xfoubounds(2,nfoumags)=xfour(1)+pmag(3,imag)
              xfoubounds(3,nfoumags)=-xfour(1)+pmag(3,imag)
              call rfft(ckoef,-mfour) !fft mit cern-routine d703
              do k=1,nint(xfoubounds(5,nfoumags))  !reelle koeffizienten
                foumags(k,nfoumags)=(-1.)**(k-1)*2.0*real(ckoef(k))
              enddo
            ENDIF !CTYP

          ENDDO   !IM

        endif !nfoumags

      ENDIF !ICAL

      if (ical.lt.10) ical=ical+1

      if (ical.eq.2) then

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     Subroutine BMAGSEQ: Magnets read from file:'
        write(lungfo,*)
        WRITE(LUNGFO,*)'     ',FILEMG(1:len_trim(filemg))
        WRITE(LUNGFO,*)'     X-shift (XSHMAGSEQ):',xshmagseq
        write(lungfo,*)
        WRITE(LUNGFO,1100)mmag,NMGSQP
1100    FORMAT('      Number of magnets: ',I6,' (Limit IS ',I6,' Magnets)')
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     Name Typ Boundaries Parameters'
        WRITE(LUNGFO,*)

        DO IM=1,mmag
          if (ctyp(im).eq.'DI') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dibounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'BEND') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dibounds(1:2,im),PMAG(1:12,IM)
            write(lungfo,*)"     vnin(1),vnin(3):",pmag(14,im),pmag(15,im)
            write(lungfo,*)"     vnout(1),vnout(3):",pmag(16,im),pmag(17,im)
          else if (ctyp(im).eq.'DIL') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dibounds(1:2,im),PMAG(1:11,IM)
          else if (ctyp(im).eq.'DCS') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dibounds(1:2,im),PMAG(1:11,IM)
          else if (ctyp(im).eq.'DQS') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dibounds(1:2,im),PMAG(1:11,IM)
          else if (ctyp(im).eq.'QF') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),qfbounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'QP') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),qfbounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'SX') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),sxbounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'DH') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dhbounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'DIF') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dibounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'DHF') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),dhbounds(1:2,im),PMAG(1:6,IM)
          else if (ctyp(im).eq.'UE') then
            WRITE(LUNGFO,1201) trim(cnam(im)),ctyp(im),uebounds(1:2,im),PMAG(1:13,IM)
          endif
1201      FORMAT('      ',a,' ',a5,15E14.6)
        enddo
        write(lungfo,*)
        WRITE(LUNGFO,*)
        ICAL=2
      ENDIF !ICAL

C--- MAGNETIC FIELD

      CALL BMAGSEQC(XIN+xshmagseq,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

      RETURN
      END
