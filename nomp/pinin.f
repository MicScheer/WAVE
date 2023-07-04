*CMZ :  4.00/15 12/02/2022  17.08.21  by  Michael Scheer
*CMZ :  4.00/13 04/12/2021  12.10.40  by  Michael Scheer
*CMZ :  4.00/01 05/04/2019  15.09.32  by  Michael Scheer
*CMZ :  3.06/00 25/02/2019  17.19.30  by  Michael Scheer
*CMZ :  3.02/08 24/06/2015  16.06.32  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.59.38  by  Michael Scheer
*CMZ :  2.67/02 28/03/2012  08.22.03  by  Michael Scheer
*CMZ :  2.66/06 22/05/2010  16.48.23  by  Michael Scheer
*CMZ :  2.64/06 15/09/2009  11.10.54  by  Michael Scheer
*CMZ :  2.64/05 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.62/03 16/07/2007  11.51.09  by  Michael Scheer
*CMZ :  2.62/02 16/07/2007  06.52.55  by  Michael Scheer
*CMZ :  2.54/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/10 05/11/2004  16.59.34  by  Michael Scheer
*CMZ :  2.51/00 17/05/2004  17.43.59  by  Michael Scheer
*CMZ :  2.47/16 11/09/2003  15.10.02  by  Michael Scheer
*CMZ :  2.34/09 18/09/2001  22.52.09  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.38.34  by  Michael Scheer
*CMZ :  2.17/00 03/11/2000  09.47.59  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  12.09.17  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.16/00 27/05/2000  14.03.56  by  Michael Scheer
*CMZ :  2.15/00 08/05/2000  13.32.10  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  2.13/02 14/12/99  16.24.13  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.47.16  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  12.37.49  by  Michael Scheer
*CMZ : 00.02/00 19/11/96  14.57.13  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  17.29.50  by  Michael Scheer
*CMZ : 00.01/06 01/02/95  16.35.43  by  Michael Scheer
*CMZ : 00.01/04 29/11/94  10.17.51  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.07.25  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.56  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.26  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PININ
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

C--- INITIALIZE GRID OF OBERSERVATION POINTS OF PINHOLE

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


      INTEGER IOB,IY,IZ,INCZ,INCY,ISOUR,INCZMX,INCYMX
      INTEGER N2POWZ,N2POWY,ICDUMZ,ICDUMY

      DOUBLE PRECISION OBSVDUM(3,1),x,y,z,xn,yn,zn,r,r0,pinwo,pinho,pinro

      INTEGER IFPHIR_A,iobsv
      DATA IFPHIR_A/0/

      pinwo=pinw
      pinho=pinh
      pinro=pinr

      IF (IPINCIRC.NE.0) THEN
        PINW=2.D0*PINR
        PINH=2.D0*PINR
      ENDIF !IPINCIRC

      IF (OBSVDZ.NE.0.D0) THEN
        MPINZ=NINT(PINW/OBSVDZ)+1
        if (mpinz.lt.3) then
          write(6,*)"*** Warning in pinin: OBSVDZ > PINW/2., will be adjusted ***"
          write(lungfo,*)"*** Warning in pinin: OBSVDZ > PINW/2., will be adjusted ***"
          mpinz=3
          obsvdz=pinw/2.0d0
        endif
        PINW=OBSVDZ*(MPINZ-1)
      ENDIF !(OBSVDZ.NE.0.D0)

      IF (OBSVDY.NE.0.D0) THEN
        MPINY=NINT(PINH/OBSVDY)+1
        if (mpiny.lt.3) then
          write(6,*)"*** Warning in pinin: OBSVDY > PINH/2., will be adjusted ***"
          write(lungfo,*)"*** Warning in pinin: OBSVDY > PINH/2., will be adjusted ***"
          mpiny=3
          obsvdy=pinh/2.0d0
        endif
        PINH=OBSVDY*(MPINY-1)
      ENDIF !(OBSVDZ.NE.0.D0)

      IF (IPINCIRC.NE.0) THEN
        PINR=min(pinw,pinh)/2.0d0
      ENDIF !IPINCIRC

      if (abs(pinwo-pinw).gt.pinw/10000.0d0) then
        print*,"*** Warning in PINI: PINW is adjusted, according to MPINZ and DOBSVDZ ***"
      endif

      if (abs(pinho-pinh).gt.pinh/10000.0d0) then
        print*,"*** Warning in PINI: PINH is adjusted, according to MPINY and DOBSVDY ***"
      endif

      if (abs(pinro-pinr).gt.pinr/10000.0d0) then
        print*,"*** Warning in PINI: PINR is adjusted, according to PINW and PINH ***"
      endif

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

C--- DATA OF PINHOLE ARE TAKEN FORM NAMELIST

      IF (IF1DIM.NE.0.AND.
     &    (MEDGEZ.NE.0.OR.MMEDGEZ.NE.0.OR.MPINZ.NE.1)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN PININ ***'
        WRITE(LUNGFO,*)'FLAG IF1DIM SET BUT'
        WRITE(LUNGFO,*)
     &    'MEDGEZ OR MMEDGEZ IN NAMELIST PINHOLE NOT ZERO OR MPINZ NOT EQUAL ONE'
        WRITE(LUNGFO,*)'ADJUSTED TO APPROPRIATE VALUES'
        WRITE(LUNGFO,*)
c        WRITE(6,*)
c        WRITE(6,*)
c        WRITE(6,*)'*** WARNING IN PININ ***'
c        WRITE(6,*)'FLAG IF1DIM SET BUT'
c        WRITE(6,*)
c     &    'MEDGEZ OR MMEDGEZ IN NAMELIST PINHOLE NOT ZERO OR MPINZ NOT EQUAL ONE'
c        WRITE(6,*)'ADJUSTED TO APPROPRIATE VALUES'
c        WRITE(6,*)
        MPINZ=1
        MEDGEZ=0
        MMEDGEZ=0
      ENDIF

      IF (IUSEM.NE.0) THEN
C MAKE MPINZ AND MPINY EVEN
        MPINZ=(MPINZ+1)/2*2
        MPINY=(MPINY+1)/2*2
      ENDIF !IUSEM

      NOBSVY=MPINY
      NOBSVZ=MPINZ
      NOBSV=NOBSVY*NOBSVZ
      IF (NOBSV.LE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'*** NEGATIVE NUMBER OF OBSERVATION POINTS!'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'*** NEGATIVE NUMBER OF OBSERVATION POINTS!'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      OBSVDUM(1,1)=PINCEN(1)

      IF (MPINZ.GT.1) THEN
        OBSVDZ=PINW/DFLOAT(MPINZ-1)
        OBSVDUM(3,1)=PINCEN(3)-PINW/2.
      ELSE
        OBSVDZ=PINW
        OBSVDUM(3,1)=PINCEN(3)
      ENDIF

      IF (MPINY.GT.1) THEN
        OBSVDY=PINH/DFLOAT(MPINY-1)
        OBSVDUM(2,1)=PINCEN(2)-PINH/2.
      ELSE
        OBSVDY=PINH
        OBSVDUM(2,1)=PINCEN(2)
      ENDIF

      IF (IF1DIM.EQ.0.AND.(MEDGEZ.LT.1.OR.MEDGEY.LT.1)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'MEDGEZ OR MEDGEY IN NAMELIST PINHOLE LOWER THAN 1'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'MEDGEZ OR MEDGEY IN NAMELIST PINHOLE LOWER THAN 1'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF

      IF (MMEDGEZ.LT.0.OR.MMEDGEY.LT.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'MMEDGEZ OR MMEDGEY IN NAMELIST PINHOLE LOWER THAN 0'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'MMEDGEZ OR MMEDGEY IN NAMELIST PINHOLE LOWER THAN 0'
        WRITE(6,*)
        WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
        STOP
      ENDIF

C- DETERMINE RMS VALUES FOR FOLDING

      IF (IFOLD.NE.0) CALL WSIGFOL

C- INCREASE PINHOLE

      MOBSVZ=NOBSVZ  !STORE VALUES
      MOBSVY=NOBSVY
      MOBSV=NOBSV

      INCZ=0
      INCY=0

      IF (IFOLD.NE.0.AND.IFOLD.NE.2) THEN

        INCZMX=-1
        INCYMX=-1
        DO ISOUR=1,NSOURCE
          INCZ=NINT(WSIGZ(ISOUR)*DGSIGZ(ISOUR)/OBSVDZ)
          IF(INCZ.GT.INCZMX) THEN
            INCZMX=INCZ
          ENDIF
          INCY=NINT(WSIGY(ISOUR)*DGSIGY(ISOUR)/OBSVDY)
          IF(INCY.GT.INCYMX) THEN
            INCYMX=INCY
          ENDIF
        ENDDO   !ISOUR

        INCZ=INCZMX
        INCY=INCYMX
        IF (IF1DIM.NE.0) INCZ=0

      ENDIF !IFOLD

C           FACTOR 2 FOR BOTH SIDES OF PINHOLE
      NOBSVZ=MOBSVZ+2*(INCZ+MEDGEZ+MMEDGEZ)
      NOBSVY=MOBSVY+2*(INCY+MEDGEY+MMEDGEY)

C--- ADJUST NUMBER OF POINTS ACCORDING TO POWERS OF IF FLAG IUSEM IS SET

      IF (IUSEM.NE.0) THEN

        N2POWZ=NINT(ALOG(FLOAT(NOBSVZ-1))/ALOG(2.))
        IF(NOBSVZ .GT. 2**N2POWZ) N2POWZ=N2POWZ+1
C150793        NOBSVZ=2**N2POWZ+1
        NOBSVZ=2**N2POWZ

        N2POWY=NINT(ALOG(FLOAT(NOBSVY-1))/ALOG(2.))
        IF(NOBSVY .GT. 2**N2POWY) N2POWY=N2POWY+1
C150793        NOBSVY=2**N2POWY+1
        NOBSVY=2**N2POWY

      ENDIF !IUSEM

      OBSVDUM(2,1)=OBSVDUM(2,1)-(NOBSVY-MOBSVY)/2*OBSVDY
      OBSVDUM(3,1)=OBSVDUM(3,1)-(NOBSVZ-MOBSVZ)/2*OBSVDZ

      NOBSV=NOBSVZ*NOBSVY
      IF (NOBSV.GT.1000000) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN PININ ***'
        WRITE(LUNGFO,*)'*** MORE THAN 1 000 0000 OBSERVATION POINTS ***'
        WRITE(LUNGFO,*)'NUMBER OF POINTS IN Y:',NOBSVY
        WRITE(LUNGFO,*)'NUMBER OF POINTS IN Z:',NOBSVZ
        IF (IFOLD.NE.0) WRITE(LUNGFO,*)'PLEASE CHECK PINHOLE SIZE AND RELATION TO BEAM EMITTANCE'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN PININ ***'
        WRITE(6,*)'*** MORE THAN 1 000 000 OBSERVATION POINTS ***'
        WRITE(6,*)'NUMBER OF POINTS IN Y:',NOBSVY
        WRITE(6,*)'NUMBER OF POINTS IN Z:',NOBSVZ
        IF (IFOLD.NE.0) WRITE(6,*)'PLEASE CHECK PINHOLE SIZE AND RELATION TO BEAM EMITTANCE'
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      IF (IOBSV_A.NE.NOBSV) THEN
        IF (IOBSV_A.NE.0) DEALLOCATE(OBSV)
        ALLOCATE(OBSV(3,NOBSV))
        IOBSV_A=NOBSV
      ENDIF !(IOBSV_A.LT.NOBSV)

      IF (IOBSVZ_A.NE.NOBSVZ) THEN
        IF (IOBSVZ_A.NE.0) DEALLOCATE(OBSVZ)
        ALLOCATE(OBSVZ(NOBSVZ))
        IOBSVZ_A=NOBSVZ
      ENDIF !(IOBSVY_A.LT.NOBSVY)

      IF (IOBSVY_A.NE.NOBSVY) THEN
        IF (IOBSVY_A.NE.0) DEALLOCATE(OBSVY)
        ALLOCATE(OBSVY(NOBSVY))
        IOBSVY_A=NOBSVY
      ENDIF !(IOBSV_A.LT.NOBSV)

      OBSV(1,1)=OBSVDUM(1,1)
      OBSV(2,1)=OBSVDUM(2,1)
      OBSV(3,1)=OBSVDUM(3,1)

      IOB=0
      DO IY=1,NOBSVY
        DO IZ=1,NOBSVZ
          IOB=IOB+1
          OBSV(3,IOB)=OBSV(3,1)+OBSVDZ*(IZ-1)
          OBSV(2,IOB)=OBSV(2,1)+OBSVDY*(IY-1)
          OBSV(1,IOB)=OBSV(1,1)
        ENDDO
      ENDDO

      DO IZ=1,NOBSVZ
        OBSVZ(IZ)=OBSV(3,IZ)
      ENDDO

      DO IY=1,NOBSVY
        IOB=(IY-1)*NOBSVZ+1
        OBSVY(IY)=OBSV(2,IOB)
      ENDDO

      IF (IPIN.EQ.0) THEN
        ICBRILL=1
      ELSE !PIN
        IF (IPBRILL.EQ.0) THEN
          ICDUMZ=NOBSVZ/2+1
          ICDUMY=NOBSVY/2+1
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.1) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+1
          ICDUMY=(NOBSVY-MOBSVY)/2+1
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.2) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+MOBSVZ
          ICDUMY=(NOBSVY-MOBSVY)/2+1
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.3) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+1
          ICDUMY=(NOBSVY-MOBSVY)/2+MOBSVY
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE IF (IPBRILL.EQ.4) THEN
          ICDUMZ=(NOBSVZ-MOBSVZ)/2+MOBSVZ
          ICDUMY=(NOBSVY-MOBSVY)/2+MOBSVY
          ICBRILL=ICDUMZ+NOBSVZ*(ICDUMY-1)
        ELSE
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN PININ: IPBRILL WRONG ***'
          WRITE(LUNGFO,*)'*** CHECK NAMELIST PINHOLE ***'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED  ***'
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN PININ: IPBRILL WRONG ***'
          WRITE(6,*)'*** CHECK NAMELIST PINHOLE ***'
          WRITE(6,*)
          WRITE(6,*)'*** PROGRAM WAVE ABORTED  ***'
          STOP
        ENDIF
      ENDIF   !IPIN

      IF (MPINR.NE.0) IRPHI=1
      IF (IPINCIRC.NE.0.AND.IFPHIR_A.NE.NOBSV) THEN
        ALLOCATE(FPHIR(NOBSV))
        IF (IRPHI.NE.0) THEN
          ALLOCATE(XC(NOBSV))
          ALLOCATE(YC(NOBSV))
        ENDIF
        IFPHIR_A=NOBSV
      ENDIF !IPINCIRC

      IF (MPINR.NE.0.AND.IPIN.EQ.0) THEN
        MPINR=0
        WRITE(6,*)
        WRITE(6,*)
     &    '*** WARNING IN PININ: IPIN.EQ.0, THEREFORE MPINR SET ZERO ***'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING IN PININ: IPIN.EQ.0, THEREFORE MPINR SET ZERO ***'
        WRITE(LUNGFO,*)
      ENDIF

      IF (MPINR.NE.0.AND.IF1DIM.NE.0) THEN
        MPINR=0
        WRITE(6,*)
        WRITE(6,*)
     &    '*** WARNING IN PININ: IF1DIM.NE.0, THEREFORE MPINR SET ZERO ***'
        WRITE(6,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '*** WARNING IN PININ: IF1DIM.EQ.0, THEREFORE MPINR SET ZERO ***'
        WRITE(LUNGFO,*)
      ENDIF

      IF (MPINR.NE.0) THEN
        CALL PININR
      ENDIF

      if (rpinsph.eq.-9999.0d0) rpinsph=pincen(1)

      if (rpinsph.ne.0) then
        r0=obsv(1,1)-rpinsph
        do iobsv=1,nobsv
          x=obsv(1,iobsv)-r0
          y=obsv(2,iobsv)
          z=obsv(3,iobsv)
          r=sqrt(x*x+y*y+z*z)
          xn=x/r
          yn=y/r
          zn=z/r
          obsv(1,iobsv)=r0+xn*rpinsph
          obsv(2,iobsv)=yn*rpinsph
          obsv(3,iobsv)=zn*rpinsph
        enddo
        if (mpinr.ne.0) then
          do iobsv=1,nobsvrphi
            x=obsvrphi(1,iobsv)-r0
            y=obsvrphi(2,iobsv)
            z=obsvrphi(3,iobsv)
            r=sqrt(x*x+y*y+z*z)
            xn=x/r
            yn=y/r
            zn=z/r
            obsvrphi(1,iobsv)=xn*rpinsph
            obsvrphi(2,iobsv)=yn*rpinsph
            obsvrphi(3,iobsv)=zn*rpinsph
          enddo
        endif

        DO IZ=1,NOBSVZ
          OBSVZ(IZ)=OBSV(3,IZ)
        ENDDO

        DO IY=1,NOBSVY
          IOB=(IY-1)*NOBSVZ+1
          OBSVY(IY)=OBSV(2,IOB)
        ENDDO

      endif

      RETURN
      END
