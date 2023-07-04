*CMZ :  3.02/00 26/08/2014  08.51.38  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.09.17  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/13 07/07/2010  11.06.21  by  Michael Scheer
*CMZ :  2.16/08 25/06/2010  12.15.46  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.36  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.03  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.15.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.21  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.09  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WDISPER
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

C--- CALCULATES LINEARE DISPERSION
C    THE INTERNAL DISPERSION OF THE WLS IS CALCULATE FROM THE TRACKING OF
C    CHROMATIC TRAJECTORIES. THE EXTERNAL DISPERSION IS CALCULATED FROM
C    THE START VALUES AND THE LINEAR TRANSFER MATRIX


      IMPLICIT NONE

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
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.


      INTEGER IE,ICOUNT,IEE,IP,IPP,IERROR,IERRFLG,iwarn

      DOUBLE PRECISION X1,Y1,Z1,X2,Y2,Z2
      DOUBLE PRECISION XF0,YF0,ZF0,VXF0,VYF0,VZF0
      DOUBLE PRECISION VX1,VY1,VZ1,VX2,VY2,VZ2

      DOUBLE PRECISION GAMMA,DS02

      DATA IERRFLG/0/

      ALLOCATE(WDIS(2,NCO))

      DS02=2.D0*DS0

C---LOOP OVER CHROMATIC CLOSED ORBITS

      IEE=0

      DO IE=-1,1,2

      IEE=IEE+1
      GAMMA=DMYGAMMA+DFLOAT(IE)*DELGAM*DMYGAMMA

C--- START OF TRAJECTORY

      X1=WSXYZ(1,1)
      Y1=WSXYZ(2,1)
      Z1=WSXYZ(3,1)

      VX1=WVXYZ(1,1)
      VY1=WVXYZ(2,1)
      VZ1=WVXYZ(3,1)

      WDIS(IEE,1)=0.0

C--- LOOP OVER POINTS

      IPP=0

      DO ICOUNT=1,NCO-1

        IPP=ICOUNT+1

        IF (abs(Y1).gt.1.0d-10.and.iwarn.eq.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*) '*** WARNING IN WDISPER ***'
          WRITE(LUNGFO,*) 'ORBIT NOT PLANAR'
          WRITE(LUNGFO,*)
          WRITE(6,*) '*** WARNING IN WDISPER ***'
          iwarn=1
C            STOP
        ENDIF

        XF0=WSXYZ(1,IPP)
        YF0=WSXYZ(2,IPP)
        ZF0=WSXYZ(3,IPP)

        VXF0=WVXYZ(1,IPP)
        VYF0=WVXYZ(2,IPP)
        VZF0=WVXYZ(3,IPP)

C- ROTATE ELECTRON TO NEXT POINT SUCH THAT ALL CORRESPONDING POINTS
C I.E. WITH THE SAME ICOUNT LAY IN ONE PLANE

        CALL WBTRACK(X1,Y1,Z1,VX1,VY1,VZ1,
     &    XF0,YF0,ZF0,VXF0,VYF0,VZF0,
     &    X2,Y2,Z2,VX2,VY2,VZ2,
     &    DTIM0,GAMMA,IERROR)

        IF (IERROR.NE.0) IERRFLG=1

        WDIS(IEE,IPP)=DSIGN
     &    (
     &    DSQRT
     &    (
     &    (X2-XF0)**2+(Z2-ZF0)**2+(Y2-YF0)**2
     &    )
     &    , (Z2-ZF0)*VXF0-(X2-XF0)*VZF0    !Y-COMPONENT OFF VECTOR-PRODUCT
     &    )
     &    /(DFLOAT(IE)*DELGAM)


        X1=X2
        Y1=Y2
        Z1=Z2

        VX1=VX2
        VY1=VY2
        VZ1=VZ2

      ENDDO !ICOUNT

      ENDDO !IE

C--- CALCULATE INTERNAL DISPERSION

      DO IP=1,NCO
        WBETA(6,IP)=(WDIS(1,IP)+WDIS(2,IP))/2.D0
      ENDDO

C--- CALCULATE DERIVATIVE OF INTERNAL DISPERSION

      DO IP=2,NCO-1
        WBETA(7,IP)=(WBETA(6,IP+1)-WBETA(6,IP-1))/DS02
      ENDDO

      WBETA(7,1)=WBETA(7,2)-(WBETA(7,3)-WBETA(7,2))
      WBETA(7,NCO)=WBETA(7,NCO-1)+(WBETA(7,NCO-1)-WBETA(7,NCO-2))


C--- TOTAL DISPERSION, I.E. TRANSFORMED EXTERNAL PLUS INTERNAL

      DO IP=1,NCO

        WBETA(6,IP)=WLTM(1,1,IP)*DISP0+WLTM(1,2,IP)*DDISP0+WBETA(6,IP)
        WBETA(7,IP)=WLTM(1,3,IP)*DISP0+WLTM(1,4,IP)*DDISP0+WBETA(7,IP)

      ENDDO !NCO

      IF (IERRFLG.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*) '*** WARNING SR WDISPER ***'
        WRITE(LUNGFO,*) 'PROBLEMS WITH SR WBTRACK'
        WRITE(LUNGFO,*) 'CHANGE STEP SIZE, WATCH RESULTS'
        WRITE(LUNGFO,*)
      ENDIF

      DEALLOCATE(WDIS)

      RETURN
      END
