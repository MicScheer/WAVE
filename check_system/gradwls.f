*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.15/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.51.58  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.24  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE GRADWLS
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

C--- CALCULATE TRANSVERSAL GRADIENT 2n(s)/rho**3=(dB/B/ds) *dz/dx/B/rho**2
C    GRADIENT IS ONLY CORRECT IF THE COORDINATE ALONG THE TRAJECTORY
C    CAN BE APPROXIMATED BY THE LONGITUDINAL COORDINATE X
c**********************************************************
C     TO MAKE SURE THAT THE ROUTINE WORKS RELIABLE, THE GRADIENT MUST
C       NOT BE TO LARGE I.E. NO STRONG EDGE FOCUSSING LIKE IN HARD-EDGE
C     MAGNETS. THEREFORE THIS ROUTINE SHOULD BE CALLED AFTER OR BY
C     SR BETAWLS, BECAUSE BETAWLS DETECTS PATHOLOGICAL CASES
c**********************************************************
      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,betawls.
      include 'betawls.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER IPOI,IP,IM,I0
      DOUBLE PRECISION BX,BY0,BZ,BYP,BYM,AX,AY,AZ,XP,XM,X0,DXP0,DXM0
      DOUBLE PRECISION DBCUT,RHO2I,BRHO2I,DBDX,DBDXP,DBDXM

      DATA DBCUT/10000./   !ALLOWED CHANGE IN dB WITHIN DXD

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*) '*** SR GRADWLS ***'
      WRITE(LUNGFO,*) 'BE AWARE THAT YOUR WLS HAS NO STRONG EDGE-FOCUSSING!'
      WRITE(LUNGFO,*) 'OTHERWISE RESULTS ARE NOT RELIABLE'
      WRITE(LUNGFO,*)

C--- LOOP OVER ALL POINTS MBETA

      DO IPOI=1,MBETA

      IF (IPOI.EQ.1) THEN
          I0=1
          IP=2
          IM=2
      ELSEIF(IPOI.EQ.MBETA) THEN
          I0=MBETA
          IP=I0-1
          IM=I0-1
      ELSE
          I0=IPOI
          IP=IPOI+1
          IM=IPOI-1
      ENDIF

      X0=XBETA(I0)
      XP=X0+(XBETA(IP)-X0)*0.5
      XM=X0+(XBETA(IM)-X0)*0.5
      DXP0=XP-X0
      DXM0=XM-X0

      CALL MYBFELD(XP,0.D0,0.D0,BX,BYP,BZ,AX,AY,AZ)   !FIELD OF WLS
      CALL MYBFELD(X0,0.D0,0.D0,BX,BY0,BZ,AX,AY,AZ)   !FIELD OF WLS
      CALL MYBFELD(XM,0.D0,0.D0,BX,BYM,BZ,AX,AY,AZ)   !FIELD OF WLS

      DBDXP=(BYP-BY0)/DXP0
      DBDXM=(BYM-BY0)/DXM0
      DBDX=0.5*(DBDXP+DBDXM)

      RHO2I=(CLIGHT1*BY0/EMOM)**2 !ACTUEL INVERSE BENDING-RADIUS
      BRHO2I=0.
      IF(RHO2I.NE.0) BRHO2I=RHO2I/BY0

      IF(DABS(BRHO2I).GT.DBCUT) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** SR GRADWLS: ERROR ***'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'FIELD GRADIENT/RHO**3 EXCEEDS CUT!'
          WRITE(LUNGFO,*)'RHO/B*dB/dZ * dZ/dX / RHO**3:',BRHO2I
          WRITE(LUNGFO,*)'CUT:',DBCUT
          WRITE(LUNGFO,*)
      ENDIF

      G4BETA(IPOI)=2.*BRHO2I*DBDX*ZBETAP(IPOI)
C6.4.92 GROSS K=-k+1/rho**2 soll berechnet werden
      BETAK(IPOI)=-DBDX/EMOM*CLIGHT1*ZBETAP(IPOI)+RHO2I

      ENDDO

      RETURN
      END
