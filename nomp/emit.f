*CMZ :  2.52/15 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  10.23.48  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.50.09  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.46  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************
      SUBROUTINE EMIT(B0,XLAM0,FASYM,
     &         E,RHODIP,TAU0E5,BETA0,DI2RING,DI5RING,
     &                  EMIRING,EMI,EMIOPT,EMITOT,EMITOPT,B3AWLS,
     &         F,F0,BETOPT,TAU,TAU1GEV,POLLEV,POLLV1G,ZMAX,
     &                  DI2,DI5,DI5OPT,LUN,
     &                  DISP0,DISPOPT,EMI2OPT,EMIT2OPT,BET2OPT,BETUNI,
     &                  BET2UNI,DI52OPT)
C****************************************************************
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

C CALCULATES APPROXIMATELY EMITTANCE EFFECTS OF AN ASYMMETRIC WLS

C INPUT: B0,XLAM0,FASYM,E,RHODIP,TAU0E5,BETA0,DI2RING,DI5RING,LUN

      IMPLICIT NONE

      DOUBLE PRECISION B0,XLAM0,FASYM                !INPUT
      DOUBLE PRECISION GAMMA,RHODIP,TAU0E5,BETA0,DI2RING,DI5RING    !INPUT
      DOUBLE PRECISION EMI,EMIOPT,EMITOT,EMITOPT,B3RING,B3RG1G      !OUTPUT

      DOUBLE PRECISION E,RHO0,RHO01,XL
      DOUBLE PRECISION F0,F,FN,BETOPT,FBET,PHI0,PHI,PHIN,PHIMAX
      DOUBLE PRECISION DI2,DI5,DI5OPT,DI2TOT,DI5TOT,DI5TOPT,EMIRING
      DOUBLE PRECISION XK,B3WLS,B3AWLS,POLFAC,POLFC1G,TAU,TAU0,TAU1GEV
      DOUBLE PRECISION POLLEV,POLLV1G,ZMAX
      DOUBLE PRECISION XHOMK,DBHOM,XKK
      DOUBLE PRECISION FBETFUN,F0FNFUN,F2,F20,F2BET,DI5T2OPT,DOPTFUN
     &                ,DISP0,DISPOPT,EMI2OPT,EMIT2OPT,BET2OPT,BETUNI,
     &                 BET2UNI,DI52OPT

      INTEGER LUN

*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DATA DBHOM/1.D-4/

      GAMMA=E/EMASSG1
      RHO01=CLIGHT1*1.D-9*B0/E
      RHO0=1./RHO01
      PHI0=XLAM0/RHO0
      XL=XLAM0/2.*(1.+FASYM)
      XKK=2.* PI1/XLAM0
      XHOMK=DSQRT(2.*DBHOM)/XKK
      PHIN=DSQRT(((FASYM+1.)**2/8./FASYM))
      FN=PHIN**3

      F=F0FNFUN(FASYM,DISP0,PHI0,RHO0)
      F0=F/FN
      FBET=FBETFUN(FASYM,DISP0,PHI0,RHO0)
      BETOPT=FBET*XLAM0
      PHI=PHIN*PHI0  !NUR FUER OUTPUT
      DI2=XLAM0/4./RHO0**2*(1.+1./FASYM)
      EMIRING=CQ1*GAMMA**2*DI5RING/DI2RING

      DISPOPT=DOPTFUN(FASYM,PHI0,RHO0)
      F2=F0FNFUN(FASYM,DISPOPT,PHI0,RHO0)
      F20=F2/FN
      F2BET=FBETFUN(FASYM,DISPOPT,PHI0,RHO0)
      BET2OPT=F2BET*XLAM0

      DI5OPT=F*PHI0**3*DI2
      EMIOPT=CQ1*GAMMA**2*DI5OPT/DI2
      EMI=EMIOPT *  0.5*(BETA0/BETOPT+BETOPT/BETA0)
      DI5=EMI/GAMMA**2/CQ1*DI2

      DI2TOT=DI2+DI2RING
      DI5TOT=DI5+DI5RING
      DI5TOPT=DI5OPT+DI5RING

      EMITOPT=CQ1*GAMMA**2*DI5TOPT/DI2TOT
      EMITOT =CQ1*GAMMA**2* DI5TOT/DI2TOT

      DI52OPT=F2*PHI0**3*DI2
      EMI2OPT=CQ1*GAMMA**2*DI52OPT/DI2
      DI5T2OPT=DI52OPT+DI5RING
      EMIT2OPT=CQ1*GAMMA**2*DI5T2OPT/DI2TOT

      BETUNI=-9999.

      IF(DI5OPT/DI2.LT.DI5RING/DI2RING)
     &   BETUNI=(BETOPT*(DSQRT(DI2**2*DI5RING**2-DI5OPT**2*DI2RING
     . **2)+DI2*DI5RING))/(DI5OPT*DI2RING)

      BET2UNI=-9999.

      IF(DI52OPT/DI2.LT.DI5RING/DI2RING)
     &   BET2UNI=(BET2OPT*(DSQRT(DI2**2*DI5RING**2-DI52OPT**2*DI2RING
     . **2)+DI2*DI5RING))/(DI52OPT*DI2RING)

C--- INFLUENCE ON BEAM POLARISATION TIME AND LEVEL

      XK=2.*PI1/XLAM0
      B3WLS=4./3.*B0**3/XK*(1-1./FASYM**2)
      B3AWLS=4./3.*DABS(B0)**3/XK*(1+1./FASYM**2)

      B3RING=2.*PI1*RHODIP*(E/RHODIP/(CLIGHT1*1.D-9))**3
      B3RG1G=2.*PI1*RHODIP*(1./RHODIP/(CLIGHT1*1.D-9))**3
      POLFAC=1.+B3AWLS/B3RING
      POLFC1G=1.+B3AWLS/B3RG1G
      TAU0=TAU0E5/E**5
      TAU=TAU0/POLFAC
      TAU1GEV=TAU0E5/POLFC1G
      POLLEV=POL1CON1*(B3WLS+B3RING)/(B3AWLS+B3RING)
      POLLV1G=POL1CON1*(B3WLS+B3RG1G)/(B3AWLS+B3RG1G)

C--- MAX. DISPLACEMENT AND DEFLECTION

      PHIMAX=PHI0/(2.*PI1)
      ZMAX=1./(16.*PI1)*XLAM0**2/RHO0*
     &      (FASYM+4./PI1) !010891

      IF(LUN.EQ.0) RETURN

C     IF(LUN.NE.6) OPEN(UNIT=LUN,FILE='EMIT.DAT',STATUS='NEW')

      WRITE(LUN,*)
      WRITE(LUN,*) '     SR EMIT:'
      WRITE(LUN,*) '     ========='
      WRITE(LUN,*)
      WRITE(LUN,*)
      WRITE(LUN,*)'     Results from analytical ansatz for asymmetric WLS:'
      WRITE(LUN,*)'     --------------------------------------------------'
      WRITE(LUN,*)
     &'     (only correct if orbit is planar, WLS has no transversal gradient,'
      WRITE(LUN,*)
     &'     WLS is symmetric with respect to origin of coordinate system,'
      WRITE(LUN,*)
     &'     centered in the straight section, (VZ/VX)**2 << 1,'
      WRITE(LUN,*)
     &'     and derivation of external dispersion vanishes outside WLS)'
      WRITE(LUN,*)
      WRITE(LUN,*)
     &'     Beta function at WLS center, I5(WLS), I5/I2 (WLS):'
      WRITE(LUN,*)'     ',SNGL(BETA0),SNGL(DI5),SNGL(DI5/DI2)
      WRITE(LUN,*)
      WRITE(LUN,*)'     Total emittance change:',SNGL(EMITOT/EMIRING)
      WRITE(LUN,*)
      WRITE(LUN,*)
     &'     optimum beta function, corresponding I5(WLS) and total emittance change'
      WRITE(LUN,*)
     &'     for actual external dispersion:'
      WRITE(LUN,*)'     ',SNGL(BETOPT),SNGL(DI5OPT),SNGL(EMITOPT/EMIRING)
      WRITE(LUN,*)
      WRITE(LUN,*)'     optimum external dispersion:',SNGL(DISPOPT)
      WRITE(LUN,*)
     &'     optimum beta function, corresponding I5(WLS) and total emittance change'
      WRITE(LUN,*)
     &'     for optimum external dispersion:'
         WRITE(LUN,*)'     ',SNGL(BET2OPT),SNGL(EMI2OPT/CQ1/GAMMA**2)
     &                     ,SNGL(EMIT2OPT/EMIRING)
      WRITE(LUN,*)
      WRITE(LUN,*)'     neutral beta function:             ',SNGL(BETUNI)
      WRITE(LUN,*)'     neutral beta function for opt. eta:',SNGL(BET2UNI)
      WRITE(LUN,*)

C     IF (LUN.NE.6) CLOSE(LUN)

      RETURN
      END
