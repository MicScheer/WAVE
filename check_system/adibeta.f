*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.24.18  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.35  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.31  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ADIBETA(GAMMA,B0,XLP,FB0N,FB0M)

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

C--- BERECHNET AUS DEN DATEN DES ADI DIE DATEN, DIE IN BETA EINGEGEBEN WERDEN
C    GERAET WIRD IN ZWEI HAELFTEN EINGEGEBEN, ALSO ZWEIMAL FUENF MAGNETE

      IMPLICIT NONE

      DOUBLE PRECISION B0,XLP,FB0N,FB0M,E,GAMMA
      DOUBLE PRECISION RHOPP,RHOPM,RHOMP1,RHOMP2,RHOMM
      DOUBLE PRECISION XLPP,XLPM,XLMP1,XLMP2,XLMM
      DOUBLE PRECISION PHIPP,PHIPM,PHIMM,PHIMP1,PHIMP2
      DOUBLE PRECISION PHI0,PHI1,PHI2,PHI3,PHI4,PHI5
      DOUBLE PRECISION PHIELMP1,PHIELMM,PHIELMP2,PHIELPM,PHIELPP
      DOUBLE PRECISION PHIERMP1,PHIERMM,PHIERMP2,PHIERPM,PHIERPP
      DOUBLE PRECISION XSTRAIGHT,XDRIFT,XLTOT,CLIGHT,EMASSGEV,PELEV

      DATA XSTRAIGHT/5.8D0/ !LAENGE DES GERADEN STUECKES IM RING
      DATA CLIGHT/2.99792458D8/
      DATA EMASSGEV/511003.3732832001D-9/

      E=GAMMA*511003.3732832001D-9  !GEV
      PELEV=DSQRT(  (E-EMASSGEV)*(E+EMASSGEV)  )*1.D9

      RHOPP = PELEV/(CLIGHT*B0)
      RHOPM =-RHOPP*FB0M
      RHOMM =-RHOPP*FB0N
      RHOMP1=-RHOMM*FB0M
      RHOMP2=-RHOMM*FB0M

      XLPP =XLP/8. !NUR HALBER WLS SOLL BETRACHTET WERDEN
      XLPM =XLPP
      XLMM =XLPP*FB0N
      XLMP1=XLMM/2.
      XLMP2=XLMM/2.
      XLTOT=XLPP+XLPM+XLMM+XLMP1+XLMP2
      XDRIFT=(XSTRAIGHT-2.*XLTOT)/2.

C     PHIMP1= 0.D0 -DASIN (DSIN( 0.D0 )-XLMP1/RHOMP1)  !1. MAGNET
C     PHIMM =PHIMP1-DASIN (DSIN(PHIMP1)- XLMM/ RHOMM)  !2. MAGNET
C     PHIMP2=PHIMM -DASIN (DSIN( PHIMM)-XLMP2/RHOMP2)  !3. MAGNET
C     PHIPM =PHIMP2-DASIN (DSIN(PHIMP2)- XLPM/ RHOPM)  !4. MAGNET
C     PHIPP =PHIPM -DASIN (DSIN( PHIPM)- XLPP/ RHOPP)  !5. MAGNET

      PHIMP1= XLMP1/RHOMP1 !1. MAGNET      !MAGNET SONST NICHT ABGEGLICHEN
      PHIMM =  XLMM/ RHOMM !2. MAGNET
      PHIMP2= XLMP2/RHOMP2 !3. MAGNET
      PHIPM = XLPM/ RHOPM  !4. MAGNET
      PHIPP = XLPP/ RHOPP  !5. MAGNET

      PHI0=0.0
      PHI1=PHI0+PHIMP1
      PHI2=PHI1+PHIMM
      PHI3=PHI2+PHIMP2
      PHI4=PHI3+PHIPM
      PHI5=PHI4+PHIPP


      PHIELMP1=DSIGN(PHI0,RHOMP1)
      PHIERMP1=PHIMP1-PHIELMP1

      PHIELMM =-PHIERMP1
      PHIERMM =PHIMM -PHIELMM
      PHIELMP2=-PHIERMM
      PHIERMP2=PHIMP2-PHIELMP2
      PHIELPM =-PHIERMP2
      PHIERPM =PHIPM -PHIELPM

      PHIELPP =-PHIERPM
      PHIERPP =PHIPP-PHIELPP

      OPEN(UNIT=10,FILE='ADI.STR',
     &       STATUS='NEW',FORM='FORMATTED')

      WRITE(10,*)'E, LAENGE DES GERADEN STUECKES:',SNGL(E),SNGL(XSTRAIGHT)
      WRITE(10,*)'B0, XLP:',SNGL(B0),SNGL(XLP)
      WRITE(10,*)'n, m:   ',SNGL(FB0N),SNGL(FB0M)

      WRITE(10,1001)PHIELMP1,RHOMP1
1001  FORMAT(' EL1  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,2001) PHIMP1,RHOMP1
2001  FORMAT(' BB1  DI',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,3001) PHIERMP1,RHOMP1
3001  FORMAT(' ER1  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')

      WRITE(10,1002)PHIELMM,RHOMM
1002  FORMAT(' EL2  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,2002) PHIMM,RHOMM
2002  FORMAT(' BB2  DI',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,3002) PHIERMM,RHOMM
3002  FORMAT(' ER2  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')

      WRITE(10,1003)PHIELMP2,RHOMP2
1003  FORMAT(' EL3  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,2003) PHIMP2,RHOMP2
2003  FORMAT(' BB3  DI',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,3003) PHIERMP2,RHOMP2
3003  FORMAT(' ER3  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')

      WRITE(10,1004)PHIELPM,RHOPM
1004  FORMAT(' EL4  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,2004) PHIPM,RHOPM
2004  FORMAT(' BB4  DI',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,3004) PHIERPM,RHOPM
3004  FORMAT(' ER4  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')

      WRITE(10,1005)PHIELPP,RHOPP
1005  FORMAT(' EL5  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,2005) PHIPP,RHOPP
2005  FORMAT(' BB5  DI',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,3005) PHIERPP,RHOPP
3005  FORMAT(' ER5  CO',E14.6,E14.6,'  0.000000E+00  0.000000E+00')
      WRITE(10,4000)XDRIFT
4000  FORMAT('   ## SD',E14.6)

      WRITE(10,*)
     &  '*** ACHTUNG: KANTENWINKEL GRAPHISCH KONTROLLIEREN ! ***'

      RETURN
      END
