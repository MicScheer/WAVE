*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.03.55  by  Michael Scheer
*CMZ :  2.54/05 20/04/2005  09.45.17  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.14/02 19/04/2000  17.02.45  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  11.45.40  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.01  by  Michael Scheer
*CMZ : 00.01/09 01/12/95  15.36.40  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.48.36  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.32  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.10  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BHELM(X,Y,Z,BX,BY,BZ)

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

C  SUBROUTINE BERECHNET B-FELD NUMERISCH FUER HELMHOLTZ-SPULEN
C     SPULEN LIEGEN PARALLEL ZUR Z,X-EBENE, PHI IST DER WINKEL
C     DES RADIUS-VEKTORS ZUR Z-ACHSE

*KEEP,bhelm.
      include 'bhelm.cmn'
*KEND.

      INTEGER MINT
      PARAMETER(MINT=1000)

      DOUBLE PRECISION DL
      DIMENSION DL(2,MINT,3)

      INTEGER ICAL,I,IC
      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,PI,PHI,CPHI,SPHI,DPHI,RX,RZ,RYP,RYM,R1P
     &  ,R1P3,R1M,R1M3,DBXP,DBYP,DBZP,DBXM,DBYM,DBZM

      DATA ICAL/0/
      DATA PI/3.141592653589793D0/


C--- ARRAYS FUER INTEGRATION AUFSETZEN

      IF (ICAL.NE.1) THEN

        IF(INTHELM.GT.MINT) STOP '*** SR BHELM: INTHELM.GT.MINT ***'

        DPHI=2.D0*PI/DFLOAT(INTHELM)
        DPHI=2.D0*PI/DFLOAT(INTHELM)

        DO IC=1,3
          P0HELM(2,IC)=R0HELM(IC)*0.5D0
          CURHELM(IC)=B0HELM(IC)*R0HELM(IC)*(DSQRT(1.25D0))**3/(4.D-7*PI)
        ENDDO

        DO I=1,INTHELM
          PHI=DFLOAT(I)*DPHI
          CPHI=DCOS(PHI)
          SPHI=DSIN(PHI)
          DO IC=1,3
            DL(1,I,IC)=R0HELM(IC)*(-SPHI)*DPHI*CURHELM(IC)*1.D-7
            DL(2,I,IC)=R0HELM(IC)*  CPHI *DPHI*CURHELM(IC)*1.D-7
          END DO
        END DO

        WRITE(16,*)
        WRITE(16,*)'      BHELM:'
        WRITE(16,*)
        WRITE(16,*)'      Radius[m], B0[T], I[A]:',
     &    (SNGL(R0HELM(1))),
     &    (SNGL(B0HELM(1))),
     &    (SNGL(CURHELM(1)))
        WRITE(16,*)'      X,Y,Z[m]:',(SNGL(P0HELM(IC,1)),IC=1,3)
        WRITE(16,*)'      Radius[m], B0[T], I[A]:',
     &    (SNGL(R0HELM(2))),
     &    (SNGL(B0HELM(2))),
     &    (SNGL(CURHELM(2)))
        WRITE(16,*)'      X,Y,Z[m]:',(SNGL(P0HELM(IC,2)),IC=1,3)
        WRITE(16,*)'      Radius[m], B0[T], I[A]:',
     &    (SNGL(R0HELM(3))),
     &    (SNGL(B0HELM(3))),
     &    (SNGL(CURHELM(3)))
        WRITE(16,*)'      X,Y,Z[m]:',(SNGL(P0HELM(IC,3)),IC=1,3)

        ICAL=1
      END IF

C--- INTEGRATION DER B-FELDER

      BX=0.D0
      BY=0.D0
      BZ=0.D0

      DO IC=1,3
        DO I=1,INTHELM
          RX=X-P0HELM(1,IC)-R0HELM(IC)*DSIN(DFLOAT(I)*DPHI)
          RZ=Z-P0HELM(3,IC)-R0HELM(IC)*DCOS(DFLOAT(I)*DPHI)
          RYP=Y-P0HELM(2,IC)
          RYM=Y+P0HELM(2,IC)

          R1P=1.D0/(DSQRT(RX*RX+RYP*RYP+RZ*RZ))
          R1P3=R1P*R1P*R1P
          R1M=1.D0/(DSQRT(RX*RX+RYM*RYM+RZ*RZ))
          R1M3=R1M*R1M*R1M

          DBXP=-RYP*DL(1,I,IC)*R1P3
          DBYP =(RX *DL(1,I,IC)-RZ*DL(2,I,IC))*R1P3
          DBZP= RYP*DL(2,I,IC)*R1P3

          DBXM=-RYM*DL(1,I,IC)*R1M3
          DBYM =(RX *DL(1,I,IC)-RZ*DL(2,I,IC))*R1M3
          DBZM= RYM*DL(2,I,IC)*R1M3


          BX=BX+DBXP+DBXM
          BY=BY+DBYP+DBYM
          BZ=BZ+DBZP+DBZM
        END DO
      END DO

      RETURN
      END
