*CMZ :  3.00/00 11/03/2013  15.09.17  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/13 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.53/05 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.47/14 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.47/12 30/06/2003  15.30.28  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  16.15.48  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.36  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.03  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  10.46.38  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.09  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.10  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WZZPYYP
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

C--- CALCULATES z,z',y,y' in phase space relative to closed orbit from
c    beta-functions
c    K. Wille: Kapitel 3.11


C    RESULTS:

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,phasetrack.
      include 'phasetrack.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEND.


      DOUBLE PRECISION ALPHA0(2),BETA0(2),X0(2),XP0(2)
      DOUBLE PRECISION ALPHA,BETA,PSI,CPSI,SPSI,RQ,RM

      INTEGER I,IC

      X0(1)=PHTRZ0
      XP0(1)=PHTRZP0
      X0(2)=PHTRY0
      XP0(2)=PHTRYP0

      ALPHA0(1)=-WBETA(3,1)/2.d0
      ALPHA0(2)=-WBETA(5,1)/2.d0
      BETA0(1)=WBETA(2,1)
      BETA0(2)=WBETA(4,1)

      DO I=1,NCO

        ALPHA=-WBETA(3,I)/2.d0
        BETA=  WBETA(2,I)
        PSI=   WBETA(8,I)
        CPSI=COS(PSI)
        SPSI=SIN(PSI)
        RQ=SQRT(BETA/BETA0(1))
        RM=SQRT(BETA*BETA0(1))
        WBZZPYYP(1,I) = RQ * (CPSI+ALPHA0(1)*SPSI) * X0(1) + RM*SPSI*XP0(1)
        WBZZPYYP(2,I) =
     &    ((ALPHA0(1)-ALPHA)*CPSI - (1.D0+ALPHA0(1)*ALPHA)*SPSI) * X0(1)/RM
     &    + (CPSI-ALPHA*SPSI) * XP0(1)/RQ

        ALPHA=-WBETA(5,I)/2.d0
        BETA=  WBETA(4,I)
        PSI=   WBETA(9,I)
        CPSI=COS(PSI)
        SPSI=SIN(PSI)
        RQ=SQRT(BETA/BETA0(2))
        RM=SQRT(BETA*BETA0(2))
        WBZZPYYP(3,I) = RQ * (CPSI+ALPHA0(2)*SPSI) * X0(2) + RM*SPSI*XP0(2)
        WBZZPYYP(4,I) =
     &    ((ALPHA0(2)-ALPHA)*CPSI - (1.D0+ALPHA0(2)*ALPHA)*SPSI) * X0(2)/RM
     &    + (CPSI-ALPHA*SPSI) * XP0(2)/RQ

      ENDDO

      OPEN(UNIT=LUNZZPYYP,FILE=FILEZZPYYP,STATUS='unknown')

        DO I=1,NCO
          WRITE(LUNZZPYYP,'(5((1PE15.5)))')WBETA(1,I),(WBZZPYYP(IC,I),IC=1,4)
        ENDDO

      CLOSE(LUNZZPYYP)

      RETURN
      END
