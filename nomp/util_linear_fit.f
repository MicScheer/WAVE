*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  1.00/00 02/08/2002  19.48.06by  Michael Scheer
*CMZ : 00.02/04 10/02/97  13.35.14  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.23.52  by  Michael Scheer
*CMZ : 00.00/01 17/01/97  16.00.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_LINEAR_FIT
     &  (IFAIL,NPAR,PARAM,NDIMPOI,NPOI,NARG,NFUN,A,T,FUNDATA,CURDAT)
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

C     DOUBLE PRECISION SUBROUTINE FOR MULTIDIMENSIONAL LINEAR FIT

C     INPUT:
C       ------
C     NPAR:     NUMBER OF PARAMETERS TO FIT
C     NPOI:     NUMBER OF DATAPOINTS (NPOI MUST BE .LE. NDIMPOI)
C     NDIMPOI:    DIMENSION FOR NUMBER OF DATAPOINTS
C     NARG:     NUMBER OF ARGUMENTS FOR FUNCTION
C     NFUN:     DIMENSION OF FUNCTION
C     A:     WORKINGSPACE A(NPAR,NPAR)
C     T:     WORKINGSPACE T(NFUN,NPAR)
C     FUNDATA:    DATA POINTS FUNDATA(NARG+NFUN,NDIMPOI)
C            IF A FUNCTION VALUE IS 9999. DATA IS SKIPPED
C            NUMBER OF DATA MUST BE AT LEAST NPAR

C     OUTPUT:
C     -------
C     PARAM:       PARAMETERS TO BE FITTED
C     IFAIL:       FAILURE FLAG, ZERO IF EVERYTHING SEEMS TO BE OK

C     EXAMPLE:
C     --------
C     NPAR=4
C     NPOI=10
C     NARG=3
C     NFUN=2
C     FUNDATA(NARG+NFUN,NPOI):
C         FUNDATA(1,I)=X(I)
C         FUNDATA(2,I)=Y(I)
C         FUNDATA(3,I)=Z(I)
C         FUNDATA(4,I)=BX(I)
C         FUNDATA(5,I)=BY(I)
C

C     FIT SUCH THAT
C
C         dCHI2/dP1=0
C         dCHI2/dP2=0
C         dCHI2/dP3=0
C         dCHI2/dP4=0
C
C     WITH
C
C         DO I=1,NPOI
C          CHI2=CHI2+
C         +((P1*X(I)*Y(I)*SIN(  K*Z(I)+)-BX(IPOI))**2
C         +((P2*X(I)*Y(I)*COS(  K*Z(I)+)-BX(IPOI))**2
C         +((P3*X(I)*Y(I)*SIN(3*K*Z(I)+)-BY(IPOI))**2
C         +((P4*X(I)*Y(I)*COS(3*K*Z(I)+)-BY(IPOI))**2
C         ENDDO   !IPOI
C
C     AND P1=PARAM(1), P2=....
C

      IMPLICIT NONE

      INTEGER IFAIL,NPAR,NDIMPOI,NPOI,NARG,NFUN
      INTEGER IPAR,IPOI,IFUN
      INTEGER JPAR,JPOI,IDAT

      DOUBLE PRECISION PARAM,A,T,FUNDATA,CURDAT
      DIMENSION PARAM(NPAR),A(NPAR,NPAR)
     &           ,T(NFUN,NPAR)
     & ,FUNDATA(NARG+NFUN,NDIMPOI),CURDAT(NARG+NFUN)

      IF (NPOI.GT.NDIMPOI) STOP
     &'*** ERROR IN UTIL_LINEAR_FIT: DIMENSION NDIMPOI EXCEEDED ***'

      IFAIL=0

      JPOI=0
      DO IPOI=1,NPOI
      DO IFUN=1,NFUN
          IF (FUNDATA(NARG+IFUN,IPOI).NE.9999.) JPOI=JPOI+1
      ENDDO
      ENDDO

C ATTENTION: SINCE T IS ALSO USED AS WORKINGSPACE FOR F010
C          ITS DIMENSION MUST BE AT LEAST 2*NPAR, I.E. NPOI.GE.NPAR
      IF (JPOI.LT.NPAR) THEN
          IFAIL=9999
          RETURN
      ENDIF

      DO JPAR=1,NPAR
          PARAM(JPAR)=0.D0
      DO IPAR=1,NPAR
          A(IPAR,JPAR)=0.D0
      ENDDO
      ENDDO


C--- SET UP EQUATION SYSTEM


C -- SET UP INHOMOGENITY OF EQUATION SYSTEM

      DO IPOI=1,NPOI
C    T(IFUN,IPAR)=dChi2/dIPAR for each point and function
      DO IDAT=1,NARG+NFUN
          CURDAT(IDAT)=FUNDATA(IDAT,IPOI)
      ENDDO
      CALL UTIL_LINEAR_FIT_USER(NARG,NFUN,NPAR,CURDAT,T)
      DO IPAR=1,NPAR
      DO IFUN=1,NFUN

          PARAM(IPAR)=PARAM(IPAR)
     &                 +T(IFUN,IPAR)*CURDAT(NARG+IFUN)

      ENDDO
      ENDDO
      ENDDO

C -- SET UP MATRIX OF EQUATION SYSTEM

      DO IPOI=1,NPOI
      DO IDAT=1,NARG+NFUN
          CURDAT(IDAT)=FUNDATA(IDAT,IPOI)
      ENDDO
      CALL UTIL_LINEAR_FIT_USER(NARG,NFUN,NPAR,CURDAT,T)
      DO IFUN=1,NFUN
      DO JPAR=1,NPAR
          DO IPAR=1,NPAR
         A(IPAR,JPAR)=A(IPAR,JPAR)
     &                      +(T(IFUN,IPAR)*T(IFUN,JPAR))
          ENDDO
      ENDDO
      ENDDO
      ENDDO

C--- SOLVE EQUATION SYSTEM WITH CERN-ROUTINE F010

      CALL DEQN(NPAR,A,NPAR,T,IFAIL,1,PARAM)

      RETURN
      END
