*CMZ :  2.44/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.09.10  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  15.40.14  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.36.40  by  Michael Scheer
*-- Author :    Michael Scheer   23/01/97

        SUBROUTINE BPHARM_USER(NARG,NFUN,NPAR,CURDAT,T)
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

*KEEP,bpharm.
      include 'bpharm.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL
      INTEGER NARG,NFUN,NPAR,IPAR
      INTEGER N,NXY
      DOUBLE PRECISION CURDAT,T
      DIMENSION CURDAT(NARG+NFUN),T(NFUN,NPAR)
      DOUBLE PRECISION TT(3,2)
      DOUBLE PRECISION COSNZ,SINNXYX,SINHY,COSNXYY,COSHX,SINNZ
     &        ,COSNXYX,COSHY
      DOUBLE PRECISION SINNXYY,SINHX,yKs
     &        ,XKS(NHARMP,NTRANSP),YKC(NHARMP,NTRANSP),XKC
      DOUBLE PRECISION ZK,X,Y,Z
      DOUBLE PRECISION ASINHY(NHARMP,NTRANSP),ACOSHX(NHARMP,NTRANSP)
      DOUBLE PRECISION ACOSHY(NHARMP,NTRANSP),ASINHX(NHARMP,NTRANSP)
      DOUBLE PRECISION ACOSNZ(NHARMP),ASINNZ(NHARMP)
      DOUBLE PRECISION ACOSNXYX(NTRANSP),ASINNXYX(NTRANSP)
      DOUBLE PRECISION ACOSNXYY(NTRANSP),ASINNXYY(NTRANSP)

      COMPLEX*16 CZKN(NHARMP),CXN(NTRANSP),CYN(NTRANSP)
      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN
         ZK=2.D0*PI1/PERLENPH
                XKC=XKCPH
                YKS=YKSPH
         ICAL=1
      ENDIF !ICAL

      x=CURdat(1)
      y=CURdat(2)
      z=CURdat(3)


      CZKN(1)=CDEXP(DCMPLX(0.D0,ZK*Z))
      DO N=2,NHARM
         CZKN(N)=CZKN(N-1)*CZKN(1)
      ENDDO !NHARM

      CXN(1)=CDEXP(DCMPLX(0.D0,XKC*X))
      CYN(1)=CDEXP(DCMPLX(0.D0,YKS*Y))
      DO NXY=2,NTRANS
         CXN(NXY)=CXN(NXY-1)*CXN(1)
         CYN(NXY)=CYN(NXY-1)*CYN(1)
      ENDDO

      T(1,1)=1.D0
      T(1,2)=0.D0
      T(1,3)=0.D0
      T(2,1)=0.D0
      T(2,2)=1.D0
      T(2,3)=0.D0
      T(3,1)=0.D0
      T(3,2)=0.D0
      T(3,3)=1.D0

      IF(NHARM0.GT.0) THEN
      DO N=NHARM0,NHARM,NHARMD
         ACOSNZ(N)=DREAL(CZKN(N))
         ASINNZ(N)=DIMAG(CZKN(N))
      ENDDO !NHARM
      ENDIF

      IF (NTRANS0.GT.0) THEN
      DO NXY=NTRANS0,NTRANS,NTRANSD
         ACOSNXYX(NXY)=DREAL(CXN(NXY))
         ASINNXYX(NXY)=DIMAG(CXN(NXY))
         ACOSNXYY(NXY)=DREAL(CYN(NXY))
         ASINNXYY(NXY)=DIMAG(CYN(NXY))
      ENDDO
      ENDIF

      IF (NPAR.GT.3) THEN

      DO NXY=NTRANS0,NTRANS,NTRANSD
      DO N=NHARM0,NHARM,NHARMD

         YKC(N,NXY)=DSQRT((ZK*N)**2+(NXY*XKC)**2)
         XKS(N,NXY)=DSQRT((ZK*N)**2+(NXY*YKS)**2)
         ASINHY(N,NXY)=DSINh(YKC(N,NXY)*y)
         ASINHX(N,NXY)=DSINh(XKS(N,NXY)*x)
         ACOSHX(N,NXY)=DSQRT(1.D0+ASINHX(N,NXY)**2)
         ACOSHY(N,NXY)=DSQRT(1.D0+ASINHY(N,NXY)**2)

      ENDDO !NHARM
      ENDDO !NTRANS

      IPAR=3
      DO NXY=NTRANS0,NTRANS,NTRANSD
      DO N=NHARM0,NHARM,NHARMD

         IPAR=IPAR+2

         SINHY=ASINHY(N,NXY)
         COSHX=ACOSHX(N,NXY)
         COSHY=ACOSHY(N,NXY)
         SINHX=ASINHX(N,NXY)
         COSNZ=ACOSNZ(N)
         SINNZ=ASINNZ(N)
         COSNXYX=ACOSNXYX(NXY)
         SINNXYX=ASINNXYX(NXY)
         COSNXYY=ACOSNXYY(NXY)
         SINNXYY=ASINNXYY(NXY)

c     INCLUDE 'RED:VPOT-HARM.FOR'

      tt(1,1)=(-COSNZ*SINNXYX*SINHY*XKC*nxy
     . )/YKC(N,nxy)
      tt(1,2)=COSNXYY*COSHX*SINNZ
      tt(2,1)=COSNXYX*COSNZ*COSHY
      tt(2,2)=(-SINNXYY*SINNZ*SINHX*YKS*nxy
     . )/XKS(N,nxy)
      tt(3,1)=(-COSNXYX*SINNZ*SINHY*ZK*n)/
     . YKC(N,nxy)
      tt(3,2)=(COSNXYY*COSNZ*SINHX*ZK*n)/
     . XKS(N,nxy)

      IF (YKS.EQ.0.D0) THEN
          TT(1,2)=0.D0
          TT(2,2)=0.D0
          TT(3,2)=0.D0
      ENDIF

         T(1,IPAR-1)=TT(1,1)
         T(1,IPAR)  =TT(1,2)
         T(2,IPAR-1)=TT(2,1)
         T(2,IPAR)  =TT(2,2)
         T(3,IPAR-1)=TT(3,1)
         T(3,IPAR)  =TT(3,2)

      ENDDO !NHARM
      ENDDO !NTRANS

      ENDIF !NPAR.GT.3

      RETURN
      END
