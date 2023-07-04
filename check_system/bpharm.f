*CMZ :  2.44/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.01.45  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.05  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  17.08.37  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.38.06  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97

      SUBROUTINE BPHARM(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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

C
C     INPUT/OUTPUT COORDINATE SYSTEM: X LONG., Y VERTICAL
C     INTERNAL COORDINATE SYSTEM: Z LONG., Y VERTICAL
C
C     UNIT: METER AND TESLA
C

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpharm.
      include 'bpharm.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL,IPAR,N,NXY

      DOUBLE PRECISION XOFF
      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION X,Y,Z,BXX,BYY,BZZ,B0X,B0Y,B0Z
      DOUBLE PRECISION ZK

        DOUBLE PRECISION COSNZ,SINNXYX,SINHY,COSNXYY,COSHX,SINNZ
     &        ,COSNXYX,COSHY
        DOUBLE PRECISION SINNXYY,SINHX,yKs,XKC
        DOUBLE PRECISION ASINHY(NHARMP,NTRANSP),ACOSHX(NHARMP,NTRANSP)
        DOUBLE PRECISION ACOSHY(NHARMP,NTRANSP),ASINHX(NHARMP,NTRANSP)
        DOUBLE PRECISION YKC(NHARMP,NTRANSP),XKS(NHARMP,NTRANSP)
        DOUBLE PRECISION ACOSNZ(NHARMP),ASINNZ(NHARMP)
        DOUBLE PRECISION ACOSNXYX(NTRANSP),ASINNXYX(NTRANSP)
        DOUBLE PRECISION ACOSNXYY(NTRANSP),ASINNXYY(NTRANSP)
      DOUBLE PRECISION B0S,B0C

        COMPLEX*16 CZKN(NHARMP),CXN(NTRANSP),CYN(NTRANSP)
       CHARACTER(64) COMMENT

      DATA ICAL/0/

C--- INITIALIZATION

      IF (ICAL.EQ.0) THEN

         OPEN(UNIT=LUNPHFIT,FILE=FILEPHFIT,STATUS='OLD')

            READ(LUNPHFIT,'(A64)')COMMENT
            READ(LUNPHFIT,*)PERLENPH,PHASEPH
            READ(LUNPHFIT,*)XLENCPH,YLENSPH
            READ(LUNPHFIT,*)XPHMIN,XPHMAX
            READ(LUNPHFIT,*)YPHMIN,YPHMAX
            READ(LUNPHFIT,*)ZPHMIN,ZPHMAX
            READ(LUNPHFIT,*)NTRANS0,NTRANS,NTRANSD
            READ(LUNPHFIT,*)NHARM0,NHARM,NHARMD
          READ(LUNPHFIT,*)NPARPH

          DO IPAR=1,NPARPH
         READ(LUNPHFIT,*)PARPH(IPAR)
          ENDDO

            WRITE(LUNGFO,*)'     SUBROUTINE BPHARM:'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     Fitparameters read from file:'
            WRITE(LUNGFO,*)'     ',FILEPHFIT
            WRITE(LUNGFO,*)

            WRITE(LUNGFO,'(''      '',A64)')COMMENT
            WRITE(LUNGFO,*)
     &'     periodlength, phase, and scaling factor:'
            WRITE(LUNGFO,*)
     &'     (PERLENPH,PHASEPH,XLENCPH,YLENSPH)'
            WRITE(LUNGFO,*)'     ',PERLENPH
            WRITE(LUNGFO,*)'     ',PHASEPH
            WRITE(LUNGFO,*)'     ',XLENCPH
            WRITE(LUNGFO,*)'     ',YLENSPH
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     XPHMIN,XPHMAX:',XPHMIN,XPHMAX
            WRITE(LUNGFO,*)'     YPHMIN,YPHMAX:',YPHMIN,YPHMAX
            WRITE(LUNGFO,*)'     ZPHMIN,ZPHMAX:',ZPHMIN,ZPHMAX
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
     &'     NTRANS0,NTRANS,NTRANSD:',NTRANS0,NTRANS,NTRANSD
            WRITE(LUNGFO,*)
     &'     NHARM0,NHARM,NHARMD:   ',NHARM0,NHARM,NHARMD
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'     NPARPH:',NPARPH
            WRITE(LUNGFO,*)

          DO IPAR=1,NPARPH
         WRITE(LUNGFO,*)'          ',PARPH(IPAR)
          ENDDO

          IF (NPARPH.GT.NPARPHP) THEN
         WRITE(LUNGFO)
     &'*** ERROR IN BPHARM: DIMENSION NPARPHP IN BPHARM.CMN EXCEEDED  ***'
         WRITE(LUNGFO)
         WRITE(6,*)
     &'*** ERROR IN BPHARM: DIMENSION NPARPHP IN BPHARM.CMN EXCEEDED  ***'
         WRITE(6,*)
         WRITE(6,*) '--- PROGRAM ABORTED DUE TO ERROR ---'
         STOP
          ENDIF

          IF (NTRANS.GT.NTRANSP) THEN
         WRITE(LUNGFO)
     &'*** ERROR IN BPHARM: DIMENSION NTRANSP IN BPHARM.CMN EXCEEDED  ***'
         WRITE(LUNGFO)
         WRITE(6,*)
     &'*** ERROR IN BPHARM: DIMENSION NTRANSP IN BPHARM.CMN EXCEEDED  ***'
         WRITE(6,*)
         WRITE(6,*) '--- PROGRAM ABORTED DUE TO ERROR ---'
         STOP
          ENDIF

         CLOSE(LUNPHFIT)

        XOFF=PHASEPH*PERLENPH
        ZK=2.D0*PI1/PERLENPH
        IF (XLENCPH.NE.0.D0) THEN
           XKCPH=2.D0*PI1/XLENCPH
        ELSE
           XKCPH=0.D0
        ENDIF
        IF (YLENSPH.NE.0.D0) THEN
           YKSPH=2.D0*PI1/YLENSPH
        ELSE
           YKSPH=0.D0
        ENDIF
        XKC=XKCPH
        YKS=YKSPH

         AXOUT=0.D0
         AYOUT=0.D0
         AZOUT=0.D0

         ICAL=1

      ENDIF !ICAL

C --- CHANGE COORDINATE SYSTEMS

      X=-ZIN
      Y=YIN
      Z=XIN+XOFF

C--- MAGNETIC FIELD

      BXOUT=0.D0
      BYOUT=0.D0
      BZOUT=0.D0
      B0X=0.0
      B0Y=0.0
      B0Z=0.0


        CZKN(1)=CDEXP(DCMPLX(0.D0,ZK*Z))
        DO N=2,NHARM
                CZKN(N)=CZKN(N-1)*CZKN(1)
        ENDDO   !NHARM

        CXN(1)=CDEXP(DCMPLX(0.D0,XKC*X))
        CYN(1)=CDEXP(DCMPLX(0.D0,YKS*Y))
        DO NXY=2,NTRANS
                CXN(NXY)=CXN(NXY-1)*CXN(1)
                CYN(NXY)=CYN(NXY-1)*CYN(1)
        ENDDO


      IF (NHARM0.GT.0) THEN
        DO N=NHARM0,NHARM,NHARMD
                ACOSNZ(N)=DREAL(CZKN(N))
                ASINNZ(N)=DIMAG(CZKN(N))
        ENDDO   !NHARM
      ENDIF

      IF (NTRANS0.GT.0) THEN
        DO NXY=NTRANS0,NTRANS,NTRANSD
                ACOSNXYX(NXY)=DREAL(CXN(NXY))
                ASINNXYX(NXY)=DIMAG(CXN(NXY))
                ACOSNXYY(NXY)=DREAL(CYN(NXY))
                ASINNXYY(NXY)=DIMAG(CYN(NXY))
        ENDDO
      ENDIF


      IF (NPARPH.GT.3) THEN

        DO NXY=NTRANS0,NTRANS,NTRANSD
        DO N=NHARM0,NHARM,NHARMD

                YKC(N,NXY)=DSQRT((ZK*N)**2+(NXY*XKC)**2)
                XKS(N,NXY)=DSQRT((ZK*N)**2+(NXY*YKS)**2)
                ASINHY(N,NXY)=DSINh(YKc(N,NXY)*y)
                ASINHX(N,NXY)=DSINh(XKs(N,NXY)*x)
                ACOSHX(N,NXY)=DSQRT(1.D0+ASINHX(N,NXY)**2)
                ACOSHY(N,NXY)=DSQRT(1.D0+ASINHY(N,NXY)**2)

        ENDDO   !NHARM
        ENDDO   !NTRANS

      IPAR=3
        DO NXY=NTRANS0,NTRANS,NTRANSD
        DO N=NHARM0,NHARM,NHARMD

                IPAR=IPAR+2

           B0C=PARPH(IPAR-1)
           B0S=PARPH(IPAR)

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


c achtung: einige Sinusterme nacheditieren   INCLUDE 'RED:VPOT-HARM_BFELD.FOR'

      bxx=(-(COSNZ*SINNXYX*SINHY*b0c*XKC*
     . nxy-YKC(N,nxy)*b0x-COSNXYY*COSHX*YKC(N,nxy
     . )*SINNZ*b0s))/YKC(N,nxy)
      byy=(XKS(N,nxy)*b0y-SINNXYY*SINNZ*SINHX*b0s*YKS*nxy
     &+COSNXYX*COSNZ*COSHY*
     . XKS(N,nxy)*b0c)/XKS(N,nxy)
      bzz=((COSNXYY*COSNZ*SINHX*b0s*ZK*n+
     . XKS(N,nxy)*b0z)*YKC(N,nxy)-COSNXYX*XKS(N,nxy)*sinnz
     . *SINHY*b0c*ZK*n)/(XKS(N,nxy)*YKC(N,nxy))


C --- CHANGE COORDINATE SYSTEMS

      BXOUT=BXOUT+BXX
      BYOUT=BYOUT+BYY
      BZOUT=BZOUT+BZZ

      ENDDO !NTRANS
      ENDDO !NHARM

      ENDIF !NPARPH.GT.3

C --- CHANGE COORDINATE SYSTEMS

      BXX= BZOUT+PARPH(3)
      BYY= BYOUT+PARPH(2)
      BZZ=-BXOUT-PARPH(1)

      BXOUT=BXX
      BYOUT=BYY
      BZOUT=BZZ

      RETURN
      END
