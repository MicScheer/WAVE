*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  3.01/00 03/07/2013  15.33.33  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.09.17  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.44  by  Michael Scheer
*CMZ :  2.66/13 20/07/2010  18.26.35  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.56/00 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.53/05 15/02/2005  12.34.54  by  Michael Scheer
*CMZ :  2.52/15 05/01/2005  16.32.31  by  Michael Scheer
*CMZ :  2.47/12 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  16.15.48  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  2.13/11 22/03/2000  14.37.49  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.25.45  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.08.11  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.07  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WI2I4I5
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

C--- CALCULATE RADIATION INTEGRALS I2,I4,I5 NUMERICALLY IN THE
C    SYSTEM OF THE LONGITUDINAL COORDINATE S ALONG THE REFERENCE ORBIT
C    AND FROM THE ANALYTICAL ANSATZ FOR THE COORDINATE X ALONG THE
C    DEVICE AXIS


      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wbetaf90.
      include 'wbetaf90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

      INTEGER IP,IC1,IC2

      DOUBLE PRECISION EC1
      DOUBLE PRECISION BY,R1,R12,R13,R13A,R
      DOUBLE PRECISION X,Z,VX,VZ,XOLD,DX
      DOUBLE PRECISION BETA,BETAP,ETA,ETAP,GRAD
      DOUBLE PRECISION WEMITOT,WNOR,HRHO3,DWI4,WI3,SIGEERING,SIGEETOT
      DOUBLE PRECISION ANBETA,ANBETAP,ANETA,ANETAP,ANHRHO3,ANWI5,ANWI52,ANWEMTOT
      DOUBLE PRECISION AB,BB,AI,BI,AII,BII,CII,ANETAI,AI2OPT,BETUNI,BET2UNI
      DOUBLE PRECISION BETOPT,BET2OPT,DISOPT,ANWI5OP,ANWI52OP,ANEMIOP,ANEMI2OP
      DOUBLE PRECISION S,HH,ANHH
      REAL*4 TEMP,XFILL,XI,XE
      REAL*8 WEIGHT

      IF(IHBETA.NE.0) THEN

        XI=WBETA(1,1)-DS0/2.
        XE=WBETA(1,NCO)+DS0/2.
        call hbook1m(IDHS,'CHROMATIC VARIABLE H(s)',NCO/IHTRSMP+1,XI,XE,VMX)
        call hbook1m(IDHSR3,'H(s)/ABS(RHO**3)',NCO/IHTRSMP+1,XI,XE,VMX)

        XI=WSXYZ(1,1)-DS0/2.
        XE=WSXYZ(1,NCO)+DS0/2.
        call hbook1m(IDHX,'CHROMATIC VARIABLE H(x)',NCO/IHTRSMP+1,XI,XE,VMX)
        call hbook1m(IDHXR3,'H(x)/ABS(RHO**3)',NCO/IHTRSMP+1,XI,XE,VMX)

      ENDIF !IHBETA

      EC1=CLIGHT1/EMOM

C--- INITIALIZE INTEGRATION

      WI2=0.0
      WI3=0.0
      WI4=0.0
      WI5=0.0
      ANWI5=0.0
      WI524=0.0

      AB=0.0
      BB=0.0
      AI=0.0
      BI=0.0
      AII=0.0
      BII=0.0
      CII=0.0

      XOLD=WSXYZ(1,1)-(WSXYZ(1,2)-WSXYZ(1,1))

C--- INTEGRATION

      DO IP=1,NCO

        X=WTRA(1,1,IP)
        Z=WTRA(3,1,IP)
        S=WBETA(1,IP)

        DX=X-XOLD

        VX=WTRA(1,2,IP)
        VZ=WTRA(3,2,IP)

        BY=WTRA(2,3,IP)

        R1=BY*EC1 !POSITIV FOR POSITIV FIELDS
        R12=R1**2
        R13=R1**3
        R13A=DABS(R13)

        IF (DABS(R1).GT.1.D-20) THEN
          R=1.D0/R1
        ELSE
          R=0
        ENDIF

        BETA=WBETA(2,IP)
        BETAP=WBETA(3,IP)
        ETA=WBETA(6,IP)
        ETAP=WBETA(7,IP)

        ANBETA=BETFUN+X**2/BETFUN
        ANBETAP=2.D0*X/BETFUN
        ANETA=-Z+DISP0    !ONLY CORRECT IF DDISP0=0
        ANETAI=-Z
        ANETAP=-VZ/VX

        GRAD=WBETAK(2,IP)  !ACTUALLY GRAD=Ky*RHO**2 (SEE FURTHER DOWN)

        HH=(ETA**2+(BETA*ETAP-0.5D0*BETAP*ETA)**2)/BETA
        HRHO3=HH*R13A
        ANHH=(ANETA**2+(ANBETA*ANETAP-0.5D0*ANBETAP*ANETA)**2)/ANBETA
        ANHRHO3=ANHH*R13A

        WI2=WI2+R12*DS0
        DWI4=(1.D0*R13-2.D0*GRAD*R1)*ETA*DS0
        WI3=WI3+ABS(R1)**3*DS0
        WI4=WI4+DWI4
        WI5=WI5+HRHO3*DS0
        ANWI5=ANWI5+ANHRHO3*DX

        AB=(ANETA-X*ANETAP)**2
        BB=ANETAP**2
        AI=AI+AB*R13A*DX
        BI=BI+BB*R13A*DX

        AII=AII+R13A*DX
        BII=BII+(ANETAI-X*ANETAP)*R13A*DX
        CII=CII+(ANETAI-X*ANETAP)**2*R13A*DX

        XFILL=S
        WEIGHT=HH
        IF (IHBETA.NE.0) CALL hfillm(IDHS,XFILL,0.,WEIGHT)
        WEIGHT=HRHO3
        IF (IHBETA.NE.0) CALL hfillm(IDHSR3,XFILL,0.,WEIGHT)
        XFILL=X
        WEIGHT=ANHH
        IF (IHBETA.NE.0) CALL hfillm(IDHX,XFILL,0.,WEIGHT)
        WEIGHT=ANHRHO3
        IF (IHBETA.NE.0) CALL hfillm(IDHXR3,XFILL,0.,WEIGHT)

        XOLD=X

      ENDDO !NCO


      IF (WI2.NE.WI4) WI524=WI5/(WI2-WI4)
      IF (WI2.NE.0.0) ANWI52=ANWI5/WI2

      WEMITOT=(DI5RING+WI5)/(DI2RING+WI2-DI4RING-WI4)
      ANWEMTOT=(DI5RING+ANWI5)/(DI2RING+WI2)
      WNOR=DI5RING/(DI2RING-DI4RING)

C--- OPTIMIZATION OF EXTERNAL DISPERSION AND BETA-FUNCTION
C     (USING ANALYTICAL ANSATZ)

      IF (WI2.GT.1D-10) THEN

        DISOPT=-BII/AII
        AI2OPT=AII*DISOPT**2+2.D0*BII*DISOPT+CII !AI FUER OPTIMALE DISPERSION

        BETOPT=DSQRT(AI/BI) !DIESER WERT DER BETATRON-FUNKTION MINIMIERT I5
        BET2OPT=DSQRT((AI2OPT)/BI) !DIESER WERT DER BETATRON-FUNKTION MINIMIERT I5
        ANWI5OP=1./BETOPT*AI+BETOPT*BI
        ANWI52OP=2.D0*DSQRT((CII-BII*BII/AII)*BI)

        BETUNI=-9999.

        IF(WI2.GT.1D-10.AND.ANWI5OP/WI2.LT.DI5RING/DI2RING) THEN
          BETUNI=DI5RING*WI2/(2.D0*DI2RING*BI)
     &      +DSQRT((DI5RING*WI2/(2.D0*DI2RING*BI))**2-AI/BI)
        ENDIF

        BET2UNI=-9999.

        IF(WI2.GT.1D-10.AND.ANWI52OP/WI2.LT.DI5RING/DI2RING) THEN
          BET2UNI=DI5RING*WI2/(2.D0*DI2RING*BI)
     &      +DSQRT((DI5RING*WI2/(2.D0*DI2RING*BI))**2-AI2OPT/BI)
        ENDIF

        ANEMIOP=(DI5RING+ANWI5OP)/(DI2RING+WI2)
        ANEMI2OP=(DI5RING+ANWI52OP)/(DI2RING+WI2)

      ENDIF !WI2

C-- I3

      DI3RING=2.0D0*PI1/RDIPOL**2
      SIGEERING=SQRT(CQ1*DMYGAMMA**2*DI3RING/(2.0D0*DI2RING+DI4RING))
      SIGEETOT=
     &  SQRT(CQ1*DMYGAMMA**2*(DI3RING+WI3)/(2.0D0*(DI2RING+WI2)+DI4RING))

C--- RESULTS


      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'RADIUS OF DIPOLMAGNETS:   ',SNGL(RDIPOL)
      WRITE(LUNGFO,*)'CIRCUMFERENCE OF THE RING:',SNGL(UMFANG)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'DISPERSION-FUNCTION AND DERIVATIVE FOR WLS:'
      WRITE(LUNGFO,*)SNGL(DISP0),SNGL(DDISP0)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'MOMENTUM-DEVIATION DISPERSION CALCULATIONS:',SNGL(DELGAM)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'INTEGRAL I2 OF RING:        ',SNGL(DI2RING)
      WRITE(LUNGFO,*)'INTEGRAL I3 OF RING:        ',SNGL(DI3RING)
      WRITE(LUNGFO,*)'INTEGRAL I4 OF RING:        ',SNGL(DI4RING)
      WRITE(LUNGFO,*)'INTEGRAL I5 OF RING:        ',SNGL(DI5RING)
      WRITE(LUNGFO,*)

      IC1=NCO/2
      IC2=IC1+2

      IF (
     &    DABS((WBETA(2,IC1)+WBETA(2,IC2))/2.D0/BETFUN-1.D0).GT.0.01
     &    .OR.
     &    WBETA(3,IC1)*WBETA(3,IC2).GT.0.0) THEN

        WRITE(LUNGFO,*)
        WRITE(6,*)'*** WARNING IN WI2I4I5 ***'
        WRITE(6,*)
        WRITE(6,*)
     &    'BETA-FUNCTION FOR ANALYTICAL APPROXIMATION AND NUMERICAL CALCULATION NOT COMPATIBLE'
        WRITE(6,*)'CHECK BETFUN, BETAH, BETAPH IN NAMELIST DEPOLA'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN WI2I4I5 ***'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    'BETA-FUNCTION FOR ANALYTICAL APPROXIMATION AND NUMERICAL CALCULATION NOT COMPATIBLE'
        WRITE(LUNGFO,*)'CHECK BETFUN, BETAH, BETAPH IN NAMELIST DEPOLA'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
      ENDIF

      IF (DDISP0.NE.0.0) THEN

        WRITE(LUNGFO,*)
        WRITE(6,*)'*** WARNING IN WI2I4I5 ***'
        WRITE(6,*)
        WRITE(6,*)
     &    'EXTERNAL DISPERSION NOT CONSTANT, I.E. DDISP0 IN NAMELIST DEPOLA NOT ZERO'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(LUNGFO,*)'*** WARNING SR WI2I4I5 ***'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    'EXTERNAL DISPERSION NOT CONSTANT, I.E. DDISP0 IN NAMELIST DEPOLA NOT ZERO'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
      ENDIF

      IF (DABS(XSTART+XSTOP).GT.1.D-10) THEN
        WRITE(LUNGFO,*)
        WRITE(6,*)'*** WARNING IN WI2I4I5 ***'
        WRITE(6,*)
        WRITE(6,*)
     &    'CONSIDERED SETUP NOT CENTERED AROUND ORIGIN'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(LUNGFO,*)'*** WARNING SR WI2I4I5 ***'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    'CONSIDERED SETUP NOT CENTERED AROUND ORIGIN'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
      ENDIF
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'EPS0H:',EPS0H
      WRITE(LUNGFO,*)'EPS0V:',EPS0V
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'RESULTS FROM ANALYTICAL ANSATZ:'
      WRITE(LUNGFO,*)'---------------------------------'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'BETFUN, BETFUNV:',SNGL(BETFUN),SNGL(BETFUNV)
      WRITE(LUNGFO,*)'BETAH, BETAPH:',SNGL(BETAH),SNGL(BETAPH)
      WRITE(LUNGFO,*)'BETAV, BETAPV:',SNGL(BETAV),SNGL(BETAPV)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '(ONLY CORRECT IF ORBIT IS PLANAR, WLS HAS NO TRANSVERSAL GRADIENT,'
      WRITE(LUNGFO,*)
     &  'WLS IS SYMMETRIC WITH RESPECT TO ORIGIN OR COORDINATE-SYSTEM'
      WRITE(LUNGFO,*)'CENTERED IN STRAIGHT SECTION, (VZ/VX)**2 << 1,'
      WRITE(LUNGFO,*)
     &  'AND DERIVATION OF EXTERNAL DISPERSION VANISHES OUTSIDE WLS)'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'BETA-FUNCTION AT WLS CENTER, I5, I5/I2 of WLS:'
      TEMP=1.E30
      IF(WI2.GT.1D-10) TEMP=SNGL(ANWI5/WI2)
      WRITE(LUNGFO,*)SNGL(BETFUN),SNGL(ANWI5),TEMP
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'TOTAL EMITTANCE CHANGE:',SNGL(ANWEMTOT/WNOR)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  'OPTIMAL BETA-FUNCTION, CORRESPONDING I5 AND TOTAL EMITTANCE CHANGE FOR ACTUAL EXTERNAL DISPERSION :'
      WRITE(LUNGFO,*)SNGL(BETOPT),SNGL(ANWI5OP),SNGL(ANEMIOP/WNOR)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'OPTIMAL EXTERNAL DISPERSION:',SNGL(DISOPT)
      WRITE(LUNGFO,*)
     &'OPTIMAL BETA-FUNCTION, CORRESPONDING I5 AND TOTAL EMITTANCE CHANGE FOR OPTIMAL EXTERNAL DISPERSION :'
         WRITE(LUNGFO,*)SNGL(BET2OPT),SNGL(ANWI52OP),SNGL(ANEMI2OP/WNOR)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'NEUTRAL BETAFUNCTION:                 ',SNGL(BETUNI)
      WRITE(LUNGFO,*)'NEUTRAL BETAFUNCTION FOR OPT. DISP.:  ',SNGL(BET2UNI)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'ENERGY-SPREAD OF RING:',SNGL(SIGEERING)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'RESULTS OF NUMERICAL CALCULATION:'
      WRITE(LUNGFO,*)'---------------------------------'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'I2, I3, I4, I5, I5/(I2-I4) of WLS:'
      WRITE(LUNGFO,*)SNGL(WI2),SNGL(WI3),SNGL(WI4),SNGL(WI5),SNGL(WI524)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'TOTAL EMITTANCE CHANGE:',SNGL(WEMITOT/WNOR)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'TOTAL ENERGY-SPREAD:',SNGL(SIGEETOT)
      WRITE(LUNGFO,*)'CHANGE TOTAL ENERGY-SPREAD:',SNGL(SIGEETOT/SIGEERING)
      WRITE(LUNGFO,*)

      RETURN
      END
