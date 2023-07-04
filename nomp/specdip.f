*CMZ :  4.01/03 12/06/2023  11.25.49  by  Michael Scheer
*CMZ :  3.02/03 07/11/2014  15.50.42  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  15.33.54  by  Michael Scheer
*CMZ :  2.57/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.57/04 01/02/2006  15.08.35  by  Michael Scheer
*CMZ :  2.34/05 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.33/03 04/05/2001  10.57.10  by  Michael Scheer
*CMZ :  2.33/02 03/05/2001  17.23.28  by  Michael Scheer
*CMZ :  2.33/01 03/05/2001  14.00.28  by  Michael Scheer
*CMZ :  2.33/00 02/05/2001  12.12.06  by  Michael Scheer
*CMZ :  2.31/00 24/04/2001  11.13.55  by  Michael Scheer
*CMZ :  2.20/01 12/12/2000  10.01.31  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  17.29.14  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.13/03 10/01/2000  17.17.14  by  Michael Scheer
*CMZ :  2.00/00 15/12/98  14.29.27  by  Michael Scheer
*CMZ : 00.02/02 15/01/97  15.18.41  by  Michael Scheer
*CMZ : 00.00/07 18/05/94  14.54.15  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.12  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.24  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SPECDIP
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- FILL ARRAY SPEC WITH DIPOL SPECTRUM FUNCTION

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,specdip.
      include 'specdip.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ISOUR,IFREQ,IOBSV,IDIP,IEPS,ical

      DOUBLE PRECISION RX,RY,RZ,RR,PSI,DFDTDP,PAR,PER,POWR
      DOUBLE PRECISION BS,BRX,BRY,BRZ,BRN,RN,BX,BY,BZ,VN
      DOUBLE PRECISION OMEGAC,Y
      DOUBLE PRECISION DX2,DZY2,EPS(6),DRRED,ANS,OM,
     &  APERANG,APERV(3),APERCORR,const

      COMPLEX*16 APOL,AFREQ(3),EXPOM
      COMPLEX*16 APOLH,APOLR,APOLL,APOL45

      REAL*8 STOK1,STOK2,STOK3,STOK4

      data ical/0/

      if (ical.eq.0) then
        if (iphase.ne.0) then
          write(lungfo,*)'*** WARNING IN SPECDIP: Calculation for B-Field-Amplitude for phase propagation not yet implement ***'
          write(6,*)'*** WARNING IN SPECDIP: Calculation for B-Field-Amplitude for phase propagation not yet implement ***'
        endif
        CONST=sqrt(3.0D0)/2.0D0/pi1*ALPHA1/ECHARGE1/emassg1/banwid*0.001d0/1.0d6
c        CONST=2.457d13
        const=const*pinw/pincen(1)*1.0d3*dmyenergy
        ical=1
      endif

      IF (APERTHICK.GT.0.0D0) THEN

        APERV(2)=SIN(APERVANG)
        APERV(3)=SIN(APERHANG)
        APERV(1)=SQRT(1.0D0-(APERV(2)**2+APERV(3)**2))

      ELSE

        APERV=0.0D0
        APERHANG=0.0D0
        APERVANG=0.0D0

      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR SPECDIP CALLED:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
     &  '     cut-off parameter (SPECCUT): ',SNGL(SPECCUT)
      WRITE(LUNGFO,*)

      if (ispecdip.eq.2) then
        WRITE(LUNGFO,*)
     &    '     Due to ISPECDIP=2, G1 is calculated instead of photon flux'
      else
        WRITE(LUNGFO,*)
     &    '     Numer of dipole, field (proj.), bending radius ,Ec [eV], cut-off,'
        WRITE(LUNGFO,*)
     &    '     X,Y,Z-Pos. [m], velocity vector (norm.), proj. mag. field vector [T]:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     thickness of aperture pinhole [m]:',sngl(aperthick)
        WRITE(LUNGFO,*)'     hor. and vert. angle of aperture [rad]:',
     &    sngl(aperhang),sngl(apervang)
        WRITE(LUNGFO,*)
      endif

      DO IDIP=1,NDIP

        VN=SQRT(VXDIP(IDIP)*VXDIP(IDIP)+VYDIP(IDIP)*VYDIP(IDIP)
     &    +VZDIP(IDIP)*VZDIP(IDIP))

        IF (VN.EQ.0.D0) THEN
          VXDIP(IDIP)=1.D0
          VN=1.D0
        ENDIF   !VN

        VXDIP(IDIP)=VXDIP(IDIP)/VN
        VYDIP(IDIP)=VYDIP(IDIP)/VN
        VZDIP(IDIP)=VZDIP(IDIP)/VN

        BS=BXDIP(IDIP)*VXDIP(IDIP)
     &    +BYDIP(IDIP)*VYDIP(IDIP)
     &    +BZDIP(IDIP)*VZDIP(IDIP)

        BXDIP(IDIP)=BXDIP(IDIP)-BS*VXDIP(IDIP)
        BYDIP(IDIP)=BYDIP(IDIP)-BS*VYDIP(IDIP)
        BZDIP(IDIP)=BZDIP(IDIP)-BS*VZDIP(IDIP)

        B0DIP(IDIP)=SQRT(
     &    BXDIP(IDIP)*BXDIP(IDIP)
     &    +BYDIP(IDIP)*BYDIP(IDIP)
     &    +BZDIP(IDIP)*BZDIP(IDIP))

        IF (B0DIP(IDIP).NE.0.D0) THEN
          RHODIP(IDIP)=EMOM/ABS(B0DIP(IDIP))/CLIGHT1
        ELSE
          if (ispecdip.eq.2) then
            b0dip(idip)=1.0d0
            RHODIP(IDIP)=EMOM/ABS(B0DIP(IDIP))/CLIGHT1
          else
            WRITE(LUNGFO,*)'*** ERROR IN SPECDIP ***'
            WRITE(LUNGFO,*)'zero magnetic field for dipole ',IDIP
            WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
            WRITE(6,*)'*** ERROR IN SPECDIP ***'
            WRITE(6,*)'zero magnetic field for dipole ',IDIP
            WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
            STOP
          endif
        ENDIF

        OMEGAC=1.5D0*DMYGAMMA**3*CLIGHT1/RHODIP(IDIP)
        ECDIP(IDIP)=OMEGAC*HBAR1/ECHARGE1

        if (ispecdip.eq.2) then
          ecdip(idip)=1.0d0
        else
          WRITE(LUNGFO,'(I5,4(1PE15.3))')
     &      IDIP,SNGL(B0DIP(IDIP)),SNGL(RHODIP(IDIP)),SNGL(ECDIP(IDIP))
     &      ,SNGL(ECDIP(IDIP)*SPECCUT)
          WRITE(LUNGFO,'(XXXXX,3(1PE15.3))')
     &      SNGL(X0DIP(IDIP)),SNGL(Y0DIP(IDIP)),SNGL(Z0DIP(IDIP))
          WRITE(LUNGFO,'(XXXXX,3(1PE15.3))')
     &      SNGL(VXDIP(IDIP)),SNGL(VYDIP(IDIP)),SNGL(VZDIP(IDIP))
          WRITE(LUNGFO,'(XXXXX,3(1PE15.3))')
     &      SNGL(BXDIP(IDIP)),SNGL(BYDIP(IDIP)),SNGL(BZDIP(IDIP))
        endif
      ENDDO !NDIP

      WRITE(LUNGFO,*)


      NSOURCE=NDIP

      DO IDIP=1,NDIP

        ISOUR=IDIP

        IF (APERTHICK.GT.0.0D0) THEN

          RX=OBSV(1,ICBRILL)-X0DIP(IDIP)
          RY=OBSV(2,ICBRILL)-Y0DIP(IDIP)
          RZ=OBSV(3,ICBRILL)-Z0DIP(IDIP)

          IF (RX.LE.0.D0) THEN
            WRITE(LUNGFO,*)'*** ERROR IN SPECDIP:'
            WRITE(LUNGFO,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
            WRITE(LUNGFO,*)'CHECK INPUT FILE'
            WRITE(6,*)'*** ERROR IN SPECDIP:'
            WRITE(6,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
            WRITE(6,*)'CHECK INPUT FILE'
            STOP '--- PROGRAM WAVE ABORTED ---'
          ENDIF

          RR=RX*RX+RY*RY+RZ*RZ
          RN=SQRT(RR)

          APERANG=ACOS((RX*APERV(1)+RY*APERV(2)+RZ*APERV(3))/RN)

          IF (IPINCIRC.EQ.0) THEN
            CALL THICKAPP(IPINCIRC,APERTHICK,PINW/2.0D0,APERANG,APERCORR)
          ELSE !IPINCIRC
            CALL THICKAPP(IPINCIRC,APERTHICK,PINR,APERANG,APERCORR)
          ENDIF

        ELSE !APERTHICK

          APERCORR=1.0D0

        ENDIF !APERTHICK

        DO IFREQ=1,NFREQ

          Y=FREQ(IFREQ)/ECDIP(IDIP)

          DO IOBSV=1,NOBSV

            RX=OBSV(1,IOBSV)-X0DIP(IDIP)
            RY=OBSV(2,IOBSV)-Y0DIP(IDIP)
            RZ=OBSV(3,IOBSV)-Z0DIP(IDIP)

            IF (RX.LE.0.D0) THEN
              WRITE(LUNGFO,*)'*** ERROR IN SPECDIP:'
              WRITE(LUNGFO,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
              WRITE(LUNGFO,*)'CHECK INPUT FILE'
              WRITE(6,*)'*** ERROR IN SPECDIP:'
              WRITE(6,*)'BAD X-DISTANCE FROM SOURCE TO OBSERVER'
              WRITE(6,*)'CHECK INPUT FILE'
              STOP '--- PROGRAM WAVE ABORTED ---'
            ENDIF

            RR=RX*RX+RY*RY+RZ*RZ
            RN=SQRT(RR)

            BS=B0DIP(IDIP)
            BX=BXDIP(IDIP)
            BY=BYDIP(IDIP)
            BZ=BZDIP(IDIP)

            BRX=BY*RZ-BZ*RY
            BRY=BZ*RX-BX*RZ
            BRZ=BX*RY-BY*RX

            BRN=SQRT(BRX*BRX+BRY*BRY+BRZ*BRZ)

            IF (BS.NE.0.D0) THEN
              PSI=ACOS(MIN(BRN/BS/RN,1.D0))
     &          *SIGN(1.D0,BX*RX+BY*RY+BZ*RZ)
            ELSE
              PSI=0.D0
            ENDIF


            ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
            IOBFR=IOBSV+NOBSV*(IFREQ-1)

            if (ispecdip.ne.2) then
              SPEC(ILIOBFR)=
     &          DFDTDP(Y,PSI,DMYGAMMA,DMYCUR,BANWID,PAR,PER,POWR)/RR
     &          *APERCORR
            else
              SPEC(ILIOBFR)=
     &          DFDTDP(Y,PSI,DMYGAMMA,DMYCUR,BANWID,PAR,PER,POWR)/RR
     &          /const
            endif

            IF (SPECCUT.NE.0.D0.AND.Y.GT.SPECCUT) THEN
              SPEC(ILIOBFR)=0.D0
              PAR=0.D0
              PER=0.D0
            ENDIF   !SPECCUT

            SPECPOW(ISOUR+NSOURCE*(IOBSV-1))=POWR/RHODIP(IDIP)/RR
     &      *APERCORR

            EXPOM=(1.D0,0.D0)

            IF (IPIN.NE.0) THEN

              DX2=RX*RX
              DZY2=RZ*RZ+RY*RY

              OM=FREQ(IFREQ)/(HBAREV1*CLIGHT1)

C       TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

              IF (DZY2.GT.0.01D0*DX2) THEN
                WRITE(LUNGFO,*)
     &            '*** ERROR IN SPECDIPA: OBSERVATION ANGLE TO LARGE ***'
                WRITE(LUNGFO,*)'DECREASE SIZE OF PINHOLE OR WGWINFC ...'
                WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
                WRITE(6,*)
     &            '*** ERROR IN SPECDIPA: OBSERVATION ANGLE TO LARGE ***'
                WRITE(6,*)'DECREASE SIZE OF PINHOLE OR WGWINFC ...'
                WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
                STOP
              ENDIF    !(DZY2.GT.0.01D0*DX2)

              EPS(1)=DZY2/DX2

              DO IEPS=2,6
                EPS(IEPS)=EPS(IEPS-1)*EPS(1)
              ENDDO !IEPS

c        TAYLOR-EXPANSION DONE WITH REDUCE
c       IN "WTAY1.RED";
c       on rounded;
c       on numval;
c       precision 13;
c       F:=SQRT(1+EPS);
c       DR:=TAY1(F,EPS,6);
c       ON FORT;
c       OUT "RED.FOR";
c       DR;
c       SHUT "RED.FOR";
C ans is actually reduce by 1.0 to avoid large overall phase

              ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &          -0.0390625D0*eps(4)+
     &          0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

              DRRED=-DABS(RX*ANS)
              EXPOM=CDEXP(DCMPLX(0.D0,DRRED*OM))

            ENDIF !IPIN

            AFREQ(1)=(0.D0,0.D0)
            AFREQ(2)=DCMPLX(0.D0,-SQRT(PER/RR/SPECNOR))*EXPOM
     &      *APERCORR
            AFREQ(3)=DCMPLX(SQRT(PAR/RR/SPECNOR),0.D0)*EXPOM*SIGN(1.D0,PSI)
     &      *APERCORR

            REAIMA(1,1,IOBFR)=DREAL(AFREQ(1))
            REAIMA(1,2,IOBFR)=DIMAG(AFREQ(1))
            REAIMA(2,1,IOBFR)=DREAL(AFREQ(2))
            REAIMA(2,2,IOBFR)=DIMAG(AFREQ(2))
            REAIMA(3,1,IOBFR)=DREAL(AFREQ(3))
            REAIMA(3,2,IOBFR)=DIMAG(AFREQ(3))

            IF (IPOLA.NE.0) THEN
              APOL=
     &          AFREQ(1)*CONJG(VPOLA(1))
     &          +AFREQ(2)*CONJG(VPOLA(2))
     &          +AFREQ(3)*CONJG(VPOLA(3))
              SPEC(ILIOBFR)=
     &          DREAL(APOL*CONJG(APOL))*SPECNOR
            ENDIF !IPOLA

            IF (ISTOKES.NE.0) THEN

              APOLH=
     &          AFREQ(1)*CONJG(VSTOKES(1,1))
     &          +AFREQ(2)*CONJG(VSTOKES(1,2))
     &          +AFREQ(3)*CONJG(VSTOKES(1,3))

              APOLR=
     &          AFREQ(1)*CONJG(VSTOKES(2,1))
     &          +AFREQ(2)*CONJG(VSTOKES(2,2))
     &          +AFREQ(3)*CONJG(VSTOKES(2,3))

              APOLL=
     &          AFREQ(1)*CONJG(VSTOKES(3,1))
     &          +AFREQ(2)*CONJG(VSTOKES(3,2))
     &          +AFREQ(3)*CONJG(VSTOKES(3,3))

              APOL45=
     &          AFREQ(1)*CONJG(VSTOKES(4,1))
     &          +AFREQ(2)*CONJG(VSTOKES(4,2))
     &          +AFREQ(3)*CONJG(VSTOKES(4,3))

              STOK1=
     &          APOLR*CONJG(APOLR)+
     &          APOLL*CONJG(APOLL)

              STOK2=-STOK1+
     &          2.0d0*APOLH*CONJG(APOLH)

              STOK3=
     &          2.0d0*APOL45*CONJG(APOL45)-
     &          STOK1

              STOK4=
     &          APOLR*CONJG(APOLR)-
     &          APOLL*CONJG(APOLL)

              if (abs(stok1)*specnor.gt.1.0D-30)
     &          STOKES(1,IOBFR)=STOKES(1,IOBFR)+
     &          STOK1*SPECNOR

              if (abs(stok2)*specnor.gt.1.0D-30)
     &          STOKES(2,IOBFR)=STOKES(2,IOBFR)+
     &          STOK2*SPECNOR

              if (abs(stok3)*specnor.gt.1.0D-30)
     &          STOKES(3,IOBFR)=STOKES(3,IOBFR)+
     &          STOK3*SPECNOR

              if (abs(stok4)*specnor.gt.1.0D-30)
     &          STOKES(4,IOBFR)=STOKES(4,IOBFR)+
     &          STOK4*SPECNOR

            ENDIF !ISTOKES

            SPECTOT(IOBFR)=SPECTOT(IOBFR)
     &        +SPEC(ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1)))

          ENDDO   !OBSV
        ENDDO  !IFREQ

      ENDDO !NDIP

      RETURN
      END
