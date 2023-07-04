*CMZ :          04/07/2023  11.46.29  by  Michael Scheer
*CMZ :  4.01/03 02/06/2023  13.01.26  by  Michael Scheer
*CMZ :  4.00/14 11/02/2022  10.26.56  by  Michael Scheer
*CMZ :  3.05/01 08/05/2018  16.08.32  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.69/02 07/11/2012  13.59.21  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  12.15.44  by  Michael Scheer
*CMZ :  2.67/00 13/02/2012  10.58.17  by  Michael Scheer
*CMZ :  2.65/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.64/04 14/09/2009  15.19.42  by  Michael Scheer
*CMZ :  2.64/03 21/08/2009  17.32.56  by  Michael Scheer
*CMZ :  2.64/02 21/08/2009  17.24.35  by  Michael Scheer
*CMZ :  2.64/01 20/08/2009  15.20.55  by  Michael Scheer
*CMZ :  2.63/05 03/08/2009  16.11.05  by  Michael Scheer
*CMZ :  2.52/12 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/11 08/12/2004  13.37.55  by  Michael Scheer
*CMZ :  2.51/02 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.50/00 29/04/2004  15.29.30  by  Michael Scheer
*CMZ :  2.20/01 24/11/2000  21.15.06  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.40.08  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  17.26.23  by  Michael Scheer
*CMZ :  2.13/07 17/02/2000  15.11.12  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  11.57.55  by  Michael Scheer
*CMZ :  2.13/03 17/01/2000  17.27.08  by  Michael Scheer
*CMZ :  2.12/02 15/06/99  15.16.33  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.08.55  by  Michael Scheer
*CMZ :  2.11/01 19/05/99  14.09.50  by  Michael Scheer
*CMZ :  2.10/01 19/03/99  14.13.05  by  Michael Scheer
*CMZ :  2.00/00 06/01/99  11.12.16  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.25.26  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.43  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.46  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE ARGSUM(ISOUR,IOBSV,IBUFF)

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
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,amplif90u.
      include 'amplif90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      INTEGER ISOUR,IOBSV,IFREQ,IPOI,ICOMP,IST,ICAL,JOBSV,IBUFF

      INTEGER NTUPP
      PARAMETER (NTUPP=22)
      REAL*8 FILLT(NTUPP)

      COMPLEX*16 APOL,EXPOM,DEXPOM,ZIOM,AX0,AY0,AZ0,AX,AY,AZ,ZI,ZONE,
     &  DMODU,DMODU0,DDMODU,daff(3),baff(3),bx0,by0,bz0,bx,by,bz
      COMPLEX*8 APOLH,APOLR,APOLL,APOL45

      DOUBLE PRECISION OM,DOM,DT,DT2,T,rnx,rny,rnz,rn
      DOUBLE PRECISION :: RSPECNOR=1.0d0 !11.2.2022

      DOUBLE PRECISION RARG(8),
     &  DXEXI,CENXEXI,DTPHASE,FREQR,GAMMA21,GAMGAM,R0,R02,H2,R00,H2R2,
     &  PHI,CORRR0,GAMGAM0,AMPDT,R2,POW

      REAL*4 STOK1,STOK2,STOK3,STOK4

      INTEGER I

      DATA ICAL/0/

      DATA ZI/(0.0D0,1.0D0)/
      DATA ZONE/(1.0D0,0.0D0)/

      IF (ICAL.EQ.0) THEN

c11.2.2022        RSPECNOR=DSQRT(SPECNOR) to be consistent with hfreq and souintana

        DO IFREQ=1,NFREQ
          DO JOBSV=1,NOBSV

            IOBFR=JOBSV+NOBSV*(IFREQ-1)

            REAIMA=0.0d0

            IF (ISTOKES.NE.0) THEN
              DO IST=1,4  !LOOP OVER COMPONENTS OF STOKES VECTOR
                STOKES(IST,IOBFR)=0.0
              ENDDO   !IST
            ENDIF

          ENDDO   !JOBSV
        ENDDO   !IFREQ

        ICAL=1

      ENDIF   !ICAL

C--- PERFORMS INTEGRATION FOR ALL FREQUENCES

C    ASSUMES FREQ(I+1)=FREQ(I)*2   FOR IFREQ2P=2
C    OR FREQ(I+1)=FREQ(I)+DELTA    FOR IFREQ2P>2

      IF (ISOUR.NE.ISOURO) THEN

        DO IFREQ=1,NFREQ
          IFROB=IFREQ+NFREQ*(IOBSV-1)
          AFREQ(1:6,IFROB)=(0.0d0,0.0d0)
        ENDDO

        SPECPOW(ISOUR+NSOURCE*(IOBSV-1))=0.D0

      ENDIF !ISOUR

      DOM=(FREQ(2)-FREQ(1))/HBAREV1

C-- LOOP OVER TIME STEPS (ACTUAL INTEGRATION)

      IF (IBUFF.EQ.1) THEN
        T=-REARGUM(6,1)
      ENDIF

      DO IPOI=1,IARGUM

        DT=REARGUM(6,IPOI)
        DT2=DT*0.5D0

        T=T+DT

        RARG(1)=REARGUM(1,IPOI)*DT
        RARG(2)=REARGUM(2,IPOI)*DT
        RARG(3)=REARGUM(3,IPOI)*DT

        RARG(4)=REARGUM(4,IPOI)
        RARG(5)=REARGUM(5,IPOI)*DT

        RARG(6:8)=REARGUM(9:11,IPOI)

        rn=norm2(rarg(6:8))
        rnx=rarg(6)/rn
        rny=rarg(7)/rn
        rnz=rarg(8)/rn

        OM=FREQ(1)/HBAREV1
        IF (IVELOFIELD.EQ.2) THEN
          ZIOM=DCMPLX(0.0D0,OM)
        ELSE
          ZIOM=DCMPLX(1.0D0,0.0D0)
        ENDIF

        EXPOM=CDEXP(DCMPLX(0.D0,RARG(4)*OM))

        IF(IFREQ2P.GT.2) DEXPOM=CDEXP(DCMPLX(0.D0,RARG(4)*DOM))

        IFROB=1+NFREQ*(IOBSV-1)
        DO ICOMP=1,3
          daff(icomp)=ZIOM*DCMPLX(RARG(ICOMP))*EXPOM*REFLEC(ICOMP)
          AFREQ(ICOMP,IFROB)=AFREQ(ICOMP,IFROB)+daff(icomp)
        ENDDO

        baff(1)=conjg(rny*daff(3)-rnz*daff(2))
        baff(2)=conjg(rnz*daff(1)-rnx*daff(3))
        baff(3)=conjg(rnx*daff(2)-rny*daff(1))

        afreq(4:6,IFROB)=afreq(4:6,IFROB)+baff(1:3)/clight1

C--- LOOP OVER ALL FREQUENCES

        DO IFREQ=2,NFREQ

            OM=FREQ(IFREQ)/HBAREV1
            IF (IVELOFIELD.EQ.2) THEN
              ZIOM=DCMPLX(0.0D0,OM)
            ELSE
              ZIOM=DCMPLX(1.0D0,0.0D0)
            ENDIF

            IF    (IFREQ2P.LT.2) THEN
              EXPOM=CDEXP(DCMPLX(0.D0,RARG(4)*OM))
            ELSE IF(IFREQ2P.EQ.2) THEN
              EXPOM=EXPOM*EXPOM
            ELSE
              EXPOM=EXPOM*DEXPOM
            ENDIF

            IFROB=IFREQ+NFREQ*(IOBSV-1)
            DO ICOMP=1,3
              daff(icomp)=ZIOM*DCMPLX(RARG(ICOMP))*EXPOM
     &          *REFLEC(ICOMP)
              AFREQ(ICOMP,IFROB)=AFREQ(ICOMP,IFROB)+daff(icomp)
            ENDDO

            baff(1)=conjg(rny*daff(3)-rnz*daff(2))
            baff(2)=conjg(rnz*daff(1)-rnx*daff(3))
            baff(3)=conjg(rnx*daff(2)-rny*daff(1))

            afreq(4:6,IFROB)=afreq(4:6,IFROB)+baff(1:3)/clight1

          ENDDO   !LOOP OVER ALL FREQUENCES

        ILIOB=ISOUR+NSOURCE*(IOBSV-1)
        SPECPOW(ILIOB)=SPECPOW(ILIOB)+RARG(5)

        IF(IWFILINT.EQ.-ISOUR.AND.IOBSV.EQ.1) THEN

          IFREQ=1
          IOBSV=1

          FILLT(1)=T
          FILLT(2)=WSOU(1,1,IPOI)
          FILLT(3)=WSOU(2,1,IPOI)
          FILLT(4)=WSOU(3,1,IPOI)
          FILLT(5)=RARG(1)
          FILLT(6)=RARG(2)
          FILLT(7)=RARG(3)
          FILLT(8)=RARG(4)
          FILLT(9)=MIN(RARG(5),1.D30)
          FILLT(10)=DREAL(EXPOM)
          FILLT(11)=DIMAG(EXPOM)
          FILLT(12)=0.0d0
          FILLT(13)=IOBSV
          FILLT(14)=IFREQ
          IF (NOBSV.EQ.1) THEN
            FILLT(15)=OBS1Y
            FILLT(16)=OBS1Z
          ELSE
            FILLT(15)=OBSV(2,IOBSV)
            FILLT(16)=OBSV(3,IOBSV)
          ENDIF
          FILLT(17)=0.0d0
          FILLT(18)=0.0d0
          FILLT(19)=0.0d0
          FILLT(20)=0.0d0
          IFROB=IFREQ+NFREQ*(IOBSV-1)
          FILLT(21)=
     &      DREAL(
     &      AFREQ(1,IFROB)*CONJG(AFREQ(1,IFROB))
     &      +AFREQ(2,IFROB)*CONJG(AFREQ(2,IFROB))
     &      +AFREQ(3,IFROB)*CONJG(AFREQ(3,IFROB))
     &      )*SPECNOR
          FILLT(22)=0.0d0

          CALL hfm(NIDSOURCE,FILLT)

        ENDIF !IWFILINT.EQ.-ISOUR

      ENDDO   !LOOP OVER TIME STEPS

      IF (NSADD.NE.0) THEN

        IF (IAMPLI.LT.0) THEN

          GAMGAM0=SOURCEG(1,1,ISOUR)**2
          GAMGAM=(SOURCEG(1,1,ISOUR)+SOURCEG(2,2,ISOUR))**2
          GAMMA21=1.0D0/GAMGAM0
          DXEXI=MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &      -MAX(SOURCEAO(1,1,ISOUR),XIANF)
          CENXEXI=(MIN(SOURCEEO(1,1,ISOUR),XIEND)
     &      +MAX(SOURCEAO(1,1,ISOUR),XIANF))/2.D0
          DTPHASE=(WTRA2IS(ISOUR)+GAMMA21*DXEXI/2.D0)/CLIGHT1
     &    *GAMGAM0/GAMGAM
          AMPDT=AMPSHIFT(1)/CLIGHT1/2.0D0/GAMGAM0
          DTPHASE=DTPHASE+AMPDT
          FREQR=2.D0*PI1/DTPHASE*HBAREV1
          POW=SPECPOW(ILIOB)
          SPECPOW(ILIOB)=0.0D0

          DO I=1,-IAMPLI
            R02=(OBSV(1,IOBSV)-CENXEXI)**2+OBSV(2,IOBSV)**2+OBSV(3,IOBSV)**2
            R2=(OBSV(1,IOBSV)-CENXEXI-DXEXI*(I-ABS(IAMPLI)/2+1))**2
     &        +OBSV(2,IOBSV)**2+OBSV(3,IOBSV)**2
            SPECPOW(ILIOB)=SPECPOW(ILIOB)+POW*R02/R2
          ENDDO

          DO IFREQ=1,NFREQ

            ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
            IFROB=IFREQ+NFREQ*(IOBSV-1)
            IOBFR=IOBSV+NOBSV*(IFREQ-1)

            OM=FREQ(IFREQ)/HBAREV1

            AX0=AFREQ(1,IFROB)
            AY0=AFREQ(2,IFROB)
            AZ0=AFREQ(3,IFROB)

            AX=AX0
            AY=AY0
            AZ=AZ0

            BX0=AFREQ(4,IFROB)
            BY0=AFREQ(5,IFROB)
            BZ0=AFREQ(6,IFROB)

            BX=AX0
            BY=AY0
            BZ=AZ0

            AFREQ(1:6,IFROB)=(0.0D0,0.0D0)

            R0=OBSV(1,NOBSV/2+1)-CENXEXI
            R02=R0*R0
            R00=R0
            H2=(OBSV(2,IOBSV))**2+(OBSV(3,IOBSV))**2
            H2R2=H2/R02

            DTPHASE=(WTRA2IS(ISOUR)+(H2R2+GAMMA21)*DXEXI/2.D0)/CLIGHT1
     &        *GAMGAM0/GAMGAM
     &        +AMPDT
            PHI=2.D0*PI1*FREQ(IFREQ)*ECHARGE1/HPLANCK1*DTPHASE

            DMODU=EXP(ZI*PHI)
            DMODU0=DMODU
            DDMODU=ZONE

            DO I=1,-IAMPLI

              R0=OBSV(1,NOBSV/2+1)+DXEXI/2.D0*(-IAMPLI-2*(I-1)-1)
              CORRR0=R00/R0
              R02=R0*R0
              H2=(OBSV(2,IOBSV))**2+(OBSV(3,IOBSV))**2
              H2R2=H2/R02

              GAMGAM=(SOURCEG(1,1,ISOUR)+(I-1)*SOURCEG(2,2,ISOUR))**2
              DTPHASE=(WTRA2IS(ISOUR)+(H2R2+GAMMA21)*DXEXI/2.D0)/CLIGHT1
     &          *GAMGAM0/GAMGAM+AMPDT
              PHI=2.D0*PI1*FREQ(IFREQ)*ECHARGE1/HPLANCK1*DTPHASE

              DMODU=EXP(ZI*PHI)
              DMODU0=DMODU
              DDMODU=ZONE

              AFREQ(1,IFROB)=AFREQ(1,IFROB)+AX
              AFREQ(2,IFROB)=AFREQ(2,IFROB)+AY
              AFREQ(3,IFROB)=AFREQ(3,IFROB)+AZ
              AFREQ(4,IFROB)=AFREQ(4,IFROB)+BX
              AFREQ(5,IFROB)=AFREQ(5,IFROB)+BY
              AFREQ(6,IFROB)=AFREQ(6,IFROB)+BZ

              IF (AMPRAN.NE.0.D0) THEN
                PHI=2.D0*PI1*XRANA(I)/FREQR*FREQ(IFREQ)
                DDMODU=EXP(ZI*PHI)
              ENDIF   !(AMPRAN.NE.0.D0)

              AX0=AX0*DMODU0
              AY0=AY0*DMODU0
              AZ0=AZ0*DMODU0

              AX=AX0*CORRR0
              AY=AY0*CORRR0
              AZ=AZ0*CORRR0

              DMODU=DMODU0*DDMODU
              AX=AX*DMODU
              AY=AY*DMODU
              AZ=AZ*DMODU

              BX0=BX0*DMODU0
              BY0=BY0*DMODU0
              BZ0=BZ0*DMODU0

              BX=BX0*CORRR0
              BY=BY0*CORRR0
              BZ=BZ0*CORRR0

              BX=BX*DMODU
              BY=BY*DMODU
              BZ=BZ*DMODU

            ENDDO !IAMPLI

          ENDDO !IFREQ

        ENDIF  !(IAMPLI.LT.0)

        DO IFREQ=1,NFREQ

          ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          IFROB=IFREQ+NFREQ*(IOBSV-1)

          IF (IPOLA.EQ.0) THEN

            SPEC(ILIOBFR)=
     &        DREAL(
     &        AFREQ(1,IFROB)*CONJG(AFREQ(1,IFROB))
     &        +AFREQ(2,IFROB)*CONJG(AFREQ(2,IFROB))
     &        +AFREQ(3,IFROB)*CONJG(AFREQ(3,IFROB))
     &        )*SPECNOR

            REAIMA(1:3,1,IOBFR)=REAIMA(1:3,1,IOBFR)+
     &        DREAL(AFREQ(1:3,IFROB))*rspecnor
            REAIMA(1:3,2,IOBFR)=REAIMA(1:3,2,IOBFR)+
     &        DIMAG(AFREQ(1:3,IFROB))*rspecnor
            REAIMA(6:8,1,IOBFR)=REAIMA(6:8,1,IOBFR)+
     &        DREAL(AFREQ(1:3,IFROB))*rspecnor
            REAIMA(6:8,2,IOBFR)=REAIMA(6:8,2,IOBFR)+
     &        DIMAG(AFREQ(1:3,IFROB))*rspecnor

          ELSE    !IPOLA

            APOL=
     &        AFREQ(1,IFROB)*CONJG(VPOLA(1))
     &        +AFREQ(2,IFROB)*CONJG(VPOLA(2))
     &        +AFREQ(3,IFROB)*CONJG(VPOLA(3))

            SPEC(ILIOBFR)=
     &        DREAL(APOL*CONJG(APOL))*SPECNOR

            REAIMA(1:3,1,IOBFR)=REAIMA(1:3,1,IOBFR)+
     &        DREAL(AFREQ(1:3,IFROB))*rspecnor
            REAIMA(1:3,2,IOBFR)=REAIMA(1:3,2,IOBFR)+
     &        DIMAG(AFREQ(1:3,IFROB))*rspecnor
            REAIMA(6:8,1,IOBFR)=REAIMA(6:8,1,IOBFR)+
     &        DREAL(AFREQ(1:3,IFROB))*rspecnor
            REAIMA(6:8,2,IOBFR)=REAIMA(6:8,2,IOBFR)+
     &        DIMAG(AFREQ(1:3,IFROB))*rspecnor

          ENDIF   !IPOLA

          IOBFR=IOBSV+NOBSV*(IFREQ-1)
          IFROB=IFREQ+NFREQ*(IOBSV-1)

          IF (ISTOKES.NE.0) THEN

            APOLH=
     &        AFREQ(1,IFROB)*CONJG(VSTOKES(1,1))
     &        +AFREQ(2,IFROB)*CONJG(VSTOKES(1,2))
     &        +AFREQ(3,IFROB)*CONJG(VSTOKES(1,3))

            APOLR=
     &        AFREQ(1,IFROB)*CONJG(VSTOKES(2,1))
     &        +AFREQ(2,IFROB)*CONJG(VSTOKES(2,2))
     &        +AFREQ(3,IFROB)*CONJG(VSTOKES(2,3))

            APOLL=
     &        AFREQ(1,IFROB)*CONJG(VSTOKES(3,1))
     &        +AFREQ(2,IFROB)*CONJG(VSTOKES(3,2))
     &        +AFREQ(3,IFROB)*CONJG(VSTOKES(3,3))

            APOL45=
     &        AFREQ(1,IFROB)*CONJG(VSTOKES(4,1))
     &        +AFREQ(2,IFROB)*CONJG(VSTOKES(4,2))
     &        +AFREQ(3,IFROB)*CONJG(VSTOKES(4,3))

            STOK1=
     &        REAL(APOLR*CONJG(APOLR))+
     &        REAL(APOLL*CONJG(APOLL))

            STOK2=-STOK1+
     &        2.*REAL(APOLH*CONJG(APOLH))

            STOK3=
     &        2.*REAL(APOL45*CONJG(APOL45))-
     &        STOK1

            STOK4=
     &        REAL(APOLR*CONJG(APOLR))-
     &        REAL(APOLL*CONJG(APOLL))


            STOKES(1,IOBFR)=STOKES(1,IOBFR)+
     &        STOK1*SPECNOR

            STOKES(2,IOBFR)=STOKES(2,IOBFR)+
     &        STOK2*SPECNOR

            STOKES(3,IOBFR)=STOKES(3,IOBFR)+
     &        STOK3*SPECNOR

            STOKES(4,IOBFR)=STOKES(4,IOBFR)+
     &        STOK4*SPECNOR
          ENDIF

        ENDDO !IFREQ

        ILIOB=ISOUR+NSOURCE*(IOBSV-1)
        SPECPOW(ILIOB)=SPECPOW(ILIOB)*pownor

      ENDIF !NSADD

      RETURN
      END
