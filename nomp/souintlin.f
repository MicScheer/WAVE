*CMZ :  4.00/14 22/12/2021  18.07.25  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  15.33.15  by  Michael Scheer
*CMZ :  3.02/03 06/11/2014  14.24.54  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  14.04.56  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.66/04 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.66/03 12/11/2009  16.27.11  by  Michael Scheer
*CMZ :  2.50/00 28/10/2009  15.52.52  by  Michael Scheer
*CMZ :  2.41/10 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.41/09 14/08/2002  17.11.32  by  Michael Scheer
*CMZ :  2.20/04 09/03/2001  14.23.38  by  Michael Scheer
*CMZ :  2.20/03 22/02/2001  18.37.06  by  Michael Scheer
*CMZ :  2.20/01 11/02/2001  19.28.38  by  Michael Scheer
*CMZ :  2.16/08 31/10/2000  14.40.08  by  Michael Scheer
*CMZ :  2.16/07 21/09/2000  11.21.00  by  Michael Scheer
*CMZ :  2.16/06 28/08/2000  14.39.42  by  Michael Scheer
*CMZ :  2.16/05 02/08/2000  13.53.32  by  Michael Scheer
*CMZ :  2.16/04 19/06/2000  14.27.19  by  Michael Scheer
*CMZ :  2.16/03 16/06/2000  14.35.02  by  Michael Scheer
*CMZ :  2.15/00 05/05/2000  19.25.24  by  Michael Scheer
*CMZ :  2.14/02 26/04/2000  16.40.34  by  Michael Scheer
*CMZ :  2.13/07 17/02/2000  15.11.13  by  Michael Scheer
*CMZ :  2.13/03 18/01/2000  17.44.41  by  Michael Scheer
*CMZ :  2.12/03 21/07/99  10.47.09  by  Michael Scheer
*CMZ :  2.12/01 10/06/99  17.49.00  by  Michael Scheer
*CMZ :  2.12/00 04/06/99  10.43.55  by  Michael Scheer
*CMZ :  2.11/01 20/05/99  17.44.32  by  Michael Scheer
*CMZ :  2.10/01 24/02/99  10.20.40  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.18.17  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.05  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.44  by  Michael Scheer
*-- Author : Michael Scheer

      SUBROUTINE SOUINTLIN(ISOUR)
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
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEND.

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE

         IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

      INTEGER NUMSTEP0,NUMSTEP

      INTEGER ISOUR,IZAEHL,ITIM0,IWARN,IOBSV,ICOMP,IFREQ,IROI
      INTEGER IEXPOMT,IEXPDOMT,IX10

      INTEGER ICAL,NTUPP,ICYCLE,IC,I
      PARAMETER (NTUPP=22)
        REAL*8 FILLT(NTUPP)
      CHARACTER(4) CTUP(NTUPP)

      DOUBLE PRECISION X1,TS,X0,X10
      DOUBLE PRECISION X2,Y2,Z2
      DOUBLE PRECISION XENDSOU,DT,DT0,DTIM01,DT2,DROIX,DTE,TE,DDTS,DT1
      DOUBLE PRECISION RARGOM,RARGDOM,RARGOMO,RARGDOMO
      DOUBLE PRECISION ARG,ARGO,ARG2,ARG4,DARG,DARG2,DARG4

      DOUBLE PRECISION BX,BY,BZ,DUM1,DUM11,DOM1,DOM2,R0,C1,PX,PY,PZ,BPX,BPY,BPZ
      DOUBLE PRECISION RNBX,RNBY,RNBZ,T,R,RNX,RNY,RNZ,R1,RX,RY,RZ,RARG(5)

      DOUBLE PRECISION CGG


      COMPLEX*16 APOL
      COMPLEX*16 EXPOM,DEXPOM,DEXPOMT,DEXPDOMT
      COMPLEX*8 APOLH,APOLR,APOLL,APOL45


      DOUBLE PRECISION OM,DOM
      REAL*4 STOK1,STOK2,STOK3,STOK4

      DATA ICAL/0/
      data ctup /'t','x','y','z','rx','ry','rz','rt','p','rea','ima','roi'
     &            ,'iob','ie','yob','zob','betn','dtom','emod','dmod'
     &            ,'spec','te'/

      DATA IWARN/0/

      IF (ICAL.EQ.0) THEN

          CGG=CLIGHT1/(DMYGAMMA*DMYGAMMA)
          DOM=(FREQ(2)-FREQ(1))/HBAREV1
          OM=FREQ(1)/HBAREV1
          C1=1.D0/CLIGHT1
          IF (IWFILINT.LT.0)
     &      CALL hbookm(NIDSOURCE,'RADIATION INTEGRAL',NTUPP
     &     ,'//WAVE',nlpoi/jwfilint+2*jwfilint,CTUP)

          DTIM01=1.D0/DTIM0

          ALLOCATE(DWTRA(3,6,NCO+1))

          DO ITIM0=1,NCO

         DWTRA(1,1,ITIM0)=WTRA(1,1,ITIM0)
         DWTRA(2,1,ITIM0)=WTRA(2,1,ITIM0)
         DWTRA(3,1,ITIM0)=WTRA(3,1,ITIM0)

         DWTRA(1,3,ITIM0)=WTRA(1,2,ITIM0)*C1
         DWTRA(2,3,ITIM0)=WTRA(2,2,ITIM0)*C1
         DWTRA(3,3,ITIM0)=WTRA(3,2,ITIM0)*C1

          ENDDO   !NCO

         DWTRA(1,5,1)=(WTRA(1,2,2)-WTRA(1,2,1))*C1*DTIM01
         DWTRA(2,5,1)=(WTRA(2,2,2)-WTRA(2,2,1))*C1*DTIM01
         DWTRA(3,5,1)=(WTRA(3,2,2)-WTRA(3,2,1))*C1*DTIM01

          DO ITIM0=2,NCO-1

         DWTRA(1,5,ITIM0)=((WTRA(1,2,ITIM0+1)-WTRA(1,2,ITIM0-1))*C1*DTIM01)/2.D0
         DWTRA(2,5,ITIM0)=((WTRA(2,2,ITIM0+1)-WTRA(2,2,ITIM0-1))*C1*DTIM01)/2.D0
         DWTRA(3,5,ITIM0)=((WTRA(3,2,ITIM0+1)-WTRA(3,2,ITIM0-1))*C1*DTIM01)/2.D0

          ENDDO   !NCO

         DWTRA(1,5,NCO)=DWTRA(1,5,NCO-1)+DWTRA(1,5,NCO-1)-DWTRA(1,5,NCO-2)
         DWTRA(2,5,NCO)=DWTRA(2,5,NCO-1)+DWTRA(2,5,NCO-1)-DWTRA(2,5,NCO-2)
         DWTRA(3,5,NCO)=DWTRA(3,5,NCO-1)+DWTRA(3,5,NCO-1)-DWTRA(3,5,NCO-2)

         DWTRA(1,5,NCO+1)=DWTRA(1,5,NCO)+DWTRA(1,5,NCO)-DWTRA(1,5,NCO-1)
         DWTRA(2,5,NCO+1)=DWTRA(2,5,NCO)+DWTRA(2,5,NCO)-DWTRA(2,5,NCO-1)
         DWTRA(3,5,NCO+1)=DWTRA(3,5,NCO)+DWTRA(3,5,NCO)-DWTRA(3,5,NCO-1)

         DWTRA(1,1,NCO+1)=DWTRA(1,1,NCO)+DWTRA(1,1,NCO)-DWTRA(1,1,NCO-1)
         DWTRA(2,1,NCO+1)=DWTRA(2,1,NCO)+DWTRA(2,1,NCO)-DWTRA(2,1,NCO-1)
         DWTRA(3,1,NCO+1)=DWTRA(3,1,NCO)+DWTRA(3,1,NCO)-DWTRA(3,1,NCO-1)

         DWTRA(1,3,NCO+1)=DWTRA(1,3,NCO)+DWTRA(1,3,NCO)-DWTRA(1,3,NCO-1)
         DWTRA(2,3,NCO+1)=DWTRA(2,3,NCO)+DWTRA(2,3,NCO)-DWTRA(2,3,NCO-1)
         DWTRA(3,3,NCO+1)=DWTRA(3,3,NCO)+DWTRA(3,3,NCO)-DWTRA(3,3,NCO-1)

          DO ITIM0=1,NCO

         DWTRA(1,2,ITIM0)=DWTRA(1,1,ITIM0+1)-DWTRA(1,1,ITIM0)
         DWTRA(2,2,ITIM0)=DWTRA(2,1,ITIM0+1)-DWTRA(2,1,ITIM0)
         DWTRA(3,2,ITIM0)=DWTRA(3,1,ITIM0+1)-DWTRA(3,1,ITIM0)

         DWTRA(1,4,ITIM0)=DWTRA(1,3,ITIM0+1)-DWTRA(1,3,ITIM0)
         DWTRA(2,4,ITIM0)=DWTRA(2,3,ITIM0+1)-DWTRA(2,3,ITIM0)
         DWTRA(3,4,ITIM0)=DWTRA(3,3,ITIM0+1)-DWTRA(3,3,ITIM0)

         DWTRA(1,6,ITIM0)=DWTRA(1,5,ITIM0+1)-DWTRA(1,5,ITIM0)
         DWTRA(2,6,ITIM0)=DWTRA(2,5,ITIM0+1)-DWTRA(2,5,ITIM0)
         DWTRA(3,6,ITIM0)=DWTRA(3,5,ITIM0+1)-DWTRA(3,5,ITIM0)

          ENDDO   !NCO

          ICAL=1

      ENDIF !ICAL

      DO IOBSV=1,NOBSV

          SPECPOW(ISOUR+NSOURCE*(IOBSV-1))=0.D0

          DARGEXPO(1,IOBSV)=(1.D0,0.D0)   !DEXPOMT
          DARGEXPO(2,IOBSV)=(1.D0,0.D0)   !DEXPDOMT
          DARGEXPO(3,IOBSV)=(0.D0,0.D0)   !(RARGOM,RARGOMO)
          DARGEXPO(4,IOBSV)=(1.D0,0.D0)   !EXPOM
          DARGEXPO(5,IOBSV)=(0.D0,0.D0)   !(ARG,ARG)
          DARGEXPO(6,IOBSV)=(1.D0,0.D0)   !DEXPOM

      DO IFREQ=1,NFREQ

            IFROB=IFREQ+NFREQ*(IOBSV-1)
            AFREQ(1,IFROB)=(0.,0.)
            AFREQ(2,IFROB)=(0.,0.)
            AFREQ(3,IFROB)=(0.,0.)

      ENDDO   !IFREQ
      ENDDO   !IOBSV

      TS=SOURCET(1,ISOUR)
      ITIM0=TS*DTIM01+1

      R0=OBSV(1,1)-SOURCEAO(1,1,ISOUR)
C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS  T=-R0*C1
      T=0.D0

      X1=SOURCEAO(1,1,ISOUR)

      IZTOT(ISOUR)=0

      XENDSOU=SOURCEEO(1,1,ISOUR)    !FINAL X

        X0=X1
        X2=X1
        X10=(XENDSOU-X0)/10.1D0

      DT0=(XENDSOU-X1)/NLPOIO/CLIGHT1
      DT=DT0
      DT2=DT/2.D0
      DT1=1.D0/DT

      IF (NROI.LT.0) THEN
          DROIX=(XENDSOU-X1)/(NROIA-1)
          DO IROI=1,NROIA
         ROIX(IROI)=X1+(IROI-1)*DROIX
         ROIP(IROI)=1.D0
          ENDDO
      ENDIF !

      ROIX(1)=ROIX(1)-1.D-6
      ROIX(NROIA)=ROIX(NROIA)+1.D-6

        IF (X1.LT.ROIX(1).OR.XENDSOU.GT.ROIX(NROIA)) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN SOUINTLIN: X OUTSIDE ROIS ***'
            WRITE(LUNGFO,*)'CHECK NAMELIST $ROIN'
            WRITE(LUNGFO,*)' *** PROGRAM WAVE ABORTED ***'
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN SOUINTLIN: X OUTSIDE ROIS ***'
            WRITE(6,*)'CHECK NAMELIST $ROIN'
            WRITE(6,*)' *** PROGRAM WAVE ABORTED ***'
            STOP
        ENDIF   !IROI

C- CHECK NUMBER OF STEPS

      IF (IWARN.EQ.0) THEN

          NUMSTEP0=NLPOIO/(SOURCEEO(1,1,NSOURCE)-SOURCEAO(1,1,1))

          DO IROI=1,NROIA-1

         NUMSTEP=(ROIX(IROI+1)-ROIX(IROI))*NUMSTEP0

         IF (NUMSTEP.LT.MYINUM) THEN

             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING IN SOUINTLIN, ROI:',IROI
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
             WRITE(LUNGFO,*)'SIZE FOR TRAJECTORY!'
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
     &              'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)'*** WARNING IN SOUINTLIN, ROI:',IROI
             WRITE(6,*)
             WRITE(6,*)'STEP SIZE FOR SOURCE POINT IS LARGER THAN STEP'
             WRITE(6,*)'SIZE FOR TRAJECTORY!'
             WRITE(6,*)
             WRITE(6,*)
     &              'CHANGE NLPOI OR ROI-PARAMETERS OR BE AWARE OF STRANGE RESULTS!'
             WRITE(6,*)

           ENDIF

          ENDDO   !IROI

         IWARN=1

      ENDIF

      IZAEHL=0 !LOOP COUNTER

C--- LOOP OVER STEPS

      DO IROI=1,NROIA
          IPOIROI(IROI)=0
      ENDDO

        IROI=1
        DO I=1,NROIA
            IF (X1.GE.ROIX(I)) THEN
                IROI=I
            ENDIF
        ENDDO   !IROI

        IX10=1

        IF (ISOUR.EQ.1) THEN
            WRITE(6,*)' '
            WRITE(6,*)
     &      '      counting from 1 to 10 for first source to show progress:'
            WRITE(6,*)' '
        ENDIF   !NSOURCE

        IZAEHL=0 !LOOP COUNTER

1000    IZAEHL=IZAEHL+1

        IF (ISOUR.EQ.1) THEN
            IF (X2.GE.X0+X10*IX10) THEN
         CALL date_and_time(dtday,dttime,dtzone,idatetime)
                WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
                IX10=IX10+1
            ENDIF   !X1
        ENDIF   !NSOURCE

      IF (IROI.LE.NROIA) THEN
          IF (X2.GE.ROIX(IROI)) THEN
         DT=DT0/ROIP(IROI)
         DT2=DT/2.D0
         DT1=1.D0/DT
         IROI=IROI+1
          ENDIF   !IROI
      ENDIF !IROI

      IPOIROI(IROI)=IPOIROI(IROI)+1

      T=T+DT

      DDTS=DTIM01*(TS+DT-(ITIM0-1)*DTIM0)

      X2=DWTRA(1,1,ITIM0)+DWTRA(1,2,ITIM0)*DDTS
      Y2=DWTRA(2,1,ITIM0)+DWTRA(2,2,ITIM0)*DDTS
      Z2=DWTRA(3,1,ITIM0)+DWTRA(3,2,ITIM0)*DDTS

      BX=DWTRA(1,3,ITIM0)+DWTRA(1,4,ITIM0)*DDTS
      BY=DWTRA(2,3,ITIM0)+DWTRA(2,4,ITIM0)*DDTS
      BZ=DWTRA(3,3,ITIM0)+DWTRA(3,4,ITIM0)*DDTS

      BPX=DWTRA(1,5,ITIM0)+DWTRA(1,6,ITIM0)*DDTS
      BPY=DWTRA(2,5,ITIM0)+DWTRA(2,6,ITIM0)*DDTS
      BPZ=DWTRA(3,5,ITIM0)+DWTRA(3,6,ITIM0)*DDTS

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION {

C REAL PART OF INTEGRAND {

      DO IOBSV=1,NOBSV

          RX=OBSV(1,IOBSV)-X2
          RY=OBSV(2,IOBSV)-Y2
          RZ=OBSV(3,IOBSV)-Z2

          R=DSQRT(RX*RX+RY*RY+RZ*RZ)
          R1=1.D0/R

          RNX=RX*R1
          RNY=RY*R1
          RNZ=RZ*R1

C--- THE DISTANCE R IS INTRODUCED HERE EXPLICITLY (S. PROGRAM OF CHAOEN WANG

          DUM1=(1.D0-BX*RNX)-BY*RNY-BZ*RNZ
          DUM11=1.D0/DUM1
          DOM1=1.D0/(R*DUM1*DUM1)

          DTE=DT*DUM1

          RNBX=RNX-BX
          RNBY=RNY-BY
          RNBZ=RNZ-BZ

          PX=(RNBY*BPZ-RNBZ*BPY)
          PY=(RNBZ*BPX-RNBX*BPZ)
          PZ=(RNBX*BPY-RNBY*BPX)

            IF (IVELOFIELD.EQ.0) THEN
              DOM2=CGG*DOM1*R1
              RARG(1)=(RNY*PZ-RNZ*PY)*DOM1+(RNX-BX)*DOM2
              RARG(2)=(RNZ*PX-RNX*PZ)*DOM1+(RNY-BY)*DOM2
              RARG(3)=(RNX*PY-RNY*PX)*DOM1+(RNZ-BZ)*DOM2
            ELSE IF (IVELOFIELD.EQ.1) THEN
              RARG(1)=(RNY*PZ-RNZ*PY)*DOM1
              RARG(2)=(RNZ*PX-RNX*PZ)*DOM1
              RARG(3)=(RNX*PY-RNY*PX)*DOM1
            ELSE IF (IVELOFIELD.LT.0) THEN
              DOM2=CGG*DOM1*R1
              RARG(1)=(RNX-BX)*DOM2
              RARG(2)=(RNY-BY)*DOM2
              RARG(3)=(RNZ-BZ)*DOM2
            ELSE  !IVELOFIELD
              WRITE(6,*)
     &          '*** ERROR IN SOUINTLIN: BAD VALUE OF IVELOFIELD  ***'
              WRITE(6,*) '*** PROGRAM WAVE ABORTED  ***'
              STOP
            ENDIF !IVELOFIELD

C DO NOT USE, RESULTS IN NUMERICAL PROBLEMS      RARG(4)=T+R*C1

          TE=T+(R-R0)*C1
          RARG(4)=TE

          RARG(5)=
     &      (RARG(1)*RARG(1)+RARG(2)*RARG(2)+RARG(3)*RARG(3))*DUM11

        IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN
            ILIOB=ISOUR+NSOURCE*(IOBSV-1)
          SPECPOW(ILIOB)=SPECPOW(ILIOB)+RARG(5)*DT
        ENDIF  !XIANF

C REAL PART OF INTEGRAND }

C COMPLEX PART OF INTEGRAND {

C    ASSUMES FREQ(I+1)=2*FREQ(I)   FOR IFREQ2P=2
C    OR FREQ(I+1)=FREQ(I)+DELTA    FOR IFREQ2P>2

        RARGOM=RARG(4)*OM


        RARGOMO=DREAL(DARGEXPO(3,IOBSV))
        ARG=RARGOM-RARGOMO
        ARGO=DREAL(DARGEXPO(5,IOBSV))
        DARG=ARG-ARGO
        DARGEXPO(5,IOBSV)=DCMPLX(ARG,DIMAG(DARGEXPO(5,IOBSV)))

        IF (DABS(ARG).LT.0.1D0) THEN
           ARG2=ARG*ARG
           ARG4=ARG2*ARG2
           DEXPOMT=DCMPLX(1.D0-ARG2/2.D0+ARG4/24.D0,
     &                    ARG*(1.D0-ARG2/6.D0+ARG4/120.D0))
                IEXPOMT=1
        ELSEIF (DABS(DARG).LT.0.1D0) THEN
           DARG2=DARG*DARG
           DARG4=DARG2*DARG2
                DEXPOMT=DARGEXPO(1,IOBSV)*
     &            DCMPLX(1.D0-DARG2/2.D0+DARG4/24.D0,
     &            DARG*(1.D0-DARG2/6.D0+DARG4/120.D0))
                IEXPOMT=2
        ELSE
                DEXPOMT=CDEXP(DCMPLX(0.D0,ARG))
                IEXPOMT=3
        ENDIF

          IF(IFREQ2P.GT.2) THEN

          EXPOM=DARGEXPO(4,IOBSV)*DEXPOMT
          DARGEXPO(4,IOBSV)=EXPOM    !STORE IT HERE, WILL CHANGE FURTHER DOWN

          RARGDOM=RARG(4)*DOM
          RARGDOMO=DIMAG(DARGEXPO(3,IOBSV))
          ARG=RARGDOM-RARGDOMO
          ARG2=ARG*ARG
          ARGO=DIMAG(DARGEXPO(5,IOBSV))
          DARG=ARG-ARGO
          DARGEXPO(5,IOBSV)=DCMPLX(DREAL(DARGEXPO(5,IOBSV)),ARG)

          IF (DABS(ARG).LT.0.1D0) THEN
           ARG2=ARG*ARG
           ARG4=ARG2*ARG2
           DEXPDOMT=DCMPLX(1.D0-ARG2/2.D0+ARG4/24.D0,
     &                     ARG*(1.D0-ARG2/6.D0+ARG4/120.D0))
                IEXPDOMT=1
          ELSEIF (DABS(DARG).LT.0.1D0) THEN
           DARG2=DARG*DARG
           DARG4=DARG2*DARG2
                DEXPDOMT=DARGEXPO(2,IOBSV)*
     &            DCMPLX(1.D0-DARG2/2.D0+DARG4/24.D0,
     &            DARG*(1.D0-DARG2/6.D0+DARG4/120.D0))
                IEXPDOMT=2
          ELSE
           DEXPDOMT=CDEXP(DCMPLX(0.D0,ARG))
                IEXPDOMT=3
          ENDIF

          DEXPOM=DARGEXPO(6,IOBSV)*DEXPDOMT

        ELSEIF (IFREQ2P.GT.0) THEN

          EXPOM=DARGEXPO(4,IOBSV)*DEXPOMT
          DARGEXPO(4,IOBSV)=EXPOM    !STORE IT HERE, WILL CHANGE FURTHER DOWN

        ELSE   !IFREQ2P

          EXPOM=CDEXP(DCMPLX(0.D0,RARG(4)*OM))
          DARGEXPO(4,IOBSV)=EXPOM    !STORE IT HERE, WILL CHANGE FURTHER DOWN

        ENDIF  !IFREQ2P


        IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN
                IFROB=1+NFREQ*(IOBSV-1)
             DO ICOMP=1,3
                AFREQ(ICOMP,IFROB)=AFREQ(ICOMP,IFROB)
     &          +DCMPLX(RARG(ICOMP))*EXPOM*REFLEC(ICOMP)*DT
             ENDDO
        ENDIF   !XIANF

      IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN
      IF (IWFILINT.NE.0) THEN
      IF (IWFILINT.EQ.-ISOUR) THEN

          IFREQ=1
          FILLT(1)=T
          FILLT(2)=X2
          FILLT(3)=Y2
          FILLT(4)=Z2
          FILLT(5)=RARG(1)
          FILLT(6)=RARG(2)
          FILLT(7)=RARG(3)
          FILLT(8)=RARG(4)
          FILLT(9)=MIN(RARG(5),1.D30)
            FILLT(10)=DREAL(EXPOM)
            FILLT(11)=DIMAG(EXPOM)
          FILLT(12)=IROI-1
          FILLT(13)=IOBSV
          FILLT(14)=IFREQ
          FILLT(15)=OBSV(2,IOBSV)
          FILLT(16)=OBSV(3,IOBSV)
          FILLT(17)=DUM1
          FILLT(18)=RARGOM-RARGOMO
          FILLT(19)=IEXPOMT
          FILLT(20)=IEXPDOMT
                IFROB=IFREQ+NFREQ*(IOBSV-1)
            FILLT(21)=
     &        DREAL(
     &          AFREQ(1,IFROB)*dCONJG(AFREQ(1,IFROB))
     &         +AFREQ(2,IFROB)*dCONJG(AFREQ(2,IFROB))
     &         +AFREQ(3,IFROB)*dCONJG(AFREQ(3,IFROB))
     &         )*SPECNOR
          FILLT(22)=TE
          CALL hfm(NIDSOURCE,FILLT)

      ELSEIF (ISOUR.EQ.IWFILINT.AND.IOBSV.EQ.1) THEN


         WRITE(LUNINT,*) X2
     &    ,(SNGL(RARG(1)),IC=1,3),SNGL(RARG(4)*OM),SNGL(RARG(5))
     &    ,SNGL(DREAL(EXPOM)),SNGL(DIMAG(EXPOM))
     &    ,SNGL(DREAL(RARG(1)*EXPOM)),SNGL(DIMAG(RARG(1)*EXPOM))
     &    ,SNGL(DREAL(RARG(2)*EXPOM)),SNGL(DIMAG(RARG(2)*EXPOM))
     &    ,SNGL(DREAL(RARG(3)*EXPOM)),SNGL(DIMAG(RARG(3)*EXPOM))

      ENDIF
      ENDIF
      ENDIF   !XIANF

C--- LOOP OVER ALL FREQUENCES

          DO IFREQ=2,NFREQ


          IF    (IFREQ2P.GT.2) THEN

         EXPOM=EXPOM*DEXPOM

          ELSEIF(IFREQ2P.EQ.2) THEN

         EXPOM=EXPOM*EXPOM

          ELSE

         EXPOM=CDEXP(DCMPLX(0.D0,RARG(4)*FREQ(IFREQ)/HBAREV1))

          ENDIF



          IF (X2.GE.XIANF.AND.X2.LE.XIEND) THEN
                IFROB=IFREQ+NFREQ*(IOBSV-1)
         DO ICOMP=1,3
                AFREQ(ICOMP,IFROB)=AFREQ(ICOMP,IFROB)
     &          +DCMPLX(RARG(ICOMP))*EXPOM*REFLEC(ICOMP)*DT
         ENDDO

         IF (IWFILINT.EQ.-ISOUR) THEN

             FILLT(1)=T
             FILLT(2)=X2
             FILLT(3)=Y2
             FILLT(4)=Z2
             FILLT(5)=RARG(1)
             FILLT(6)=RARG(2)
             FILLT(7)=RARG(3)
             FILLT(8)=RARG(4)
             FILLT(9)=MIN(RARG(5),1.D30)
                  FILLT(10)=DREAL(EXPOM)
                  FILLT(11)=DIMAG(EXPOM)
             FILLT(12)=IROI-1
             FILLT(13)=IOBSV
             FILLT(14)=IFREQ
             FILLT(15)=OBSV(2,IOBSV)
             FILLT(16)=OBSV(3,IOBSV)
             FILLT(17)=DUM1
             FILLT(18)=RARGOM-RARGOMO
             FILLT(19)=IEXPOMT
             FILLT(20)=IEXPDOMT
                IFROB=IFREQ+NFREQ*(IOBSV-1)
            FILLT(21)=
     &        DREAL(
     &          AFREQ(1,IFROB)*dCONJG(AFREQ(1,IFROB))
     &         +AFREQ(2,IFROB)*dCONJG(AFREQ(2,IFROB))
     &         +AFREQ(3,IFROB)*dCONJG(AFREQ(3,IFROB))
     &         )*SPECNOR
             FILLT(22)=TE
             CALL hfm(NIDSOURCE,FILLT)

         ENDIF
          ENDIF   !XIANF

          ENDDO   !LOOP OVER ALL FREQUENCES


C COMPLEX PART OF INTEGRAND }

          DARGEXPO(1,IOBSV)=DEXPOMT
          DARGEXPO(2,IOBSV)=DEXPDOMT
          DARGEXPO(3,IOBSV)=DCMPLX(RARGOM,RARGDOM)
          DARGEXPO(6,IOBSV)=DEXPOM

      ENDDO !IOBSV

C CONTRIBUTION OF TIME STEP TO SYNCHROTRON RADIATION }

      TS=TS+DT
      ITIM0=TS*DTIM01+1

C--- END OF LOOP OVER TIME STEPS

      IF (X2.LT.XENDSOU)  GOTO 1000

C- STORE NUMBER OF POINTS FOR INTEGRATION

      IPOISOU(ISOUR)=IZAEHL
      IZTOT(ISOUR)=IZTOT(ISOUR)+IZAEHL

      DO IOBSV=1,NOBSV
      DO IFREQ=1,NFREQ


          ILIOBFR=ISOUR+NSOURCE*(IOBSV-1+NOBSV*(IFREQ-1))
          IFROB=IFREQ+NFREQ*(IOBSV-1)
          IOBFR=IOBSV+NOBSV*(IFREQ-1)

          IF (IPOLA.EQ.0) THEN

            SPEC(ILIOBFR)=
     &        DREAL(
     &          AFREQ(1,IFROB)*CONJG(AFREQ(1,IFROB))
     &         +AFREQ(2,IFROB)*CONJG(AFREQ(2,IFROB))
     &         +AFREQ(3,IFROB)*CONJG(AFREQ(3,IFROB))
     &         )*SPECNOR

            REAIMA(1,1,IOBFR)=REAIMA(1,1,IOBFR)+
     &        DREAL(AFREQ(1,IFROB))
            REAIMA(2,1,IOBFR)=REAIMA(2,1,IOBFR)+
     &        DREAL(AFREQ(2,IFROB))
            REAIMA(3,1,IOBFR)=REAIMA(3,1,IOBFR)+
     &        DREAL(AFREQ(3,IFROB))

            REAIMA(1,2,IOBFR)=REAIMA(1,2,IOBFR)+
     &        DIMAG(AFREQ(1,IFROB))
            REAIMA(2,2,IOBFR)=REAIMA(2,2,IOBFR)+
     &        DIMAG(AFREQ(2,IFROB))
            REAIMA(3,2,IOBFR)=REAIMA(3,2,IOBFR)+
     &        DIMAG(AFREQ(3,IFROB))

          ELSE    !IPOLA

          APOL=
     &          AFREQ(1,IFROB)*CONJG(VPOLA(1))
     &         +AFREQ(2,IFROB)*CONJG(VPOLA(2))
     &         +AFREQ(3,IFROB)*CONJG(VPOLA(3))

            SPEC(ILIOBFR)=
     &        DREAL(APOL*CONJG(APOL))*SPECNOR
            REAIMA(1,1,IOBFR)=REAIMA(1,1,IOBFR)+
     &        DREAL(AFREQ(1,IFROB))
            REAIMA(2,1,IOBFR)=REAIMA(2,1,IOBFR)+
     &        DREAL(AFREQ(2,IFROB))
            REAIMA(3,1,IOBFR)=REAIMA(3,1,IOBFR)+
     &        DREAL(AFREQ(3,IFROB))

            REAIMA(1,2,IOBFR)=REAIMA(1,2,IOBFR)+
     &        DIMAG(AFREQ(1,IFROB))
            REAIMA(2,2,IOBFR)=REAIMA(2,2,IOBFR)+
     &        DIMAG(AFREQ(2,IFROB))
            REAIMA(3,2,IOBFR)=REAIMA(3,2,IOBFR)+
     &        DIMAG(AFREQ(3,IFROB))

          ENDIF   !IPOLA

      IF (ISTOKES.NE.0) THEN

          APOLH=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(1,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(1,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(1,3))

          APOLR=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(2,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(2,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(2,3))

          APOLL=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(3,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(3,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(3,3))

          APOL45=
     &          AFREQ(1,IFROB)*CONJG(VSTOKES(4,1))
     &         +AFREQ(2,IFROB)*CONJG(VSTOKES(4,2))
     &         +AFREQ(3,IFROB)*CONJG(VSTOKES(4,3))

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
     &                          STOK1*SPECNOR

            STOKES(2,IOBFR)=STOKES(2,IOBFR)+
     &                          STOK2*SPECNOR

            STOKES(3,IOBFR)=STOKES(3,IOBFR)+
     &                          STOK3*SPECNOR

            STOKES(4,IOBFR)=STOKES(4,IOBFR)+
     &                          STOK4*SPECNOR

      ENDIF !ISTOKES

        ENDDO !IFREQ

        ILIOB=ISOUR+NSOURCE*(IOBSV-1)
        SPECPOW(ILIOB)=SPECPOW(ILIOB)
     &                *ECHARGE1/16.D0/PI1/PI1/EPS01/CLIGHT1
     &                *DMYCUR     !NUMBER OF e-

        ENDDO !IOBSV

      IF (IWFILINT.LT.0.AND.ISOUR.EQ.NSOURCE) THEN
                CALL MHROUT(NIDSOURCE,ICYCLE,' ')
                CALL hdeletm(NIDSOURCE)
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       SR SOUINTLIN, SOURCE:', ISOUR
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       ROIs:'
      WRITE(LUNGFO,*)

      DO IROI=1,NROIA
               WRITE(LUNGFO,*)
     &           IROI,SNGL(ROIX(IROI)),SNGL(ROIP(IROI)),IPOIROI(IROI+1)
      ENDDO
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'       TOTAL NUMBER OF STEPS:',IZAEHL
      WRITE(LUNGFO,*)'       (controlled by NLPOI and namelist $ROIN)'

        IF (ISOUR.EQ.1.and.nsource.gt.1) THEN
                WRITE(6,*)' '
                WRITE(6,*)' '
                WRITE(6,*)' sources treated so far:'
                WRITE(6,*)' '
      ENDIF

      IF (ISOUR.EQ.NSOURCE) THEN
          DEALLOCATE(DWTRA)
      ENDIF

      RETURN
      END
