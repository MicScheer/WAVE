*CMZ :  4.00/15 17/05/2022  08.32.28  by  Michael Scheer
*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 07/12/2021  18.47.10  by  Michael Scheer
*CMZ :  3.02/04 18/03/2015  09.49.00  by  Michael Scheer
*CMZ :  3.02/03 04/11/2014  16.10.50  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.67/02 26/04/2012  14.52.47  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  10.38.44  by  Michael Scheer
*CMZ :  2.64/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.63/05 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.62/04 03/01/2008  17.15.05  by  Michael Scheer
*CMZ :  2.47/12 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.36/00 08/11/2001  11.40.47  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  12.21.53  by  Michael Scheer
*CMZ :  2.15/00 16/05/2000  13.53.47  by  Michael Scheer
*CMZ :  2.13/11 20/03/2000  13.17.03  by  Michael Scheer
*CMZ :  2.13/02 09/12/99  10.41.02  by  Michael Scheer
*CMZ :  2.12/00 02/06/99  12.03.53  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.01/09 19/10/95  15.35.33  by  Michael Scheer
*CMZ : 00.01/07 09/03/95  12.30.05  by  Michael Scheer
*CMZ : 00.01/06 13/02/95  14.15.56  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.22.28  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.30  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HTRACK
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

C--- HISTOGRAM FOR REFERENCE ORBIT AND MAGNETIC FIELD

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEND.

*KEEP,track.
      include 'track.cmn'
*KEEP,wbtab.
      include 'wbtab.cmn'
*KEEP,berror.
      include 'berror.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: XPHANA,PHPHANA,RESPHANA,
     &  pherr
      DOUBLE PRECISION XLENDEV2,A0,A1,CHI2
      DOUBLE PRECISION XPAR(3),YPAR(3),A(3),YPPAR(3),XOPT,AOPT,XG,XG2
      DOUBLE PRECISION PHASEMEAN,OMEGAR,FREQR,DXG,
     &  erra,errb

      INTEGER IFAIL,IXG,J,I1,IXG1,IXG2,IDIXG,nallo

      INTEGER I,NBIN
      INTEGER ICYCLE,IFILL
      REAL*4 XI,XE,X
      REAL*8 Z,Y,BX,BY,BZ,AX,AY,AZ

      INTEGER NTUPP,NPHANA
      PARAMETER (NTUPP=15,NPHANA=3)
      REAL*8 TUP(NTUPP),TUPPHANA(NPHANA)

      CHARACTER(2) CHTAGS(NTUPP)
      data chtags/'x','y','z','bx','by','bz','vx','vy','vz'
     &  ,'ax','ay','az','t','ph','pr'/

      CHARACTER(2) CHPHANA(NPHANA)
      DATA CHPHANA/'x','ph','dp'/

      IF (XSTARTH.EQ.9999.) THEN
         XI=WTRA(1,1,1)-DS0/2.
      ELSE
         XI=XSTARTH-DS0/2.
      ENDIF !XSTARTH

      IF (XSTOPH.EQ.9999.) THEN
        XE=WTRA(1,1,NCO)+DS0/2.
      ELSE
        XE=XSTOPH+DS0/2.
      ENDIF !XSTARTH

      NBIN=(XE-XI)/DS0
      XE=DS0*NBIN+XI

      IF (IHTRACK.GT.0) THEN

        call hbook1m(IDTRACKZ,'Z OF REF. ORBIT',NBIN,XI,XE,VMX)
        call hbook1m(IDTRACKY,'Y OF REF. ORBIT',NBIN,XI,XE,VMX)

        call hbook1m(IDBX,'Bx',NBIN,XI,XE,VMX)
        call hbook1m(IDBY,'By',NBIN,XI,XE,VMX)
        call hbook1m(IDBZ,'Bz',NBIN,XI,XE,VMX)

        call hbook1m(IDAX,'Ax',NBIN,XI,XE,VMX)
        call hbook1m(IDAY,'Ay',NBIN,XI,XE,VMX)
        call hbook1m(IDAZ,'Az',NBIN,XI,XE,VMX)

      ENDIF !IHTRACK

      IF (IHTRSMP.EQ.0) IHTRSMP=1

      nallo=nco/ihtrsmp+1024
      CALL hbookm(NIDTRACK,'TRAJECTORY',NTUPP,'//WAVE',nallo,CHTAGS)
      if (abs(ihtrsmp).eq.9999) then
        nallo=NPWBTAB+1024
      else
        nallo=nint(nco*ds0/(abs(ihtrack)*0.0005))+1024
      endif
      CALL hbookm(NIDTRCKG,'TRAJECTORY ON GRID',NTUPP,'//WAVE',nallo,
     &  CHTAGS)
      IF (NCO.GE.3.AND.IBERROR.NE.0) THEN
        CALL hbookm(NIDPHANA,'PHASE ADVANCE',NPHANA,'//WAVE',nco+1024,CHPHANA)
      endif
      nallo=max((XSTOP-XSTART)/ZLENERR*2.,dble(nberror))+10

      ALLOCATE(XPHANA(nallo))
      ALLOCATE(PHPHANA(nallo))
      ALLOCATE(RESPHANA(nallo))
      ALLOCATE(pherr(nallo))

      xphana=0.0d0
      phphana=0.0d0
      resphana=0.0d0
      pherr=0.0d0

      XI=WTRA(1,1,1)-DS0
      IF (IHTRSMP.EQ.0) IHTRSMP=1

      DO I=1,NCO,IHTRSMP

        X=DS0*I+XI
        Y=WTRA(2,1,I)
        Z=WTRA(3,1,I)
        BX=WTRA(1,3,I)
        BY=WTRA(2,3,I)
        BZ=WTRA(3,3,I)
        AX=WTRA(1,4,I)
        AY=WTRA(2,4,I)
        AZ=WTRA(3,4,I)

        IF (IHTRACK.GT.0) THEN
          CALL hfillm(IDTRACKZ,X,0.,Z)
          CALL hfillm(IDTRACKY,X,0.,Y)
          CALL hfillm(IDBX,X,0.,BX)
          CALL hfillm(IDBY,X,0.,BY)
          CALL hfillm(IDBZ,X,0.,BZ)
          CALL hfillm(IDAX,X,0.,AX)
          CALL hfillm(IDAY,X,0.,AY)
          CALL hfillm(IDAZ,X,0.,AZ)
        ENDIF  !IHTRACK

        TUP(1)=WTRA(1,1,I)
        TUP(2)=WTRA(2,1,I)
        TUP(3)=WTRA(3,1,I)
        TUP(4)=WTRA(1,3,I)
        TUP(5)=WTRA(2,3,I)
        TUP(6)=WTRA(3,3,I)
        TUP(7)=WTRA(1,2,I)
        TUP(8)=WTRA(2,2,I)
        TUP(9)=WTRA(3,2,I)
        TUP(10)=WTRA(1,4,I)
        TUP(11)=WTRA(2,4,I)
        TUP(12)=WTRA(3,4,I)
        TUP(13)=WTIM0(I)
        TUP(14)=HTRA2(I)
        TUP(15)=WTRA2(I)
        CALL hfm(NIDTRACK,TUP)

      ENDDO !NCO

      IF (NCO.GE.3.AND.IBERROR.NE.0) THEN

        I1=1
        IFILL=0

        XLENDEV2=ZLENERR/2.D0*(NBERROR-1)/2.D0
        DO IXG=1,NBERROR+1

          XG=XCENERR-XLENDEV2+ZLENERR/2.D0*(IXG-1)

          XG2=XG*XG

          DO I=I1,NCO
            IF (WTRA(1,1,I).GE.XG) THEN
              I1=I
              GOTO 20
            ENDIF
          ENDDO

20        IF (I.EQ.1) THEN
            XPAR(1)=WTRA(1,1,1)
            XPAR(2)=WTRA(1,1,2)
            XPAR(3)=WTRA(1,1,3)
          ELSE IF (I.LT.NCO) THEN
            XPAR(1)=WTRA(1,1,I-1)
            XPAR(2)=WTRA(1,1,I)
            XPAR(3)=WTRA(1,1,I+1)
          ELSE
            XPAR(1)=WTRA(1,1,NCO-2)
            XPAR(2)=WTRA(1,1,NCO-1)
            XPAR(3)=WTRA(1,1,NCO)
          ENDIF

          IF (I.EQ.1) THEN
            YPAR(1)=HTRA2(1)
            YPAR(2)=HTRA2(2)
            YPAR(3)=HTRA2(3)
          ELSE IF (I.LT.NCO) THEN
            YPAR(1)=HTRA2(I-1)
            YPAR(2)=HTRA2(I)
            YPAR(3)=HTRA2(I+1)
          ELSE
            YPAR(1)=HTRA2(NCO-2)
            YPAR(2)=HTRA2(NCO-1)
            YPAR(3)=HTRA2(NCO)
          ENDIF
          CALL UTIL_PARABEL(XPAR,YPAR,A,YPPAR,XOPT,AOPT,IFAIL)

          IFILL=IFILL+1
          if (ifill.le.nallo) then
            XPHANA(IFILL)=XG
            PHPHANA(IFILL)=A(1)+A(2)*XG+A(3)*XG2
            !PHPHANA(IFILL)=phphana(ifill)-phphana(1)
          else
            write(lungfo,*)'*** ERROR IN HTRACK: ARRAY SIZE TO SMALL'
            write(lungfo,*)'*** Check namelist BERRORN'
            write(6,*)'*** WARNING IN HTRACK: ARRAY SIZE TO SMALL'
            write(6,*)'*** Check namelist BERRORN'
            stop '*** Program WAVE aborted ***'
          endif

        ENDDO  !NCO

C FIT A STRAIGHT LINE

cmsh 20150318 99      CALL DLSQP1(IFILL,XPHANA,PHPHANA,A0,A1,CHI2,IFAIL)
99      CALL util_straight_line_fit(ifill,xphana,phphana,pherr,
     &    a1,a0,chi2,erra,errb,ifail)

        IF (IFAIL.NE.0) THEN
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN HTRACK: STRAIGHT '//
     &      'LINE FIT OF PHASE ADVANCE FAILED ***'
          WRITE(6,*)
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN HTRACK: STRAIGHT '//
     &      'LINE FIT OF PHASE ADVANCE FAILED ***'
          WRITE(LUNGFO,*)
        ENDIF  !IFAIL

        PHASEMEAN=A1*ZLENERR/CLIGHT1
        OMEGAR=2.D0*PI1/PHASEMEAN
        FREQR=OMEGAR*HBAREV1

        RESRMS=0.0D0
        DO I=1,IFILL
          RESPHANA(I)=
     &      (PHPHANA(I)-(A0+A1*XPHANA(I)))/CLIGHT1*OMEGAR/PI1/2.0d0 !
          RESRMS=RESRMS+RESPHANA(I)**2
          TUPPHANA(1)=XPHANA(I)
          TUPPHANA(2)=PHPHANA(I)/CLIGHT1*OMEGAR/PI1/2.0d0
          TUPPHANA(3)=RESPHANA(I)
          CALL hfm(NIDPHANA,TUPPHANA)
        ENDDO  !IFILL
        RESRMS=SQRT(RESRMS/IFILL)

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'      Phase error analysis of magnetic field'
        WRITE(LUNGFO,*)'      (phase is sampled according to ZLENERR/2'
        WRITE(LUNGFO,*)'       in namelist BERRORN):'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      first harmonic (eV) corresponding to mean phase advance:'
        WRITE(LUNGFO,*)
     &    '      ',FREQR
        WRITE(LUNGFO,*)
     &    '      rms phase error for 1. harmonical (rad, deg.):'
        WRITE(LUNGFO,*)
     &    '      ',SNGL(RESRMS*PI1*2.0d0),' ',SNGL(RESRMS*360.D0)
        WRITE(LUNGFO,*)
     &    '      rel. rms phase error for 1. harmonical (per pole, per period):'
        WRITE(LUNGFO,*)
     &    '      ',SNGL(RESRMS/2.0d0),' ',SNGL(RESRMS)
        WRITE(LUNGFO,*)

      ENDIF !NCO.GE.3.AND.IBERROR.NE.0

      IF (NCO.GE.3) THEN

        I1=1

        IF (IABS(IHTRACK).EQ.9999) THEN
          IDIXG=1
          IF (BTABS.EQ.9999.) BTABS=WSXYZ(1,1)
          IF (BTABE.EQ.9999.) BTABE=WSXYZ(1,NCO) !DO NOT USE XSTOP
          IF (NPWBTAB.EQ.9999.) NPWBTAB=NCO
          DXG=(BTABE-BTABS)/(NPWBTAB-1)
          IXG1=NINT(BTABS/DXG)
          IXG2=NINT(BTABE/DXG)
        ELSE
          DXG=1.D0/2000.D0
          IXG1=NINT(XSTART/DXG)
          IXG2=NINT(XSTOP/DXG)
          IDIXG=IABS(IHTRACK)
        ENDIF

        DO IXG=IXG1,IXG2,IDIXG

          XG=FLOAT(IXG)*DXG
          XG2=XG*XG
          TUP(1)=XG

          DO I=I1,NCO
            IF (WTRA(1,1,I).GE.XG) THEN
              I1=I
              GOTO 10
            ENDIF
          ENDDO
10        IF (I.EQ.1) THEN
            XPAR(1)=WTRA(1,1,1)
            XPAR(2)=WTRA(1,1,2)
           ELSEIF (I.GT.NCO) THEN
            XPAR(1)=WTRA(1,1,NCO-1)
            XPAR(2)=WTRA(1,1,NCO)
          ELSE
            XPAR(1)=WTRA(1,1,I-1)
            XPAR(2)=WTRA(1,1,I)
          ENDIF

          DO J=2,3
            IF (I.EQ.1) THEN
              YPAR(1)=WTRA(J,1,1)
              YPAR(2)=WTRA(J,1,2)
           ELSEIF (I.GT.NCO) THEN
              YPAR(1)=WTRA(J,1,NCO-1)
              YPAR(2)=WTRA(J,1,NCO)
          ELSE
              YPAR(1)=WTRA(J,1,I-1)
              YPAR(2)=WTRA(J,1,I)
          ENDIF
            TUP(J)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))
          ENDDO   !J

          DO J=1,3
            IF (I.EQ.1) THEN
              YPAR(1)=WTRA(J,3,1)
              YPAR(2)=WTRA(J,3,2)
           ELSEIF (I.GT.NCO) THEN
              YPAR(1)=WTRA(J,3,NCO-1)
              YPAR(2)=WTRA(J,3,NCO)
            ELSE
              YPAR(1)=WTRA(J,3,I-1)
              YPAR(2)=WTRA(J,3,I)
            ENDIF
            TUP(3+J)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))
          ENDDO   !J

          DO J=1,3
            IF (I.EQ.1) THEN
              YPAR(1)=WTRA(J,2,1)
              YPAR(2)=WTRA(J,2,2)
           ELSEIF (I.GT.NCO) THEN
              YPAR(1)=WTRA(J,2,NCO-1)
              YPAR(2)=WTRA(J,2,NCO)
            ELSE
              YPAR(1)=WTRA(J,2,I-1)
              YPAR(2)=WTRA(J,2,I)
            ENDIF
            TUP(6+J)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))
          ENDDO   !J

          DO J=1,3
            IF (I.EQ.1) THEN
              YPAR(1)=WTRA(J,4,1)
              YPAR(2)=WTRA(J,4,2)
           ELSEIF (I.GT.NCO) THEN
              YPAR(1)=WTRA(J,4,NCO-1)
              YPAR(2)=WTRA(J,4,NCO)
            ELSE
              YPAR(1)=WTRA(J,4,I-1)
              YPAR(2)=WTRA(J,4,I)
            ENDIF
            TUP(9+J)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))
          ENDDO   !J

          IF (I.EQ.1) THEN
            YPAR(1)=WTIM0(1)
            YPAR(2)=WTIM0(2)
           ELSEIF (I.GT.NCO) THEN
              YPAR(1)=WTIM0(NCO-1)
              YPAR(2)=WTIM0(NCO)
          ELSE
            YPAR(1)=WTIM0(I-1)
            YPAR(2)=WTIM0(I)
          ENDIF
          TUP(13)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))

          IF (I.EQ.1) THEN
            YPAR(1)=HTRA2(1)
            YPAR(2)=HTRA2(2)
          ELSE IF (I.GT.NCO) THEN
            YPAR(1)=HTRA2(NCO-1)
            YPAR(2)=HTRA2(NCO)
          ELSE
            YPAR(1)=HTRA2(I-1)
            YPAR(2)=HTRA2(I)
          ENDIF
          TUP(14)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))

          IF (I.EQ.1) THEN
            YPAR(1)=WTRA2(1)
            YPAR(2)=WTRA2(2)
          ELSE IF (I.GT.NCO) THEN
            YPAR(1)=WTRA2(NCO-1)
            YPAR(2)=WTRA2(NCO)
          ELSE
            YPAR(1)=WTRA2(I-1)
            YPAR(2)=WTRA2(I)
          ENDIF
          TUP(15)=YPAR(1)+(YPAR(2)-YPAR(1))/(XPAR(2)-XPAR(1))*(XG-XPAR(1))
          CALL hfm(NIDTRCKG,TUP)

        ENDDO  !IXG

      ELSE  !NCO.GE.3
        WRITE(6,*)
     &    'WARNING IN HTRACK: NCO.LT.3, NO NTUPLE OF TRACKJETROY ON GRID'
      ENDIF !NCO.GE.3

      IF (IHTRACK.GT.0) THEN
        CALL MHROUT(IDTRACKZ,ICYCLE,' ')
        CALL MHROUT(IDTRACKY,ICYCLE,' ')
        CALL MHROUT(IDBX,ICYCLE,' ')
        CALL MHROUT(IDBY,ICYCLE,' ')
        CALL MHROUT(IDBZ,ICYCLE,' ')
        CALL MHROUT(IDAX,ICYCLE,' ')
        CALL MHROUT(IDAY,ICYCLE,' ')
        CALL MHROUT(IDAZ,ICYCLE,' ')

        CALL hdeletm(IDTRACKZ)
        CALL hdeletm(IDTRACKY)
        CALL hdeletm(IDBX)
        CALL hdeletm(IDBY)
        CALL hdeletm(IDBZ)
        CALL hdeletm(IDAX)
        CALL hdeletm(IDAY)
        CALL hdeletm(IDAZ)

      ENDIF

      DEALLOCATE(XPHANA)
      DEALLOCATE(PHPHANA)
      DEALLOCATE(RESPHANA)

      RETURN
      END
