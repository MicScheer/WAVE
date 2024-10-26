*CMZ :  4.01/04 06/12/2023  10.37.05  by  Michael Scheer
*CMZ :  4.01/02 26/04/2023  07.35.28  by  Michael Scheer
*CMZ :  4.00/15 28/04/2022  09.27.00  by  Michael Scheer
*CMZ :  4.00/07 09/07/2020  12.37.21  by  Michael Scheer
*CMZ :  3.06/00 26/02/2019  14.12.24  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.03/04 03/08/2017  15.33.14  by  Michael Scheer
*CMZ :  3.03/02 31/08/2016  15.06.48  by  Michael Scheer
*CMZ :  3.01/02 02/09/2013  15.00.48  by  Michael Scheer
*CMZ :  3.01/01 19/07/2013  16.20.01  by  Michael Scheer
*CMZ :  3.00/02 10/04/2013  09.27.41  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.67/02 18/04/2012  14.13.08  by  Michael Scheer
*CMZ :  2.66/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.66/07 02/03/2010  09.40.31  by  Michael Scheer
*CMZ :  2.62/02 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.16/08 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  16.58.05  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.44.33  by  Michael Scheer
*CMZ : 00.01/07 28/02/95  14.22.48  by  Michael Scheer
*CMZ : 00.01/06 17/02/95  12.25.05  by  Michael Scheer
*CMZ : 00.01/04 19/01/95  09.36.21  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.45.15  by  Michael Scheer
*CMZ : 00.01/01 21/09/94  17.30.43  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.46.58  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.41  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BELLIP(X,Y,Z,BX,BY,BZ,AX,AY,AZ)

      use f1k

*KEEP,GPLHINT.
*KEND.

      IMPLICIT NONE

      DOUBLE PRECISION BX,BY,BZ,AX,AY,AZ,X,Y,Z
      DOUBLE PRECISION WLEN1,PARK,B0EFF,EHARM1,parkv,parkh
      DOUBLE PRECISION RHV,pbue(13),f1

      INTEGER,save :: ICAL=0
      INTEGER i

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,ellip.
      include 'ellip.cmn'
*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEND.


      IF (ICAL.EQ.0) THEN

        park=parkell

        if (nharmell.ne.0.and.harmell.ne.0.0d0) then
          if (harmell.eq.-9999.0d0) then
c            if (ifreq2p.eq.1) then
c              harmell=freqlow
c            else
              harmell=(freqlow+freqhig)/2.0d0
c            endif
          endif
          if (harmell.lt.0.0d0) then
            harmell=-wtoe1/harmell
          endif
          WLEN1=wtoe1/abs(harmell/nharmell)
          park=2.0d0*(wlen1/(xlellip*1.0D9/2.0d0/DMYGAMMA**2)-1.0d0)
          if (park.lt.0.0d0) then
            write(6,*)
     &        '*** Error in BELLIP:'
            write(6,*)
     &        'Inconsistent values of NHARMELL, HARMELL, and XLELLIP'
            write(6,*)' '
            write(lungfo,*)
     &        '*** Error in BELLIP:'
            write(lungfo,*)
     &        'Inconsistent values of NHARMELL, HARMELL, and XELLIP'
            write(lungfo,*)' '
            stop
          endif
          park=sqrt(park)
          parkell=park
        endif

        IF (parkell.NE.0.0) THEN

          B0EFF=parkell/(echarge1*XLELLIP/(2.*PI1*EMASSKG1*CLIGHT1))

          if (b0elliph.eq.0.0d0.and.b0ellipv.ne.0d0) then
            b0ellipv=b0ellipv/abs(b0ellipv)*b0eff
          else if (b0ellipv.eq.0.0d0.and.b0elliph.ne.0d0) then
            b0elliph=b0elliph/abs(b0elliph)*b0eff
          else

            rhv=b0elliph/b0ellipv

            b0elliph=b0eff/sqrt(1.0d0+1.0d0/rhv**2)*b0elliph/abs(b0elliph)
            b0ellipv=b0elliph/rhv

          endif

        ENDIF

        IF (B0ELLIPV.NE.0.0) THEN
          RHV=B0ELLIPH/B0ELLIPV
        ELSE
          RHV=0.0D0
        ENDIF

        XLENELL=DABS((PERELLIP+ELLSHFT)*XLELLIP)
        IF (XLELLIP.NE.0.0) then
          XKELLIP=2.D0*PI1/XLELLIP
          zampell=b0ellipv*clight1/emom/xkellip**2
          yampell=b0elliph*clight1/emom/xkellip**2
        endif

        B0EFF=DSQRT(B0ELLIPH**2+B0ELLIPV**2)
        PARKv=ECHARGE1*DABS(B0elliph)*XLELLIP/(2.*PI1*EMASSKG1*CLIGHT1)
        PARKh=ECHARGE1*DABS(B0ellipv)*XLELLIP/(2.*PI1*EMASSKG1*CLIGHT1)
        PARK=ECHARGE1*DABS(B0EFF)*XLELLIP/(2.*PI1*EMASSKG1*CLIGHT1)
        WLEN1=(1.0d0+PARK**2/2.0d0)/2.0d0/DMYGAMMA**2*XLELLIP*1.0D9
        IF (WLEN1.NE.0.0) EHARM1=WTOE1/WLEN1

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     SR BELLIP, Parameters of elliptical undulator:'
        WRITE(LUNGFO,*)
     &    '     horizontal and vertical peak field [T]:'
     &    ,SNGL(B0ELLIPH),SNGL(B0ELLIPV)
        WRITE(LUNGFO,*)
     &    '     number of periods, period and device length [m]:'
        WRITE(LUNGFO,*)
     &    '     ',SNGL(PERELLIP),SNGL(XLELLIP),SNGL(XLENELL)
        WRITE(LUNGFO,*)
     &    '     shift parameter [periods], B_h/B_v:    '
     &    ,SNGL(ELLSHFT),SNGL(RHV)
        WRITE(LUNGFO,*)
     &    '     B0_eff [T], deflection parameter K:'
     &    ,SNGL(B0EFF),SNGL(PARK)
        WRITE(LUNGFO,*)
     &    '     Approximated path hori. and. vert. amplitudes [m]:',
     &    SNGL(zampell),sngl(yampell)
        WRITE(LUNGFO,*)
     &    '     Hori. and. vert. deflection angle [rad]:',
     &    SNGL(parkv/dmygamma),sngl(parkh/dmygamma)
        WRITE(LUNGFO,*)
     &    '     1. harmonical [nm] and [eV], omega [1/s]:',
     &    sngl(wlen1),SNGL(EHARM1),sngl(eharm1/hbarev1)
        WRITE(LUNGFO,*)
     &    '     critical energy [eV]:',SNGL(ecdipev1*DABS(B0EFF)*DMYENERGY**2)
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     longitudinal position of device center [m]:',xcenell
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     taper factor to compensate energy-loss:',elltap
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     Estimate of 1. harm. on-axis flux-density in 10m / mm**2 / 0.1%BW'
        do i=1,1000
          if (park.ge.f1kk(i)) cycle
          f1=f1kf(i)
          exit
        enddo
        WRITE(LUNGFO,*)
     &    '     ',
     &    sngl(1.744e14*perellip**2*dmyenergy**2*dmycur*f1/100.0d0)
        WRITE(LUNGFO,*)

        ICAL=1

      ENDIF

      pbue=PMAG(1:13,NMGSQP)

      PMAG(1,NMGSQP)=park
      PMAG(2,NMGSQP)=b0ellipv
      PMAG(3,NMGSQP)=b0elliph
cerror 31.8.2016      PMAG(4,NMGSQP)=ellshft*xlellip
      PMAG(4,NMGSQP)=ellshft
      PMAG(5,NMGSQP)=xcenell
      PMAG(6,NMGSQP)=xlellip
      PMAG(7,NMGSQP)=perellip
      PMAG(8,NMGSQP)=0.0d0
      pmag(11,nmgsqp)=elltap

      call bue(x,y,z,bx,by,bz,ax,ay,az,nmgsqp)

      PMAG(1:13,NMGSQP)=pbue

      RETURN
      END
