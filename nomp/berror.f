*CMZ :  4.00/15 06/06/2022  08.45.26  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.63/02 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.03.48  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.45.16  by  Michael Scheer
*CMZ : 00.01/12 27/09/96  15.28.21  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.23  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.59  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BERROR(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

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

C CALCULATES MAGNETIC FIELD ERRORS; FOR FORMULAS SEE BHALBASY

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,ellana.
      include 'ellana.cmn'
*KEEP,berror.
      include 'berror.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,halbach.
      include 'halbach.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL,NPERMXP,IX,NTOTIN,NTOT2IN,kfibona,LFIBO,MODANA,ixz,lunerr
c      PARAMETER (NPERMXP=1024)

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  XKX,YKY,ZKZ,DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &  ,BXH,BYH,BZH,AXH,AYH,AZH,X,TOTLEN,TOTLEN2,SUM,RMS,POLWID,
     &  xxran,xxranz,byellmx,bzellmx

      real rr(2)
      REAL, dimension (:), allocatable :: xran,yran,fibo

      DATA ICAL/0/

      save

      BxOUT=0.0d0
      ByOUT=0.0d0
      BzOUT=0.0d0
      axOUT=0.0d0
      ayOUT=0.0d0
      azOUT=0.0d0

      IF (ICAL.EQ.0) THEN

        IF (XLENERR.LT.0.0D0) THEN
          MODANA=1
          call bellana(0.0d0,0.0d0,0.0d0,bxh,byh,bzh,axh,ayh,azh)
          byellmx=abs(byh)
          call bellana(zlellana/4.0d0,0.0d0,0.0d0,bxh,byh,bzh,axh,ayh,azh)
          bzellmx=abs(bzh)
        ELSE
          MODANA=0
        ENDIF

        IF (XSTART.EQ.9999.0d0) RETURN
        IF (XSTOP.EQ.9999.0d0) RETURN

        if (nberror.eq.-1) then

          kfibona=1

          if (modana.eq.0) then
            NBERROR=(XSTOP-XSTART-ZLENERR)/(ZLENERR/2.0d0)
          else
            NBERROR=(NPERELLA-2)*2
          endif

          NBERROR=(NBERROR/2)*2-1

          allocate(XRAN(nberror),YRAN(nberror),FIBO(nberror))

c          IF (NBERROR.GT.NPERMXP) THEN
c            WRITE(6,*)'*** ERROR IN BERROR: DIMENSION NPERMXP EXCEEDED'
c            WRITE(6,*)
c     &        '*** USE FEWER PERIODS IN NAMELIST BERROR OR INCREASE PARAMETER'
c            STOP '*** PROGRAM WAVE ABORTED ***'
c          ENDIF

          call fibonacci(nberror,xran,etafibo,icutfibo)

          DO IX=1,NBERROR
            IF (XRAN(IX).GT.0.) THEN
              YRAN(IX)=1.
              FIBO(IX)=1.
            ELSE IF (XRAN(IX).LT.0.) THEN
              YRAN(IX)=-1.
              FIBO(IX)=-1.
            ELSE
              FIBO(IX)=0.
              YRAN(IX)=0.
            ENDIF
            IF (XRAN(IX).NE.0) LFIBO=IX
          ENDDO

          IF (NBERRMOD.EQ.1) THEN
            SUM=0.0D0
            DO IX=2,LFIBO
              FIBO(IX-1)=YRAN(IX-1)-YRAN(IX)
              XRAN(IX-1)=IX-1
              SUM=SUM+FIBO(IX-1)
            ENDDO
            SUM=SUM+FIBO(LFIBO)
            FIBO(1)=FIBO(1)-SUM
          ELSE IF (NBERRMOD.NE.2) THEN
            if (lfibo.ge.nberror) then
              lfibo=nberror-1
              XRAN(nberror)=0.
              YRAN(nberror)=0.
              fibo(nberror)=0.
            endif
            SUM=0.0D0
            DO IX=2,LFIBO-1
              FIBO(IX)=YRAN(IX)-YRAN(IX-1)/2.-YRAN(IX+1)/2.
              XRAN(IX-1)=IX-1
              XRAN(IX+1)=IX+1
              SUM=SUM+FIBO(IX)
            ENDDO
            FIBO(1)=-YRAN(2)/2.
            FIBO(NBERROR)=-SUM
          ENDIF   !NBERRMOD

        else

          kfibona=0

        endif

        IF (NBERROR.EQ.-9999) THEN
c          if (khalbasy.ne.0) then
c            nberror=int(ahwpol-2)
c          else if (khalba.ne.0) then
c            nberror=int(perhal*2.0d0)-2
c          else
            NBERROR=(XSTOP-XSTART-ZLENERR)/(ZLENERR/2.0d0)-2
c          endif
c          if (nberrmod.ne.0) then
c            NBERROR=(NBERROR/2)*2-1
c          endif
        ENDIF

        if (kfibona.eq.0) allocate(XRAN(nberror),YRAN(nberror))

c        IF (NBERROR.GT.NPERMXP) THEN
c          WRITE(6,*)'*** ERROR IN BERROR: DIMENSION NPERMXP EXCEEDED'
c          WRITE(6,*)
c     &      '*** USE FEWER PERIODS IN NAMELIST BERROR OR INCREASE PARAMETER'
c          STOP '*** PROGRAM WAVE ABORTED ***'
c        ENDIF

c        IF (NBERROR.LT.4) THEN
c          WRITE(LUNGFO,*)
c     &      '*** WARNING IN BERROR: LESS THAN THREE POLES FOR FIELDERRORS SPECIFIED'
c          WRITE(LUNGFO,*)
c     &      '*** NBERROR SET TO THREE'
c          WRITE(LUNGFO,*)
c          WRITE(6,*)
c     &      '*** WARNING IN BERROR: LESS THAN THREE POLES FOR FIELDERRORS SPECIFIED'
c          WRITE(6,*)
c     &      '*** NBERROR SET TO THREE'
c          WRITE(6,*)
c          NBERROR=3
c        ENDIF

        POLWID=ZLENERR/2.0D0
        TOTLEN=NBERROR*POLWID
        TOTLEN2=TOTLEN/2.0d0

C--- K-VALUES

        XKERROR=0.0d0
        YKERROR=0.0d0
        ZKERROR=0.0d0

        IF (ZLENERR.NE.0.0d0) ZKERROR=2.0d0*PI1/ZLENERR
        IF (YLENERR.NE.0.0d0) YKERROR=2.0d0*PI1/YLENERR
        IF (XLENERR.NE.0.0d0) XKERROR=2.0d0*PI1/XLENERR

C--- ADJUST K-VALUES

        YKERROR=DSQRT(ZKERROR**2+XKERROR**2)
        YLENERR=2.0d0*PI1/YKERROR

        IF (KFIBONA.EQ.0) THEN

          if (nberrmod.lt.10) then
            IF (IBERRSEED.NE.0) THEN
              CALL RMARIN(IBERRSEED,NTOTIN,NTOT2IN) !CERN V113
            ENDIF
            CALL RNORML(XRAN,NBERROR,rr)
          else
            nberrmod=nberrmod-1
            open(newunit=lunerr,file='wave_berror.dat')
            close(lunerr)
          endif

          YRAN(1:nberror)=XRAN(1:nberror)

          iF (NBERRMOD.eq.0) THEN
c            nberrmod=nberrmod/3*3
            nberror=nberror/3*3
            DO IX=2,NBERROR-1,3
              XRAN(IX-1)=yRAN(IX)/2.
              XRAN(IX+1)=YRAN(IX)/2.
            ENDDO
          else if (NBERRMOD.eq.3) THEN
c            nberrmod=nberrmod/3*3
            nberror=nberror/4*4
            DO IX=2,NBERROR-2,4
              XRAN(IX-1)=yRAN(IX)
              XRAN(IX+1)=yRAN(IX)
              XRAN(IX+2)=-yRAN(IX)
            ENDDO
          else IF (NBERRMOD.EQ.1) THEN
            YRAN(NBERROR)=0.
            DO IX=2,NBERROR,2
              XRAN(IX-1)=YRAN(IX)
            ENDDO
          else IF (NBERRMOD.ne.2) THEN
            stop "*** Error in BERROR: Bad NBERRMOD ***"

          ENDIF   !NBERRMOD

          SUM=0.0d0
          DO IX=1,NBERROR-1,2
            SUM=SUM+XRAN(IX)
          ENDDO
          DO IX=2,NBERROR,2
            SUM=SUM+XRAN(IX)
          ENDDO
          RMS=0.0d0
          SUM=SUM/NBERROR
          DO IX=1,NBERROR
            RMS=RMS+(SUM-XRAN(IX))**2
          ENDDO
          RMS=SQRT(RMS/NBERROR)*B0ERROR

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      SUBROUTINE BERROR:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      NBERROR, IBERRSEED,NBERRMOD:'
     &      ,NBERROR, IBERRSEED,NBERRMOD
          WRITE(LUNGFO,*)'      B0ERROR,RMS:',SNGL(B0ERROR),SNGL(RMS)
          WRITE(LUNGFO,*)'      ZLENERR,XLENERR,XCENERR:'
     &      ,SNGL(ZLENERR),SNGL(XLENERR),SNGL(XCENERR)
          WRITE(LUNGFO,*)'      MODANA:',MODANA
          WRITE(LUNGFO,*)

        ELSE !KFIBONA

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      SUBROUTINE BERROR:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      Fibonacci-Sequence applied'
          WRITE(LUNGFO,*)'      ETAFIBO,ICUTFIBO:',ETAFIBO,ICUTFIBO
          WRITE(LUNGFO,*)'      see fibonacci.dat for sequence'
          WRITE(LUNGFO,*)'      NBERROR, NBERRMOD:',NBERROR, NBERRMOD
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'      B0ERROR:',SNGL(B0ERROR)
          WRITE(LUNGFO,*)'      ZLENERR,XLENERR,XCENERR:'
          WRITE(LUNGFO,*)'      ',ZLENERR,XLENERR,XCENERR
          WRITE(LUNGFO,*)'      MODANA:',MODANA
          WRITE(LUNGFO,*)

        ENDIF !KFIBONA

        ICAL=1

      ENDIF

      X=XIN-XCENERR+totlen2
      IX=X/POLWID+1

      IF (IX.LT.1.OR.IX.GT.NBERROR) THEN
        return
      else
        XXRAN=XRAN(IX)
      endif

      IF (KFIBONA.EQ.0) THEN

        IF (IX.LT.1.OR.IX.GT.NBERROR) THEN
          return
        else
          XXRAN=XRAN(IX)
        endif

      ELSE !KFIBONA

        IF (IX.LT.1.OR.IX.GT.NBERROR) THEN
          xxran=0.0d0
        else
          IF (FIBO(IX).NE.0) THEN
            IX=ABS(XRAN(IX))
            XXRAN=FIBO(IX)
          ELSE
            XXRAN=0.0
          ENDIF
        ENDIF

        IXz=(X-zlellana/4.0d0)/POLWID+1

        IF (IXz.LT.1.OR.IXz.GT.NBERROR) THEN
          xxranz=.0d0
        else
          IF (FIBO(IXz).NE.0) THEN
            IXz=ABS(XRAN(IXz))
            XXRANz=FIBO(IXz)
          ELSE
            XXRANz=0.0
          ENDIF
        ENDIF

      ENDIF !KFIBONA

      IF (MODANA.EQ.0) THEN

        XKX=XKERROR*(-ZIN)
        YKY=YKERROR*YIN
        ZKZ=ZKERROR*(XIN-XCENERR)

        DSNXKX=DSIN(XKX)
        DCSXKX=DCOS(XKX)
        DSHYKY=DSINH(YKY)
        DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
        DSNZKZ=DSIN(ZKZ)
        DCSZKZ=DCOS(ZKZ)

        BXH=-XKERROR/YKERROR*B0ERROR*xxran*DSNXKX*DSHYKY*DCSZKZ
        BYH=                 B0ERROR*xxran*DCSXKX*DCHYKY*DCSZKZ
        BZH=-ZKERROR/YKERROR*B0ERROR*xxran*DCSXKX*DSHYKY*DSNZKZ

        AXH=B0ERROR*xxran/ZKERROR*                DCSXKX*DCHYKY*DSNZKZ
        AYH=B0ERROR*xxran/ZKERROR*XKERROR/YKERROR*DSNXKX*DSHYKY*DSNZKZ
        AZH=0.0d0

        BZOUT=-BXH
        BYOUT=BYH
        BXOUT=BZH

        AZOUT=-AXH
        AYOUT=AYH
        AXOUT=AZH

      ELSE !MODANA

        if (xxran.ne.0.0d0) then

          XKX=XKERROR*(-ZIN)
          YKY=YKERROR*YIN
          ZKZ=ZKERROR*X

          DSNXKX=DSIN(XKX)
          DCSXKX=DCOS(XKX)
          DSHYKY=DSINH(YKY)
          DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
          DSNZKZ=DSIN(ZKZ)
          DCSZKZ=DCOS(ZKZ)

          BXH=-XKERROR/YKERROR*B0ERROR*xxran*DSNXKX*DSHYKY*DCSZKZ
          BYH=                 B0ERROR*xxran*DCSXKX*DCHYKY*DCSZKZ
          BZH=-ZKERROR/YKERROR*B0ERROR*xxran*DCSXKX*DSHYKY*DSNZKZ

          AXH=B0ERROR*xxran/ZKERROR*                DCSXKX*DCHYKY*DSNZKZ
          AYH=B0ERROR*xxran/ZKERROR*XKERROR/YKERROR*DSNXKX*DSHYKY*DSNZKZ
          AZH=0.0d0

          BZOUT=-BXH
          BYOUT=BYH
          BXOUT=BZH

          AZOUT=-AXH
          AYOUT=AYH
          AXOUT=AZH

          BxOUT=bxout*byellmx
          ByOUT=byout*byellmx
          BzOUT=bzout*byellmx

          axOUT=axout*byellmx
          ayOUT=ayout*byellmx
          azOUT=azout*byellmx

        endif

        if (xxranz.ne.0.0d0) then

          X=XIN-XCENERR+TOTLEN2
          X=X-(IXz+0.5)*POLWID

          XKX=XKERROR*(-ZIN)
          YKY=YKERROR*YIN
c          ZKZ=ZKERROR*X
          ZKZ=ZKERROR*(X+polwid/2.0d0)

          DSNXKX=DSIN(XKX)
          DCSXKX=DCOS(XKX)
          DSHYKY=DSINH(YKY)
          DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
          DSNZKZ=DSIN(ZKZ)
          DCSZKZ=DCOS(ZKZ)

          BXH=-XKERROR/YKERROR*B0ERROR*xxranz*DSNXKX*DSHYKY*DCSZKZ
          BYH=                 B0ERROR*xxranz*DCSXKX*DCHYKY*DCSZKZ
          BZH=-ZKERROR/YKERROR*B0ERROR*xxranz*DCSXKX*DSHYKY*DSNZKZ

          AXH=B0ERROR*xxranz/ZKERROR*                DCSXKX*DCHYKY*DSNZKZ
          AYH=B0ERROR*xxranz/ZKERROR*XKERROR/YKERROR*DSNXKX*DSHYKY*DSNZKZ
          AZH=0.0d0

          BxOUT=bxout+bzh*bzellmx
          ByOUT=byout-bxh*bzellmx
          BzOUT=bzout-byh*bzellmx

          AxOUT=axout+azh*bzellmx
          AyOUT=ayout-axh*bzellmx
          AzOUT=azout-ayh*bzellmx

        endif !xxranz

      ENDIF !MODANA

9999  continue

      RETURN
      END
