*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.05/02 09/05/2018  11.13.36  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.67/04 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.63/05 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.63/02 24/01/2008  15.09.06  by  Michael Scheer
*CMZ :  2.52/09 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/01 30/06/2004  16.42.15  by  Michael Scheer
*CMZ :  2.51/00 25/05/2004  18.17.34  by  Michael Scheer
*CMZ :  2.50/00 29/04/2004  18.11.18  by  Michael Scheer
*CMZ :  2.47/05 16/04/2004  09.24.47  by  Michael Scheer
*CMZ :  2.16/08 25/10/2000  12.21.53  by  Michael Scheer
*CMZ :  2.16/07 14/09/2000  17.12.07  by  Michael Scheer
*CMZ :  2.16/05 28/07/2000  16.14.29  by  Michael Scheer
*CMZ :  2.16/04 20/07/2000  12.23.59  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.12/03 17/06/99  13.51.34  by  Michael Scheer
*CMZ :  2.10/01 18/03/99  11.30.11  by  Michael Scheer
*CMZ :  2.02/00 15/02/99  10.19.37  by  Michael Scheer
*CMZ :  2.00/00 05/01/99  14.45.11  by  Michael Scheer
*-- Author :    Michael Scheer   15/12/98
      SUBROUTINE SETNLPOI
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
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEND.

C--- SETNLPOI ESTIMATES VALUES FOR NLPOI
C--- FORMULA FROM X-RAY DATA BOOKLET

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER NLPOIMX,NLPERMN,NLPOIMN,ISOUR,NLSOUMN
      DOUBLE PRECISION B0DIP,YLOW,YHIGH,WGLOW,WGHIGH,ECDIP,PHIDEFL,
     &  DEFLEC,DLENG,SOULEN,DNPER,HARM,HARMWIDTH

      DATA NLPOIMX/2000000/,NLPERMN/100/,NLSOUMN/10/

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR SETNLPOI CALLED:'
      WRITE(LUNGFO,*)

        PHIDEFL=ANGRMS*SQRT(2.)
        B0DIP=BRMS*SQRT(2.)
        DEFLEC=PHIDEFL*DMYGAMMA

        IF (B0DIP.NE.0.D0) THEN
          DLENG=DEFLEC/93.4/B0DIP
        ELSE
          DLENG=0.D0
        ENDIF

        IF (DLENG.EQ.0.0) DLENG=(XMX-XMN)

        HARM=0.95D3*DMYENERGY**2/(1.D0+DEFLEC**2/2.D0)/DLENG/100.d0

        IF (DLENG.EQ.0.0) DLENG=(XMX-XMN)

        DNPER=(XMX-XMN)/DLENG

        IF (DNPER.LT.1.D0) DNPER=1.D0

        NLPOIMN=NLPERMN*DNPER*DEFLEC
        HARMWIDTH=HARM/DNPER

        IF (B0DIP.EQ.0.D0) THEN
          WRITE(LUNGFO,*)'*** WARNING IN SETNLPOI:  ***'
          WRITE(LUNGFO,*)'ZERO MAG. FIELD'
          WRITE(6,*)'*** WARNING IN SETNLPOI:  ***'
          WRITE(6,*)'ZERO MAG. FIELD'
          NLPOI=1000
          GOTO 9999
        ENDIF

        ECDIP=ecdipev1*DMYENERGY**2*B0DIP

        WRITE(LUNGFO,*)
     &    '     mag. field [T], Ec [eV]:',SNGL(B0DIP),SNGL(ECDIP)
        WRITE(LUNGFO,*)
     &    '     first harmonical [eV] (estimate):'
        WRITE(LUNGFO,*)
     &    '     ',SNGL(HARM)
     &    ,SNGL(HARM-HARMWIDTH),' -> ',SNGL(HARM+HARMWIDTH)
        WRITE(LUNGFO,*)
     &    '     K, lx, N of effective periodical device:'
        WRITE(LUNGFO,*)
     &    '     ',SNGL(DEFLEC),SNGL(DLENG),SNGL(DNPER)
        WRITE(LUNGFO,*)

        YLOW=FREQ(1)/ECDIP
        YHIGH=FREQ(NFREQ)/ECDIP

        IF (YLOW.LT.1.) THEN
          WGLOW=0.408/DMYENERGY*YLOW**(-0.354)/1000.
        ELSE
          WGLOW=0.408/DMYENERGY*YLOW**(-0.549)/1000.
        ENDIF

        IF (YHIGH.LT.1.) THEN
          WGHIGH=0.408/DMYENERGY*YHIGH**(-0.354)/1000.
        ELSE
          WGHIGH=0.408/DMYENERGY*YHIGH**(-0.549)/1000.
        ENDIF

        WRITE(LUNGFO,*)
     &    '      low and high size of radiaton cone [1/gamma]: '
        WRITE(LUNGFO,*)'     '
     &    ,SNGL(WGLOW*DMYGAMMA),SNGL(WGHIGH*DMYGAMMA)
        WRITE(LUNGFO,*)

        IF (WGLOW.GT.PHIDEFL) THEN
          WGLOW=PHIDEFL
        ENDIF

        IF (WGHIGH.GT.PHIDEFL) THEN
          WGHIGH=PHIDEFL
        ENDIF

        IF (NLPOI.EQ.-9999) THEN

          NLPOI=100*WGWINFC*WGLOW/WGHIGH
          IF (NLPOI.LT.NLPOIMN) THEN
            NLPOI=NLPOIMN
          ENDIF

          SOULEN=-1.D30
          DO ISOUR=1,NSOURCE
            IF (SOULEN.LT.SOURCEE(1,1,ISOUR)-SOURCEA(1,1,ISOUR)) THEN
              SOULEN=SOURCEE(1,1,ISOUR)-SOURCEA(1,1,ISOUR)
            ENDIF
          ENDDO

          IF (SOULEN.LE.0.D0) THEN
            WRITE(LUNGFO,*)'*** ERROR IN SETNLPOI: ZERO SOURCE LENGTH  ***'
            WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
            WRITE(6,*)'*** ERROR IN SETNLPOI: ZERO SOURCE LENGTH  ***'
            WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
            STOP
          ENDIF

          IF (NLPOI/SOULEN.LT.MYINUM) THEN
            NLPOI=NLSOUMN*MYINUM*SOULEN
          ENDIF

          IF (IUNDULATOR.NE.0) THEN
            NLPOI=NLPOIMN
            IF (NLPOI/SOULEN.LT.MYINUM) THEN
              NLPOI=MYINUM*SOULEN*1.1
            ENDIF
          ENDIF

          IF (NLPOI.GT.NLPOIMX) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** WARNING IN SETNLPOI: NLPOI VERY LARGE'
            WRITE(LUNGFO,*)'    NLPOI: ',NLPOI
            WRITE(LUNGFO,*)'    NLPOI limited to: ',NLPOIMX
            WRITE(6,*)
            WRITE(6,*)
            WRITE(6,*)'*** WARNING IN SETNLPOI: NLPOI VERY LARGE'
            WRITE(6,*)'    NLPOI: ',NLPOI
            WRITE(6,*)'    NLPOI limited to: ',NLPOIMX
            WRITE(6,*)
            NLPOI=NLPOIMX
          ENDIF

        ENDIF

9999  WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'      NLPOI set to: ',NLPOI
      WRITE(LUNGFO,*)

      RETURN
      END
