*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  3.01/00 17/06/2013  09.14.35  by  Michael Scheer
*CMZ :  2.67/04 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.16/08 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.16/04 19/07/2000  11.26.53  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  2.00/00 05/01/99  14.47.39  by  Michael Scheer
*-- Author :    Michael Scheer   15/12/98
      SUBROUTINE SETWGWIN
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

C--- ESTIMATES VALUE FOR WGWINFC
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
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION B0DIP,YLOW,WGLOW,ECDIP,PHIDEFL,WGEXPAND

      DATA WGEXPAND/20.D0/

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR SETWGWIN CALLED:'
      WRITE(LUNGFO,*)

C        PHIDEFL=(PHIMX-PHIMN)/2.
C     B0DIP=DSQRT(DABS(BMAXGL2))
        PHIDEFL=ANGRMS*SQRT(2.)
      B0DIP=BRMS*SQRT(2.)

      IF (B0DIP.EQ.0.D0) THEN
          WRITE(LUNGFO,*)'*** ERROR IN SETWGWIN: B0DIP.EQ.0  ***'
          WRITE(LUNGFO,*)'CHECK MAG. FIELD'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(6,*)'*** ERROR IN SETWGWIN: B0DIP.EQ.0  ***'
          WRITE(6,*)'CHECK MAG. FIELD'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
      ENDIF

      ECDIP=ecdipev1*DMYENERGY**2*B0DIP

      WRITE(LUNGFO,*)
     &'     mag. field [T], Ec [eV]:',SNGL(B0DIP),SNGL(ECDIP)
      WRITE(LUNGFO,*)

      YLOW=FREQ(1)/ECDIP

      IF (YLOW.LT.1.) THEN
          WGLOW=0.408/DMYENERGY*YLOW**(-0.354)/1000.
      ELSE
          WGLOW=0.408/DMYENERGY*YLOW**(-0.549)/1000.
      ENDIF

      IF (WGLOW*WGEXPAND.GT.PHIDEFL) THEN
          WGLOW=PHIDEFL
            WGWINFC=PHIDEFL*DMYGAMMA*1.2  !1.2 TO HAVE SOME OVERLAP
            WRITE(LUNGFO,*)
     &'    WGWINFC corresponds to angle greater than max. deflection'
      ENDIF

      IF (WGWINFC.EQ.9999.) THEN
          WGWINFC=WGLOW*WGEXPAND*DMYGAMMA
      ENDIF

      WRITE(LUNGFO,*)'      WGWINFC set to: ',WGWINFC
      WRITE(LUNGFO,*)

      RETURN
      END
