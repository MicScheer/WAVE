*CMZ :  4.00/04 17/05/2019  11.46.54  by  Michael Scheer
*CMZ :  2.66/09 08/04/2010  12.10.57  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.48/03 03/03/2004  12.49.38  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 11/06/98  10.26.40  by  Michael Scheer
*CMZ : 00.01/04 16/01/95  13.52.05  by  Michael Scheer
*CMZ : 00.01/03 28/11/94  12.04.43  by  Michael Scheer
*CMZ :  0.00/03 15/11/94  16.48.07  by  Michael Scheer
*CMZ :  0.00/02 02/11/94  16.15.54  by  Michael Scheer
*CMZ :  0.00/01 31/10/94  10.05.04  by  Michael Scheer
*CMZ :  0.00/00 28/10/94  16.14.51  by  Michael Scheer
*-- Author :    Michael Scheer   28/10/94
      SUBROUTINE BPOLYHARM(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)
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

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      INTEGER NORDP
      PARAMETER (NORDP=20)

      CHARACTER(50) FILECOE

      INTEGER IPOLY
      INTEGER ICAL
      INTEGER LUNCOE
      INTEGER NFIRSTX,NORDX,NSTEPX
      INTEGER NFIRSTY,NORDY,NSTEPY
      INTEGER IFHALBA

      DATA ICAL/0/

      DATA LUNCOE/10/
      DATA FILECOE/'wave_bpolyharm_coef.dat'/

      DOUBLE PRECISION X,Y,Z
      DOUBLE PRECISION XIN,YIN,ZIN
      DOUBLE PRECISION Q(NORDP,NORDP),QA0(NORDP),QA(NORDP,NORDP)
      DOUBLE PRECISION BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION BX,BY,BZ

      DOUBLE PRECISION XLX,YLY,ZLZ,XKX,YKY,ZKZ
      DOUBLE PRECISION GAP2PI,WIDTH

      IF (ICAL.EQ.0) THEN


       OPEN(UNIT=LUNCOE,FILE=FILECOE,STATUS='OLD')

         READ(LUNCOE,*)IFHALBA,IPOLY
       CLOSE(LUNCOE)

       CALL BRCOEF(LUNCOE,FILECOE
     &            ,XLX,YLY,ZLZ
     &            ,XKX,YKY,ZKZ
     &            ,NFIRSTX,NORDX,NSTEPX
     &            ,NFIRSTY,NORDY,NSTEPY
     &            ,NORDP,Q,QA0,QA,IFHALBA,GAP2PI,WIDTH)

      IF (XSTART.EQ.9999.) XSTART=-ZLZ/2.
      IF (XSTOP.EQ.9999.) XSTOP=+ZLZ/2.

       WRITE(LUNGFO,*)
       WRITE(LUNGFO,*)
     &'*** WARNING SR BPOLYHARM: VECTOR POTENTIAL NOT AVAILABLE (SET TO ZERO) ***'
       WRITE(LUNGFO,*)
       WRITE(6,*)
       WRITE(6,*)
     &'*** WARNING SR BPOLYHARM: VECTOR POTENTIAL NOT AVAILABLE (SET TO ZERO) ***'
       WRITE(6,*)
       WRITE(LUNGFO,*)
       WRITE(LUNGFO,*)'     SR BPOLYHARM:'
       WRITE(LUNGFO,*)
     &'     COEFFICIENTS READ FROM FILE wave_bpolyharm_coef.dat'
       WRITE(LUNGFO,*)
     &'     POLYHARM RUN NUMBER ON FILE:',IPOLY
       WRITE(LUNGFO,*)
       CALL wave_print_file(lungfo,'wave_bpolyharm_coef.dat')
       WRITE(LUNGFO,*)
          ICAL=1
      ENDIF !ICAL

      AXOUT=0.0
      AYOUT=0.0
      AZOUT=0.0

      X=-ZIN
      Y=YIN
        Z=XIN

            CALL BHARM(X,Y,Z,BX,BY,BZ
     &                ,NFIRSTX,NORDX,NSTEPX
     &                ,NFIRSTY,NORDY,NSTEPY
     &                ,Q,QA0,QA,NORDP,XKX,YKY,ZKZ,IFHALBA,GAP2PI,WIDTH)

      BXOUT=BZ
      BYOUT=BY
      BZOUT=-BX

      RETURN
      END
