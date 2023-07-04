*CMZ :  3.05/01 04/05/2018  19.58.00  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.50/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.20/01 29/11/2000  18.31.06  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.15/00 04/05/2000  17.02.01  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.26.26  by  Michael Scheer
*CMZ :  2.10/01 24/02/99  10.20.40  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.17.32  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.26  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.48  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WFILINT
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
*KEND.

C--- DUMP INTEGRAND

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEND.

      INTEGER IPOI,IC,ICAL


      COMPLEX*16 AX,AY,AZ,EXPOM
      DOUBLE PRECISION OM

      DATA ICAL/0/

      OM=FREQ(1)/HBAREV1

C-- LOOP OVER TIME STEPS (ACTUAL INTEGRATION)

      IF (ICAL.EQ.0) THEN  !CV2
C         OPEN(UNIT=LUNINT,FILE=FILEINT,STATUS='NEW')
        ICAL=1
      ENDIF !ICAL     CV2

      DO IPOI=1,IARGUM

        IF (WSOU(1,1,IPOI).GE.XIANF.AND.WSOU(1,1,IPOI).LE.XIEND) THEN

          EXPOM=CDEXP(DCMPLX(0.D0,REARGUM(4,IPOI)*OM))

          AX=DCMPLX(REARGUM(1,IPOI))*EXPOM
          AY=DCMPLX(REARGUM(2,IPOI))*EXPOM
          AZ=DCMPLX(REARGUM(3,IPOI))*EXPOM

          WRITE(LUNINT,*) WSOU(1,1,IPOI)
     &      ,(SNGL(REARGUM(IC,IPOI)),IC=1,3)
     &      ,SNGL(REARGUM(4,IPOI)*OM),SNGL(REARGUM(5,IPOI))
     &      ,SNGL(DREAL(EXPOM)),SNGL(DIMAG(EXPOM))
     &      ,SNGL(DREAL(AX)),SNGL(DIMAG(AX))
     &      ,SNGL(DREAL(AY)),SNGL(DIMAG(AY))
     &      ,SNGL(DREAL(AZ)),SNGL(DIMAG(AZ))

        ENDIF   !XIANF

      ENDDO   !LOOP OVER TIME STEPS

C      IF (NSADD.NE.0) CLOSE(LUNINT)

CV2   CLOSE(LUNINT)

      RETURN
      END
